"""Download, verify and build the GDAL C library from source.

This is the last resort of :mod:`TopoAnalysis.gdal_setup`, for machines with
no conda, no system package manager and no administrator rights.  Everything
here installs into a single prefix that the caller chooses -- normally the
active virtualenv -- so nothing outside it is touched and removing the prefix
undoes the whole thing.

Three components, built in order, because each needs the one before it:

    sqlite3   only when the system copy is unusable.  PROJ needs both the
              library *and* the ``sqlite3`` command, which it runs during its
              own build to assemble ``proj.db``.
    PROJ      mandatory for GDAL 3.x.  Built without TIFF and curl support,
              which are only needed for downloading and reading the optional
              high-accuracy datum grids.
    GDAL      built against its own bundled copies of libtiff, libgeotiff,
              libpng, libjpeg and zlib (``GDAL_USE_INTERNAL_LIBS``), so the
              only external dependency left is PROJ.  That is what makes an
              unattended build realistic: no hunting for system ``-dev``
              packages.  The result reads GeoTIFF, ArcInfo ASCII, ENVI, HFA
              and the rest of the built-in drivers, which is everything
              TopoAnalysis uses.

The Python bindings are deliberately *not* built here.  ``cmake`` would build
them for whichever interpreter it happens to find; :mod:`TopoAnalysis.gdal_setup`
installs them with pip instead, against the interpreter we actually care
about.

Every tarball is pinned by SHA-256.  The pins below were taken from bytes
checked against what the projects themselves publish -- GDAL's and PROJ's
``.md5`` sidecars, SQLite's ``PRODUCT`` line -- so a download that does not
match them is a download that should not be built.
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

__all__ = [
    "SOURCES",
    "apply_environment_overrides",
    "build_one",
    "components_needed",
    "default_cache_dir",
    "missing_build_tools",
    "system_sqlite_is_usable",
]


class BuildError(RuntimeError):
    """A download, configure, compile or install step failed."""


# ---------------------------------------------------------------------------
# What we build, and exactly which bytes
# ---------------------------------------------------------------------------

SOURCES = {
    "sqlite": {
        "version": "3.53.4",
        "url": "https://www.sqlite.org/2026/sqlite-autoconf-3530400.tar.gz",
        "sha256": "0e9483900e92cd5de8fd48d16bf9200145a61f7fd5be542a5ac81d8a9516eb9c",
        "minutes": 2,
    },
    "proj": {
        "version": "9.4.1",
        "url": "https://download.osgeo.org/proj/proj-9.4.1.tar.gz",
        "sha256": "ffe20170ee2b952207adf8a195e2141eab12cda181e49fdeb54425d98c7171d7",
        "minutes": 6,
    },
    "gdal": {
        # Not merely "a recent GDAL": the version matters for this build.
        # GDAL bundles its own zlib, libpng, libtiff and friends, and the
        # copies shipped up to about 3.9 do not compile against a current
        # Apple SDK -- both zlib's ``zutil.h`` and libpng's ``pngpriv.h``
        # take their "Mac OS classic" branch, because that is selected on
        # ``TARGET_OS_MAC``, which modern SDKs define on every Apple
        # platform.  zlib then redefines ``fdopen`` to NULL and the system
        # ``_stdio.h`` fails to parse.  Later GDAL releases carry repaired
        # copies, so the pin has to stay reasonably current.
        "version": "3.13.3",
        "url": "https://github.com/OSGeo/gdal/releases/download/v3.13.3/gdal-3.13.3.tar.gz",
        "sha256": "5e0c388d83da2d686cc00a40272882432cdb54edff43d4af173e532844a0a0ea",
        "minutes": 20,
    },
}

#: Order matters: each entry links against the ones before it.
BUILD_ORDER = ("sqlite", "proj", "gdal")

#: Seconds to wait for any single read from a download before giving up.
DOWNLOAD_TIMEOUT = 60.0

#: The floor declared by both pinned projects: PROJ 9.4's CMakeLists and
#: GDAL 3.13's ``GdalCMakeMinimumRequired.cmake`` each say 3.16.
MIN_CMAKE_VERSION = (3, 16)


def _cmake_version(cmake):
    """``(major, minor)`` reported by ``cmake --version``, or None."""
    try:
        result = subprocess.run(
            [cmake, "--version"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    text = result.stdout.decode("utf8", "replace")
    match = re.search(r"(\d+)\.(\d+)", text)
    return (int(match.group(1)), int(match.group(2))) if match else None


def apply_environment_overrides():
    """Let ``TOPOANALYSIS_GDAL_SOURCE_*`` environment variables move a pin.

    A pinned version is a dead end the day it stops being suitable -- a
    platform that needs a newer GDAL, a site that mirrors its own tarballs.
    Overriding one means supplying its checksum too, so the integrity check
    cannot be turned off by accident::

        TOPOANALYSIS_GDAL_SOURCE_GDAL_VERSION=3.12.4
        TOPOANALYSIS_GDAL_SOURCE_GDAL_URL=https://example.org/gdal-3.12.4.tar.gz
        TOPOANALYSIS_GDAL_SOURCE_GDAL_SHA256=<64 hex digits>
    """
    for component in SOURCES:
        prefix = "TOPOANALYSIS_GDAL_SOURCE_%s_" % component.upper()
        url = os.environ.get(prefix + "URL")
        sha256 = os.environ.get(prefix + "SHA256")
        version = os.environ.get(prefix + "VERSION")
        if url or sha256:
            if not (url and sha256):
                raise BuildError(
                    "%sURL and %sSHA256 must be set together" % (prefix, prefix)
                )
            SOURCES[component]["url"] = url
            SOURCES[component]["sha256"] = sha256
        if version:
            SOURCES[component]["version"] = version


def default_cache_dir():
    """Where downloaded tarballs are kept between runs."""
    root = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache"
    )
    return os.path.join(root, "TopoAnalysis", "downloads")


# ---------------------------------------------------------------------------
# Prerequisites
# ---------------------------------------------------------------------------


def missing_build_tools():
    """Names of build tools that are needed for a source build and absent.

    Returns an empty list when the machine can build.  ``ninja`` is not
    required -- CMake falls back to makefiles -- so it is not listed.
    """
    missing = []
    cmake = shutil.which("cmake")
    if cmake is None:
        missing.append("cmake")
    else:
        version = _cmake_version(cmake)
        # Both PROJ 9 and GDAL 3.13 declare this floor.  Checking it here
        # turns a puzzling failure part-way through `cmake -S` into a
        # sentence naming the version that is needed.
        if version is not None and version < MIN_CMAKE_VERSION:
            missing.append(
                "cmake %d.%d or newer (found %s)"
                % (MIN_CMAKE_VERSION[0], MIN_CMAKE_VERSION[1],
                   ".".join(str(part) for part in version))
            )
    # ninja is not an alternative here: no -G is passed, so CMake generates
    # makefiles, and SQLite's tarball is a plain autoconf build that calls
    # make directly.
    if shutil.which("make") is None:
        missing.append("make")
    if _c_compiler() is None:
        missing.append("a C compiler")
    # PROJ and GDAL are C++17 projects, and a machine can genuinely have cc
    # without a C++ compiler -- a bare gcc package without g++ is the usual
    # way.  Checking only for a C compiler would let such a machine start a
    # build that CMake aborts a few seconds later.
    if _cxx_compiler() is None:
        missing.append("a C++ compiler")
    return missing


def _find_compiler(variable, candidates):
    """The compiler named by ``variable``, else the first of ``candidates``."""
    explicit = os.environ.get(variable)
    if explicit:
        # The variable may carry arguments ("ccache gcc"); only the program
        # matters here.
        return shutil.which(explicit.split()[0])
    for name in candidates:
        found = shutil.which(name)
        if found is not None:
            return found
    return None


def _c_compiler():
    """The C compiler to probe with, honouring ``CC`` as the toolchain does."""
    return _find_compiler("CC", ("cc", "gcc", "clang"))


def _cxx_compiler():
    """The C++ compiler, honouring ``CXX``."""
    return _find_compiler("CXX", ("c++", "g++", "clang++"))


def system_sqlite_is_usable():
    """True when the system SQLite can satisfy PROJ's build.

    PROJ needs three things, and a machine can easily have one or two of
    them: the shared library, the ``sqlite3.h`` header (often packaged
    separately as ``libsqlite3-dev``), and the ``sqlite3`` command-line shell,
    which PROJ runs to build ``proj.db``.  Checking the header by compiling
    against it is the only honest test -- the file can be present but
    unreachable from the compiler's default include path, which is the normal
    situation inside a conda environment.
    """
    if shutil.which("sqlite3") is None:
        return False
    compiler = _c_compiler()
    if compiler is None:
        return False

    program = (
        "#include <sqlite3.h>\n"
        "#if SQLITE_VERSION_NUMBER < 3011000\n"
        '#error "PROJ needs SQLite 3.11 or newer"\n'
        "#endif\n"
        "int main(void) { return sqlite3_libversion_number() > 0 ? 0 : 1; }\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        source = os.path.join(tmp, "probe.c")
        with open(source, "w") as handle:
            handle.write(program)
        # Compile only: linking would additionally need -lsqlite3 to resolve,
        # and a missing library is a different (rarer) problem than a missing
        # header.  PROJ's CMake will find the library if the header is there.
        result = subprocess.run(
            [compiler, "-c", source, "-o", os.path.join(tmp, "probe.o")],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    return result.returncode == 0


def components_needed():
    """Which of :data:`BUILD_ORDER` this machine actually has to build."""
    needed = list(BUILD_ORDER)
    if system_sqlite_is_usable():
        needed.remove("sqlite")
    return needed


# ---------------------------------------------------------------------------
# Download and unpack
# ---------------------------------------------------------------------------


def _download(url, sha256, cache_dir, log):
    """Fetch ``url`` into ``cache_dir``, verifying its SHA-256.

    A cached file that fails the check is deleted and fetched again: that is
    far more often a truncated download than anything sinister, and leaving it
    in place would make every later run fail the same way.
    """
    os.makedirs(cache_dir, exist_ok=True)
    target = os.path.join(cache_dir, url.rsplit("/", 1)[-1])

    if os.path.exists(target):
        if _sha256(target) == sha256:
            log("    using cached %s" % os.path.basename(target))
            return target
        log("    cached %s failed its checksum; re-downloading" % os.path.basename(target))
        os.remove(target)

    log("    downloading %s" % url)
    partial = target + ".part"
    try:
        # With no timeout a stalled connection blocks the install for ever,
        # with nothing on screen to say why.  This bounds each individual
        # read, not the whole transfer, so a slow-but-alive mirror still
        # works.
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT) as response:
            with open(partial, "wb") as handle:
                shutil.copyfileobj(response, handle)
    except Exception as exc:  # network, DNS, TLS, HTTP -- all equally fatal here
        if os.path.exists(partial):
            os.remove(partial)
        raise BuildError("could not download %s: %s" % (url, exc))

    got = _sha256(partial)
    if got != sha256:
        os.remove(partial)
        raise BuildError(
            "checksum mismatch for %s\n  expected %s\n  got      %s\n"
            "Refusing to build code that is not the code we pinned." % (url, sha256, got)
        )
    os.replace(partial, target)
    return target


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract(tarball, into):
    """Unpack ``tarball`` and return the directory it created."""
    with tarfile.open(tarball) as archive:
        names = archive.getnames()
        # ``filter="data"`` rejects absolute paths, parent-directory escapes,
        # links out of the tree and device nodes.  It is the default from
        # Python 3.14; ask for it explicitly wherever it exists.
        try:
            archive.extractall(into, filter="data")
        except TypeError:  # pragma: no cover - Python < 3.11.4
            _check_members(archive, names)
            archive.extractall(into)
    roots = {name.split("/", 1)[0] for name in names if name and not name.startswith("/")}
    if len(roots) != 1:
        raise BuildError("unexpected archive layout in %s" % tarball)
    return os.path.join(into, roots.pop())


def _check_members(archive, names):  # pragma: no cover - old Python only
    """The safety half of ``filter='data'`` for interpreters that lack it."""
    for name in names:
        if name.startswith("/") or os.path.isabs(name) or ".." in name.split("/"):
            raise BuildError("refusing to extract unsafe path %r" % name)
    for member in archive.getmembers():
        if member.isdev():
            raise BuildError("refusing to extract device node %r" % member.name)


# ---------------------------------------------------------------------------
# Running the build
# ---------------------------------------------------------------------------


def _build_env(prefix):
    """Environment that makes each component find the ones built before it."""
    env = dict(os.environ)
    bindir = os.path.join(prefix, "bin")
    env["PATH"] = bindir + os.pathsep + env.get("PATH", "")

    libdir = os.path.join(prefix, "lib")
    pkgconfig = os.pathsep.join(
        [os.path.join(libdir, "pkgconfig"), os.path.join(prefix, "share", "pkgconfig")]
    )
    existing = env.get("PKG_CONFIG_PATH")
    env["PKG_CONFIG_PATH"] = pkgconfig + (os.pathsep + existing if existing else "")

    # So the freshly built tools (PROJ's own, and sqlite3) can run during the
    # build before anything has been relinked.
    loader = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    env[loader] = libdir + (os.pathsep + env[loader] if env.get(loader) else "")
    return env


def _run(command, cwd, env, log, logfile):
    """Run a build command, tee-ing its output to ``logfile``."""
    log("    $ %s" % " ".join(command))
    with open(logfile, "a") as handle:
        handle.write("\n$ %s\n" % " ".join(command))
        handle.flush()
        result = subprocess.run(
            command, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT
        )
    if result.returncode != 0:
        raise BuildError(
            "%s failed (exit %d).\nThe full build log is at:\n  %s"
            % (command[0], result.returncode, logfile)
        )


def _jobs(requested):
    if requested:
        return max(1, int(requested))
    return max(1, (os.cpu_count() or 2))


def _cmake_common(prefix):
    """CMake settings shared by PROJ and GDAL."""
    settings = [
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX=" + prefix,
        "-DCMAKE_PREFIX_PATH=" + prefix,
        "-DBUILD_TESTING=OFF",
        "-DBUILD_SHARED_LIBS=ON",
        # Pin the library directory instead of letting GNUInstallDirs choose.
        # On Fedora, RHEL and openSUSE it defaults to lib64 on 64-bit hosts,
        # and everything else here -- the RPATH below, the loader path in
        # _build_env, the -L handed to the bindings build -- assumes one
        # known location.  A private prefix has no distro convention to
        # honour, so "lib" everywhere is both simpler and correct.
        "-DCMAKE_INSTALL_LIBDIR=lib",
        # Without an install RPATH the shared libraries only resolve each
        # other when LD_LIBRARY_PATH happens to be set, which it will not be
        # in the user's next shell.  Baking the path in makes the prefix work
        # on its own.
        "-DCMAKE_INSTALL_RPATH=" + os.path.join(prefix, "lib"),
        "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON",
        "-DCMAKE_INSTALL_NAME_DIR=" + os.path.join(prefix, "lib"),
    ]
    if sys.platform == "darwin":
        # Build for the architecture the *interpreter* runs as.  On an Apple
        # Silicon machine running an x86_64 Python under Rosetta, the default
        # would be arm64 and the bindings would refuse to link.
        settings.append("-DCMAKE_OSX_ARCHITECTURES=" + platform.machine())
    return settings


def build_one(component, prefix, jobs=None, log=None, cache_dir=None, keep_build=False):
    """Download, build and install one component into ``prefix``.

    Parameters
    ----------
    component : {'sqlite', 'proj', 'gdal'}
    prefix : str
        Installation prefix.  Created if it does not exist.
    jobs : int, optional
        Parallel compile jobs; defaults to the CPU count.
    log : callable, optional
        Receives progress lines.  Defaults to printing.
    cache_dir : str, optional
        Where tarballs are cached; defaults to :func:`default_cache_dir`.
    keep_build : bool
        Leave the unpacked source tree behind, for debugging a failure.

    Returns
    -------
    str
        The path of the build log, which is worth showing the user whether or
        not the build succeeded.
    """
    if component not in SOURCES:
        raise ValueError("unknown component %r" % component)
    log = log or (lambda line: print(line, flush=True))
    cache_dir = cache_dir or default_cache_dir()
    apply_environment_overrides()
    spec = SOURCES[component]
    prefix = os.path.abspath(prefix)
    os.makedirs(prefix, exist_ok=True)

    logdir = os.path.join(prefix, "var", "log", "topoanalysis")
    os.makedirs(logdir, exist_ok=True)
    logfile = os.path.join(logdir, "build-%s.log" % component)
    # Start each attempt with a clean log so a failure message points at the
    # failure, not at the tail of some earlier successful run.
    open(logfile, "w").close()

    log("  %s %s" % (component, spec["version"]))
    tarball = _download(spec["url"], spec["sha256"], cache_dir, log)

    workdir = tempfile.mkdtemp(prefix="topoanalysis-%s-" % component)
    try:
        source = _extract(tarball, workdir)
        env = _build_env(prefix)
        builder = {
            "sqlite": _build_sqlite,
            "proj": _build_proj,
            "gdal": _build_gdal,
        }[component]
        builder(source, prefix, env, _jobs(jobs), log, logfile)
        log("  %s installed into %s" % (component, prefix))
    finally:
        if keep_build:
            log("    build tree kept at %s" % workdir)
        else:
            shutil.rmtree(workdir, ignore_errors=True)
    return logfile


def _build_sqlite(source, prefix, env, jobs, log, logfile):
    """SQLite's amalgamation tarball: plain configure/make/install.

    Only flags that both of SQLite's build systems understand are used here.
    Releases up to 3.48 shipped classic autoconf; 3.49 and later ship
    autosetup, which has no ``--enable-shared`` at all (shared is the default
    and only the negative form exists).  Passing the positive form would make
    the newer releases fail on an unknown option.
    """
    _run(
        [
            os.path.join(source, "configure"),
            "--prefix=" + prefix,
            "--disable-static",
            # PROJ needs the core library and the shell, nothing else.
            "--disable-readline",
        ],
        source,
        env,
        log,
        logfile,
    )
    _run(["make", "-j%d" % jobs], source, env, log, logfile)
    _run(["make", "install"], source, env, log, logfile)


def _build_proj(source, prefix, env, jobs, log, logfile):
    build = os.path.join(source, "build")
    os.makedirs(build, exist_ok=True)
    settings = _cmake_common(prefix) + [
        # Both are only used for the optional high-accuracy grid files, which
        # need a network at run time anyway.  Leaving them out removes the
        # need for libtiff and libcurl development packages.
        "-DENABLE_TIFF=OFF",
        "-DENABLE_CURL=OFF",
        # The command-line tools are not needed by GDAL, but projinfo is
        # useful enough when diagnosing a projection problem to keep.
        "-DBUILD_PROJSYNC=OFF",
        "-DBUILD_TESTING=OFF",
    ]
    _run(["cmake", "-S", source, "-B", build] + settings, source, env, log, logfile)
    _run(["cmake", "--build", build, "-j", str(jobs)], source, env, log, logfile)
    _run(["cmake", "--install", build], source, env, log, logfile)


def _build_gdal(source, prefix, env, jobs, log, logfile):
    build = os.path.join(source, "build")
    os.makedirs(build, exist_ok=True)
    settings = _cmake_common(prefix) + [
        # The point of the whole exercise: use GDAL's vendored copies of
        # libtiff, libgeotiff, libpng, libjpeg, libz and friends rather than
        # requiring system development packages we have no way to install.
        "-DGDAL_USE_INTERNAL_LIBS=ON",
        "-DGDAL_USE_EXTERNAL_LIBS=OFF",
        # PROJ is the exception: it is mandatory and we just built it.
        "-DGDAL_USE_PROJ=ON",
        "-DPROJ_ROOT=" + prefix,
        # Bindings are pip's job, against the interpreter we mean to use.
        "-DBUILD_PYTHON_BINDINGS=OFF",
        "-DBUILD_JAVA_BINDINGS=OFF",
        "-DBUILD_CSHARP_BINDINGS=OFF",
        "-DGDAL_BUILD_OPTIONAL_DRIVERS=ON",
        "-DOGR_BUILD_OPTIONAL_DRIVERS=ON",
    ]
    _run(["cmake", "-S", source, "-B", build] + settings, source, env, log, logfile)
    _run(["cmake", "--build", build, "-j", str(jobs)], source, env, log, logfile)
    _run(["cmake", "--install", build], source, env, log, logfile)
