"""Find, diagnose and install GDAL -- both the C library and the Python bindings.

Why this module exists
----------------------
GDAL is the one dependency ``pip install`` cannot supply by itself.  There is
no GDAL wheel on PyPI: at *every* version, for *every* platform, the only
artifact published is an sdist whose ``setup.py`` shells out to
``gdal-config`` and compiles the SWIG bindings against a GDAL C library that
has to be on the machine already.  So "install GDAL" is really two jobs that
must agree with each other:

    1. the C library, ``libgdal``, which pip knows nothing about, and
    2. the Python bindings, ``osgeo``, whose version must match it.

This module does both.  It probes what is present, works out the cheapest
route to a working install, and takes it.  In rough order of preference:

    ``bindings``  the C library is already here -- just build ``osgeo``
                  against it.  Seconds.
    ``conda``     a conda environment is active; conda-forge ships both
                  halves, already matched.  A minute or two.
    ``wheel``     Windows only, where prebuilt wheels that bundle the C
                  library are the only practical option.
    ``system``    apt/dnf/pacman/zypper/apk/brew supply the C library, then
                  fall through to ``bindings``.  Needs administrator rights
                  except under Homebrew.
    ``source``    download and build PROJ and GDAL from pinned, checksummed
                  tarballs into a prefix of our own, then fall through to
                  ``bindings``.  Half an hour, but needs nothing but a
                  compiler and CMake -- no root, no package manager.

Two traps this module exists to avoid
-------------------------------------
**The silent ``gdal_array`` trap.**  GDAL sdists only grew a
``pyproject.toml`` in 3.9.  Before that they declare no build dependencies at
all, so under pip's default build isolation numpy is absent when the bindings
compile -- and the numpy bridge, ``osgeo.gdal_array``, is quietly skipped.
Everything imports; ``Band.ReadAsArray()`` then fails at run time with an
unhelpful error.  Ubuntu 22.04 ships GDAL 3.4 and 24.04 ships 3.8, so this is
the common case, not the exotic one.  We install the bindings with
``--no-build-isolation`` against an interpreter that already has numpy, which
both fixes the trap and builds the bridge against the numpy the user will
actually run.

**The architecture trap.**  A Homebrew GDAL on Apple Silicon is arm64; a
Python installed under Rosetta is x86_64.  The bindings cannot link the two,
and the error names a missing symbol rather than the real cause.  The probe
checks for this before anything is built.

Command line
------------
::

    topoanalysis-gdal             # what is installed, and what is wrong
    topoanalysis-gdal plan        # what an install would do, step by step
    topoanalysis-gdal install     # do it
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import _gdal_source
except ImportError:  # pragma: no cover
    import _gdal_source

__all__ = [
    "Plan",
    "Probe",
    "Step",
    "activate_environment",
    "apply_plan",
    "ensure",
    "main",
    "make_plan",
    "probe",
]


#: The oldest GDAL whose C API TopoAnalysis is known to work against.
MIN_LIBRARY_VERSION = (3, 0)

#: GDAL sdists from this version on declare numpy in ``[build-system]
#: requires``.  Below it, build isolation leaves numpy out and the
#: ``gdal_array`` bridge is silently dropped.  See the module docstring.
NUMPY_BUILD_REQUIRES_FROM = (3, 9)

#: Fallback when the GitHub API cannot be reached; see :func:`_windows_find_links`.
CGOHLKE_FALLBACK_TAG = "v2026.8.20"


# ---------------------------------------------------------------------------
# Probing
# ---------------------------------------------------------------------------

# Run in a subprocess, never in-process.  A bindings/library version mismatch
# or an architecture mismatch can abort the interpreter outright rather than
# raising ImportError, and a diagnostic tool that dies while diagnosing is
# worse than useless.
_PROBE_SCRIPT = r"""
import json, os, platform, sys
out = {"executable": sys.executable, "python": sys.version.split()[0]}
# The *target's* architecture, not ours.  An x86_64 Python under Rosetta on
# an arm64 machine is the whole reason the architecture check exists, and
# asking the wrong process gives exactly the wrong answer.
out["machine"] = platform.machine()
# Asked here rather than inferred from the parent process: --python can point
# at a different interpreter, and it is *its* environment we install into.
out["sys_prefix"] = sys.prefix
out["base_prefix"] = getattr(sys, "base_prefix", sys.prefix)

# The same data-file fixup TopoAnalysis applies at import time, repeated here
# rather than imported: the interpreter being probed very often does not have
# TopoAnalysis installed yet.  Without it the probe would report a broken
# GDAL that works perfectly well in use.
def _activate():
    if os.environ.get("TOPOANALYSIS_NO_GDAL_ENV"):
        return
    root = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache")
    for path in (os.path.join(sys.prefix, "share", "topoanalysis", "gdal-prefix.json"),
                 os.path.join(root, "TopoAnalysis", "gdal-prefix.json")):
        try:
            with open(path) as handle:
                manifest = json.load(handle)
        except Exception:
            continue
        recorded = manifest.get("python_prefix")
        if recorded is not None and recorded != sys.prefix:
            return
        proj_lib = manifest.get("proj_lib")
        if proj_lib and os.path.isfile(os.path.join(proj_lib, "proj.db")):
            os.environ["PROJ_LIB"] = proj_lib
            os.environ["PROJ_DATA"] = proj_lib
        gdal_data = manifest.get("gdal_data")
        if gdal_data and os.path.isdir(gdal_data):
            os.environ["GDAL_DATA"] = gdal_data
        return

try:
    _activate()
except Exception:
    pass
try:
    import numpy
    out["numpy"] = numpy.__version__
except Exception:
    out["numpy"] = None
try:
    from osgeo import gdal
    out["bindings_version"] = gdal.__version__
    out["osgeo_file"] = getattr(gdal, "__file__", None)
except Exception as exc:
    out["bindings_error"] = "%s: %s" % (type(exc).__name__, exc)
try:
    from osgeo import gdal_array  # noqa: F401
    out["has_gdal_array"] = True
except Exception as exc:
    out["has_gdal_array"] = False
    out["gdal_array_error"] = "%s: %s" % (type(exc).__name__, exc)
out["proj_env"] = {
    name: os.environ.get(name)
    for name in ("PROJ_LIB", "PROJ_DATA", "GDAL_DATA")
    if os.environ.get(name)
}
try:
    from osgeo import osr
    try:
        out["proj_search_paths"] = list(osr.GetPROJSearchPaths() or [])
    except Exception:
        pass
    # Without this a failed lookup returns an error code rather than raising,
    # and the diagnosis loses the one message that explains it -- PROJ names
    # the offending proj.db and why it was rejected.
    osr.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    out["projection_works"] = bool(srs.ExportToWkt())
    if not out["projection_works"]:
        out["projection_error"] = "EPSG:4326 produced an empty definition"
except Exception as exc:
    out["projection_works"] = False
    out["projection_error"] = "%s: %s" % (type(exc).__name__, exc)
sys.stdout.write(json.dumps(out))
"""


class Probe(object):
    """A snapshot of the GDAL situation on one interpreter."""

    def __init__(self, **fields):
        self.python = fields.get("python_executable") or sys.executable
        self.python_version = fields.get("python")
        self.numpy_version = fields.get("numpy")
        self.bindings_version = fields.get("bindings_version")
        self.bindings_error = fields.get("bindings_error")
        self.osgeo_file = fields.get("osgeo_file")
        self.has_gdal_array = bool(fields.get("has_gdal_array"))
        self.gdal_array_error = fields.get("gdal_array_error")
        self.projection_works = bool(fields.get("projection_works"))
        self.projection_error = fields.get("projection_error")
        self.library_version = fields.get("library_version")
        self.gdal_config = fields.get("gdal_config")
        self.prefix = fields.get("prefix")
        self.crashed = bool(fields.get("crashed"))
        self.problems = tuple(fields.get("problems", ()))
        self.sys_prefix = fields.get("sys_prefix")
        self.base_prefix = fields.get("base_prefix")
        self.machine = fields.get("machine") or platform.machine()
        self.proj_env = dict(fields.get("proj_env") or {})
        self.proj_search_paths = tuple(fields.get("proj_search_paths") or ())

    # -- derived -----------------------------------------------------------

    @property
    def is_virtualenv(self):
        return bool(self.sys_prefix) and self.sys_prefix != self.base_prefix

    @property
    def is_conda(self):
        return bool(self.sys_prefix) and os.path.isdir(
            os.path.join(self.sys_prefix, "conda-meta")
        )

    @property
    def target_prefix(self):
        """The environment an install should go into.

        A virtualenv or conda environment is self-contained, so a C library
        built into it is removed with it and its ``bin`` is already on PATH
        when the environment is active -- which is exactly what building the
        bindings needs.  Anything else falls back to a per-user directory
        rather than writing into a system prefix.
        """
        if self.is_virtualenv or self.is_conda:
            return self.sys_prefix
        return default_source_prefix()

    @property
    def ok(self):
        """True when nothing needs doing."""
        return not self.problems

    @property
    def codes(self):
        return tuple(code for code, _ in self.problems)

    def has(self, code):
        return code in self.codes

    def describe(self):
        """A human-readable report, as a list of lines."""
        lines = [
            "python            : %s (%s)" % (self.python_version or "?", self.python),
            "numpy             : %s" % (self.numpy_version or "not installed"),
            "GDAL bindings     : %s" % (self.bindings_version or "not installed"),
            "GDAL C library    : %s" % (self.library_version or "not found"),
        ]
        if self.gdal_config:
            lines.append("gdal-config       : %s" % self.gdal_config)
        lines.append(
            "numpy bridge      : %s"
            % ("present" if self.has_gdal_array else "MISSING (ReadAsArray will fail)")
        )
        lines.append(
            "projection lookup : %s" % ("works" if self.projection_works else "FAILS")
        )
        if self.problems:
            lines.append("")
            lines.append("problems:")
            for code, message in self.problems:
                lines.append("  [%s] %s" % (code, message))
        else:
            lines.append("")
            lines.append("GDAL is installed and working.")
        return lines


def probe(python=None):
    """Inspect ``python`` (default: this interpreter) for a working GDAL."""
    python = python or sys.executable
    fields = {"python_executable": python}

    result = subprocess.run(
        [python, "-c", _PROBE_SCRIPT],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        fields.update(json.loads(result.stdout.decode("utf8", "replace")))
    except ValueError:
        # The probe did not even get to print: importing osgeo took the
        # interpreter down with it.  That is itself the finding.
        fields["crashed"] = True
        fields["bindings_error"] = (
            "the interpreter died while importing osgeo (exit %s). %s"
            % (result.returncode, result.stderr.decode("utf8", "replace").strip()[:400])
        )

    config = _gdal_config_info(sys_prefix=fields.get("sys_prefix"))
    if config:
        fields["gdal_config"] = config["path"]
        fields["library_version"] = config["version"]
        fields["prefix"] = config["prefix"]

    found = Probe(**fields)
    found.problems = tuple(_diagnose(found, config))
    return found


def _diagnose(found, config):
    """Turn a raw probe into a list of ``(code, message)`` problems."""
    problems = []

    if found.bindings_version is None:
        problems.append(
            (
                "no-bindings",
                "the osgeo Python bindings are not installed (%s)"
                % (found.bindings_error or "import failed"),
            )
        )

    if config is None:
        # Not a problem in itself when the bindings work -- a wheel or conda
        # package can carry the library without shipping gdal-config -- but
        # it is what blocks building bindings from an sdist.
        if found.bindings_version is None:
            problems.append(
                (
                    "no-library",
                    "no GDAL C library found: gdal-config is not on PATH",
                )
            )
    else:
        library = _parse_version(config["version"])
        if library and library[:2] < MIN_LIBRARY_VERSION:
            problems.append(
                (
                    "old-library",
                    "GDAL %s is older than the %d.%d TopoAnalysis needs"
                    % (config["version"], MIN_LIBRARY_VERSION[0], MIN_LIBRARY_VERSION[1]),
                )
            )

        bindings = _parse_version(found.bindings_version)
        if bindings and library and bindings[:2] != library[:2]:
            problems.append(
                (
                    "version-mismatch",
                    "the bindings are %s but the C library is %s; they must "
                    "agree to at least major.minor"
                    % (found.bindings_version, config["version"]),
                )
            )

        mismatch = _architecture_mismatch(config["prefix"], found.machine)
        if mismatch:
            problems.append(("arch-mismatch", mismatch))

    if found.bindings_version is not None and not found.has_gdal_array:
        problems.append(
            (
                "no-gdal-array",
                "osgeo.gdal_array is missing, so Band.ReadAsArray() will fail "
                "(%s). The bindings were built without numpy present."
                % (found.gdal_array_error or "import failed"),
            )
        )

    if found.bindings_version is not None and not found.projection_works:
        message = "looking up EPSG:4326 failed (%s)" % (
            found.projection_error or "unknown error"
        )
        if found.proj_env:
            message += ". These are set in the environment and override the "
            message += "paths compiled into the library: " + ", ".join(
                "%s=%s" % item for item in sorted(found.proj_env.items())
            )
        if found.proj_search_paths:
            message += ". PROJ searched: " + ", ".join(found.proj_search_paths)
        problems.append(("no-proj-data", message))

    return problems


def _manifest_for(python_prefix):
    """The GDAL manifest recorded inside ``python_prefix``, or None."""
    try:
        with open(manifest_path(python_prefix)) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _config_search_order(sys_prefix=None):
    """Where to look for ``gdal-config``, most specific first.

    PATH is the *last* resort, not the first.  ``--python`` can name an
    interpreter in a different environment, and what matters is the library
    that interpreter will load -- not whichever ``gdal-config`` happens to be
    on the PATH of the shell this command was typed into.  Getting that wrong
    produces a confident and completely wrong version-mismatch diagnosis.
    """
    candidates = []
    explicit = os.environ.get("GDAL_CONFIG")
    if explicit:
        candidates.append(explicit)
    if sys_prefix:
        manifest = _manifest_for(sys_prefix)
        if manifest and manifest.get("prefix"):
            candidates.append(
                os.path.join(manifest["prefix"], "bin", "gdal-config")
            )
        candidates.append(os.path.join(sys_prefix, "bin", "gdal-config"))
    candidates.append("gdal-config")
    return candidates


def _gdal_config_info(gdal_config=None, sys_prefix=None):
    """Ask ``gdal-config`` about the installed C library, or return None."""
    if gdal_config:
        candidates = [gdal_config]
    else:
        candidates = _config_search_order(sys_prefix)

    for program in candidates:
        path = program if os.sep in program else shutil.which(program)
        if not path or not os.path.exists(path):
            continue
        info = _interrogate_gdal_config(path)
        if info is not None:
            return info
    return None


def _interrogate_gdal_config(path):
    info = {"path": path}
    for option in ("version", "prefix", "cflags", "libs"):
        try:
            result = subprocess.run(
                [path, "--" + option], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
            )
        except OSError:
            return None
        if result.returncode != 0:
            return None
        info[option] = result.stdout.decode("utf8", "replace").strip()
    return info if info.get("version") else None


def _parse_version(text):
    """``'3.9.3'`` -> ``(3, 9, 3)``; None when it does not look like one."""
    if not text:
        return None
    match = re.match(r"(\d+)\.(\d+)(?:\.(\d+))?", str(text))
    if not match:
        return None
    return tuple(int(part) for part in match.groups(default="0"))


def _architecture_mismatch(prefix, machine=None):
    """Describe an architecture clash between libgdal and the target Python.

    ``machine`` is the target interpreter's ``platform.machine()``, which is
    not necessarily ours.  Returns None when there is no clash, or when we
    cannot tell -- an unknown answer must not block an install that would
    have worked.
    """
    if sys.platform not in ("darwin", "linux"):
        return None
    library = _find_library(prefix)
    if library is None:
        return None

    tool = shutil.which("file")
    if tool is None:
        return None
    try:
        result = subprocess.run(
            [tool, "-L", library], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
    except OSError:
        return None
    description = result.stdout.decode("utf8", "replace").lower()

    machine = (machine or platform.machine()).lower()
    aliases = {
        "x86_64": ("x86_64", "x86-64", "amd64"),
        "amd64": ("x86_64", "x86-64", "amd64"),
        "arm64": ("arm64", "aarch64"),
        "aarch64": ("arm64", "aarch64"),
    }.get(machine)
    if aliases is None:
        return None
    if any(alias in description for alias in aliases):
        return None

    # The library is a real binary of some architecture, just not ours.
    for name, names in (("x86_64", ("x86_64", "x86-64")), ("arm64", ("arm64", "aarch64"))):
        if any(alias in description for alias in names):
            return (
                "the GDAL C library at %s is %s but this Python is %s; the "
                "bindings cannot link against it. Install a GDAL built for "
                "%s, or use a Python that matches the library."
                % (library, name, machine, machine)
            )
    return None


def _find_library(prefix):
    if not prefix:
        return None
    suffixes = (".dylib",) if sys.platform == "darwin" else (".so",)
    for libdir in ("lib", "lib64"):
        directory = os.path.join(prefix, libdir)
        if not os.path.isdir(directory):
            continue
        for entry in sorted(os.listdir(directory)):
            if entry.startswith("libgdal") and any(
                suffix in entry for suffix in suffixes
            ):
                return os.path.join(directory, entry)
    return None


# ---------------------------------------------------------------------------
# Environments we might install into
# ---------------------------------------------------------------------------


def _conda_prefix():
    """The active conda environment's prefix, or None.

    ``conda-meta`` is the reliable marker: ``CONDA_PREFIX`` can be left over
    from a shell that has since switched interpreters.
    """
    if os.path.isdir(os.path.join(sys.prefix, "conda-meta")):
        return sys.prefix
    return None


def _conda_executable():
    """A conda-family program that can install into the current environment.

    ``CONDA_EXE`` comes first because ``conda`` itself is frequently only a
    shell function, invisible to ``shutil.which``.  mamba is preferred when
    present purely because it is very much faster.
    """
    for name in ("mamba", "micromamba"):
        found = shutil.which(name)
        if found:
            return found
    from_env = os.environ.get("CONDA_EXE")
    if from_env and os.path.exists(from_env):
        return from_env
    return shutil.which("conda")


def default_source_prefix():
    """Per-user fallback prefix, used when there is no environment to use.

    :attr:`Probe.target_prefix` prefers the interpreter's own virtualenv or
    conda environment; this is where it lands when there is neither.
    """
    root = os.environ.get("XDG_DATA_HOME") or os.path.join(
        os.path.expanduser("~"), ".local", "share"
    )
    return os.path.join(root, "TopoAnalysis", "gdal")


#: System package managers, in the order we look for them.  Each entry is
#: (program, install arguments, whether it needs administrator rights).
_SYSTEM_MANAGERS = (
    ("apt-get", ["install", "-y", "gdal-bin", "libgdal-dev"], True),
    ("dnf", ["install", "-y", "gdal", "gdal-devel"], True),
    ("yum", ["install", "-y", "gdal", "gdal-devel"], True),
    ("zypper", ["install", "-y", "gdal", "libgdal-devel"], True),
    ("pacman", ["-S", "--noconfirm", "gdal"], True),
    ("apk", ["add", "--no-cache", "gdal", "gdal-dev"], True),
    # Homebrew refuses to run as root, and does not need to.
    ("brew", ["install", "gdal"], False),
)


def _system_manager():
    for program, arguments, needs_root in _SYSTEM_MANAGERS:
        found = shutil.which(program)
        if found:
            return found, arguments, needs_root
    return None


def _is_root():
    return hasattr(os, "geteuid") and os.geteuid() == 0


# ---------------------------------------------------------------------------
# Plans
# ---------------------------------------------------------------------------


class Step(object):
    """One thing an install does, describable without doing it."""

    def __init__(
        self,
        summary,
        command=None,
        action=None,
        env=None,
        fallbacks=(),
        needs_root=False,
        minutes=1,
    ):
        self.summary = summary
        self.command = command
        self.action = action
        self.env = env or {}
        self.fallbacks = tuple(fallbacks)
        self.needs_root = needs_root
        self.minutes = minutes

    def __repr__(self):  # pragma: no cover - debugging aid
        return "<Step %s>" % self.summary


class Plan(object):
    """An ordered list of :class:`Step` plus why it was chosen.

    ``blockers`` are reasons the plan cannot run at all.  They are kept
    separate from the steps so that ``plan`` can still show what *would*
    happen -- which is the most useful thing to print when the answer is
    "install these two tools first".
    """

    def __init__(self, strategy, steps=(), notes=(), blockers=()):
        self.strategy = strategy
        self.steps = tuple(steps)
        self.notes = tuple(notes)
        self.blockers = tuple(blockers)

    @property
    def minutes(self):
        return sum(step.minutes for step in self.steps)

    @property
    def needs_root(self):
        return any(step.needs_root for step in self.steps)

    def describe(self):
        lines = ["strategy: %s" % self.strategy]
        for note in self.notes:
            lines.append("  note: %s" % note)
        for blocker in self.blockers:
            lines.append("  CANNOT RUN: %s" % blocker)
        if not self.steps:
            lines.append("  (nothing to do)")
            return lines
        lines.append("")
        for number, step in enumerate(self.steps, 1):
            lines.append("  %d. %s" % (number, step.summary))
            if step.command:
                lines.append("     $ %s" % " ".join(step.command))
        lines.append("")
        lines.append("  estimated time: about %d minute(s)" % max(1, self.minutes))
        return lines


def make_plan(
    found=None,
    strategy="auto",
    python=None,
    prefix=None,
    allow_root=False,
    jobs=None,
):
    """Decide how to get a working GDAL, without doing anything yet.

    Parameters
    ----------
    found : Probe, optional
        A previous probe; taken fresh when omitted.
    strategy : str
        ``'auto'``, or one of ``'bindings'``, ``'conda'``, ``'wheel'``,
        ``'system'``, ``'source'`` to force a route.
    python : str, optional
        Interpreter to install for; defaults to the current one.
    prefix : str, optional
        Where a source build should install; see
        :func:`default_source_prefix`.
    allow_root : bool
        Permit steps that need administrator rights.  Without it they are
        still planned, so the user can see the command, but refused at
        execution time.
    jobs : int, optional
        Parallel compile jobs for a source build.
    """
    python = python or sys.executable
    found = found if found is not None else probe(python)
    notes = []

    if strategy == "auto":
        strategy, auto_notes = _choose_strategy(found, allow_root=allow_root)
        notes.extend(auto_notes)

    if strategy == "none":
        return Plan("none", (), ("GDAL is already installed and working.",))

    if strategy == "proj-data":
        return Plan(
            "proj-data",
            (),
            (
                "GDAL and its bindings are fine; PROJ cannot read its "
                "database. Nothing to install -- this is an environment "
                "problem.",
            ),
            _proj_data_advice(found),
        )

    builder = {
        "bindings": _plan_bindings,
        "conda": _plan_conda,
        "wheel": _plan_wheel,
        "system": _plan_system,
        "source": _plan_source,
    }.get(strategy)
    if builder is None:
        raise ValueError("unknown strategy %r" % strategy)

    produced = builder(
        found, python=python, prefix=prefix, allow_root=allow_root, jobs=jobs
    )
    # Builders return (steps, notes), or (steps, notes, blockers) when they
    # can tell in advance that this route cannot work on this machine.
    steps, more_notes = produced[0], produced[1]
    blockers = produced[2] if len(produced) > 2 else ()
    return Plan(strategy, steps, notes + list(more_notes), blockers)


def _proj_data_advice(found):
    """What to actually do about a PROJ database that will not load."""
    advice = []
    if found is not None and found.proj_env:
        advice.append(
            "%s set in the environment, overriding the paths compiled into "
            "the library. Activating a conda environment sets these, and they "
            "then follow you into every other GDAL on the machine. Unset them "
            "and try again:  unset %s"
            % (
                " and ".join(sorted(found.proj_env)),
                " ".join(sorted(found.proj_env)),
            )
        )
    else:
        advice.append(
            "PROJ's data files are missing from the installation itself. "
            "Reinstall PROJ -- `conda install -c conda-forge proj` in a conda "
            "environment, the system's proj-data package otherwise -- or "
            "build GDAL with `topoanalysis-gdal install --strategy source`, "
            "which installs a matching PROJ alongside it."
        )
    return advice


def _choose_strategy(found, allow_root=False):
    """Pick the cheapest route that stands a chance of working.

    ``allow_root`` matters here and not only at execution time: asking for it
    is what makes the system package manager a usable route at all, and
    without it a machine with apt would silently fall through to a half-hour
    source build.
    """
    notes = []

    if found.ok:
        return "none", notes

    # A pure numpy-bridge or PROJ-data failure is repaired by rebuilding the
    # bindings against the library that is already here; no need to touch the
    # C library at all.
    # A PROJ data problem on its own is not fixed by installing anything:
    # the library and bindings are fine, and something in the environment is
    # pointing them at the wrong database.  Reinstalling GDAL would leave
    # that variable exactly where it was.
    if set(found.codes) == {"no-proj-data"}:
        return "proj-data", notes

    # ``no-proj-data`` *alongside* other problems is different: it usually
    # just follows from bindings that are broken or mismatched, and goes away
    # when they are rebuilt.
    repairable = {"no-bindings", "no-gdal-array", "version-mismatch", "no-proj-data"}
    wrong_architecture = found.has("arch-mismatch")
    if wrong_architecture:
        # Nothing installed on top of this library will fix it; say so plainly.
        notes.append(
            "the C library's architecture does not match this interpreter, so "
            "the existing library cannot be reused"
        )
    elif found.gdal_config and not found.has("old-library"):
        if set(found.codes) <= repairable:
            return "bindings", notes

    if found.is_conda and _conda_executable():
        notes.append("a conda environment is active, which is the most reliable route")
        return "conda", notes

    if sys.platform == "win32":
        notes.append(
            "on Windows the only practical source of a GDAL C library is a "
            "prebuilt wheel or conda"
        )
        return "wheel", notes

    manager = _system_manager()
    if manager is not None and found.has("old-library"):
        # The too-old library almost certainly came from this very package
        # manager -- a distro pinned to GDAL 2.x -- so asking it again just
        # reinstalls the same version, and the bindings step would then
        # faithfully build bindings for a library we have already rejected.
        notes.append(
            "%s is available but is where the too-old library came from; "
            "a newer GDAL has to come from somewhere else"
            % os.path.basename(manager[0])
        )
    elif manager is not None and wrong_architecture:
        # The package manager is very often what produced the mismatch in the
        # first place -- a Homebrew arm64 GDAL under a Rosetta x86_64 Python --
        # so asking it again would install the same wrong thing and loop.
        notes.append(
            "%s is available but would supply the same architecture again; "
            "building from source is what actually resolves this"
            % os.path.basename(manager[0])
        )
    elif manager is not None:
        program, _, needs_root = manager
        if not needs_root or _is_root() or allow_root:
            return "system", notes
        notes.append(
            "%s could install the C library but needs administrator rights; "
            "pass --allow-root to use it, or fall back to a source build"
            % os.path.basename(program)
        )

    notes.append(
        "no package manager can be used without help, so GDAL will be built "
        "from source -- this needs only a compiler and CMake"
    )
    return "source", notes


def _pip(python, *arguments):
    return [python, "-m", "pip", "install"] + list(arguments)


def _bindings_requirements(library_version):
    """Version specifiers to try for the bindings, best first.

    The bindings and the C library have to agree.  An exact match is ideal;
    failing that any patch release of the same major.minor will do, because
    GDAL does not break its C API within a minor version.
    """
    parsed = _parse_version(library_version)
    if not parsed:
        return ["gdal"]
    major, minor, patch = parsed
    return [
        "gdal==%d.%d.%d" % (major, minor, patch),
        "gdal==%d.%d.*" % (major, minor),
    ]


def _plan_bindings(
    found, python, prefix=None, allow_root=False, jobs=None, config=None, force=None
):
    """Build ``osgeo`` against a C library that is already installed.

    ``force`` re-installs even when the requested version is already
    present.  It has to default to on whenever we are *repairing* bindings,
    because the commonest repair -- a missing ``gdal_array`` -- leaves an
    ``osgeo`` of exactly the right version in place.  Plain ``pip install
    gdal==X.Y.Z`` would report "already satisfied" and change nothing, and
    even with a reinstall pip would happily serve the broken wheel back out
    of its own cache, since that wheel's name and version are correct.  Only
    ``--force-reinstall --no-cache-dir`` actually rebuilds.
    """
    config = config or _gdal_config_info(
        sys_prefix=found.sys_prefix if found is not None else None
    )
    steps = []
    notes = []
    if config is None:
        raise ValueError(
            "the 'bindings' strategy needs a GDAL C library, but gdal-config "
            "is not on PATH"
        )

    version = config["version"]
    env = {}
    # Make sure the gdal-config we chose is the one the sdist's setup.py
    # finds, whatever else is on PATH.
    env["PATH"] = os.path.dirname(config["path"]) + os.pathsep + os.environ.get("PATH", "")
    libdir = _library_directory(config)
    if libdir:
        # Bake the library location into the extension so importing osgeo
        # works without the user setting LD_LIBRARY_PATH first.
        env["LDFLAGS"] = (
            "-L%s -Wl,-rpath,%s %s" % (libdir, libdir, os.environ.get("LDFLAGS", ""))
        ).strip()

    steps.append(
        Step(
            "make sure the build tools and numpy are present",
            command=_pip(python, "--upgrade", "setuptools", "wheel", "numpy"),
            minutes=1,
        )
    )

    parsed = _parse_version(version)
    if parsed and parsed[:2] < NUMPY_BUILD_REQUIRES_FROM:
        notes.append(
            "GDAL %s predates the sdist declaring numpy as a build dependency, "
            "so the bindings are installed without build isolation -- otherwise "
            "osgeo.gdal_array would be silently left out and ReadAsArray() "
            "would fail at run time" % version
        )

    if force is None:
        force = found is not None and found.bindings_version is not None
    extra = ["--no-build-isolation"]
    if force:
        # --no-deps because step 1 has already put numpy in place, and
        # --force-reinstall would otherwise tear down and rebuild the user's
        # numpy as a side effect of repairing GDAL.
        extra += ["--force-reinstall", "--no-cache-dir", "--no-deps"]
        notes.append(
            "the bindings are being rebuilt rather than merely requested: "
            "this version pin can be satisfied by an install already present, "
            "or by a wheel in pip's cache, that was built without numpy or "
            "against a different C library"
        )

    requirements = _bindings_requirements(version)
    steps.append(
        Step(
            "build the GDAL %s Python bindings against the installed C library"
            % version,
            command=_pip(python, *(extra + [requirements[0]])),
            fallbacks=[
                _pip(python, *(extra + [requirement]))
                for requirement in requirements[1:]
            ],
            env=env,
            minutes=3,
        )
    )
    return steps, notes


def _library_directory(config):
    """Where libgdal actually lives, according to gdal-config.

    Taken from ``--libs`` rather than assumed to be ``<prefix>/lib``: on
    Fedora, RHEL and openSUSE a 64-bit library lands in ``lib64``, and an
    RPATH pointing at an empty ``lib`` would leave the bindings unable to
    find the library they were just built against.
    """
    for token in (config.get("libs") or "").split():
        if token.startswith("-L"):
            candidate = token[2:]
            if candidate and os.path.isdir(candidate):
                return candidate
    # No usable -L: find the library itself rather than guess between lib
    # and lib64, both of which commonly exist.
    library = _find_library(config.get("prefix"))
    return os.path.dirname(library) if library else None


def _plan_conda(found, python, prefix=None, allow_root=False, jobs=None):
    executable = _conda_executable()
    if executable is None:
        raise ValueError("no conda, mamba or micromamba found")
    target = (found.sys_prefix if found is not None else None) or _conda_prefix() or sys.prefix
    command = [
        executable,
        "install",
        "-y",
        "-p",
        target,
        "-c",
        "conda-forge",
        "gdal>=%d.%d" % MIN_LIBRARY_VERSION,
    ]
    notes = [
        "conda-forge ships the C library and the bindings as one package, "
        "already matched, which removes both of the traps this module exists "
        "to avoid"
    ]
    return [
        Step(
            "install GDAL from conda-forge into %s" % target,
            command=command,
            minutes=3,
        )
    ], notes


def _plan_wheel(found, python, prefix=None, allow_root=False, jobs=None):
    """Windows: install a wheel that bundles the C library."""
    find_links = _windows_find_links()
    notes = [
        "these wheels bundle the GDAL C library, so nothing else has to be "
        "installed; they are built by Christoph Gohlke, who has maintained "
        "the de-facto Windows builds of the scientific stack for years"
    ]
    return [
        Step(
            "install prebuilt GDAL wheels for Windows",
            command=_pip(
                python, "--only-binary", ":all:", "--find-links", find_links, "gdal"
            ),
            minutes=2,
        )
    ], notes


def _windows_find_links():
    """The newest cgohlke/geospatial-wheels release page, or a pinned one."""
    try:
        import urllib.request

        with urllib.request.urlopen(
            "https://api.github.com/repos/cgohlke/geospatial-wheels/releases/latest",
            timeout=15,
        ) as response:
            tag = json.load(response)["tag_name"]
    except Exception:
        tag = CGOHLKE_FALLBACK_TAG
    return (
        "https://github.com/cgohlke/geospatial-wheels/releases/expanded_assets/%s" % tag
    )


def _plan_system(found, python, prefix=None, allow_root=False, jobs=None):
    manager = _system_manager()
    if manager is None:
        raise ValueError("no supported system package manager found")
    program, arguments, needs_root = manager

    elevate = needs_root and not _is_root()

    def _command(arguments):
        return (["sudo"] if elevate else []) + [program] + arguments

    steps = []
    notes = []
    if os.path.basename(program) == "apt-get":
        # On a container image whose package lists were stripped -- which is
        # every official Debian and Ubuntu image -- `apt-get install` fails
        # outright with "Unable to locate package" until the lists are
        # refreshed.
        steps.append(
            Step(
                "refresh the apt package lists",
                command=_command(["update"]),
                needs_root=elevate,
                minutes=1,
            )
        )

    steps.append(
        Step(
            "install the GDAL C library with %s" % os.path.basename(program),
            command=_command(arguments),
            needs_root=elevate,
            minutes=3,
        )
    )

    # The bindings step cannot be planned in detail yet -- gdal-config does
    # not exist until the step above has run -- so it is deferred and built
    # when we get there.
    steps.append(
        Step(
            "build the Python bindings against the newly installed C library",
            action=_DeferredBindings(python),
            minutes=3,
        )
    )
    return steps, notes


def _plan_source(found, python, prefix=None, allow_root=False, jobs=None):
    prefix = prefix or (
        found.target_prefix if found is not None else default_source_prefix()
    )
    # So the plan quotes the versions the build will actually fetch.
    _gdal_source.apply_environment_overrides()
    missing = _gdal_source.missing_build_tools()
    notes = [
        "GDAL will be built from pinned, checksummed tarballs and installed "
        "into %s; nothing outside that prefix is modified" % prefix
    ]
    blockers = []
    if missing:
        blockers.append(
            "a source build needs %s, which %s not installed. %s"
            % (
                " and ".join(missing),
                "are" if len(missing) > 1 else "is",
                _build_tool_hint(),
            )
        )

    steps = []
    for component in _gdal_source.components_needed():
        spec = _gdal_source.SOURCES[component]
        steps.append(
            Step(
                "download and build %s %s" % (component, spec["version"]),
                action=_SourceBuild(component, prefix, jobs),
                minutes=spec["minutes"],
            )
        )
    if "sqlite" not in _gdal_source.components_needed():
        notes.append("the system SQLite is usable, so it will not be rebuilt")

    steps.append(
        Step(
            "record the prefix so TopoAnalysis can find GDAL's data files",
            action=_RecordPrefix(
                prefix, found.sys_prefix if found is not None else None
            ),
            minutes=1,
        )
    )
    steps.append(
        Step(
            "build the Python bindings against the GDAL we just built",
            action=_DeferredBindings(python, prefix=prefix),
            minutes=3,
        )
    )
    return steps, notes, blockers


def _build_tool_hint():
    """How to get a compiler and CMake on this platform."""
    if sys.platform == "darwin":
        return "Install them with: xcode-select --install && brew install cmake"
    if sys.platform == "win32":
        return (
            "Install the Visual Studio C++ build tools and CMake, or use "
            "--strategy conda instead."
        )
    manager = _system_manager()
    if manager and os.path.basename(manager[0]) == "apt-get":
        return "Install them with: sudo apt-get install build-essential cmake"
    if manager and os.path.basename(manager[0]) in ("dnf", "yum"):
        return "Install them with: sudo dnf install gcc-c++ make cmake"
    return "Install a C/C++ compiler and CMake, or use --strategy conda."


class _SourceBuild(object):
    """Deferred call into :mod:`TopoAnalysis._gdal_source`."""

    def __init__(self, component, prefix, jobs):
        self.component = component
        self.prefix = prefix
        self.jobs = jobs

    def __call__(self, log):
        _gdal_source.build_one(self.component, self.prefix, jobs=self.jobs, log=log)


class _RecordPrefix(object):
    """Note which environment's ``osgeo`` belongs to this GDAL prefix."""

    def __init__(self, prefix, python_prefix=None):
        self.prefix = prefix
        self.python_prefix = python_prefix

    def __call__(self, log):
        path = write_manifest(self.prefix, self.python_prefix)
        log("    recorded %s" % path)


class _DeferredBindings(object):
    """Plan and run the bindings step once the C library actually exists.

    ``system`` and ``source`` both install a C library as an earlier step, so
    ``gdal-config`` cannot be interrogated when the plan is drawn up.
    """

    def __init__(self, python, prefix=None):
        self.python = python
        self.prefix = prefix

    def __call__(self, log):
        search = None
        if self.prefix:
            search = os.path.join(self.prefix, "bin", "gdal-config")
            if not os.path.exists(search):
                raise RuntimeError(
                    "expected gdal-config at %s after the build, but it is not "
                    "there" % search
                )
        config = _gdal_config_info(search)
        if config is None:
            raise RuntimeError(
                "the C library was installed but gdal-config is still not on "
                "PATH, so the bindings cannot be built"
            )
        # The C library changed underneath whatever bindings may already be
        # installed, so they must be rebuilt even if their version still
        # nominally matches.
        steps, notes = _plan_bindings(None, self.python, config=config, force=True)
        for note in notes:
            log("    note: %s" % note)
        for step in steps:
            log("    %s" % step.summary)
            _run_step(step, log)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def apply_plan(plan, allow_root=False, log=None, dry_run=False):
    """Carry out ``plan``.  Returns True when every step succeeded."""
    log = log or (lambda line: print(line, flush=True))

    if dry_run:
        for line in plan.describe():
            log(line)
        return True

    if plan.blockers:
        # Stopping here turns "cmake: command not found" three minutes into a
        # build into a sentence that says what to install.
        log("")
        log("This route cannot run on this machine:")
        for blocker in plan.blockers:
            log("    %s" % blocker)
        return False

    if plan.needs_root and not allow_root and not _is_root():
        log("")
        log("This route needs administrator rights. Either run:")
        for step in plan.steps:
            if step.needs_root and step.command:
                log("    %s" % " ".join(step.command))
        log("")
        log("and then re-run this command, or pass --allow-root to let it use")
        log("sudo directly, or choose --strategy source to build GDAL without")
        log("any administrator rights at all.")
        return False

    for number, step in enumerate(plan.steps, 1):
        log("[%d/%d] %s" % (number, len(plan.steps), step.summary))
        _run_step(step, log)
    return True


def _run_step(step, log):
    if step.action is not None:
        step.action(log)
        return

    attempts = [step.command] + list(step.fallbacks)
    environment = dict(os.environ)
    environment.update(step.env)

    last = None
    for index, command in enumerate(attempts):
        log("    $ %s" % " ".join(command))
        result = subprocess.run(command, env=environment)
        if result.returncode == 0:
            return
        last = result.returncode
        if index + 1 < len(attempts):
            log("    that failed (exit %d); trying %s" % (last, attempts[index + 1][-1]))
    raise RuntimeError("%s failed (exit %s)" % (" ".join(attempts[0]), last))


# ---------------------------------------------------------------------------
# Remembering a source-built prefix
# ---------------------------------------------------------------------------


MANIFEST_NAME = "gdal-prefix.json"


def manifest_path(python_prefix=None):
    """Where the record of a TopoAnalysis-built GDAL lives.

    Preferably inside the environment the bindings were installed into, so
    that two virtualenvs can each have their own build.  The per-user cache
    is the fallback for the case where that environment is not writable.
    """
    if python_prefix:
        return os.path.join(python_prefix, "share", "topoanalysis", MANIFEST_NAME)
    root = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache"
    )
    return os.path.join(root, "TopoAnalysis", MANIFEST_NAME)


def write_manifest(prefix, python_prefix=None):
    """Record a GDAL prefix that TopoAnalysis built, and return the path.

    ``python_prefix`` is the environment whose ``osgeo`` was built against
    ``prefix``.  It is what makes :func:`activate_environment` safe: the data
    files are only forced on an interpreter that is actually running the
    library they belong to.
    """
    python_prefix = python_prefix or sys.prefix
    payload = {
        "prefix": prefix,
        "gdal_data": os.path.join(prefix, "share", "gdal"),
        "proj_lib": os.path.join(prefix, "share", "proj"),
        "gdal_version": _gdal_source.SOURCES["gdal"]["version"],
        "python_prefix": python_prefix,
    }

    # Inside the environment if we can, and only otherwise in the per-user
    # cache.  Writing both would leave a cache copy describing an environment
    # it does not belong to, which is exactly the confusion the
    # ``python_prefix`` field exists to prevent.
    for candidate in (manifest_path(python_prefix), manifest_path()):
        try:
            os.makedirs(os.path.dirname(candidate), exist_ok=True)
            with open(candidate, "w") as handle:
                json.dump(payload, handle, indent=2)
        except OSError:
            continue
        return candidate
    raise RuntimeError(  # pragma: no cover - both locations unwritable
        "could not record the GDAL prefix in %s or %s"
        % (manifest_path(python_prefix), manifest_path())
    )


def read_manifest():
    """The manifest for this interpreter, or None."""
    for candidate in (manifest_path(sys.prefix), manifest_path()):
        try:
            with open(candidate) as handle:
                return json.load(handle)
        except (OSError, ValueError):
            continue
    return None


def activate_environment():
    """Point GDAL and PROJ at the data files belonging to the library in use.

    This exists because of a failure that is easy to hit and very hard to
    read.  Activating a conda environment exports ``GDAL_DATA`` and
    ``PROJ_LIB``, and those variables outrank the paths compiled into any
    *other* GDAL on the machine.  A PROJ 9 library handed PROJ 8's
    ``proj.db`` does not fall back gracefully; it refuses every coordinate
    lookup with a message about ``DATABASE.LAYOUT.VERSION`` metadata that
    names neither the environment variable nor the installation that set it.

    So these variables are *overridden*, not merely filled in when unset --
    a stale value is precisely the problem, and an unset one is comparatively
    harmless.  Two things keep that from being reckless:

    * it only happens when the manifest was written for this very
      interpreter, which means its ``osgeo`` was built against the prefix
      being pointed at, so these really are the matching data files; and
    * ``TOPOANALYSIS_NO_GDAL_ENV=1`` turns it off entirely.

    Any failure here is swallowed: this must never be the reason an import
    breaks.
    """
    try:
        if os.environ.get("TOPOANALYSIS_NO_GDAL_ENV"):
            return
        manifest = read_manifest()
        if not manifest:
            return
        # A manifest from a different environment describes a library this
        # interpreter is not running; its data files would be as wrong as
        # the ones we are trying to displace.
        recorded = manifest.get("python_prefix")
        if recorded is not None and recorded != sys.prefix:
            return

        proj_lib = manifest.get("proj_lib")
        if proj_lib and os.path.isfile(os.path.join(proj_lib, "proj.db")):
            # PROJ 9.1 renamed PROJ_LIB to PROJ_DATA and still honours both.
            os.environ["PROJ_LIB"] = proj_lib
            os.environ["PROJ_DATA"] = proj_lib

        gdal_data = manifest.get("gdal_data")
        if gdal_data and os.path.isdir(gdal_data):
            os.environ["GDAL_DATA"] = gdal_data
    except Exception:  # pragma: no cover - defensive by design
        pass


# ---------------------------------------------------------------------------
# Convenience entry point
# ---------------------------------------------------------------------------


def ensure(
    strategy="auto",
    python=None,
    prefix=None,
    allow_root=False,
    jobs=None,
    dry_run=False,
    log=None,
):
    """Make GDAL work, and return the resulting :class:`Probe`.

    This is what :file:`install.py` calls.  It is a no-op when GDAL already
    works, so it is safe to call on every install.
    """
    log = log or (lambda line: print(line, flush=True))
    python = python or sys.executable

    before = probe(python)
    if before.ok:
        log("GDAL %s is already installed and working." % before.bindings_version)
        return before

    plan = make_plan(
        before,
        strategy=strategy,
        python=python,
        prefix=prefix,
        allow_root=allow_root,
        jobs=jobs,
    )
    # apply_plan prints the plan itself on a dry run, so describing it here
    # too would say everything twice.
    if not dry_run:
        for line in plan.describe():
            log(line)
        log("")

    if not apply_plan(plan, allow_root=allow_root, log=log, dry_run=dry_run):
        return before
    if dry_run:
        return before

    after = probe(python)
    if after.ok:
        log("")
        log("GDAL %s is installed and working." % after.bindings_version)
    return after


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="topoanalysis-gdal",
        description=(
            "Find, diagnose and install GDAL -- both the C library and the "
            "Python bindings."
        ),
    )
    parser.add_argument(
        "command",
        nargs="?",
        default="doctor",
        choices=("doctor", "plan", "install"),
        help=(
            "doctor: report what is installed and what is wrong (default). "
            "plan: show what an install would do. install: do it."
        ),
    )
    parser.add_argument(
        "--strategy",
        default="auto",
        choices=("auto", "bindings", "conda", "wheel", "system", "source"),
        help="force a particular route instead of choosing automatically",
    )
    parser.add_argument(
        "--prefix",
        help="where a source build installs (default: the active virtualenv)",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="install for this interpreter instead of the current one",
    )
    parser.add_argument(
        "--allow-root",
        action="store_true",
        help="permit steps that need sudo, instead of only printing them",
    )
    parser.add_argument(
        "--jobs", type=int, help="parallel compile jobs for a source build"
    )
    arguments = parser.parse_args(argv)

    found = probe(arguments.python)

    if arguments.command == "doctor":
        for line in found.describe():
            print(line)
        if not found.ok:
            print("")
            print("Run 'topoanalysis-gdal plan' to see how this would be fixed.")
        return 0 if found.ok else 1

    try:
        plan = make_plan(
            found,
            strategy=arguments.strategy,
            python=arguments.python,
            prefix=arguments.prefix,
            allow_root=arguments.allow_root,
            jobs=arguments.jobs,
        )
    except ValueError as exc:
        print("cannot plan an install: %s" % exc, file=sys.stderr)
        return 2

    if arguments.command == "plan":
        for line in plan.describe():
            print(line)
        return 0

    try:
        ok = apply_plan(plan, allow_root=arguments.allow_root)
    except Exception as exc:
        print("")
        print("install failed: %s" % exc, file=sys.stderr)
        return 1
    if not ok:
        return 1

    after = probe(arguments.python)
    print("")
    for line in after.describe():
        print(line)
    return 0 if after.ok else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
