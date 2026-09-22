#!/usr/bin/env python3
"""Install TopoAnalysis and everything it needs, from a fresh clone.

    git clone <this repository>
    cd TopoAnalysis
    python3 install.py

That is the whole procedure.  The script creates an environment, installs the
dependencies, compiles the Priority-Flood kernels, installs GDAL -- the C
library *and* the Python bindings -- and then checks that the result actually
works by reading and writing a GeoTIFF.

It uses nothing but the standard library, so it runs on a machine where
nothing has been set up yet.

Why it makes a virtualenv by default
------------------------------------
Because installing into a system Python usually fails, and when it does not
fail it is a bad idea.  Debian, Ubuntu, Fedora and Homebrew all now mark their
Python installations "externally managed" (PEP 668), so ``pip install``
refuses to touch them.  Creating ``.venv`` sidesteps that, keeps this project
from disturbing anything else, and gives the from-source GDAL build somewhere
self-contained to live.  If you are already inside a virtualenv or a conda
environment, that one is used instead and nothing new is created.

Useful options
--------------
``--no-gdal``            skip GDAL entirely; everything except file I/O works
``--gdal-strategy``      force a route: conda, system, source, bindings, wheel
``--allow-root``         let the system package manager use sudo
``--system``             install into the current interpreter, no virtualenv
``--jobs N``             parallel compile jobs
``--dry-run``            print what would happen and stop
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_PYTHON = (3, 8)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

_STEP = [0]


def heading(text):
    _STEP[0] += 1
    print("")
    print("=" * 72)
    print("%d. %s" % (_STEP[0], text))
    print("=" * 72)
    sys.stdout.flush()


def say(text=""):
    print(text)
    sys.stdout.flush()


def die(text):
    print("")
    print("ERROR: %s" % text, file=sys.stderr)
    sys.exit(1)


def run(command, env=None, check=True):
    say("    $ %s" % " ".join(command))
    result = subprocess.run(command, env=env)
    if check and result.returncode != 0:
        die("command failed (exit %d):\n    %s" % (result.returncode, " ".join(command)))
    return result.returncode


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def in_virtualenv():
    return sys.prefix != getattr(sys, "base_prefix", sys.prefix)


def in_conda():
    return os.path.isdir(os.path.join(sys.prefix, "conda-meta"))


def venv_python(directory):
    """The interpreter inside a virtualenv directory."""
    if sys.platform == "win32":
        return os.path.join(directory, "Scripts", "python.exe")
    return os.path.join(directory, "bin", "python")


def venv_is_usable(python):
    """True when a virtualenv's interpreter exists and has a working pip.

    Checking pip matters because a *partial* virtualenv is a real and common
    thing.  ``venv`` creates the interpreter before it bootstraps pip, and on
    Debian and Ubuntu -- where pip's bootstrap lives in a separate
    ``python3-venv`` package -- it fails in between, leaving a directory that
    has ``bin/python`` and no way to install anything into it.
    """
    if not os.path.exists(python):
        return False
    result = subprocess.run(
        [python, "-c", "import pip"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    return result.returncode == 0


def create_venv(directory):
    """Make a virtualenv at ``directory``, reusing a usable one already there."""
    python = venv_python(directory)
    if venv_is_usable(python):
        say("    reusing the existing environment at %s" % directory)
        return python

    existed = os.path.exists(directory)
    if existed and os.path.exists(python):
        die(
            "there is a virtual environment in %s but it has no working pip,\n"
            "so nothing can be installed into it. Remove it and re-run:\n"
            "    rm -rf %s\n"
            "On Debian and Ubuntu the usual cause is a missing python3-venv:\n"
            "    sudo apt install python3-venv" % (directory, directory)
        )

    say("    creating a virtual environment in %s" % directory)
    result = subprocess.run([sys.executable, "-m", "venv", directory])
    if result.returncode != 0 or not venv_is_usable(python):
        # We made this directory a moment ago and it is not usable, so it is
        # ours to clean up.  Leaving the half-built shell behind would make
        # the next run take the "reuse" path into the same dead end.
        if not existed:
            shutil.rmtree(directory, ignore_errors=True)
        die(
            "could not create a working virtual environment in %s.\n"
            "On Debian and Ubuntu this is almost always a missing package:\n"
            "    sudo apt install python3-venv\n"
            "Alternatively re-run with --system to install into this\n"
            "interpreter instead." % directory
        )
    return python


def externally_managed():
    """The PEP 668 marker path if this interpreter forbids pip installs.

    Debian, Ubuntu, Fedora and Homebrew all ship one.  pip refuses to install
    anything into such an interpreter, with an error that does not mention
    this script at all -- so it is worth saying plainly up front.
    """
    for key in ("stdlib", "platstdlib"):
        directory = sysconfig.get_path(key)
        if not directory:
            continue
        marker = os.path.join(directory, "EXTERNALLY-MANAGED")
        if os.path.exists(marker):
            return marker
    return None


def activation_hint(directory):
    if sys.platform == "win32":
        return "%s\\Scripts\\activate" % directory
    return "source %s/bin/activate" % directory


# ---------------------------------------------------------------------------
# The steps
# ---------------------------------------------------------------------------


def check_python():
    if sys.version_info[:2] < MIN_PYTHON:
        die(
            "TopoAnalysis needs Python %d.%d or newer; this is %s"
            % (MIN_PYTHON[0], MIN_PYTHON[1], sys.version.split()[0])
        )
    say("    Python %s at %s" % (sys.version.split()[0], sys.executable))


#: pip needs to understand PEP 660 to do an editable install of this project.
MIN_PIP = (21, 3)


def install_build_tools(python):
    """setuptools, wheel, pybind11 and numpy, before anything else.

    pybind11 has to be present for the C++ kernels to compile, and numpy has
    to be present before the GDAL bindings are built -- otherwise the numpy
    bridge is silently skipped.  Both are build-time requirements of things
    installed later, so they go in first.

    pip itself is upgraded only when it is too old to do the job.  This
    script may be pointed at an environment the user uses for other work, and
    replacing their pip as a side effect of installing one package is not
    ours to do.
    """
    version = pip_version(python)
    if version is None or version < MIN_PIP:
        say(
            "    pip %s is older than %d.%d; upgrading it"
            % (".".join(map(str, version)) if version else "(unknown)", *MIN_PIP)
        )
        run([python, "-m", "pip", "install", "--upgrade", "pip"])
    else:
        say("    pip %s is recent enough; leaving it alone" % ".".join(map(str, version)))

    run([python, "-m", "pip", "install", "--upgrade", "setuptools", "wheel"])
    run([python, "-m", "pip", "install", "pybind11>=2.6", "numpy>=1.20"])


def pip_version(python):
    """``(major, minor)`` of the target interpreter's pip, or None."""
    result = subprocess.run(
        [python, "-c", "import pip, sys; sys.stdout.write(pip.__version__)"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        return None
    parts = result.stdout.decode("utf8", "replace").strip().split(".")
    try:
        return int(parts[0]), int(parts[1])
    except (IndexError, ValueError):
        return None


def install_package(python, extras, editable):
    target = "." if not extras else ".[%s]" % extras
    command = [python, "-m", "pip", "install"]
    if editable:
        command.append("-e")
    command.append(target)
    run(command)


def install_gdal(python, strategy, prefix, allow_root, jobs, dry_run):
    """Hand off to the GDAL bootstrapper, which lives with the package."""
    sys.path.insert(0, HERE)
    try:
        import gdal_setup
    except ImportError as exc:  # pragma: no cover - only a broken checkout
        die("could not load the GDAL bootstrapper from %s: %s" % (HERE, exc))

    try:
        return gdal_setup.ensure(
            strategy=strategy,
            python=python,
            prefix=prefix,
            allow_root=allow_root,
            jobs=jobs,
            dry_run=dry_run,
            log=say,
        )
    except Exception as exc:
        # A failed GDAL install is a partial failure, not a crash: everything
        # except file I/O still works, and a traceback here would obscure the
        # fact that the rest of the installation succeeded.
        say("")
        say("GDAL could not be installed: %s" % exc)
        return None


_VERIFY = r"""
import sys, tempfile, os
import numpy as np
import TopoAnalysis
from TopoAnalysis import Elevation, FilledElevation, FlowDirectionD8, Area

print("    TopoAnalysis %s" % TopoAnalysis.__version__)
print("    compute backend : %s" % TopoAnalysis.backend())

# A small DEM with a pit, so the fill has something to do.
z = np.array([[5.0, 4.0, 3.0, 4.0, 5.0],
              [4.0, 3.0, 2.0, 3.0, 4.0],
              [3.0, 2.0, 0.5, 2.0, 3.0],
              [2.0, 1.0, 1.5, 1.0, 2.0],
              [1.0, 0.0, 1.0, 0.0, 1.0]])
dem = Elevation(dx=30.0, grid=z)
filled = FilledElevation(elevation=dem)
d8 = FlowDirectionD8(flooded_dem=filled)
area = Area(flow_direction=d8)
assert np.isfinite(area._griddata).any(), "flow accumulation produced nothing"
print("    in-memory pipeline : ok (max area %.0f m^2)" % np.nanmax(area._griddata))

available = TopoAnalysis.gdal_is_available()
print("    GDAL available  : %s" % available)
if not available:
    print("    (file I/O is unavailable; everything else works)")
    sys.exit(0)

from osgeo import gdal, osr
print("    GDAL version    : %s" % gdal.__version__)

# The real test: write a georeferenced GeoTIFF and read it back through
# TopoAnalysis.  This exercises the C library, the numpy bridge and PROJ.
tmp = tempfile.mkdtemp()
path = os.path.join(tmp, "probe.tif")
try:
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(path, 5, 5, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((500000.0, 30.0, 0.0, 4000000.0, 0.0, -30.0))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32611)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(z.astype("float32"))
    ds = None

    back = Elevation(gdal_filename=path)
    assert back._georef_info.nx == 5 and back._georef_info.ny == 5, "wrong shape"
    assert np.allclose(back._griddata, z), "values did not survive the round trip"
    print("    GeoTIFF round trip : ok (EPSG:32611, %dx%d)"
          % (back._georef_info.nx, back._georef_info.ny))
finally:
    import shutil as _sh
    _sh.rmtree(tmp, ignore_errors=True)
"""


def verify(python):
    """Prove the install works rather than assuming it does.

    Run from outside the source tree: importing from the checkout directory
    would succeed even if the install had not, because the package directory
    *is* the repository.  Running elsewhere makes the check mean something.
    """
    result = subprocess.run([python, "-c", _VERIFY], cwd=tempfile.gettempdir())
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Install TopoAnalysis and all of its dependencies, GDAL included.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--venv",
        default=None,
        help="virtualenv directory to create or reuse (default: .venv). "
        "Giving this forces a virtualenv even from inside a conda "
        "environment, which is otherwise used as-is.",
    )
    parser.add_argument(
        "--system",
        action="store_true",
        help="install into the current interpreter instead of a virtualenv",
    )
    parser.add_argument(
        "--break-system-packages",
        action="store_true",
        help="with --system, install into an externally managed (PEP 668) "
        "Python anyway",
    )
    parser.add_argument(
        "--extras",
        default="test",
        help="optional dependency groups to install (default: test; 'all' "
        "adds statsmodels, dask and rasterio; '' for none). GDAL is not an "
        "extra -- this script installs it separately, C library and all.",
    )
    parser.add_argument(
        "--no-editable",
        dest="editable",
        action="store_false",
        help="install a copy instead of linking to this checkout",
    )
    parser.add_argument("--no-gdal", action="store_true", help="do not install GDAL")
    parser.add_argument(
        "--gdal-strategy",
        default="auto",
        choices=("auto", "bindings", "conda", "wheel", "system", "source"),
        help="force a particular route for GDAL (default: auto)",
    )
    parser.add_argument(
        "--gdal-prefix", help="where a from-source GDAL build should install"
    )
    parser.add_argument(
        "--allow-root",
        action="store_true",
        help="permit system package manager steps that need sudo",
    )
    parser.add_argument("--jobs", type=int, help="parallel compile jobs")
    parser.add_argument(
        "--dry-run", action="store_true", help="show what would happen, then stop"
    )
    arguments = parser.parse_args(argv)

    say("TopoAnalysis installer")
    say("running in %s" % HERE)

    if arguments.break_system_packages:
        # Set in the environment rather than passed as a flag so that it also
        # reaches the pip calls made later by the GDAL bootstrapper.
        os.environ["PIP_BREAK_SYSTEM_PACKAGES"] = "1"

    # -- 1 ----------------------------------------------------------------
    heading("Checking the interpreter")
    check_python()

    # -- 2 ----------------------------------------------------------------
    heading("Choosing an environment")
    if arguments.system:
        python = sys.executable
        marker = externally_managed()
        if marker and not arguments.break_system_packages:
            die(
                "this Python is marked externally managed (PEP 668), so pip\n"
                "will refuse to install into it:\n"
                "    %s\n"
                "It is managed by the operating system's package manager, and\n"
                "installing over it can break system tools. Drop --system to\n"
                "use a virtual environment instead, or pass\n"
                "--break-system-packages if you really mean it." % marker
            )
        say("    --system given: installing into %s" % python)
        activate = None
    elif not arguments.venv and (in_virtualenv() or in_conda()):
        # An environment the user already chose and activated. Respect it --
        # but only when they did not name one explicitly, because --venv is
        # how you ask for a separate environment from inside conda.
        python = sys.executable
        say(
            "    already inside %s; installing into it"
            % ("a conda environment" if in_conda() else "a virtualenv")
        )
        say("    %s" % sys.prefix)
        activate = None
    else:
        directory = arguments.venv or os.path.join(HERE, ".venv")
        if arguments.dry_run:
            # A real run would reuse a usable environment that is already
            # there, so a dry run has to report on that one -- describing the
            # outer interpreter instead would be describing a different
            # machine state than the one an install would act on.
            existing = venv_python(directory)
            if venv_is_usable(existing):
                say("    would reuse the existing environment at %s" % directory)
                python = existing
            else:
                say("    would create a virtualenv in %s" % directory)
                python = sys.executable
        else:
            python = create_venv(directory)
        activate = directory

    if arguments.dry_run:
        say("")
        say("--dry-run: stopping before anything is installed.")
        if activate and python == sys.executable:
            say("Note: that environment does not exist yet, so the GDAL report")
            say("below describes %s instead. A fresh" % sys.executable)
            say("virtualenv starts with no GDAL bindings at all, whatever it")
            say("says here.")
            say("")
        say("The GDAL plan for this machine would be:")
        sys.path.insert(0, HERE)
        import gdal_setup

        found = gdal_setup.probe(python)
        for line in found.describe():
            say("    " + line)
        if not found.ok:
            say("")
            plan = gdal_setup.make_plan(
                found,
                strategy=arguments.gdal_strategy,
                python=python,
                prefix=arguments.gdal_prefix,
                allow_root=arguments.allow_root,
                jobs=arguments.jobs,
            )
            for line in plan.describe():
                say("    " + line)
        return 0

    # -- 3 ----------------------------------------------------------------
    heading("Installing build tools")
    install_build_tools(python)

    # -- 4 ----------------------------------------------------------------
    heading("Installing TopoAnalysis and its dependencies")
    os.chdir(HERE)
    install_package(python, arguments.extras, arguments.editable)

    # -- 5 ----------------------------------------------------------------
    gdal_ok = True
    if arguments.no_gdal:
        heading("Skipping GDAL (--no-gdal)")
        say("    everything except reading and writing raster files will work")
    else:
        heading("Installing GDAL (C library and Python bindings)")
        found = install_gdal(
            python,
            arguments.gdal_strategy,
            arguments.gdal_prefix,
            arguments.allow_root,
            arguments.jobs,
            dry_run=False,
        )
        gdal_ok = found is not None and found.ok
        if not gdal_ok:
            say("")
            if found is not None:
                say("GDAL is still not working:")
                for line in found.describe():
                    say("    " + line)
                say("")
            say("TopoAnalysis itself is installed and everything except file")
            say("I/O will work. To retry just the GDAL part:")
            say("    %s -m TopoAnalysis.gdal_setup install --strategy source" % python)

    # -- 6 ----------------------------------------------------------------
    heading("Checking the installation")
    ok = verify(python)

    say("")
    say("=" * 72)
    if ok and gdal_ok:
        say("TopoAnalysis is installed.")
    elif ok:
        say("TopoAnalysis is installed, but GDAL is not -- see above.")
    else:
        say("TopoAnalysis is installed, but the check did not pass -- see above.")
    say("=" * 72)
    if activate:
        say("")
        say("Activate the environment with:")
        say("    %s" % activation_hint(activate))
    say("")
    say("Then try:")
    say("    topoanalysis-process --info")
    say("    topoanalysis-gdal doctor")

    # A GDAL that was asked for and did not arrive is a partial failure, and
    # a script that called this needs to be able to tell.  The message above
    # makes clear that the rest of the installation is fine; --no-gdal is how
    # you say you did not want it in the first place.
    return 0 if (ok and gdal_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
