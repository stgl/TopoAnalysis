# Installation

```bash
git clone <this repository>
cd TopoAnalysis
python3 install.py
```

That is the whole procedure. `install.py` creates a virtual environment,
installs the dependencies, compiles the Priority-Flood kernels, installs GDAL
— the C library *and* the Python bindings — and then proves the result works
by writing a georeferenced GeoTIFF and reading it back.

It uses nothing but the standard library, so it runs on a machine where
nothing has been set up yet. When it finishes it tells you how to activate
the environment.

```
========================================================================
6. Checking the installation
========================================================================
    TopoAnalysis 1.0.0
    compute backend : c++
    in-memory pipeline : ok (max area 8100 m^2)
    GDAL available  : True
    GDAL version    : 3.13.3
    GeoTIFF round trip : ok (EPSG:32611, 5x5)
```

## What it installs, and what is optional

| Component | Needed for | If it is missing |
|---|---|---|
| NumPy, SciPy, matplotlib | everything | hard requirement |
| a C++17 compiler | the fast kernels | falls back to NumPy, same numbers, ~50× slower |
| GDAL | reading and writing raster files | everything except file I/O still works |

Neither the compiler nor GDAL is required, and `install.py` will finish
without them rather than failing.

## Options

```bash
python3 install.py --no-gdal              # skip GDAL entirely
python3 install.py --venv ~/envs/topo     # a virtualenv somewhere specific
python3 install.py --system               # no virtualenv; use this interpreter
python3 install.py --gdal-strategy source # force a from-source GDAL build
python3 install.py --allow-root           # let apt/dnf use sudo
python3 install.py --jobs 8               # parallel compile jobs
python3 install.py --dry-run              # show the plan, change nothing
```

`--dry-run` is worth running first on an unfamiliar machine. It reports what
GDAL is present, what is wrong with it, and exactly which commands an install
would run.

`--venv` forces a virtual environment even from inside a conda environment,
which is otherwise used as it stands.

**Exit status.** `0` means everything asked for is installed and verified.
`1` means something is missing — including the case where TopoAnalysis
installed correctly but GDAL did not, which the output states explicitly. Use
`--no-gdal` if you do not want GDAL and do not want its absence to be an
error.

**`--system` and PEP 668.** Debian, Ubuntu, Fedora and Homebrew mark their
Python installations as externally managed, and pip refuses to install into
them. `install.py --system` stops with an explanation rather than letting pip
fail obscurely. Inside a container, where there is no host environment to
protect, add `--break-system-packages`.

### Why it makes a virtualenv by default

Because installing into a system Python usually fails, and when it does not
fail it is a bad idea. Debian, Ubuntu, Fedora and Homebrew all mark their
Python installations "externally managed" (PEP 668), so `pip install` refuses
to touch them. Creating `.venv` sidesteps that and gives a from-source GDAL
build somewhere self-contained to live. If you are already inside a
virtualenv or a conda environment, that one is used and nothing new is
created.

---

# GDAL

GDAL is the one dependency pip cannot supply on its own, so it gets its own
tool. If you only want the summary: run `topoanalysis-gdal install` and it
will work out what to do.

## Why it is difficult

There is no GDAL wheel on PyPI. Not for your platform — for *any* platform,
at *any* version. Every release is an sdist whose `setup.py` shells out to
`gdal-config` and compiles the SWIG bindings against a GDAL C library that
must already be on the machine. So installing GDAL is two jobs that have to
agree with each other:

1. the **C library**, `libgdal`, which pip knows nothing about, and
2. the **Python bindings**, `osgeo`, whose version must match it.

`pip install GDAL` does only the second half, and fails confusingly when the
first half is absent or a different version.

## The tool

```bash
topoanalysis-gdal            # what is installed, and what is wrong with it
topoanalysis-gdal plan       # what an install would do, step by step
topoanalysis-gdal install    # do it
```

`doctor` is the one to reach for when something is wrong:

```
python            : 3.11.9 (/home/you/TopoAnalysis/.venv/bin/python)
numpy             : 2.1.0
GDAL bindings     : 3.4.1
GDAL C library    : 3.8.4
gdal-config       : /usr/bin/gdal-config
numpy bridge      : MISSING (ReadAsArray will fail)
projection lookup : works

problems:
  [version-mismatch] the bindings are 3.4.1 but the C library is 3.8.4; they
                     must agree to at least major.minor
  [no-gdal-array]    osgeo.gdal_array is missing, so Band.ReadAsArray() will
                     fail. The bindings were built without numpy present.
```

## The routes it can take

`topoanalysis-gdal install` picks the cheapest one that stands a chance.
`--strategy` forces a particular choice.

| Strategy | When it is chosen | Cost |
|---|---|---|
| `bindings` | the C library is already here; just build `osgeo` against it | seconds |
| `conda` | a conda environment is active | a minute or two |
| `wheel` | Windows, where prebuilt wheels bundling the C library are the only practical option | a minute |
| `system` | apt/dnf/pacman/zypper/apk/brew can supply the C library | a few minutes, needs root except under Homebrew |
| `source` | nothing else is available | about half an hour |

### `source` — building GDAL without root

This is the route that makes an unattended install possible on a machine with
no conda, no package manager and no administrator rights. It downloads,
verifies and builds:

```
sqlite3   only when the system copy is unusable. PROJ needs the library, the
          header, and the sqlite3 command, which it runs to assemble proj.db.
PROJ      mandatory for GDAL 3.x. Built without TIFF and curl, which are only
          needed for the optional high-accuracy datum grids.
GDAL      built against its own bundled libtiff, libgeotiff, libpng, libjpeg
          and zlib, so the only external dependency left is PROJ.
```

Using GDAL's bundled libraries is what makes this realistic: there are no
system `-dev` packages to hunt for. The result reads GeoTIFF, ArcInfo ASCII,
ENVI, HFA and the other built-in drivers — everything TopoAnalysis uses.

Everything installs into one prefix, by default the active virtualenv, so
nothing outside it is touched and deleting the environment removes GDAL with
it. The only prerequisites are a C/C++ compiler and CMake.

Every tarball is pinned by SHA-256, checked against what the projects
themselves publish. A download that does not match the pin is not built. To
move a pin — a newer GDAL, or a local mirror — set all three variables
together, so the integrity check cannot be disabled by accident:

```bash
export TOPOANALYSIS_GDAL_SOURCE_GDAL_VERSION=3.12.4
export TOPOANALYSIS_GDAL_SOURCE_GDAL_URL=https://example.org/gdal-3.12.4.tar.gz
export TOPOANALYSIS_GDAL_SOURCE_GDAL_SHA256=<64 hex digits>
```

The pinned GDAL is deliberately recent. GDAL bundles its own copies of zlib
and libpng, and the copies shipped up to about 3.9 do not compile against a
current Apple SDK: both take their "Mac OS classic" branch, because that is
selected on `TARGET_OS_MAC`, which modern SDKs define on every Apple
platform. zlib then redefines `fdopen` to `NULL` and the system `_stdio.h`
fails to parse.

## Platform notes

**conda (any platform)** — the most reliable route, and the only one that is
equally good on Windows:

```bash
conda env create -f environment.yml
conda activate topoanalysis
python -m pip install -e .
```

`environment.yml` puts `conda-forge` first, because mixing GDAL from
`defaults` with compilers from `conda-forge` produces a C++ ABI mismatch that
only appears at import time.

**Ubuntu / Debian** — `sudo apt install gdal-bin libgdal-dev`, then
`topoanalysis-gdal install`. Note that 22.04 ships GDAL 3.4 and 24.04 ships
3.8; both are below the version where the sdist started declaring numpy as a
build dependency, so install the bindings with `topoanalysis-gdal` rather than
plain pip — see the trap below.

**macOS** — Homebrew's GDAL is built for the Mac's native architecture. If
your Python is x86_64 running under Rosetta on an Apple Silicon machine, the
two cannot link. `topoanalysis-gdal doctor` detects this and will not try to
"fix" it by reinstalling the same thing; use conda, a matching Python, or
`--strategy source`, which builds for whichever architecture the interpreter
is.

**Windows** — use conda, or let `topoanalysis-gdal install` fetch a prebuilt
wheel that bundles the C library. Building GDAL from source on Windows is not
attempted.

**Docker** — the `osgeo/gdal` images already have both halves installed and
matched; `topoanalysis-gdal doctor` will report nothing to do and
`pip install -e .` is enough.

## Two traps worth knowing about

**The silent `gdal_array` trap.** GDAL sdists only grew a `pyproject.toml` in
version 3.9. Before that they declare no build dependencies at all, so under
pip's default build isolation numpy is absent when the bindings compile — and
the numpy bridge, `osgeo.gdal_array`, is quietly skipped. Everything imports
cleanly; `Band.ReadAsArray()` then fails at run time. Ubuntu 22.04 and 24.04
both ship GDAL versions on the wrong side of that line, so this is the common
case, not an exotic one. `topoanalysis-gdal` installs the bindings with
`--no-build-isolation` against an interpreter that already has numpy, which
both avoids the trap and builds the bridge against the numpy you will
actually run.

**The architecture trap.** Described under macOS above. The symptom is a
missing-symbol error at link or import time that says nothing about
architectures.

---

## Checking the installation

```bash
topoanalysis-process --info
topoanalysis-gdal doctor
```

```
TopoAnalysis 1.0.0
  compute backend : c++
  GDAL available  : True
  fill algorithm  : Barnes, R., Lehman, C., Mulla, D. (2014). ...
```

Then run the tests:

```bash
python -m pytest -q                              # a few seconds
TOPOANALYSIS_PURE_PYTHON=1 python -m pytest -q   # same tests, NumPy kernels
```

The second run is worth doing after any change to the C++. It re-runs every
test against the NumPy kernels, and the suite contains explicit
backend-comparison tests (`test_kernels.py::test_fill_backends_agree` and
`test_routing_backends_agree`) that call both implementations on the same
inputs and require identical output — so a divergence shows up as a failure
rather than as a quietly different answer.

## When the compiler is missing

The extension is declared `optional`, so a failed compile is a warning, not an
installation failure:

```
TopoAnalysis: pybind11 not found; building without the compiled kernels.
```

The library still installs and works. To get the speed back:

| Platform | What to install |
|---|---|
| Ubuntu / Debian | `sudo apt install build-essential python3-dev` |
| RHEL / Fedora | `sudo dnf install gcc-c++ python3-devel` |
| macOS | `xcode-select --install` |
| Windows | [Build Tools for Visual Studio](https://visualstudio.microsoft.com/visual-cpp-build-tools/), "Desktop development with C++" |
| conda (any) | `conda install -c conda-forge compilers pybind11` |

Then reinstall:

```bash
python -m pip install --force-reinstall --no-deps -e .
python -c "import TopoAnalysis; print(TopoAnalysis.backend())"   # -> c++
```

Editing a header under `_cpp/` does trigger a rebuild — the extension declares
them in `depends` — but a stale `build/` directory can still confuse
setuptools. `make distclean && make build` clears it.

## Installing without a network

`pip install` needs `pybind11` at build time. If the machine is offline, fetch
the wheel somewhere with a network and carry it over:

```bash
# on a connected machine
pip download pybind11 setuptools wheel -d wheelhouse

# on the target
python -m pip install --no-index --find-links wheelhouse pybind11 setuptools wheel
python -m pip install --no-build-isolation --no-index -e .
```

`--no-build-isolation` is the important flag: without it pip builds in a fresh
environment and tries to download `pybind11` again.

A from-source GDAL build also needs a network, but only to fetch its three
tarballs. They are cached under `~/.cache/TopoAnalysis/downloads`, so copying
that directory across is enough to make the build work offline.

## Building from an sdist

```bash
python -m build --sdist
python -m pip install dist/TopoAnalysis-1.0.0.tar.gz
```

The sdist carries `_cpp/*.hpp` and `_cpp/*.cpp`, so the kernels are rebuilt on
the target machine for its own Python and architecture.

## Both import styles

The repository directory *is* the package, so which style works depends on
what is on `sys.path`:

```python
from TopoAnalysis import dem          # parent directory on the path (installed)
import dem                            # the package directory itself on the path
```

Every module supports both. The first is preferred; the second is how the
original scripts in this repository loaded it, and is kept working
deliberately.

## Optional extras

```bash
pip install -e ".[stats]"      # statsmodels, for the steepness regressions
pip install -e ".[parallel]"   # dask, for multiscale curvature
pip install -e ".[quadrats]"   # rasterio, for analysis.Quadrats
pip install -e ".[all]"        # all of the above
pip install -e ".[docs]"       # sphinx, if you want to build a site
```

`[all]` deliberately excludes GDAL. Listing it there would make
`pip install .[all]` fail on any machine without a GDAL C library already
present, taking the rest of the package down with it.

Each is imported only when the code that needs it runs, and the error names
the missing package and how to get it.

There is a `[gdal]` extra, but it will not work on its own: it installs the
bindings, which cannot build without the C library. Use `topoanalysis-gdal
install` instead, which installs both.

## Troubleshooting

**`ImportError: TopoAnalysis needs the GDAL Python bindings ...`** — you called
something that reads or writes a file. Run `topoanalysis-gdal install`, or
build grids in memory with `dx=..., grid=...`.

**`AttributeError` or a numpy error from `ReadAsArray`** — the `gdal_array`
trap above. `topoanalysis-gdal doctor` confirms it; `topoanalysis-gdal
install` rebuilds the bindings correctly.

**`symbol not found` or a segfault when importing `osgeo`** — either a C++ ABI
mismatch from mixing conda channels, or the architecture trap. Rebuild the
environment from `environment.yml` with `conda-forge` first, or run
`topoanalysis-gdal doctor`, which checks the architecture explicitly.

**`ERROR 1: PROJ: proj_create_from_database: Cannot find proj.db`** — PROJ's
data files are missing or `PROJ_LIB` points somewhere wrong.
`topoanalysis-gdal doctor` reports this as `projection lookup : FAILS`.

**`backend()` says `python` after installing a compiler** — the wheel was
cached. `pip install --force-reinstall --no-deps --no-cache-dir -e .`

**`RuntimeWarning: TopoAnalysis is using its pure-Python kernels on a grid of N
cells`** — exactly what it says; see the compiler table above.

**Tests fail with `ModuleNotFoundError: No module named 'TopoAnalysis'`** — run
them from the package directory (`cd TopoAnalysis && pytest`) or from its
parent (`pytest TopoAnalysis/tests`). Both work, and so does a checkout under
any other name (a GitHub zip unpacks to `TopoAnalysis-main`);
`tests/conftest.py` sorts out the import either way.
