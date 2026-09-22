# TopoAnalysis

Digital elevation model analysis for tectonic geomorphology: depression
filling, D8 flow routing, drainage area, chi and channel steepness, and the
terrain derivatives that go with them.

Depression filling uses the Priority-Flood algorithms of **Barnes, Lehman and
Mulla (2014)**, compiled to C++ — roughly **50× faster** than the equivalent
NumPy, and 100–500× faster than the pure-Python priority queue this library
used before. If the extension cannot be built, an equivalent NumPy
implementation takes over automatically and produces identical numbers.

A companion module, [`TopoAnalysis.topotoolbox`](docs/topotoolbox_parity.md),
reproduces the results of the MATLAB package
[TopoToolbox 2](https://topotoolbox.wordpress.com/) function by function, so
results can be cross-checked between the two.

---

## Install

```bash
git clone <this repository>
cd TopoAnalysis
python3 install.py
```

That is the whole procedure. It creates a virtual environment, installs the
dependencies, compiles the kernels, installs **GDAL — both the C library and
the Python bindings** — and then checks the result by writing a georeferenced
GeoTIFF and reading it back.

GDAL needs its own installer because pip cannot install it: there is no GDAL
wheel on PyPI for any platform at any version, only an sdist that compiles
against a C library which must already exist. `install.py` supplies both
halves, using conda, the system package manager, a prebuilt wheel on Windows,
or — when none of those is available — by building PROJ and GDAL from pinned,
checksummed source into the virtual environment, which needs no
administrator rights at all.

```bash
topoanalysis-gdal            # what GDAL is installed, and what is wrong with it
topoanalysis-gdal install    # install or repair it
python3 install.py --dry-run # show the plan for this machine, change nothing
```

Neither GDAL nor a C++ compiler is required: without a compiler the library
falls back to NumPy kernels that produce identical numbers, and without GDAL
everything except file I/O works. See
[docs/installation.md](docs/installation.md) for platform notes, offline
installs, and troubleshooting.

## A first analysis

```python
from TopoAnalysis import Elevation, FilledElevation, FlowDirectionD8, Area, FlowLength

dem    = Elevation(gdal_filename='srtm.tif')
filled = FilledElevation(elevation=dem)          # Priority-Flood
d8     = FlowDirectionD8(flooded_dem=filled)     # steepest-descent D8
area   = Area(flow_direction=d8)                 # contributing area, m^2
length = FlowLength(flow_direction=d8)           # longest upstream path

area.save('area.tif')
```

Grids are built **by keyword**, and which keywords you pass decides how the
grid is made:

| What you pass | What you get |
|---|---|
| `gdal_filename='dem.tif'` | read from a file |
| `ai_ascii_filename=..., EPSGprojectionCode=32611` | read an ArcInfo ASCII grid |
| `dx=30.0, grid=array` | wrap an array already in memory |
| `nx=..., ny=..., projection=..., geo_transform=...` | an empty georeferenced grid |
| `flow_direction=d8`, `elevation=dem`, … | derive from another grid |

Every class lists the combinations it accepts in
`required_inputs_and_actions`; passing something else raises an error that
tells you what the class does accept.

## Chi and channel steepness

```python
from TopoAnalysis import Chi, Ksi, ScaledRelief

outlets = area.areas_greater_than(1e8)           # big basins

chi = Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e6, outlets=outlets)
ksi = Ksi(area=area, flow_direction=d8, flow_length=length, theta=0.45, Ao=1e6)
relief = ScaledRelief(flow_direction=d8, elevation=dem, flow_length=length,
                      Ao=1e6, theta=0.45, area=area)
```

`Chi` integrates upstream from the outlets you name; `Ksi` integrates along
longest-flow-length paths across the whole grid. Plotting `ScaledRelief`
against `Ksi` gives a chi-elevation plot whose slope is the steepness index.

## Cross-checking against TopoToolbox

```python
from TopoAnalysis import topotoolbox as tt

zf = tt.fillsinks(dem)             # flat fill, == GRIDobj/fillsinks
fd = tt.flowobj(dem)               # == FLOWobj(DEM,'preprocess','fill')
a  = tt.flowacc(fd)                # == flowacc(FD) -- in CELLS, not m^2
c  = tt.chitransform(fd, a)        # == chitransform(S,A), mn=0.45, a0=1e6
k  = tt.ksn(fd, dem, a)            # == ksn(S,DEM,A)

tt.gradient8(dem)                  # == GRIDobj/gradient8
tt.curvature(dem, 'profc')         # == GRIDobj/curvature
tt.hillshade(dem)                  # == GRIDobj/hillshade  (cosine, not 0-255)
tt.aspect(dem)                     # == GRIDobj/aspect
```

The conventions differ between the two packages in a few places that matter —
drainage area in cells versus square metres, chi integrated with the
trapezoidal rule versus the left-endpoint rule, hillshade as a cosine versus
ESRI's 0–255. [docs/topotoolbox_parity.md](docs/topotoolbox_parity.md)
tabulates all of them and says where exact agreement is and is not
achievable.

## Command line

```bash
topoanalysis-process srtm.tif -o results --products filled,flowdir,area,length
topoanalysis-process --info       # which backend is active, is GDAL present
```

## Performance

`python -m TopoAnalysis.benchmark 1024` — a 1024 × 1024 DEM, about a million
cells, on a 2019 laptop:

| Stage | C++ | NumPy |
|---|---:|---:|
| Priority-Flood fill | 0.10 s | 5.5 s |
| D8 flow directions | 0.02 s | 0.04 s |
| Flow accumulation | 0.01 s | 0.41 s |
| Flow length | 0.02 s | 0.83 s |
| Chi | 0.04 s | 2.6 s |
| **whole pipeline** | **0.20 s** | **9.4 s** |

The compiled kernels release the GIL, so several grids can be processed in
parallel from threads.

## What changed in 1.0

Version 1.0 repairs about 70 defects — several of which stopped the library
importing at all on current NumPy and SciPy — and replaces the recursive
basin traversals, which needed a million-frame recursion limit and still
overflowed on real DEMs, with iterative kernels.

**The public interface is unchanged**: every class, method and keyword still
exists and still means what it did. A handful of results are *different
because they were wrong*; they are listed in [CHANGELOG.md](CHANGELOG.md)
under "Results that change". Read that section before comparing new output
with old.

## Documentation

- [Installation](docs/installation.md) — including offline and no-compiler cases
- [User guide](docs/guide.md) — the grid model, the analysis chain, worked examples
- [API reference](docs/api.md) — every public class and function
- [Algorithms](docs/algorithms.md) — what each computation actually does, and why
- [TopoToolbox parity](docs/topotoolbox_parity.md) — the correspondence table
- [Changelog](CHANGELOG.md) — fixed defects and behavioural changes

## References

Barnes, R., Lehman, C. and Mulla, D. (2014). Priority-flood: An optimal
depression-filling and watershed-labeling algorithm for digital elevation
models. *Computers & Geosciences* 62, 117–127.

Perron, J.T. and Royden, L. (2013). An integral approach to bedrock river
profile analysis. *Earth Surface Processes and Landforms* 38, 570–576.

Schwanghart, W. and Scherler, D. (2014). TopoToolbox 2 — MATLAB-based
software for topographic analysis and modeling in Earth surface sciences.
*Earth Surface Dynamics* 2, 1–7.

Willett, S.D., McCoy, S.W., Perron, J.T., Goren, L. and Chen, C.-Y. (2014).
Dynamic reorganization of river basins. *Science* 343, 1248765.

## Credits

Originally written by Sam Johnstone (January 2015); maintained since by the
Stanford geomorphology group.
