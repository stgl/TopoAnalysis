# TopoAnalysis documentation

Digital elevation model analysis for tectonic geomorphology.

```bash
git clone <this repository>
cd TopoAnalysis
python3 install.py
```

That installs everything, GDAL included — both its C library and its Python
bindings, which pip cannot supply on its own.

## Start here

- **[Installation](installation.md)** — the one-command install, how GDAL is
  obtained on each platform, and the offline and no-compiler cases.
- **[User guide](guide.md)** — the grid model, the analysis chain from DEM
  to channel steepness, and worked examples.
- **[Algorithms](algorithms.md)** — what each computation actually does,
  the conventions it assumes, and where it will disagree with another
  package.
- **[API reference](api.md)** — every public class and function, generated
  from the docstrings.
- **[TopoToolbox parity](topotoolbox_parity.md)** — the function-by-function
  correspondence with the MATLAB package, and the places where exact
  agreement is not achievable.
- **[Changelog](../CHANGELOG.md)** — the defects repaired in 1.0, and the
  handful of results that change because they were wrong.

## Thirty seconds

```python
from TopoAnalysis import Elevation, FilledElevation, FlowDirectionD8, Area

dem    = Elevation(gdal_filename='srtm.tif')
filled = FilledElevation(elevation=dem)       # Priority-Flood, Barnes et al. 2014
d8     = FlowDirectionD8(flooded_dem=filled)  # steepest-descent D8
area   = Area(flow_direction=d8)              # contributing area, m^2
area.save('area.tif')
```

## Two things to know before you start

**Row 0 is the north edge.** Grids are stored the way a GeoTIFF is, so the
first row of `_griddata` is the top of the map and row index increases
southwards. Every convention in the library follows from that: the D8 code
for "south" is 4, `y` decreases as the row index grows, and
`geoTransform[5]` is negative.

**Drainage area is in square metres, not cells.** `Area` accumulates
`dx**2` per cell. TopoToolbox's `flowacc` counts cells. When comparing the
two, or reading a paper that reports one, check which is meant —
`TopoAnalysis.topotoolbox.flowacc` gives the cell count.

## Which backend am I on?

```python
import TopoAnalysis
TopoAnalysis.backend()      # 'c++' or 'python'
```

or from a shell:

```bash
topoanalysis-process --info
```

Both backends produce identical numbers; the C++ one is about 50× faster.
`TOPOANALYSIS_PURE_PYTHON=1` forces the NumPy path.

## Regenerating these pages

`api.md` is generated from the docstrings:

```bash
python docs/generate_api.py     # or: make docs
```

The rest are written by hand. They are plain Markdown and are meant to be
read as they are, on the repository page or in an editor -- there is no
Sphinx build. The `[docs]` extra exists for anyone who wants to wire one up.
