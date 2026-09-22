# User guide

## The grid model

Everything in TopoAnalysis is a grid. A grid holds a 2-D NumPy array and the
georeferencing that places it on the ground:

```python
grid._griddata            # the array; row 0 is the NORTH edge
grid._georef_info.dx      # cell size, in map units
grid._georef_info.nx, .ny # columns, rows
grid._georef_info.xllcenter, .yllcenter   # centre of the SOUTH-WEST cell
grid._georef_info.geoTransform            # the GDAL 6-tuple
grid._georef_info.projection              # WKT, or 0 if unset
```

`NaN` is no-data in floating-point grids.

Subclasses add meaning rather than structure. `Elevation` is a grid of
heights, `FlowDirectionD8` a grid of direction codes, `Area` a grid of
contributing areas — all of them are `BaseSpatialGrid` underneath and share
its clipping, tiling, resampling and plotting.

## Constructing grids

Grids are built by keyword, and the combination you pass selects the
constructor:

```python
from TopoAnalysis import Elevation

Elevation(gdal_filename='dem.tif')                      # from a file
Elevation(ai_ascii_filename='dem.asc', EPSGprojectionCode=32611)
Elevation(dx=30.0, grid=array)                          # from an array
Elevation(nx=100, ny=80, projection=wkt, geo_transform=gt)   # empty
Elevation(nx=100, ny=80, dx=30.0)                       # random values
```

Derived grids name the grid they come from:

```python
FilledElevation(elevation=dem)
FlowDirectionD8(flooded_dem=filled)
Area(flow_direction=d8)
Ksi(area=area, flow_direction=d8, flow_length=length, theta=0.45, Ao=1e6)
```

Each class lists what it accepts in `required_inputs_and_actions`, and an
unrecognised combination raises an error that prints the list.

## The analysis chain

```python
from TopoAnalysis import (Elevation, FilledElevation, FlowDirectionD8,
                          Area, FlowLength, Chi, Ksi, ScaledRelief)

dem    = Elevation(gdal_filename='srtm.tif')
filled = FilledElevation(elevation=dem)
d8     = FlowDirectionD8(flooded_dem=filled)
area   = Area(flow_direction=d8)
length = FlowLength(flow_direction=d8)
```

Each step depends on the one before. Filling removes closed depressions so
every cell has a downhill path out; routing turns that surface into a
network; accumulation and length are properties of the network.

### Filling

```python
filled = FilledElevation(elevation=dem)                        # default
filled = FilledElevation(elevation=dem, aggradation_slope=0.0) # flat
filled = FilledElevation(elevation=dem, maximum_pit_depth=20.0)
filled = FilledElevation(elevation=dem, outlets=[(x, y)])
print(filled.fill_report)   # {'cells_filled': ..., 'max_fill_depth': ...}
```

By default the filled surface carries a very slight gradient
(`aggradation_slope = 1e-12` per map unit) across what were depressions, so
D8 routing has a direction to follow. `aggradation_slope=0` gives a flat
fill, which is what TopoToolbox's `fillsinks` returns — but then the flats
have no defined flow direction.

`maximum_pit_depth` leaves genuinely closed basins — playas, craters —
unfilled, while still routing the terrain beyond them.

### Flow direction

D8 codes follow the ArcGIS convention, with row 0 at the north:

```
| 32  64 128 |
| 16   X   1 |
|  8   4   2 |
```

0 means the cell has no downstream neighbour: a pit, an outlet at the grid
edge, or no-data.

The grid also exposes the network directly, which is what makes the
downstream computations linear-time:

```python
d8.receivers            # flat index of each cell's receiver, or -1
d8.topological_order    # donors before receivers
d8.cycle_count          # 0 for a valid network
d8.step_lengths()       # distance from each cell to its receiver
d8.basin_mask([(x, y)]) # boolean mask of everything draining to a point
```

### Drainage area

```python
area = Area(flow_direction=d8)              # square map units
```

Each cell contributes its own area plus everything upstream. For a
geographic (lat/lon) grid use `GeographicArea`, which uses the true spherical
cell area — cells shrink towards the poles.

To count *cells* instead, which is what TopoToolbox reports:

```python
from TopoAnalysis import topotoolbox as tt
cells = tt.flowacc(d8)
```

## Chi and channel steepness

Chi is the integral

    chi(x) = integral from the outlet to x of (A0 / A)^theta dx'

and a channel in steady state plots as a straight line in chi–elevation
space, with slope k_sn. There are two ways to compute it here.

**`Chi`** integrates upstream from outlets you name, over the whole network
above each one:

```python
outlets = area.areas_greater_than(1e8)
chi = Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e6, outlets=outlets)
```

**`Ksi`** integrates along longest-flow-length paths across the whole grid,
with no outlets needed:

```python
ksi = Ksi(area=area, flow_direction=d8, flow_length=length,
          theta=0.45, Ao=1e6)
```

Pair either with the matching relief grid and the slope of the resulting
plot is the steepness index:

```python
relief = ScaledRelief(flow_direction=d8, elevation=dem, flow_length=length,
                      Ao=1e6, theta=0.45, area=area)

from TopoAnalysis import plot
plot(ksi, relief, xlabel='chi (m)', ylabel='scaled relief (m)')
```

`Chi` takes `trapezoid=True` to average the integrand across each step
instead of taking it at the upstream cell. That is the rule TopoToolbox uses,
and the one to pick when comparing against it.

### Fitting steepness cell by cell

`KsFromChiWithSmoothing` fits chi–elevation over a moving window along each
flow path, giving a map of k_sn rather than one number per basin:

```python
ks = KsFromChiWithSmoothing(elevation=dem, area=area, flow_direction=d8,
                            theta=0.45, vertical_interval=50.0,
                            area_threshold=1e6)
ks.save('ksn.tif')     # 7 bands: ks, n, mse, ss, r2, p, n_regression
```

Use `horizontal_interval` instead of `vertical_interval` to set the window
by distance rather than by drop. `ThetaFromChiWithSmoothing` fits the
concavity instead of holding it fixed. Both need `statsmodels`.

### Profiles

For individual profiles rather than maps, extract the basin as a nested
dictionary and work on that:

```python
from TopoAnalysis import demRecursionTools, plotting

tree = d8.map_values_to_recursive_list(outlet, elevation=dem, area=area)
e, c = demRecursionTools.chi_elevation(tree, area._mean_pixel_dimension(),
                                       [0.45], xo=1000.0)
ks, r2 = demRecursionTools.best_ks_with_r2_list(tree, de, [0.45])

plotting.plot_profiles(dem, d8, area, outlet, 'k-')
plotting.plot_chi_profiles(dem, d8, area, outlet, 'k-', theta=0.45)
```

Each node of the tree carries `index`, `distance_scale` (the step length
from its parent, 1 or sqrt(2)), `next` (its upstream children), and one entry
per grid you passed as a keyword.

## Terrain attributes

```python
from TopoAnalysis import Hillshade, MaxSlope, Laplacian, LocalRelief

Hillshade(elevation=dem, azimuth=315, inclination=45)   # ESRI, 0-255
MaxSlope(elevation=dem)                                 # gradient magnitude
Laplacian(elevation=dem)                                # + in valleys
LocalRelief(elevation=dem, pixel_radius=25)             # range in a disc

dem.principal_curvatures()                              # (k1, k2)
dem.calculate_gradient_over_length_scale(500.0)         # (dz/dx, dz/dy)
dem.calculate_laplacian_over_length_scale(500.0)
```

Measuring over an explicit length scale is usually what you want on a noisy
DEM: a 3×3 stencil on 1 m lidar measures the noise, not the landscape.

## Working with large DEMs

The compiled kernels handle a few million cells comfortably. Beyond that,
process in tiles:

```python
tiles = dem.tile(tile_xdim=2000, tile_ydim=2000,
                 tile_xpadding=100, tile_ypadding=100)

results = []
for tile in tiles:
    filled = FilledElevation(elevation=tile)
    d8 = FlowDirectionD8(flooded_dem=filled)
    results.append(Area(flow_direction=d8).remove_padding(100, 100))

mosaic = Area.mosaic(results)
```

The padding gives each tile enough context that the interior is correct; it
is stripped before the tiles are reassembled. Pick padding larger than the
longest flow path you need to resolve within a tile — for drainage area,
which integrates over the whole basin, tiling is only valid when basins fit
inside a tile.

## Choosing a backend

```python
import TopoAnalysis
TopoAnalysis.backend()      # 'c++' or 'python'
```

Set `TOPOANALYSIS_PURE_PYTHON=1` to force the NumPy kernels. The two produce
identical results — the test suite asserts it — so this is for debugging and
for comparing, not for changing answers.

## Reading and writing

```python
grid.save('out.tif')                  # LZW-compressed GeoTIFF, NaN as nodata
grid.write_to_ai('out.asc')           # ArcInfo ASCII, -9999 as nodata
grid.vectorize('polygons')            # -> polygons.shp
Elevation.load('out.tif')
```

`FlowLength.save` also writes a `*_directions` side-file holding its
main-stem direction codes; `FlowLength.load` expects to find it.

## Cross-checking against TopoToolbox

See [topotoolbox_parity.md](topotoolbox_parity.md). The short version:

```python
from TopoAnalysis import topotoolbox as tt

tt.fillsinks(dem)         # flat fill
tt.flowobj(dem)           # FLOWobj(DEM,'preprocess','fill')
tt.flowacc(fd)            # in CELLS, not m^2
tt.chitransform(fd, a)    # trapezoidal, mn=0.45, a0=1e6
tt.ksn(fd, dem, a)        # carves at 1e-5 first, no a0
```
