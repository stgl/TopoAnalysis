# Parity with TopoToolbox

[TopoToolbox 2](https://topotoolbox.wordpress.com/) (Schwanghart & Scherler
2014) is the MATLAB package that covers much of the same ground as
TopoAnalysis. This page says, function by function, where the two agree,
where they differ, and what to do about it.

`TopoAnalysis.topotoolbox` holds functions written to reproduce
TopoToolbox's numbers exactly. The rest of TopoAnalysis keeps its own
conventions — area in square metres, ESRI hillshade, left-endpoint chi —
because those are what the rest of the library is built on.

```python
from TopoAnalysis import dem, topotoolbox as tt

z  = dem.Elevation(gdal_filename='srtm.tif')
fd = tt.flowobj(z)          # == FLOWobj(DEM,'preprocess','fill')
a  = tt.flowacc(fd)         # == flowacc(FD)
c  = tt.chitransform(fd, a) # == chitransform(S,A)
```

---

## Correspondence table

**Agreement** is one of:

- **exact** — the same arithmetic on the same inputs;
- **exact¹** — exact up to float64-vs-float32 rounding (TopoToolbox works in
  single precision when the DEM is single);
- **equivalent** — the same quantity, related by a stated conversion;
- **differs** — genuinely different results; the reason is given.

### Depressions and flats

| TopoToolbox | TopoAnalysis | Agreement | Notes |
|---|---|---|---|
| `fillsinks(DEM)` | `tt.fillsinks(dem)` | exact | Flat fill; seeds are the perimeter plus every cell touching NaN; 8-connected. |
| | `FilledElevation(elevation=dem)` | differs | Adds `aggradation_slope` (default 1e-12/m) across filled areas. Pass `aggradation_slope=0` for the flat fill. |
| `fillsinks(DEM, maxdepth)` | `tt.fillsinks(dem, maxdepth=...)` | exact | All-or-nothing per connected depression, iterated as MATLAB does — not a per-cell depth test. |
| `identifyflats(DEM)` | `tt.identifyflats(dem)` | exact | Returns `(flats, sills)`; `output='all'` adds the closed basins. Includes MATLAB's cleared rim and no-data halo, and its treatment of an isolated pit as a flat. |
| `imposemin(FD, DEM, sl)` | `tt.imposemin(fd, dem, sl)` | exact | Carves rather than fills: only lowers cells. |
| `elevateminima` | — | — | Not implemented. |

### Flow

| TopoToolbox | TopoAnalysis | Agreement | Notes |
|---|---|---|---|
| `FLOWobj(DEM,'preprocess','fill')` | `tt.flowobj(dem)` | **differs inside flats** — see below | Identical elsewhere. |
| `FLOWobj(DEM,'preprocess','carve')` | — | — | Not implemented; `tt.flowobj` raises rather than guessing. |
| `FLOWobj2GRIDobj(FD)` | `d8._griddata` | exact | Same ArcGIS codes, same row-0-is-north orientation. |
| `flowacc(FD)` | `tt.flowacc(fd)` | exact | **Number of cells**, minimum 1. |
| | `Area(flow_direction=d8)` | equivalent | `Area = flowacc × dx²`. |
| `flowacc(FD, W0)` | `tt.flowacc(fd, weights)` | exact | |
| `flowdistance(FD,'upstream')` | `tt.flowdistance(fd,'upstream')` | exact¹ | Distance to the outlet. TopoToolbox accumulates in float32. |
| `flowdistance(FD,'downstream')` | `tt.flowdistance(fd,'downstream')` | exact¹ | Longest path reaching each cell. |
| | `FlowLength(flow_direction=d8)` | exact¹ | The same quantity as `'downstream'`. |
| `drainagebasins(FD)` | `tt.drainagebasins(fd)` | **partition** exact, labels differ | Both number basins in edge-list order; the orders differ. Compare partitions. |
| `drainagebasins(FD,IX)` | `tt.drainagebasins(fd, outlets=...)` | exact | Labels are positions in the outlet list. |
| `streamorder(S,'strahler')` | `tt.streamorder(fd, mask)` | exact | |
| `streamorder(S,'shreve')` | `tt.streamorder(fd, mask, 'shreve')` | exact | |

### Terrain attributes

| TopoToolbox | TopoAnalysis | Agreement | Notes |
|---|---|---|---|
| `gradient8(DEM)` | `tt.gradient8(dem)` | exact | Steepest descent over 8 neighbours, clamped at 0. Diagonals use `√2·dx`. |
| `gradient8(DEM,unit)` | `tt.gradient8(dem, unit)` | exact | `tangent`, `degree`, `radian`, `sine`, `percent`. |
| `arcslope(DEM)` | `tt.arcslope(dem)` | exact² | Horn 3×3, `8·dx` denominator. ²NaN infill uses SciPy's distance transform; ties next to a NaN hole can break differently from MATLAB's `bwdist`. |
| | `MaxSlope(elevation=dem)` | differs | Centred-difference gradient magnitude — a different estimator. |
| `hillshade(DEM)` | `tt.hillshade(dem)` | exact | Defaults azimuth 315, altitude 60 (the *code's* defaults; the MATLAB docstring disagrees with itself). Returns a cosine in [−1, 1]. |
| | `Hillshade(elevation=dem, azimuth=, inclination=)` | differs | ESRI formulation: 0–255, and an aspect of `atan2(dz/dy, −dz/dx)` from centred differences. Strongly correlated with the TopoToolbox version but not equal. |
| `hillshade(...,'method','mdow')` | — | — | Not implemented. |
| `aspect(DEM)` | `tt.aspect(dem)` | exact | Degrees clockwise from north; flat cells report 90. |
| `aspect(DEM,true)` | `tt.aspect(dem, classify=True)` | exact | The 8-class Gómez-Plaza ranking. |
| `curvature(DEM,type)` | `tt.curvature(dem, ctype)` | exact | `profc`, `planc`, `tangc`, `meanc`, `total`. |
| | `Laplacian(elevation=dem)` | differs | `∇²z` from a 5-point stencil — a different quantity from any of the five. |
| `cellarea(DEM)` | `tt.cellarea(grid)` | exact³ | ³Assumes a sphere of total area 510 072 000 km². `GeographicArea` uses R = 6 371 000 m, which is 1.5 × 10⁻⁵ smaller. |
| `localtopography(DEM,r,'range')` | `tt.localtopography(dem, r, 'range')` | **≈** | TopoToolbox's default disc is a periodic-line approximation; this uses an exact disc (TopoToolbox's `N=0`). `'mean'` averages over that binary disc where MATLAB's `fspecial('disk')` weights the rim fractionally; `'std'` uses the N−1 sample deviation, as `stdfilt` does. |
| | `LocalRelief(elevation=dem, pixel_radius=)` | exact⁵ | Identical to `tt.localtopography(..., 'range')` with the matching radius, on grids with no no-data. ⁵`LocalRelief` does not mask `NaN`, so around a hole the two differ; use `tt.localtopography` there. |
| `roughness`, `evansslope` | — | — | Not implemented. |

### Chi and steepness

| TopoToolbox | TopoAnalysis | Agreement | Notes |
|---|---|---|---|
| `chitransform(S,A)` | `tt.chitransform(fd, a)` | exact | `mn=0.45`, `a0=1e6`, trapezoidal, area in **cells**. |
| | `Chi(..., trapezoid=True)` | exact | Given the same outlets and area in m². |
| | `Chi(...)` (default) | differs | Left-endpoint rule instead of trapezoidal. |
| | `Ksi(...)` | differs | Integrates along longest-flow-length paths, not upstream from outlets. |
| `cumtrapz(S,g)` | — | equivalent | The integration inside `tt.chitransform`. |
| `ksn(S,DEM,A)` | `tt.ksn(fd, dem, a)` | exact | Carves at `sl=1e-5` first; no `a0` (reference area 1 m²). |
| `mchi(S,DEM,chi)` | `tt.mchi(fd, dem, chi)` | exact | `dz/dχ`, no carving. |
| `chiplot`, `slopearea` | — | — | Not implemented; `KsFromChiWithSmoothing` is the nearest equivalent. |

### Referencing

| TopoToolbox | TopoAnalysis | Agreement | Notes |
|---|---|---|---|
| `getcoordinates(DEM)` | `tt.getcoordinates(grid)` | exact | `y` runs **north to south**, matching row order. |
| | `grid.coordinate_vectors()` | differs in order | `y` runs south to north. |
| `getextent(DEM)` | `tt.getextent(grid)` | exact | Outermost **cell centres**. |
| | `grid.extent()` | differs | Outer pixel **edges** — what `imshow` wants; half a cell larger each way. |
| `coord2sub`, `coord2ind` | `tt.coord2sub(grid, x, y)` | exact⁴ | ⁴0-based here, 1-based in MATLAB. |
| | `grid._xy_to_rowscols(v)` | **≈** | Rounds halves to even (Python) rather than away from zero (MATLAB); differs only for points landing exactly on a cell boundary. |
| `sub2coord`, `ind2coord` | `tt.sub2coord(grid, r, c)` | exact⁴ | |
| `GRIDobj2ascii` | `grid.write_to_ai(path)` | equivalent | Both centre-registered; TopoAnalysis writes `-9999` for no-data. |
| `GRIDobj2geotiff` | `grid.save(path)` | equivalent | TopoAnalysis adds LZW compression and writes NaN as the no-data value. |

---

## Where exact agreement is not achievable

### Flow directions inside flats

This is the one difference that changes a drainage network.

After filling, a depression becomes a plateau with no downhill direction.
The two packages resolve it differently:

- **TopoToolbox** fills flat and then builds an *auxiliary topography* — a
  grey-weighted geodesic distance transform inside each flat — so cells
  drain towards the nearest sill. Flow directions inside a flat come from
  the ranking that surface induces.
- **TopoAnalysis** adds a tiny gradient while filling
  (`aggradation_slope`), so each cell's direction points back along the
  path the flood arrived by.

Both give a valid, connected network, and both drain each flat to the same
sill. **Which cells they route through inside the flat is different.**
Drainage area therefore differs inside flats and immediately downstream of
them. Everything outside flats is identical.

There is no way to reconcile this short of implementing the geodesic
transform. If exact agreement inside flats matters, compute the flow
directions in MATLAB and load them:

```python
d8 = FlowDirectionD8(gdal_filename='flowdir_from_matlab.tif')
```

The code convention is identical, so the grid transfers without conversion.

TopoToolbox's default is actually `'carve'`, not `'fill'` — it *lowers* the
DEM through depressions rather than raising it. That produces a different
surface again, and is not implemented here.

### Tie-breaking on equal gradients

Where two neighbours are equally steep, TopoToolbox prefers the diagonal
(its test is `G1 <= G2`) and then the larger linear index; TopoAnalysis
takes the first in the order E, SE, S, SW, W, NW, N, NE. Exact ties are rare
on a real DEM with an ε-filled surface, but common on synthetic or
integer-valued ones.

### Border cells of `hillshade` and `aspect`

None — these agree exactly. MATLAB's `surfnorm` expands the grid with a
*quadratic* ghost cell (`3·z₁ − 3·z₂ + z₃`), and a centred difference across
that ghost is algebraically the second-order one-sided difference, which is
what `numpy.gradient(..., edge_order=2)` computes.
`test_surfnorm_matches_matlab_including_the_border` checks it against a
literal port of `surfnorm.m`.

### Precision

TopoToolbox works in whatever class the DEM is, usually float32 from a
GeoTIFF; TopoAnalysis works in float64 throughout. Accumulated quantities —
flow length especially, which sums thousands of terms — drift apart in the
last few digits. Compare with a relative tolerance of about 1e-6, not
exactly.

### Basin labels

Both packages number basins in the order their edge list happens to meet
them, and the two orders differ. The *partition* is identical. Compare
partitions:

```python
import numpy as np

def same_partition(a, b):
    """True when two labellings group the cells identically."""
    pairs = np.unique(np.stack([a.ravel(), b.ravel()]), axis=1)
    return (len(np.unique(pairs[0])) == pairs.shape[1]
            and len(np.unique(pairs[1])) == pairs.shape[1])
```

### Earth radius

`tt.cellarea` uses TopoToolbox's sphere (total area 510 072 000 km², so
R = 6 371 047 m); `GeographicArea` uses R = 6 371 000 m. Areas differ by
1.5 × 10⁻⁵ — negligible for most work, but it will show up in an exact
comparison. Setting the radius brings them to within 3 × 10⁻¹²:

```python
grid.earth_radius = (tt.EARTH_SURFACE_AREA_KM2 / (4 * np.pi)) ** 0.5 * 1000.0
grid._GeographicGridMixin__app = None     # the area grid is cached
```

The residue is because `GeographicArea` evaluates the longitude difference
per cell, so its areas vary in the last few bits along a row, whereas
`cellarea` computes one value per row and broadcasts it.

---

## Cross-checking in practice

A recipe for comparing a result against MATLAB:

```python
import numpy as np
from TopoAnalysis import dem, topotoolbox as tt

z  = dem.Elevation(gdal_filename='srtm.tif')
fd = tt.flowobj(z)
a  = tt.flowacc(fd)                       # cells, like flowacc(FD)
c  = tt.chitransform(fd, a, mn=0.45, a0=1e6)

matlab = dem.BaseSpatialGrid(gdal_filename='chi_from_matlab.tif')
difference = c._griddata - matlab._griddata
print(np.nanmax(np.abs(difference)),
      np.nanpercentile(np.abs(difference), 99))
```

Expect agreement to ~1e-6 relative away from flats, and visible differences
inside them. If the difference is a constant factor, check the units: area
in cells versus m², or `a0` versus `a0^mn`.

## How this was verified

Without MATLAB available, parity was established three ways:

1. **From the source.** Every function above was transcribed from the
   TopoToolbox MATLAB source in this repository — the exact stencils,
   constants, default values, rounding rules and NaN handling — and each
   transcription was independently re-checked against the source.
2. **Against closed forms.** Planes, quadratics and spherical cells have
   answers that can be written down: a 1:1 east-facing plane has gradient
   exactly 1, aspect exactly 90°, curvature exactly 0. `test_topotoolbox.py`
   checks those.
3. **Against each other.** Where TopoAnalysis and the parity module compute
   the same quantity by different routes, the tests assert the documented
   conversion holds — `flowacc × dx² == Area`, `chitransform ==
   Chi(trapezoid=True)`, `localtopography('range') == LocalRelief`.

What this does **not** establish is bit-for-bit agreement with a specific
MATLAB release. If you have MATLAB, running the same DEM through both and
comparing is worth doing, and the differences above are what to expect.
