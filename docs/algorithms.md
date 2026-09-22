# Algorithms

What each computation actually does, the conventions it assumes, and where
it will disagree with another package.

---

## Depression filling — Priority-Flood

A DEM straight off a sensor is full of closed depressions: real ones
(playas, craters, karst) and spurious ones (noise, bridges, vegetation
returns). Flow routing needs every cell to have a downhill path to the edge,
so the depressions are filled first.

TopoAnalysis uses the Priority-Flood family of Barnes, Lehman and Mulla
(2014). The idea is simple: start from the cells water can leave through,
and flood inwards, always processing the lowest cell you have reached so far.

```
push every seed cell onto a min-heap, mark it closed
while the heap is not empty:
    c = pop the lowest cell
    for each of c's eight neighbours n that is not closed:
        close n
        raise n to at least c's elevation
        push n
```

Because cells are processed in nondecreasing elevation order, the first time
you reach a cell you have reached it by the lowest possible path, so the
elevation you raise it to is the lowest that lets water out. The filled
surface is the pointwise-minimal surface with no depressions.

### Seeds

Water leaves the DEM at the grid perimeter and at no-data holes, so the seed
set is:

- every cell on the outer ring, and
- every valid cell that touches a `NaN`.

This is the same seed set TopoToolbox's `fillsinks` uses. A depression that
touches a no-data hole is therefore *not* filled — it drains into the hole.

You can supply your own seeds with `outlets=[(x, y), ...]`, which floods
inwards from those points instead.

### The three variants

`TopoAnalysis.kernels.priority_flood` implements all three algorithms from
the paper:

| Variant | Paper | When it is used |
|---|---|---|
| Priority-Flood | Algorithm 1 | `reference_algorithm=True`, as a test oracle |
| Improved Priority-Flood | Algorithm 2 | the default for a flat fill |
| Priority-Flood + ε | Algorithm 4 | whenever `aggradation_slope > 0` |

Algorithm 2 adds a FIFO queue alongside the heap: a cell that needs no
filling sits at the current flood level and cannot be reordered, so it can
skip the heap entirely. That is what makes the practical cost O(n) instead
of O(n log n). The result is identical to Algorithm 1 — the test suite
asserts it on random surfaces.

The FIFO shortcut is only sound while every cell it holds is at exactly the
current level. Two options break that, and both fall back to Algorithm 1:

- **ε filling**, where each filled cell is raised by its own increment; and
- **`maximum_pit_depth`**, which re-queues an unfilled cell *below* the
  current level.

### Flat fill or ε fill?

A flat fill produces exactly level plateaus, which have no downhill
direction — D8 routing has nothing to follow across them. TopoAnalysis
therefore defaults to raising each filled cell by

    parent_elevation + aggradation_slope × step_length

with `aggradation_slope = 1e-12` per map unit. That is far below any
meaningful elevation difference but strictly monotonic, so every filled cell
has a defined downslope neighbour.

Set `aggradation_slope=0` for the flat fill TopoToolbox produces. TopoToolbox
then handles the flats separately, by building an auxiliary "geodesic"
surface inside them so cells drain towards the nearest sill. The two
approaches give *different but equally valid* flow directions inside flats;
see [topotoolbox_parity.md](topotoolbox_parity.md).

### Leaving real basins alone

`maximum_pit_depth=20.0` leaves any depression deeper than 20 m unfilled.
The flood still traverses it — so the terrain upstream is routed — but the
cells keep their elevations, and the basin becomes a local base level.

> The original implementation marked such cells closed without queueing
> them, which orphaned everything upstream of a deep pit.

### Cost

O(n) time and O(n) memory in practice. About 100 ms per million cells in
C++; 5 s in NumPy.

---

## Flow direction — D8

Each cell drains to exactly one of its eight neighbours: the one reached by
the steepest *gradient*,

    slope = (z_centre − z_neighbour) / distance

with `distance = dx` for the four cardinal neighbours and `√2 · dx` for the
four diagonals. Dividing by distance matters: a diagonal step covers 41%
more ground, so an equal drop over a diagonal is a shallower slope. Choosing
by raw drop biases the network towards diagonals and inflates flow lengths.

Ties go to the first direction in the order E, SE, S, SW, W, NW, N, NE. Any
rule is arbitrary; what matters is that it is deterministic, so two runs
agree.

Codes follow the ArcGIS convention, with row 0 at the north edge:

```
| 32  64 128 |          | (i-1,j-1) (i-1,j) (i-1,j+1) |
| 16   X   1 |    ==    | (i  ,j-1)    X    (i  ,j+1) |
|  8   4   2 |          | (i+1,j-1) (i+1,j) (i+1,j+1) |
```

0 means the cell has no downstream neighbour: a pit, an outlet at the grid
edge, or no-data. On a filled DEM the only zeros are edge outlets.

For a geographic (lat/lon) grid the distance uses the local cell size, which
varies with latitude.

---

## The flow network

Once directions are known, the grid is a forest: every cell has one
*receiver* and any number of *donors*, and following receivers always
terminates.

```python
d8.receivers          # flat index of each cell's receiver, or -1
d8.topological_order  # every donor appears before its receiver
```

The order comes from Kahn's algorithm on the donor counts, **not** from
sorting by elevation. That distinction matters: a flow-direction grid loaded
from a file carries no elevations, so an elevation sort is unavailable —
and the original code silently sorted such grids by their *direction codes*,
producing wrong accumulations.

Kahn's algorithm also detects cycles, which a hand-edited or externally
produced flow-direction grid can contain. `d8.cycle_count` reports them and
a warning is raised; a valid network has none.

Given the order, everything downstream is one linear pass:

| Quantity | Pass | Recurrence |
|---|---|---|
| accumulation | forward | `A[recv] += A[k]` |
| flow length | forward | `L[recv] = max(L[recv], L[k] + step)` |
| distance to outlet | backward | `D[k] = D[recv] + step` |
| chi | backward | `χ[k] = χ[recv] + integrand · step` |
| basin labels | backward | `label[k] = label[recv]` |
| minima imposition | forward | `z[recv] = min(z[recv], z[k] − sl·step)` |
| Strahler order | forward | the order-increment rule |

---

## Drainage area

`Area` accumulates the per-cell area downstream, so each cell holds its own
area plus every cell upstream of it, in **square map units**.

TopoToolbox's `flowacc` accumulates unit weights instead and reports the
**number of cells**. `TopoAnalysis.topotoolbox.flowacc` does the same; the
two differ by exactly `dx²` on a projected grid.

For a geographic grid, `GeographicArea` uses the true spherical cell area

    A(lat) = R² · Δλ · (sin(lat + Δφ/2) − sin(lat − Δφ/2))

with `R = 6 371 000 m`. Cells shrink towards the poles: at 60° latitude they
are half the area of equatorial ones, and using a constant `dx²` would
overestimate drainage area by that factor.

### Gating

`evaluate_at` or `mask` restrict accumulation to a subset of cells. A gated-
out cell contributes nothing and passes nothing downstream — the first cell
*outside* the gate still receives what the gate delivers, and the flow
stops there. `ValleyArea` uses this to accumulate only over valley floors.

---

## Flow length and the main stem

`FlowLength` records, for each cell, the length of the **longest** flow path
reaching it — the standard "flow length" of a basin, and the distance used
to define its relief.

Alongside the distances it records which upstream neighbour supplies each
cell's longest path. Chaining those donors from an outlet traces the main
stem. Ties are broken towards the donor with the lower flat index, so the
main stem does not depend on which valid topological order was used.

`Relief`, `ScaledRelief` and `Ksi` all propagate along that main-stem
network, which is why they need a `FlowLength` as well as a flow direction.

> Before version 1.0 these codes used a row-flipped convention, internally
> consistent but incompatible with `FlowDirectionD8`. `*_directions` files
> written by older versions must be regenerated.

---

## Chi

For a river in steady state under the stream-power law, elevation is linear
in the integral quantity

    χ(x) = ∫[outlet → x] (A₀ / A(x'))^θ dx'

(Perron & Royden 2013). The slope of elevation against χ is the steepness
index k_sn, and A₀ makes χ dimensionally a length — usually 1 km² = 10⁶ m².

Two implementations, because they answer different questions.

**`Chi`** integrates upstream from outlets you name, over the whole network
above each one. Use it when the basins are the unit of analysis.

**`Ksi`** integrates along longest-flow-length paths across the whole grid,
without outlets. Use it for a regional map.

### The integration rule

Per step from a cell to its receiver:

- *left endpoint* (the default): `χ[k] = χ[recv] + (A₀/A[k])^θ · dl`
- *trapezoid* (`trapezoid=True`): averages the integrand across the step,
  `χ[k] = χ[recv] + ½((A₀/A[k])^θ + (A₀/A[recv])^θ) · dl`

TopoToolbox's `chitransform` uses the trapezoid rule. Area grows downstream
so the integrand shrinks downstream, which makes the left-endpoint rule the
larger of the two. Use `trapezoid=True` when comparing.

χ is 0 at each outlet and grows upstream.

> The `Ksi` integrand used to be `(A₀ / (A − A₀))^θ`, which diverges as area
> approaches the threshold. It is now `(A₀ / A)^θ`.

---

## Terrain derivatives

### Slope

`MaxSlope` is the magnitude of the centred-difference gradient,
`√((dz/dx)² + (dz/dy)²)`, with the boundary padded by replicating the edge.

`topotoolbox.gradient8` is the steepest *downward* gradient over the eight
neighbours, clamped at 0 — a different quantity, always non-negative, 0 in
pits and flats.

`topotoolbox.arcslope` is the ArcGIS/Horn 3×3 estimator: a 1-2-1 weighted
stencil over an `8·dx` denominator, with no √2 anywhere. It is less
sensitive to single-cell noise than a plain centred difference.

### Measuring over a length scale

On a 1 m lidar DEM a 3×3 stencil measures the noise, not the landscape.
`calculate_gradient_over_length_scale(L)` and
`calculate_laplacian_over_length_scale(L)` use a stencil reaching
`N = ⌈L/dx⌉` cells either side, so the sampled cells are `2N·dx` apart. A
band `N` cells wide around the edge comes back as `NaN`.

### Curvature

`Laplacian` is `∇²z` from SciPy's 5-point stencil, scaled by `1/dx²`.
Positive in valleys, negative on ridges, and the basis of the
valley-floor masks that `ValleyArea` uses.

`principal_curvatures()` returns the two principal curvatures of the surface
from the second fundamental form. `topotoolbox.curvature` offers the five
TopoToolbox flavours — profile, plan, tangential, mean and total — from the
Evans 3×3 quadratic fit.

### Hillshade

`Hillshade` follows the ESRI definition and returns 0–255:

    H = 255 · (cos(zenith)·cos(slope) + sin(zenith)·sin(slope)·cos(azimuth − aspect))

with `zenith = 90° − altitude`, `azimuth_math = 360° − azimuth + 90°`, and
`aspect = atan2(dz/dy, −dz/dx)`. Values below 0 — slopes facing away from
the light — are clamped to 0.

> The `−dz/dx` is easy to miss and easy to get wrong; without it the
> illumination is mirrored east–west, which is what the original did.

`topotoolbox.hillshade` returns the raw cosine in [−1, 1] from a
central-difference surface normal, which is what TopoToolbox does (it only
scales to 0–255 when plotting).

---

## Valley extraction

`ValleyArea` builds a valley-floor mask in three steps:

1. **Convergence.** Cells whose Laplacian exceeds `valley_laplace_value`.
2. **Connectivity.** `PriorityFillGrid` floods outwards from cells with more
   than `min_area_value` of drainage area through the convergent cells, so
   only the parts that connect to a real channel survive. Isolated
   convergent patches on hillslopes are dropped.
3. **Closing** (optional, `iterations=n`). Dilate then erode, to bridge
   one-cell gaps without growing the network.

Drainage area is then accumulated over that mask only.
`MainstemValleyArea` restricts it further, to the longest flow path.

`MultiscaleCurvatureValleyWidth` takes a different route: it fits a local
quadratic at a range of window sizes and records, at each cell, the window
that minimises the smaller principal curvature. That scale is the valley
width, because a window matched to the valley resolves it best.

---

## Restored topography

`RestoredElevation` asks what the landscape would look like if channel
steepness were uniform. Starting from the outlets it integrates

    dz/dx = k_s · A^(−θ)

upstream, then lets each divide migrate towards whichever side ended up
lower, recomputes drainage area, and repeats. Comparing the result with the
real DEM shows where the landscape is out of steady state.

`randomize=True` randomises elevations inside the basin first and re-routes,
so the result does not inherit the real network's planform.

---

## Flexural isostasy

`Deflection` solves the thin elastic plate equation spectrally:

    w(k) = −ρ_c g h(k) / [(ρ_m − ρ_c) g + D (2π|k|)⁴]

Long wavelengths approach the Airy limit `−h ρ_c/(ρ_m − ρ_c)`; short ones
are supported by plate strength. Because it is spectral the load is
implicitly periodic, so pad the DEM if edge effects matter.

Pass `restored_elevation` as well to get the *change* in deflection between
two surfaces — the isostatic response to the erosion between them.

---

## Numerical conventions

| | |
|---|---|
| Row order | row 0 is the **north** edge |
| No-data | `NaN` in floating-point grids; 0 in flow-direction grids |
| Cell registration | coordinates refer to cell **centres** |
| `geoTransform` | GDAL's 6-tuple, origin at the **outer corner** of the NW cell |
| Diagonal distance | `√2 · dx`, never 1.41 |
| Precision | float64 throughout (TopoToolbox often works in float32) |
| Cells | assumed square; `dx` is used for both axes |

---

## References

Barnes, R., Lehman, C. and Mulla, D. (2014). Priority-flood: An optimal
depression-filling and watershed-labeling algorithm for digital elevation
models. *Computers & Geosciences* 62, 117–127.

Horn, B.K.P. (1981). Hill shading and the reflectance map. *Proceedings of
the IEEE* 69, 14–47.

O'Callaghan, J.F. and Mark, D.M. (1984). The extraction of drainage networks
from digital elevation data. *Computer Vision, Graphics and Image
Processing* 28, 323–344.

Perron, J.T. and Royden, L. (2013). An integral approach to bedrock river
profile analysis. *Earth Surface Processes and Landforms* 38, 570–576.

Schwanghart, W. and Scherler, D. (2014). TopoToolbox 2 — MATLAB-based
software for topographic analysis and modeling in Earth surface sciences.
*Earth Surface Dynamics* 2, 1–7.

Willett, S.D., McCoy, S.W., Perron, J.T., Goren, L. and Chen, C.-Y. (2014).
Dynamic reorganization of river basins. *Science* 343, 1248765.
