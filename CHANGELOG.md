# Changelog

## 1.0.0

The first release with a test suite, an install story, and compiled kernels.

The **public interface is unchanged**: every class, method and keyword that
existed before still exists and still means the same thing. What changed is
that a lot of it now works.

---

## Results that change

Read this section before comparing new output against old. Each item is a
place where the previous answer was wrong, so the new number differs on
purpose.

### Chi (`Ksi`)

The integrand was `(Ao / (A - Ao))**theta`; it is now `(Ao / A)**theta`. The
old form diverges as drainage area approaches the threshold `Ao`, so chi was
inflated near every channel head and the error propagated downstream.
Chi is defined as ∫(A₀/A)^θ dx (Perron & Royden 2013).

### Chi (`Chi`)

Chi at the outlet is now exactly 0. It used to be one cell-step of the
integrand, so every chi value in a basin was offset by that amount.

`Chi` also accepts `trapezoid=True`, which averages the integrand across each
step instead of taking it at the upstream cell. That is the rule TopoToolbox
uses; the default remains the original left-endpoint rule.

### Hillshade

The aspect term was `atan2(dz/dy, dz/dx)` where the ESRI definition the code
claimed to implement is `atan2(dz/dy, -dz/dx)`. Illumination was mirrored
east–west: slopes that should have been lit were shadowed and vice versa.
Hillshades are now clipped to 0–255 and stored as `uint8` rather than
wrapping negative values around the byte into bright speckle.

### D8 flow directions

Three separate problems, all fixed by moving to a single shared kernel:

- The windowed comparison stopped two columns and two rows short, so the
  last row and the last column of every grid received **no flow direction at
  all**.
- The diagonal distance factor was `1.41`, not `sqrt(2)`.
- Ties between equally steep neighbours resolved towards the north-east
  because of the order the codes were written; they now resolve
  deterministically in the order E, SE, S, SW, W, NW, N, NE.

`FlowDirectionD8.update_flow_codes_in_mask`, the per-cell re-router, had the
scaling **inverted** — it divided the cardinal drops by 1.41 and left the
diagonals alone — which biased routing towards the diagonals.

### Flow accumulation (`Area`)

The accumulation order came from `flow_direction.sort()`. That returns an
elevation sort only when the grid was built from a filled DEM in the same
session; for a flow-direction grid **loaded from a file** it sorted the cells
by their direction *code*, which is meaningless, and the resulting areas were
wrong. Accumulation now derives its own topological order from the flow
network, so it is correct however the grid was obtained.

The `mask` keyword had no effect: the test was `mask[i, j] is not None`, and
`__getitem__` returns the value, so it was true even for a masked-out cell.
Masks now actually mask.

### `MainstemValleyArea`

Its `__calcD8Area` override was name-mangled to
`_MainstemValleyArea__calcD8Area` while the inherited caller looked for
`_Area__calcD8Area`, so it never ran: the class silently behaved exactly
like `ValleyArea`. It now accumulates along the longest-flow-length path as
intended. (The method also compared a whole grid object against a scalar,
which would have raised the moment it did run.)

### `FlowLength` direction codes

`FlowLength` stored its main-stem directions in a **row-flipped** convention:
internally consistent, but contradicting the ArcGIS convention used by
`FlowDirectionD8` and by the rest of the library. It now uses the same
convention as everything else.

**`*_directions` side-files written by earlier versions must be
regenerated.** `FlowLength.load` cannot detect which convention a file uses.

### Recursive-profile `distance_scale`

In the nested dictionaries from `map_values_to_recursive_list`,
`distance_scale` was written onto the *parent* node and overwritten by each
child in turn, so it ended up describing whichever child was visited last.
Chi profiles built from these dictionaries were wrong at every confluence
where a diagonal and a cardinal tributary met. Each node now carries the
step length from its own parent.

The root node's `index` is now a plain `(row, col)` pair like every other
node's, instead of the nested `((row, col),)` the caller happened to pass.

### Flexural deflection

`(2 * 3.141)**4` is 0.2% low; it is now `(2 * np.pi)**4`. The wavenumber
axes were built from half-width ranges whose length did not have to match
the grid, so the spectral multiply could fail outright; they now come from
`np.fft.fftfreq`.

### Georeferencing of `nx/ny/projection/geo_transform` grids

`__populate_georef_info_using_geoTransform` read `ny` two lines before it
assigned it, so a grid built from an explicit geotransform picked up the
default of 0 and its `yllcenter` came out `ny * dx` too far north. Every
coordinate derived from it — `_rowscols_to_xy`, `_xy_to_rowscols`,
`extent()`, the ArcInfo ASCII header, outlets given in map coordinates —
inherited the error, and the same grid round-tripped through GDAL disagreed
with itself. The other four construction paths were always correct.

### Geographic curvature

`GeographicLaplacian` divided by `dx`, which on a lat/lon grid is in
**degrees**: curvatures came out about 10¹⁰ too large. It now uses the true
metric cell size, like every other geographic class.

### Strahler stream order

The order was applied one donor at a time, which makes a junction of orders
{2, 2, 3} come out as 4 or 3 depending on which donor arrives first — and
the two backends emit different (equally valid) topological orders, so they
disagreed on about 3% of networks. Junctions are now resolved as a whole:
the largest incoming order, plus one if two or more donors carry it.

Note this also differs from TopoToolbox, whose incremental form has the same
order-dependence.

### Smaller behavioural changes

**`search_down_flow_direction_with_length`** now includes the terminal cell
of the path, and follows a path that merely *starts* on the grid perimeter.
The stop condition is "no downstream neighbour inside the grid" rather than
"this cell is on the perimeter", so a cell on the western edge draining east
is now followed instead of returning an empty path.

**`FlowDirectionD8.divides()`** now evaluates the perimeter as well as the
interior. Border cells with no upstream neighbour are divides, and used to
be reported as 0 regardless.

**`Chi(mask=...)`** treats the mask as a boolean gate on the integration
rather than multiplying the finished chi grid by the mask's *values*. The
old behaviour scaled chi by whatever the mask happened to contain, and
carried a hard-coded `area > 1e4` escape hatch.

**`_xy_to_rowscols`** rejects a point exactly one cell beyond the eastern or
southern edge. It used to return the subscript `nx` or `ny`, which raised an
`IndexError` somewhere else entirely.

**`FlowDirectionD8` caches its flow network** (`receivers`,
`topological_order`) because every downstream computation now reads it.
Writing a flow code with `fd[i, j] = ...` or `set_value_at_rowscols` drops
the cache, so hand edits take effect as they always did; call
`fd.invalidate_network()` after modifying `_griddata` directly.

**`tile(..., tile_xpadding=0)`** no longer duplicates the edge tiles' data.
`griddata[:, -0:]` selects the whole array rather than nothing, so a zero
padding produced double-width edge tiles.

**`resample`** warns when SciPy cannot provide the requested interpolation
order and it falls back to linear, instead of quietly changing the answer.

**`Elevation.track_flow_downhill`** and
**`DiscreteFlowAccumulation`** stop on a `NaN` neighbour as before, but the
test is now a real `NaN` check: `None not in numpy_array` never matched, so
off-grid neighbours did not stop the walk the way the code intended.

### Other corrected numerics

| What | Was | Now |
|---|---|---|
| `assignBCs` south-east corner | copied the south-west corner | copies the south-east corner |
| `get_XY_matricies` | returned `n-1` coordinates per axis | returns `n` |
| `resample` | dropped `xllcenter` from the axis stop, and flipped the grid left-right unconditionally | keeps the south-west cell centre fixed; uses `RegularGridInterpolator` |
| `sort(mask=...)` | indexed the mask with positions in the sorted array | selects the sorted positions inside the mask |
| `average_over_distance` | kernel centred half a cell off; `ifftshift` where `fftshift` was needed | kernel centred on cell (0,0), no shift needed |
| `extent_of_data` | scanned rows for the column limits | scans the right axis for each |
| gradient / Laplacian over a length scale | divided by `(2N+1)dx`; float slice indices | divides by `2N·dx`; integer indices |
| `LocalRelief` window | `np.ogrid[-r:r]`, one cell short on the far side | symmetric about the centre |
| `LogArea` | declared `uint8`, truncating log₁₀(area) to integers | `float64` |
| `tile` geotransform | north edge one cell too far north | consistent with the reader |
| `CircularMovingWindow` | squared the row offset twice, giving a vertical band | a disc |
| `datasets` `slope_y` | multiplied the relief by a ramp, zeroing row 0 | adds a ramp |
| geographic cell area | `np.float64(range(n))` — a `TypeError` on modern NumPy | `np.arange` |
| ArcInfo ASCII reader | ignored `NODATA_value`, turning −9999 flags into elevations | maps it to `NaN` |
| ArcInfo ASCII writer | wrote the literal text `nan` | writes −9999 by default |

---

## Crashes fixed

These raised immediately on any current scientific-Python stack.

**Removed upstream APIs.** `np.NAN`, `np.NaN`, `np.float`,
`scipy.ndimage.morphology`, `scipy.ndimage.filters`,
`scipy.interpolate.interp2d`, `np.int`, and `matplotlib.pyplot.hold` were all
in use. A regression test now fails if any of them reappears.

**Python 3 iterators.** `zip()` objects were used as NumPy subscripts and
subscripted directly in `ValueGrid.set_value_at_indexes`,
`CrossDivideDChi`, `NormalizedCrossDivideDChi` and
`FlowDirectionD8.divides_for_outlets`. Several functions returned a `zip`
that callers iterated twice, getting nothing the second time.

**Comparisons with `None`.** `bounds_of_basin_for_outlet` evaluated
`lon > bounds[0][1] or bounds[0][1] is None` — the comparison runs first, so
it raised on its own first iteration. `areas_between` compared `None` with a
float when a cell had no upstream neighbour.

**Names that did not exist.**

| Call | Problem |
|---|---|
| `PriorityQueueMixIn.randomize_subbasins_with_mask` | called `self.__flood`, mangled to `_PriorityQueueMixIn__flood`; the method is `_flood` |
| `MovingWindow.apply_moving_window` | called `self.__adjustKernel`; the method is `__adjust_kernel` |
| `ScaledRelief(flooded_dem=...)` | called `_create_from_flow_direction_flooded_dem_and_elevation`, never defined |
| `ChannelSlopeWithSmoothing(horizontal_interval=...)` | dispatched to `_create_from_elevation_flow_direction`, never defined |
| `BaseSpatialGrid.apply_moving_window` | `out_grid = None` then `out_grid._georef_info = ...` |
| `convert_rivertools_directions_to_arc` | `int()` of an array; referenced an undefined `self.noData` |
| `dem.plot(xlabel=...)` | the local name was only bound when the argument was *absent* |
| `utils.calc_ks_for_outlet` | `kwargs.pop('xo')` with no default |
| `MovingWindow` subclasses | read `self.window_radius`, which the constructor never set |
| `DiscreteFlowAccumulation._mean_pixel_dimension` | required a `kwargs['elevation']` its callers did not pass |
| `GeographicThetaFromChiWithSmoothing` | inherited `KsFromChiWithSmoothing`, so it fitted ks instead of theta |

**Stack overflows.** `get_indexes_of_upstream_cells`,
`map_values_to_recursive_list`, `Chi.__recurse_chi`,
`RestoredElevation.__fill_upstream_points` and
`FlowLength.__get_upstream_indexes` were all recursive, which is why the
module raised the recursion limit to 1 000 000 at import — a limit CPython
segfaults long before reaching. All of them are now iterative, and the
`setrecursionlimit` call is gone.

**Undercounting.** `self._n[rows, cols] += 1` with repeated subscripts
applies each increment only once; the profile-crossing counts in
`KsFromChiWithSmoothing` and `ThetaFromChiWithSmoothing` now use
`np.add.at`.

---

## Thread safety

The compiled kernels release the GIL, but eleven of the thirteen bindings
resolved their NumPy buffers *inside* the released region.
`py::array::request()` calls `PyObject_GetBuffer`, which touches reference
counts, so running a kernel from two threads at once leaked references and
could abort the interpreter. Every binding now resolves and validates its
buffers while the GIL is held and runs nothing but plain C++ afterwards.
`test_kernels_are_thread_safe` pins it down.

The bindings also validate array sizes, so a mismatched `order` array raises
instead of reading past the end of it.

## Added

- **`install.py`** — one command that takes a fresh clone to a working
  installation: it creates a virtual environment, installs the dependencies,
  compiles the kernels, installs GDAL, and then proves the result works by
  round-tripping a georeferenced GeoTIFF. Standard library only, so it runs
  on a machine where nothing is set up yet. `--dry-run` reports what it would
  do without doing it.
- **`TopoAnalysis.gdal_setup`** and the `topoanalysis-gdal` command — find,
  diagnose and install GDAL, meaning *both* the C library and the Python
  bindings. There is no GDAL wheel on PyPI at any version for any platform,
  so pip alone can never finish the job. `doctor` reports what is installed
  and what is wrong with it; `plan` shows what an install would do; `install`
  does it, choosing between an existing C library, conda, a prebuilt Windows
  wheel, the system package manager, or a from-source build.

  Two failure modes it exists to catch:

  - **The silent `gdal_array` trap.** GDAL sdists only declared numpy as a
    build dependency from 3.9 on. Below that, pip's build isolation leaves
    numpy out and the numpy bridge is quietly omitted — `osgeo` imports
    cleanly and `Band.ReadAsArray()` fails at run time. Ubuntu 22.04 (GDAL
    3.4) and 24.04 (3.8) are both on the wrong side of that line. The
    bindings are now installed with `--no-build-isolation` against an
    interpreter that already has numpy.
  - **Architecture mismatch**, such as a Homebrew arm64 GDAL under an x86_64
    Python running via Rosetta. This is detected before anything is built,
    and the package manager that caused it is not asked to fix it.
  - **A stale `PROJ_LIB` or `GDAL_DATA`.** Activating a conda environment
    exports both, and they then outrank the paths compiled into every other
    GDAL on the machine; a PROJ 9 library handed PROJ 8's `proj.db` refuses
    every coordinate lookup. `doctor` names the offending variable and its
    value instead of proposing a reinstall, which would not have helped.

  Repairing broken bindings rebuilds them rather than re-requesting them:
  the usual fault leaves `osgeo` at exactly the right version, so plain
  `pip install gdal==X.Y.Z` reports "already satisfied", and pip's wheel
  cache holds the very build that is broken.

  `install.py` exits non-zero when GDAL was asked for and did not arrive,
  while saying plainly that the rest of the installation is fine. `--no-gdal`
  is how to say it was never wanted.
- **`TopoAnalysis._gdal_source`** — builds SQLite, PROJ and GDAL from pinned,
  SHA-256-checked tarballs into a single prefix, by default the active
  virtual environment. It needs only a compiler and CMake: no root, no
  package manager, and no system `-dev` packages, because GDAL is built
  against its own bundled libtiff, libgeotiff, libpng, libjpeg and zlib.
  Pins can be moved with `TOPOANALYSIS_GDAL_SOURCE_*` environment variables,
  which require a checksum alongside any new URL.
- **`TopoAnalysis.kernels`** — compute kernels in C++ (`TopoAnalysis._cpp`)
  with an equivalent NumPy implementation in `TopoAnalysis.fastops`. The
  test suite asserts the two agree exactly. Set
  `TOPOANALYSIS_PURE_PYTHON=1` to force the NumPy path.
- **`TopoAnalysis.topotoolbox`** — functions that reproduce TopoToolbox 2
  results: `fillsinks`, `identifyflats`, `imposemin`, `flowobj`, `flowacc`,
  `flowdistance`, `drainagebasins`, `streamorder`, `gradient8`, `arcslope`,
  `hillshade`, `aspect`, `curvature`, `cellarea`, `localtopography`,
  `chitransform`, `ksn`, `mchi`, and the referencing helpers.
- **`topoanalysis-process`** — a command-line pipeline. `--info` reports the
  active backend.
- **`python -m TopoAnalysis.benchmark`** — times both backends.
- `FilledElevation(aggradation_slope=0.0)` for a flat fill, and a
  `fill_report` attribute recording how much was filled.
- `FlowDirectionD8.receivers`, `.topological_order`, `.cycle_count`,
  `.basin_mask()` and `.step_lengths()`, which expose the flow network
  directly.
- `FlowDirectionD8(elevation=...)` as an alias for `flooded_dem=`, for DEMs
  conditioned elsewhere.
- Every grid class now accepts the in-memory constructors `dx=..., grid=...`
  and `nx=..., ny=..., dx=...`, not just `BaseSpatialGrid`.
- `BaseSpatialGrid.extent()`, `.coordinate_vectors()` and the correctly
  spelled `.get_XY_matrices()` alias.
- A circular flow path in a hand-edited flow-direction grid is now reported
  (`cycle_count`, and a warning) instead of producing silently wrong results.
- `identifyflats(..., output='all')` also returns the closed basins.
- `FilledElevation`'s flood reports which cells a `maximum_pit_depth` left
  unfilled, so `clip_to_fill` still excises them — the idiom for masking out
  quarries and lakes, which the traversal fix would otherwise have broken.
- `FlowDirectionD8.invalidate_network()`.

## Changed

- **GDAL is imported lazily.** The package imports, and most of it runs,
  without GDAL installed. `gdal_is_available()` reports whether it is there.
  GDAL exceptions are enabled, so a bad path raises instead of returning
  `None` and failing three frames later.
- **statsmodels is imported lazily** and only by the classes that fit
  regressions; it used to print a warning at import time.
- **rasterio is imported lazily** by `analysis.Quadrats.load_data`.
- Both import styles work: `from TopoAnalysis import dem` and, with the
  package directory on `sys.path`, `import dem`.
- `save()` writes `NaN` as the GeoTIFF no-data value and flushes before
  returning.
- `mosaicFolder` uses the GDAL Python API instead of shelling out to
  `gdal_merge.py`, which is not on most installations' `PATH`.
- Progress reporting is quieter and accepts `verbose=False`.
- Bare `except:` clauses — which swallowed `KeyboardInterrupt` and real
  programming errors, and were used for bounds checking in the hot loops —
  are gone.
- The NumPy kernels reject a non-C-contiguous grid rather than filling a
  throwaway copy and reporting success. The compiled kernels always did.
- `from TopoAnalysis import *` no longer names classes that were never
  imported into the package namespace.
- The test suite runs from a checkout under any name — a GitHub zip unpacks
  to `TopoAnalysis-main` — as well as from an installed copy.
