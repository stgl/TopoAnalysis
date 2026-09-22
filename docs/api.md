# API reference

Generated from the docstrings by `docs/generate_api.py`, so
`help(TopoAnalysis.Area)` in a session shows the same text. For how the
pieces fit together start with [guide.md](guide.md); for what the
computations actually do see [algorithms.md](algorithms.md).

Attributes prefixed with a single underscore (`_griddata`,
`_georef_info`) are not private in practice -- they are how you reach
the data -- but they are not listed here; see
[guide.md](guide.md#the-grid-model).

## `TopoAnalysis.dem` -- the grid classes

Contents:

- **Core grid types**: `BaseSpatialGrid`, `Georef_info`, `BaseSpatialShape`, `ValueGrid`, `Elevation`, `GeographicElevation`
- **Depression filling**: `PriorityQueueMixIn`, `FilledElevation`, `PriorityFillGrid`
- **Flow routing**: `FlowDirection`, `FlowDirectionD8`, `FlowLength`, `GeographicFlowLength`, `Mask`, `DiscreteFlowAccumulation`, `GeographicDiscreteFlowAccumulation`
- **Drainage area**: `Area`, `GeographicArea`, `LogArea`, `ValleyArea`, `GeographicValleyArea`, `MainstemValleyArea`, `GeographicMainstemValleyArea`
- **Chi and channel steepness**: `Chi`, `GeographicChi`, `Ksi`, `GeographicKsi`, `Relief`, `ScaledRelief`, `ChiScaledRelief`, `ChannelSlope`, `CrossDivideDChi`, `NormalizedCrossDivideDChi`, `RestoredElevation`, `GeographicRestoredElevation`
- **Along-flow smoothing**: `AlongFlowSmoothing`, `KsFromChiWithSmoothing`, `GeographicKsFromChiWithSmoothing`, `ThetaFromChiWithSmoothing`, `GeographicThetaFromChiWithSmoothing`, `ChannelSlopeWithSmoothing`, `ChannelDownSlopeWithSmoothing`, `ChannelUpSlopeWithSmoothing`
- **Terrain attributes**: `CalculationMixin`, `Hillshade`, `GeographicHillshade`, `MaxSlope`, `GeographicMaxSlope`, `Laplacian`, `GeographicLaplacian`, `Gradient`, `LocalRelief`, `MultiscaleCurvatureValleyWidth`, `ScarpWavelet`
- **Isostasy**: `Deflection`, `GeographicDeflection`
- **Mixins**: `GDALMixin`, `GeographicGridMixin`, `MaxFlowLengthTrackingMixin`

### Core grid types

#### `BaseSpatialGrid`

*inherits* `GDALMixin`

A georeferenced raster, and the base class of every grid here.

Construct one by keyword; the combination you pass selects how it is
built (see :attr:`required_inputs_and_actions`)::

    BaseSpatialGrid(gdal_filename='grid.tif')
    BaseSpatialGrid(dx=30.0, grid=array)
    BaseSpatialGrid(nx=100, ny=80, dx=30.0)          # random values
    BaseSpatialGrid(nx=100, ny=80, projection=wkt, geo_transform=gt)

The data live in ``_griddata``, a 2-D NumPy array with row 0 at the
**north** edge and ``NaN`` for no-data; the georeferencing lives in
``_georef_info`` (a :class:`Georef_info`).

Indexing is by ``(row, col)`` and is forgiving: reading outside the grid
returns ``None`` and writing outside it is ignored, which is what the
neighbourhood traversals in this module rely on.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `nx`, `ny`, `dx`
- `dx`, `grid`

*Stored as* `float64`

| Member | Signature | Summary |
|---|---|---|
| `apply_moving_window` | `(moving_window)` | Apply a :class:`~TopoAnalysis.MovingWindow.MovingWindow` to the grid. |
| `average_over_distance` | `(distance, grid=None)` | Mean of the grid within ``distance`` of each cell. |
| `calculate_gradient_over_length_scale` | `(length_scale)` | Centred finite-difference gradient measured over ``length_scale``. |
| `calculate_laplacian_over_length_scale` | `(length_scale)` | Laplacian measured over ``length_scale`` with a centred stencil. |
| `clip_to_bounds` | `(bounds)` | Clip to ``((xmin, xmax), (ymin, ymax))``. |
| `clip_to_extent` | `(extent)` | Clip to ``(xmin, xmax, ymin, ymax)``, keeping every cell inside it. |
| `clip_to_mask_grid` | `(mask_grid)` | Reproject and resample this grid onto ``mask_grid``'s frame, in place. |
| `clip_to_shapefile` | `(shapefile)` | Cut the grid to a vector boundary, in place. |
| `coordinate_vectors` | `()` | ``(x, y)`` coordinates of the cell centres, west-to-east and |
| `export_arc_e00_grid` | `(filename)` | Write an ArcInfo E00 interchange grid. |
| `extent` | `()` | ``(xmin, xmax, ymin, ymax)`` of the grid's outer pixel edges. |
| `extent_of_data` | `()` | ``(xmin, xmax, ymin, ymax)`` of the smallest box holding all valid data. |
| `find_nearest_cell_with_greatest_value` | `(index, pixel_radius=5)` | Index of the highest-valued cell within ``pixel_radius``. |
| `find_nearest_cell_with_value` | `(index, value, pixel_radius)` | Index of the cell within ``pixel_radius`` whose value is closest to ``value``. |
| `find_nearest_cell_with_value_greater_than` | `(index, value, pixel_radius)` | Index of the cell within ``pixel_radius`` whose value is the smallest |
| `get_XY_matrices` | `()` | Meshgrid of cell-centre coordinates. |
| `get_XY_matricies` | `()` | Meshgrid of cell-centre coordinates. |
| `load` | `(filename)` | Read a grid previously written by :meth:`save`. *(classmethod)* |
| `location_in_grid` | `(xo)` | True when ``xo = (x, y)`` lands on a cell holding real data. |
| `mosaic` | `(tiles)` | Reassemble tiles produced by :meth:`tile` into a single grid. *(classmethod)* |
| `plot` | `(**kwargs)` | Display the grid with ``imshow``, georeferenced in map coordinates. |
| `principal_curvatures` | `()` | Maximum and minimum principal curvature of the surface. |
| `remove_padding` | `(xpadding=10, ypadding=10)` | Strip the overlap that :meth:`tile` added, ready for :meth:`mosaic`. |
| `resample` | `(de, interpolation='quintic')` | Resample onto a grid of spacing ``de``. |
| `save` | `(filename)` | Write the grid to an LZW-compressed GeoTIFF. |
| `set_value_at_rowscols` | `(value, rowscols)` | Set ``value`` at each ``(row, col)`` in ``rowscols``. |
| `snap_locations_to_closest_value` | `(v, value, pixel_radius=5)` | Move each ``(x, y)`` to the nearby cell whose value best matches. |
| `snap_locations_to_greatest_value` | `(v, pixel_radius=5)` | Move each ``(x, y)`` to the highest-valued cell within ``pixel_radius``. |
| `sort` | `(reverse=True, force=False, mask=None)` | Flat indices of the cells, ordered by value. |
| `tile` | `(tile_xdim=400, tile_ydim=400, tile_xpadding=10, tile_ypadding=10)` | Split the grid into overlapping tiles. |
| `vectorize` | `(filename)` | Polygonize the grid into an ESRI shapefile called ``filename.shp``. |
| `write_to_ai` | `(filename, nodata_value=-9999.0)` | Write the grid as an ArcInfo ASCII grid. |

#### `Georef_info`

Where a grid sits on the ground.

Attributes
----------
geoTransform : tuple
    The GDAL 6-tuple ``(x_origin, dx, 0, y_origin, 0, -dy)``, whose
    origin is the *outer corner* of the north-west cell.
projection : str or int
    WKT, or the integer 0 when no projection is set.
dx : float
    Cell size in map units.  Cells are assumed square.
xllcenter, yllcenter : float
    Centre of the south-west cell.
nx, ny : int
    Columns and rows.

#### `BaseSpatialShape`

A vector layer, used to cut rasters to a boundary.

::

    shape = BaseSpatialShape(shapefile_name='basin.shp')
    mask = shape.createMaskFromShape(dem._georef_info,
                                     dem._georef_info.projection,
                                     gdal.GDT_Byte)

| Member | Signature | Summary |
|---|---|---|
| `createMaskFromShape` | `(geoRefInfo, projection, dtype, noDataValue=0)` | Rasterise the shape onto ``geoRefInfo``'s frame, 1 inside and 0 out. |

#### `ValueGrid`

*inherits* `BaseSpatialGrid`

A plain raster of values, with a helper for scattered writes.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `dx`, `grid`

| Member | Signature | Summary |
|---|---|---|
| `set_value_at_indexes` | `(indexes, value)` | Set ``value`` at each ``(row, col)`` pair in ``indexes``. |

#### `Elevation`

*inherits* `CalculationMixin`, `BaseSpatialGrid`

A DEM.

Adds the terrain-derivative helpers of :class:`CalculationMixin` to the
base raster, plus the edge detection that depression filling starts from.

| Member | Signature | Summary |
|---|---|---|
| `findDEMedge` | `()` | Cells where water can leave the DEM. |
| `outlets_at_coastlines` | `(iterations=3)` | Map coordinates of the cells forming the land/sea boundary. |
| `track_flow_downhill` | `(starting_point, maximum_pit_depth=20)` | Walk downhill from ``starting_point``, hopping over shallow pits. |

#### `GeographicElevation`

*inherits* `GeographicGridMixin`, `Elevation`

A DEM on a latitude/longitude grid.

### Depression filling

#### `PriorityQueueMixIn`

Depression filling by Priority-Flood.

Mix into a grid class to give it :meth:`_flood`, which removes closed
depressions so that every cell has a downhill path to an outlet.

The algorithms are those of

    Barnes, R., Lehman, C., Mulla, D. (2014). "Priority-flood: An optimal
    depression-filling and watershed-labeling algorithm for digital
    elevation models." *Computers & Geosciences* 62, 117-127.

and run in :mod:`TopoAnalysis.kernels`, in C++ where the extension is
built and in NumPy otherwise.

Attributes
----------
aggradation_slope : float
    Gradient of the increment added while filling.  The default of 1e-12
    makes filled surfaces very slightly convergent, which gives D8 routing
    something to follow across what would otherwise be a dead flat
    plateau.  Set it to 0 for the flat fill that TopoToolbox's
    ``fillsinks`` produces.

| Member | Signature | Summary |
|---|---|---|
| `randomize_subbasins_with_mask` | `(*args, **kwargs)` | Replace elevations inside a mask with random values and re-flood. |

#### `FilledElevation`

*inherits* `PriorityQueueMixIn`, `Elevation`

A DEM with its closed depressions removed.

::

    filled = FilledElevation(elevation=dem)

Filling uses Priority-Flood (Barnes et al. 2014).  By default the filled
surface carries a gradient of
:attr:`PriorityQueueMixIn.aggradation_slope` across former depressions so
that D8 routing has a direction to follow; pass
``aggradation_slope=0`` -- or set the class attribute -- for the flat fill
that TopoToolbox's ``fillsinks`` produces.

Other keywords (``mask``, ``outlets``, ``maximum_pit_depth``,
``clip_to_fill``) are passed to :meth:`PriorityQueueMixIn._flood`.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`

#### `PriorityFillGrid`

*inherits* `PriorityQueueMixIn`, `BaseSpatialGrid`

Cells of a mask that are connected to a set of outlets.

Floods outwards from ``outlets`` through the cells where ``mask`` is 1
and returns a 0/1 grid of everything reached.  Used to grow a valley
network from a curvature mask, keeping only the parts that connect to a
real channel.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `outlets`, `mask`

*Stored as* `uint8`

### Flow routing

#### `FlowDirection`

*inherits* `BaseSpatialGrid`

Base class for flow-direction grids.

#### `FlowDirectionD8`

*inherits* `FlowDirection`

Single-direction (D8) flow directions.

Values are ArcGIS direction codes; 0 means the cell has no downstream
neighbour (a pit, an outlet at the grid edge, or no-data)::

    | 32  64 128 |
    | 16   X   1 |
    |  8   4   2 |

Build one from a filled DEM::

    d8 = FlowDirectionD8(flooded_dem=FilledElevation(elevation=dem))

Direction is chosen by steepest *gradient*, so the drop to a diagonal
neighbour is divided by ``sqrt(2) * dx`` before being compared with the
cardinal ones.  Ties go to the first direction in the order E, SE, S, SW,
W, NW, N, NE.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flooded_dem`
- `elevation`

*Stored as* `uint8`

| Member | Signature | Summary |
|---|---|---|
| `basin_mask` | `(outlets)` | Boolean array marking every cell draining to any of ``outlets``. |
| `bounds_of_basin_for_outlet` | `(outlet)` | ``((xmin, xmax), (ymin, ymax))`` of the basin draining to ``outlet``. |
| `convert_rivertools_directions_to_arc` | `(no_data=0)` | Convert RiverTools flow codes in place to the ArcGIS convention. |
| `cycle_count` | `` | How many cells lie on circular flow paths (0 for a valid network). |
| `divides` | `()` | Grid marking drainage divides: cells with no upstream neighbour. |
| `divides_for_outlets` | `(outlet1, outlet2)` | Cells on each side of the divide shared by two basins. |
| `get_flow_to_cell` | `(i, j)` | Index of the cell that ``(i, j)`` drains into. |
| `get_indexes_of_upstream_cells` | `(i, j)` | Tuple of ``(row, col)`` for every cell upstream of ``(i, j)``. |
| `get_indexes_of_upstream_cells_for_location` | `(x, y)` | As :meth:`get_indexes_of_upstream_cells`, from a map coordinate. |
| `get_upstream_cell_indexes` | `(i, j)` | Indices of the immediate neighbours that drain into ``(i, j)``. |
| `invalidate_network` | `()` | Drop the cached receivers and ordering. |
| `locations_of_paired_hollows` | `(outlet1, outlet2, area, Ao=100000.0)` | Channel heads reached by flowing down from a shared divide. |
| `map_values_to_recursive_list` | `(outlet, **kwargs)` | Nested dictionary of the basin upstream of ``outlet``. |
| `paired_divides` | `(mask=None)` | Pair each divide cell with the cell directly across the divide. |
| `pixel_scale` | `(dtype=<class 'numpy.float32'>)` | Step-length multiplier for each cell: 1 cardinal, sqrt(2) diagonal. |
| `receivers` | `` | Flat index of the cell each cell drains to, or -1 for none. |
| `search_down_flow_direction` | `(start, search_length=inf)` | The cells along the flow path from ``start``, without distances. |
| `search_down_flow_direction_from_rowscols_location` | `(start, return_rowscols=False, search_length=inf)` | Follow flow downstream from a ``(row, col)`` start. |
| `search_down_flow_direction_from_xy_location` | `(start, return_rowscols=False, search_length=inf)` | Follow flow downstream from an ``(x, y)`` start. |
| `search_down_flow_direction_with_length` | `(start, search_length=inf)` | Follow flow downstream from ``start``, recording distance travelled. |
| `set_value_at_rowscols` | `(value, rowscols)` | Set ``value`` at each ``(row, col)`` in ``rowscols``. |
| `step_lengths` | `(pixel_dimension=None)` | Distance from each cell to its receiver, in map units. |
| `topological_order` | `` | Flat cell indices ordered so every donor precedes its receiver. |
| `update_flow_codes_in_mask` | `(*args, **kwargs)` | Re-derive flow codes inside a mask after the DEM has changed. |

#### `FlowLength`

*inherits* `BaseSpatialGrid`

Longest upstream flow distance reaching each cell.

::

    length = FlowLength(flow_direction=d8)

Alongside the distances, the grid records which upstream neighbour
supplies each cell's longest path.  That "main stem" network is what
:class:`Relief`, :class:`ScaledRelief` and :class:`Ksi` follow, and what
:meth:`locations_along_flow_path_from_outlet` walks.

.. note::
   The main-stem codes use the same ArcGIS convention as
   :class:`FlowDirectionD8`, pointing *from* a cell *at* its upstream
   donor.  Before version 1.0 they used a row-flipped convention that
   was internally consistent but disagreed with every other
   flow-direction grid in the library, so ``*_directions`` side-files
   written by older versions must be regenerated.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`

| Member | Signature | Summary |
|---|---|---|
| `clip_to_extent` | `(extent)` | Clip both the lengths and the main-stem codes to ``extent``. |
| `indexes_along_flow_path_from_outlet` | `(outlet)` | ``(row, col)`` along the longest flow path upstream of ``outlet``. |
| `is_along_flow_length` | `(from_index, to_index)` | True when ``from_index`` supplies ``to_index``'s longest path. |
| `load` | `(filename)` | Read back a grid written by :meth:`save`. *(classmethod)* |
| `locations_along_flow_path_from_outlet` | `(outlet)` | As above, in map coordinates. |
| `locations_along_flow_path_from_outlets` | `(outlets)` | Concatenated main-stem paths for several outlets. |
| `main_stem_directions` | `` | Codes pointing from each cell at its longest-path donor. |
| `map_values_to_recursive_list` | `(outlet, **kwargs)` | Nested main-stem dictionary upstream of ``outlet``. |
| `points_with_length` | `(length, fd)` | Map coordinates where the flow path first exceeds ``length``. |
| `save` | `(filename)` | Write the lengths, plus the main-stem codes to ``filename_directions``. |

#### `GeographicFlowLength`

*inherits* `GeographicGridMixin`, `FlowLength`

Flow length on a latitude/longitude grid.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`

#### `Mask`

*inherits* `BaseSpatialGrid`

0/1 grid marking the cells that drain to a set of outlets.

::

    mask = Mask(flow_direction=d8, outlets=[(x, y), ...])

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `outlets`

*Stored as* `uint8`

| Member | Signature | Summary |
|---|---|---|
| `perform_erosion` | `(structure=None, iterations=1)` | Morphological erosion: shrink the mask by ``iterations`` cells. |
| `perform_opening` | `(structure=None, iterations=1)` | Morphological opening: erode then dilate, removing specks. |

#### `DiscreteFlowAccumulation`

*inherits* `BaseSpatialGrid`

Area drained by a single steepest-descent walk from each outlet.

Unlike :class:`Area`, which accumulates over the whole D8 network, this
follows one path per outlet, stepping to the lowest unvisited neighbour
each time, and records the running area along it.  It is used to trace
individual flow paths across surfaces (debris flows, lava) where the
network abstraction does not apply.

::

    dfa = DiscreteFlowAccumulation(elevation=dem, outlets=[(x, y)])

Keyword arguments
-----------------
mask : BaseSpatialGrid
    Confine the walk to cells where the mask is 1.
terminations_only : bool
    Record the total area only at the end of each path.
display_output : bool
    Report progress outlet by outlet.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `outlets`
- `elevation`

#### `GeographicDiscreteFlowAccumulation`

*inherits* `GeographicGridMixin`, `DiscreteFlowAccumulation`

Discrete flow accumulation on a latitude/longitude grid.

### Drainage area

#### `Area`

*inherits* `BaseSpatialGrid`

D8 contributing (drainage) area.

::

    area = Area(flow_direction=d8)

Values are in squared map units: each cell contributes its own area
(``dx**2``, or the true spherical cell area for a
:class:`GeographicArea`) plus everything upstream of it.

To count contributing *cells* instead -- which is what TopoToolbox's
``flowacc`` returns -- divide by the per-cell area, or pass
``weights=np.ones(...)``.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`
- `dx`, `grid`

| Member | Signature | Summary |
|---|---|---|
| `areas_between` | `(fd, min_area, max_area)` | Channel heads: cells in an area band whose donors are all below it. |
| `areas_greater_than` | `(min_area)` | Map coordinates of every cell whose area is at least ``min_area``. |

#### `GeographicArea`

*inherits* `GeographicGridMixin`, `Area`

Drainage area on a latitude/longitude grid, using true cell areas.

#### `LogArea`

*inherits* `BaseSpatialGrid`

Base-10 logarithm of drainage area.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `area`

*Stored as* `float64`

#### `ValleyArea`

*inherits* `_ValleyMaskMixin`, `Area`

Drainage area accumulated only over valley-floor cells.

::

    va = ValleyArea(flow_direction=d8, area=area, laplace=laplacian,
                    valley_laplace_value=0.01, min_area_value=1e6)

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `area`, `laplace`, `valley_laplace_value`, `min_area_value`

#### `GeographicValleyArea`

*inherits* `GeographicGridMixin`, `ValleyArea`

Valley-floor area on a latitude/longitude grid.

#### `MainstemValleyArea`

*inherits* `_ValleyMaskMixin`, `Area`

Valley-floor area accumulated along the main stem only.

Where :class:`ValleyArea` sums every tributary, this follows only the
longest flow path into each cell, so the result is the valley area of the
trunk stream rather than of the whole basin.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `area`, `laplace`, `valley_laplace_value`, `min_area_value`, `flow_direction`

#### `GeographicMainstemValleyArea`

*inherits* `GeographicGridMixin`, `MainstemValleyArea`

Main-stem valley area on a latitude/longitude grid.

### Chi and channel steepness

#### `Chi`

*inherits* `BaseSpatialGrid`

The chi coordinate, integrated upstream from chosen outlets.

.. math::
    \chi(x) = \int_{x_{\mathrm{outlet}}}^{x}
              \left(\frac{A_0}{A(x')}\right)^{\theta}\,\mathrm{d}x'

::

    chi = Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e6,
              outlets=[(x, y)])

Chi is 0 at each outlet and increases upstream.  Cells outside the
outlets' drainage areas are left as 0.

Keyword arguments
-----------------
trapezoid : bool
    Average the integrand across each step instead of taking it at the
    upstream cell.  This is the rule TopoToolbox's ``chitransform`` uses;
    the default (``False``) reproduces the original TopoAnalysis result.
maximum_length : float
    Stop integrating beyond this flow distance from the outlet.
mask : BaseSpatialGrid
    Restrict the integration to cells where the mask is non-zero.

.. note::
   Before version 1.0, chi at the outlet itself was one cell-step of the
   integrand rather than 0.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `area`, `flow_direction`, `theta`, `Ao`, `outlets`
- `area`, `flow_direction`, `flow_length`, `theta`, `Ao`, `basin_length`

#### `GeographicChi`

*inherits* `GeographicGridMixin`, `Chi`

Chi on a latitude/longitude grid.

#### `Ksi`

*inherits* `BaseSpatialGrid`, `MaxFlowLengthTrackingMixin`

The chi coordinate integrated along longest-flow-length paths.

.. math::
    \chi = \int \left(\frac{A_0}{A}\right)^{\theta} \mathrm{d}x

::

    ksi = Ksi(area=area, flow_direction=d8, flow_length=length,
              theta=0.45, Ao=1e6)

Cells with less than ``Ao`` of drainage area are left out of the
integral, so chi is defined only on the channel network.

.. note::
   Before version 1.0 the integrand was ``(Ao / (A - Ao))**theta``, which
   diverges as ``A`` approaches ``Ao`` and does not match the definition
   of chi in the literature. Values from older runs are not comparable
   with these.

See also :class:`Chi`, which integrates upstream from chosen outlets over
the whole network rather than along main stems.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `area`, `flow_direction`, `theta`, `Ao`, `flow_length`

#### `GeographicKsi`

*inherits* `GeographicGridMixin`, `Ksi`

Chi along longest-flow-length paths on a latitude/longitude grid.

#### `Relief`

*inherits* `BaseSpatialGrid`, `MaxFlowLengthTrackingMixin`

Height of the head of the longest flow path above each cell.

::

    relief = Relief(flow_direction=d8, elevation=dem, flow_length=length)

Pass ``area`` and ``Ao`` together to stop the relief signal being carried
down from cells with less than ``Ao`` of drainage area, which keeps
hillslope noise out of the channel network.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `elevation`, `flow_length`

#### `ScaledRelief`

*inherits* `Relief`

:class:`Relief` scaled by ``Ao**theta``.

That scaling puts relief in the same units as :class:`Ksi`, so the ratio
of the two is a channel steepness.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `flooded_dem`, `elevation`, `flow_length`, `Ao`, `theta`
- `flow_direction`, `elevation`, `flow_length`, `Ao`, `theta`

#### `ChiScaledRelief`

*inherits* `BaseSpatialGrid`

Height above the basin outlet, scaled by ``Ao**theta``.

Plotted against :class:`Chi`, the slope of this quantity is the channel
steepness index.

::

    relief = ChiScaledRelief(elevation=dem, flow_direction=d8,
                             theta=0.45, Ao=1e6, outlets=[(x, y)])

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `flow_direction`, `theta`, `Ao`, `outlets`
- `elevation`, `flow_direction`, `flow_length`, `theta`, `Ao`, `basin_length`

#### `ChannelSlope`

*inherits* `BaseSpatialGrid`

Downstream gradient of each cell along its D8 flow direction.

::

    slope = ChannelSlope(flow_direction=d8, elevation=dem)

Positive values are downhill.  Cells with no downstream neighbour stay 0.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `elevation`

#### `CrossDivideDChi`

*inherits* `_CrossDivideMixin`, `BaseSpatialGrid`

Difference in chi across each drainage divide.

A large value means the two sides of the divide are far from
steady state, so the divide is expected to migrate towards the
higher-chi side (Willett et al. 2014).

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `chi`, `flow_direction`

#### `NormalizedCrossDivideDChi`

*inherits* `_CrossDivideMixin`, `BaseSpatialGrid`

Chi contrast across each divide, divided by the mean chi there.

Normalising makes divides in high- and low-chi parts of a landscape
comparable.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `chi`, `flow_direction`

#### `RestoredElevation`

*inherits* `BaseSpatialGrid`

Topography implied by a uniform channel steepness.

Starting from the given outlets, elevations are rebuilt upstream from

.. math::  \frac{\mathrm{d}z}{\mathrm{d}x} = k_s A^{-\theta}

and the divides are then allowed to migrate towards whichever side ends
up lower.  Repeating that (``iterations`` times) converges on the
landscape a steady, spatially uniform ``ks`` would produce, which can be
compared against the real DEM to find where it is out of steady state.

::

    restored = RestoredElevation(flow_direction=d8, elevation=dem,
                                 area=area, theta=0.45, ks=50,
                                 outlets=[(x, y)], iterations=5)

Keyword arguments
-----------------
randomize : bool
    Randomise elevations inside the basin first, then re-route, so the
    result does not inherit the real network's planform.
fix_external_outlets : bool
    Hold the outer boundary of the basin fixed while divides migrate.
verbose : bool
    Report convergence each iteration (default ``True``).

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `flow_direction`, `elevation`, `area`, `theta`, `ks`, `outlets`, `iterations`

#### `GeographicRestoredElevation`

*inherits* `GeographicGridMixin`, `RestoredElevation`

Restored topography on a latitude/longitude grid.

### Along-flow smoothing

#### `AlongFlowSmoothing`

Collect the cells within a window along each flow path.

Mixed into the ``*WithSmoothing`` classes, which fit a model over that
window at every cell.  The window is set either by drop
(``vertical_interval``) or by distance (``horizontal_interval``), and
extends both upstream and downstream of the cell being evaluated.

#### `KsFromChiWithSmoothing`

*inherits* `BaseSpatialGrid`, `AlongFlowSmoothing`

Channel steepness fitted to chi--elevation in a moving window.

::

    ks = KsFromChiWithSmoothing(elevation=dem, area=area,
                                flow_direction=d8, theta=0.45,
                                vertical_interval=50.0)

At every channel cell, chi and elevation are extracted over a window
along the flow path and a line through the origin is fitted; its slope
is k_s.  ``save()`` writes seven bands: ks, the number of profiles
crossing each cell, mean squared error, sum of squares, r-squared,
p-value, and the number of points in the regression.

Needs ``statsmodels``.  Pass ``verbose=False`` to silence the progress
ticker.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `flow_direction`, `theta`, `vertical_interval`
- `elevation`, `area`, `flow_direction`, `theta`, `horizontal_interval`

| Member | Signature | Summary |
|---|---|---|
| `load` | `(filename)` | Read back the multi-band grid written by :meth:`save`. *(classmethod)* |
| `save` | `(filename)` | Write the grid to an LZW-compressed GeoTIFF. |

#### `GeographicKsFromChiWithSmoothing`

*inherits* `GeographicGridMixin`, `KsFromChiWithSmoothing`

Windowed steepness fit on a latitude/longitude grid.

#### `ThetaFromChiWithSmoothing`

*inherits* `BaseSpatialGrid`, `AlongFlowSmoothing`

Concavity fitted to chi--elevation in a moving window.

As :class:`KsFromChiWithSmoothing`, but ``theta`` is optimised at each
cell rather than held fixed: the value returned is the one that makes
the chi--elevation relation most nearly linear.  The search is bounded
at +/-10.

Needs ``statsmodels``.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `flow_direction`, `vertical_interval`, `min_area`
- `elevation`, `area`, `flow_direction`, `horizontal_interval`, `min_area`

| Member | Signature | Summary |
|---|---|---|
| `load` | `(filename)` | Read back the multi-band grid written by :meth:`save`. *(classmethod)* |
| `save` | `(filename)` | Write the grid to an LZW-compressed GeoTIFF. |

#### `GeographicThetaFromChiWithSmoothing`

*inherits* `GeographicGridMixin`, `ThetaFromChiWithSmoothing`

Windowed concavity fit on a latitude/longitude grid.

#### `ChannelSlopeWithSmoothing`

*inherits* `BaseSpatialGrid`, `AlongFlowSmoothing`

Channel gradient measured over a window along the flow path.

Far less noisy than a cell-to-cell slope, which on a real DEM is mostly
vertical quantisation.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `flow_direction`, `vertical_interval`
- `elevation`, `area`, `flow_direction`, `horizontal_interval`

| Member | Signature | Summary |
|---|---|---|
| `calc_channel_slope` | `(i, j, elevation, de, find_points_along_path)` | Gradient across the window centred on ``(i, j)``. |
| `points_along_path` | `(points, i, j)` | Which part of the window to measure over; the whole of it here. |

#### `ChannelDownSlopeWithSmoothing`

*inherits* `ChannelSlopeWithSmoothing`

Channel gradient over the downstream half of the window.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `flow_direction`, `vertical_interval`
- `elevation`, `area`, `flow_direction`, `horizontal_interval`

| Member | Signature | Summary |
|---|---|---|
| `points_along_path` | `(points, i, j)` | Only the downstream half of the window. |

#### `ChannelUpSlopeWithSmoothing`

*inherits* `ChannelSlopeWithSmoothing`

Channel gradient over the upstream half of the window.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `flow_direction`, `vertical_interval`
- `elevation`, `area`, `flow_direction`, `horizontal_interval`

| Member | Signature | Summary |
|---|---|---|
| `points_along_path` | `(points, i, j)` | Only the upstream half of the window. |

### Terrain attributes

#### `CalculationMixin`

Finite-difference terrain derivatives.

Supplies the slope, curvature and boundary-padding helpers that
:class:`Elevation`, :class:`Hillshade` and :class:`MaxSlope` share.

| Member | Signature | Summary |
|---|---|---|
| `assignBCs` | `(grid, nx, ny)` | Pad a grid by one cell, replicating the edge values. |
| `calcAverageSlopeOfGridSubset` | `(gridSubset, dx)` | Least-squares plane through a block of cells; returns its slopes. |
| `calcContourCurvature` | `(grid, dx)` | Contour (plan) curvature: the curvature of the contour lines. |
| `calcFiniteCurv` | `(grid, dx)` | Laplacian from a 5-point centred stencil, same shape as the input. |

#### `Hillshade`

*inherits* `CalculationMixin`, `BaseSpatialGrid`

Shaded relief, following the ESRI hillshade definition.

::

    hs = Hillshade(elevation=dem, azimuth=315, inclination=45)

``azimuth`` is the compass bearing of the light source in degrees
(0 = north, 90 = east) and ``inclination`` its height above the horizon.
Values run 0-255.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `azimuth`, `inclination`

*Stored as* `uint8`

| Member | Signature | Summary |
|---|---|---|
| `calcHillshade` | `(az, elev, z_factor=1.0)` | Fill the grid with an ESRI-style hillshade. |

#### `GeographicHillshade`

*inherits* `GeographicGridMixin`, `Hillshade`

Hillshade on a latitude/longitude grid.

#### `MaxSlope`

*inherits* `CalculationMixin`, `BaseSpatialGrid`

Magnitude of the topographic gradient, from centred differences.

::

    slope = MaxSlope(elevation=dem)

Values are dimensionless (rise over run); take ``arctan`` for degrees.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`

*Stored as* `float64`

| Member | Signature | Summary |
|---|---|---|
| `calcSlope` | `()` | Replace the grid with the gradient magnitude of its own values. |

#### `GeographicMaxSlope`

*inherits* `GeographicGridMixin`, `MaxSlope`

Gradient magnitude on a latitude/longitude grid.

#### `Laplacian`

*inherits* `CalculationMixin`, `BaseSpatialGrid`

Laplacian of elevation -- positive in valleys, negative on ridges.

::

    curvature = Laplacian(elevation=dem)

Uses SciPy's 5-point stencil, scaled by ``1 / dx**2``, so the units are
inverse length.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`

*Stored as* `float64`

#### `GeographicLaplacian`

*inherits* `GeographicGridMixin`, `Laplacian`

Laplacian of elevation on a latitude/longitude grid.

#### `Gradient`

*inherits* `BaseSpatialGrid`

The two components of the topographic gradient.

::

    gradient = Gradient(elevation=dem)
    gradient._gx, gradient._gy

Stored as two bands rather than one, so direction survives ``save()``.
:meth:`average_gradient` smooths both components over a length scale,
which is how a regional aspect is obtained.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`

| Member | Signature | Summary |
|---|---|---|
| `average_gradient` | `(distance)` | Both gradient components averaged over a disc of radius ``distance``. |
| `load` | `(filename)` | Read back the multi-band grid written by :meth:`save`. *(classmethod)* |
| `plot` | `(**kwargs)` | Display the grid with ``imshow``, georeferenced in map coordinates. |
| `save` | `(filename)` | Write the grid to an LZW-compressed GeoTIFF. |

#### `LocalRelief`

*inherits* `BaseSpatialGrid`

Elevation range within a circular window.

::

    relief = LocalRelief(elevation=dem, pixel_radius=25)

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `pixel_radius`

#### `MultiscaleCurvatureValleyWidth`

*inherits* `BaseSpatialGrid`

Valley width from the scale at which curvature is most negative.

A local quadratic is fitted at a range of window sizes; the width
recorded at each cell is the window that minimises the smaller principal
curvature, i.e. the scale at which the valley is best resolved.  The
curvature itself is kept in ``_minC``.

::

    width = MultiscaleCurvatureValleyWidth(
        elevation=dem, area=area, area_cutoff=1e6,
        min_width=30.0, max_width=600.0)

Pass ``use_dask=True`` to evaluate the scales in parallel.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `area`, `area_cutoff`, `max_width`, `min_width`

| Member | Signature | Summary |
|---|---|---|
| `load` | `(filename)` | Read back the multi-band grid written by :meth:`save`. *(classmethod)* |
| `mosaic` | `(tiles)` | Reassemble tiles, carrying the curvature band along with the widths. *(classmethod)* |
| `remove_padding` | `(xpadding=10, ypadding=10)` | Strip tile overlap from both the widths and the curvatures. |
| `save` | `(filename)` | Write the grid to an LZW-compressed GeoTIFF. |

#### `ScarpWavelet`

*inherits* `BaseSpatialGrid`

Template-matched scarp amplitude, age, orientation and signal-to-noise.

Holds the four-band output of a scarp template search.  Load one with
:meth:`load`; the bands become ``_A``, ``_kt``, ``_orientation`` and
``_SNR``.  :meth:`template_window` needs the separate ``scarplet``
package.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`

| Member | Signature | Summary |
|---|---|---|
| `calculate_valid_orientations` | `(window_size, age, num=46, discard_max_min=True)` | Range of template orientations that stay inside the valid data. |
| `load` | `(filename)` | Read back the multi-band grid written by :meth:`save`. *(classmethod)* |
| `load_elevation` | `(filename)` | Attach the DEM the templates were matched against. |
| `plot_orientations` | `(*args, **kwargs)` | Draw scarp orientations over a hillshade, faded by signal-to-noise. |
| `template_window` | `(window_size, age, orientation, use_pixels=False)` | The scarp template mask at one size, age and orientation. |
| `valid_data` | `()` | 1 where the attached DEM has usable, above-sea-level data; NaN elsewhere. |

### Isostasy

#### `Deflection`

*inherits* `BaseSpatialGrid`

Flexural deflection of an elastic plate under a topographic load.

Solves the thin-plate equation in the Fourier domain:

.. math::
    w(k) = \frac{-\rho_c\,g\,h(k)}
                {(\rho_m - \rho_c)\,g + D\,(2\pi k)^4}

::

    w = Deflection(elevation=dem, D=1e23, rho_m=3300, rho_c=2700, g=9.81)

Pass ``restored_elevation`` as well to get the *change* in deflection
between two surfaces, which is the isostatic response to the erosion
between them.

Because the solution is spectral, the load is implicitly periodic; pad
the DEM if edge effects matter.

*Constructor keywords* (any one combination)

- `nx`, `ny`, `projection`, `geo_transform`
- `ai_ascii_filename`, `EPSGprojectionCode`
- `gdal_filename`
- `elevation`, `D`, `rho_m`, `rho_c`, `g`

#### `GeographicDeflection`

*inherits* `GeographicGridMixin`, `Deflection`

Flexural deflection on a latitude/longitude grid.

### Mixins

#### `GDALMixin`

Reading and writing rasters through GDAL.

Mixed into :class:`BaseSpatialGrid`; nothing here is meant to be called
directly.  Every method reaches GDAL through the module-level lazy
proxies, so the mixin can be present without GDAL being installed.

| Member | Signature | Summary |
|---|---|---|
| `getDEMcoords` | `(GdalData, dx)` | Cell-centre x and y coordinate lists for an open GDAL dataset. |

#### `GeographicGridMixin`

Spherical geometry for grids in degrees of latitude and longitude.

Mix it in *before* the grid class to override the cell area and cell
dimension with their true, latitude-dependent values::

    class GeographicArea(GeographicGridMixin, Area):
        pass

Every ``Geographic*`` class in this module is exactly that.

#### `MaxFlowLengthTrackingMixin`

Propagate a value downstream along longest-flow-length paths only.

Subclasses set :attr:`_propagate_mode` to ``'carry'`` (the receiver takes
the donor's value) or ``'sum'`` (the receiver adds it) and may define
:meth:`_propagation_gate` to suppress the transfer for some cells.

### Module-level functions in `TopoAnalysis.dem`

| Function | Summary |
|---|---|
| `gdal_is_available()` | True when the GDAL Python bindings can be imported. |
| `mosaicFolder(folderPath, fileSuffix, outfile)` | Merge every ``*fileSuffix`` raster in a folder into ``outfile``. |
| `plot(*args, **kwargs)` | Scatter one grid's values against another's, cell by cell. |

## `TopoAnalysis.topotoolbox` -- TopoToolbox-equivalent functions

TopoToolbox-equivalent terrain analysis.

Every function here reproduces the numerical result of the identically named
function in TopoToolbox 2 (Schwanghart & Scherler 2014), the MATLAB package,
so that a result computed in TopoAnalysis can be checked against one computed
in MATLAB.

The rest of TopoAnalysis keeps its own idioms -- drainage area in square
metres, ESRI-style hillshade, chi integrated with a left-endpoint rule.  This
module exists for the cases where you need the *same number*, not just the
same quantity.  :doc:`docs/topotoolbox_parity` tabulates the correspondence
and lists the places where exact agreement is not achievable.

See [topotoolbox_parity.md](topotoolbox_parity.md) for the full
correspondence table and the places where exact agreement is not
achievable.

| Function | Summary |
|---|---|
| `arcslope(elevation, unit='tangent')` | ArcGIS/Horn 3x3 slope magnitude. |
| `aspect(elevation, classify=False)` | Down-slope direction in degrees clockwise from north. |
| `cellarea(grid, unit='m')` | Area of each cell of a geographic (lat/lon) grid. |
| `chitransform(flow_direction, area, mn=0.45, a0=1000000.0, outlets=None, correctcellsize=True, stream_mask=None)` | The chi integral transform. |
| `coord2sub(grid, x, y)` | Map coordinates to 0-based ``(row, col)`` subscripts. |
| `curvature(elevation, ctype='profc')` | Surface curvature from the Evans 3x3 quadratic fit. |
| `drainagebasins(flow_direction, outlets=None, valid=None)` | Label drainage basins. |
| `fillsinks(elevation, maxdepth=None)` | Fill topographic depressions to a flat surface. |
| `flowacc(flow_direction, weights=None)` | Flow accumulation in **number of cells**. |
| `flowdistance(flow_direction, direction='upstream')` | Along-flow distance in map units. |
| `flowobj(elevation, preprocess='fill')` | Build a D8 flow-direction grid the way ``FLOWobj`` does. |
| `getcoordinates(grid)` | ``(x, y)`` cell-centre coordinate vectors. |
| `getextent(grid)` | ``(xmin, xmax, ymin, ymax)`` of the outermost **cell centres**. |
| `gradient8(elevation, unit='tangent')` | Steepest downward gradient over the 8 neighbours. |
| `hillshade(elevation, azimuth=315.0, altitude=60.0, exaggerate=1.0)` | Shaded relief as the cosine of the illumination incidence angle. |
| `identifyflats(elevation, output='both')` | Locate flat areas, their sills and the closed basins. |
| `imposemin(flow_direction, elevation, sl=0.0)` | Carve elevations so every flow path descends at gradient >= ``sl``. |
| `ksn(flow_direction, elevation, area, theta=0.45, correctcellsize=True, min_gradient=1e-05)` | Normalized channel steepness index, :math:`k_{sn} = S A^{\theta}`. |
| `localtopography(elevation, radius=5000.0, type='range')` | Statistic of elevation within a circular window. |
| `mchi(flow_direction, elevation, chi)` | Gradient of elevation in chi space, :math:`M_{\chi}`. |
| `streamorder(flow_direction, stream_mask, order_type='strahler')` | Strahler or Shreve stream order. |
| `sub2coord(grid, row, col)` | 0-based ``(row, col)`` subscripts to cell-centre map coordinates. |

### `fillsinks(elevation, maxdepth=None)`

Fill topographic depressions to a flat surface.

Equivalent to ``GRIDobj/fillsinks(DEM)`` and, with ``maxdepth``, to
``fillsinks(DEM, maxdepth)``.

Unlike :class:`TopoAnalysis.dem.FilledElevation`, whose default adds a
tiny gradient across filled areas so that D8 routing has a direction to
follow, this fills to an exactly flat plateau -- which is what
TopoToolbox does, and what ``identifyflats`` then relies on.

Parameters
----------
elevation : BaseSpatialGrid
    The DEM.  ``NaN`` is no-data and, as in TopoToolbox, acts as an
    outlet: water reaching a NaN cell leaves the grid.
maxdepth : float, optional
    Leave depressions deeper than this unfilled.

Returns
-------
Elevation
    A new grid; the input is not modified.

### `identifyflats(elevation, output='both')`

Locate flat areas, their sills and the closed basins.

Equivalent to ``GRIDobj/identifyflats``.

* A **flat** is a cell with no strictly lower neighbour -- the MATLAB test
  is ``imerode(dem, ones(3)) == dem``.  An isolated pit therefore counts
  as a flat.  The grid's outer ring is cleared, as is the one-cell halo
  around any no-data.
* A **sill** is a non-flat cell at the same elevation as an adjacent
  flat: the point a flat drains through.
* A **closed basin** is a regional minimum that does not touch the grid
  edge.

Parameters
----------
elevation : BaseSpatialGrid
    Normally the output of :func:`fillsinks`.
output : {'both', 'flats', 'sills', 'closedbasins', 'all'}
    ``'both'`` returns ``(flats, sills)``; ``'all'`` adds the closed
    basins.

Returns
-------
BaseSpatialGrid or tuple
    Boolean grids.

### `imposemin(flow_direction, elevation, sl=0.0)`

Carve elevations so every flow path descends at gradient >= ``sl``.

Equivalent to ``FLOWobj/imposemin(FD, DEM, sl)``.  Cells are only ever
lowered, so channels are cut through obstructions instead of being
drowned by them.

Returns
-------
Elevation
    A new grid.

### `flowobj(elevation, preprocess='fill')`

Build a D8 flow-direction grid the way ``FLOWobj`` does.

Parameters
----------
elevation : BaseSpatialGrid
preprocess : {'fill', 'none'}
    ``'fill'`` removes depressions first (TopoToolbox's
    ``'preprocess','fill'``).  ``'carve'`` is not implemented; see the
    note below.

Returns
-------
FlowDirectionD8

Notes
-----
TopoToolbox routes flow across the flat areas left by filling using an
auxiliary "geodesic" topography, so that flat cells drain towards the
nearest sill.  TopoAnalysis instead adds a small gradient while filling
(:attr:`~TopoAnalysis.dem.PriorityQueueMixIn.aggradation_slope`).  Both
produce a valid, connected drainage network with identical basin
boundaries wherever the terrain is not flat, but **inside flats the two
networks differ**.  TopoToolbox's ``'carve'`` preprocessing lowers the
DEM instead of raising it and is a different surface again.

### `flowacc(flow_direction, weights=None)`

Flow accumulation in **number of cells**.

Equivalent to ``FLOWobj/flowacc(FD)``: every cell counts itself, so the
minimum value is 1.

Parameters
----------
flow_direction : FlowDirectionD8
weights : ndarray or BaseSpatialGrid, optional
    Per-cell weight; defaults to 1.

Returns
-------
BaseSpatialGrid

See Also
--------
TopoAnalysis.dem.Area : the same accumulation weighted by cell area.

### `flowdistance(flow_direction, direction='upstream')`

Along-flow distance in map units.

Equivalent to ``FLOWobj/flowdistance(FD, direction)`` with no seed
points.

Parameters
----------
direction : {'upstream', 'downstream', 'maxdownstream'}
    ``'upstream'`` gives the distance from each cell down to its outlet.
    ``'downstream'`` (a synonym for ``'maxdownstream'``) gives the length
    of the longest flow path reaching each cell -- the quantity
    :class:`TopoAnalysis.dem.FlowLength` computes.

### `drainagebasins(flow_direction, outlets=None, valid=None)`

Label drainage basins.

Equivalent to ``FLOWobj/drainagebasins``.  Without ``outlets`` every
terminal cell that something drains into seeds its own basin; with them,
each cell is labelled by the 1-based position of the outlet it drains to.

Parameters
----------
outlets : sequence of (x, y), optional
valid : ndarray or BaseSpatialGrid, optional
    Restrict labelling to these cells (e.g. ``~np.isnan(dem._griddata)``).

Returns
-------
BaseSpatialGrid
    Integer labels; 0 means "reaches none of the outlets".

Notes
-----
TopoToolbox numbers basins in the order its edge list happens to meet
them, so the *partition* matches but the integer labels generally do
not.  Compare partitions, not label values.

### `gradient8(elevation, unit='tangent')`

Steepest downward gradient over the 8 neighbours.

Equivalent to ``GRIDobj/gradient8``.  The cardinal and diagonal drops are
divided by ``cellsize`` and ``sqrt(2)*cellsize`` respectively, and the
result is clamped at 0, so pits and flats give exactly 0 rather than a
negative gradient.

Parameters
----------
unit : {'tangent', 'degree', 'radian', 'sine', 'percent'}

### `arcslope(elevation, unit='tangent')`

ArcGIS/Horn 3x3 slope magnitude.

Equivalent to ``GRIDobj/arcslope``.  Uses the 1-2-1 weighted 3x3 stencil
with an ``8 * cellsize`` denominator -- there is no ``sqrt(2)`` anywhere;
the diagonals enter through the weights.

``NaN`` holes are filled from their nearest valid neighbour before
differencing, then restored, which is how TopoToolbox handles them here
(and is the opposite of what :func:`gradient8` does).

### `hillshade(elevation, azimuth=315.0, altitude=60.0, exaggerate=1.0)`

Shaded relief as the cosine of the illumination incidence angle.

Equivalent to ``GRIDobj/hillshade`` with its default ``'surfnorm'``
method, including its defaults of ``azimuth=315`` and ``altitude=60``.

Returns
-------
BaseSpatialGrid
    Values in ``[-1, 1]``.  TopoToolbox only scales to 0-255 when it is
    plotting, so neither does this.

See Also
--------
TopoAnalysis.dem.Hillshade : the ESRI formulation, returning 0-255.

### `aspect(elevation, classify=False)`

Down-slope direction in degrees clockwise from north.

Equivalent to ``GRIDobj/aspect``.  0 = north, 90 = east, 180 = south,
270 = west.  Perfectly flat cells come out as 90, exactly as they do in
TopoToolbox.

Parameters
----------
classify : bool
    Return the 8-class Gomez-Plaza wetness ranking (1 = N, 3 = NE,
    5 = E, 7 = SE, 8 = S, 6 = SW, 4 = W, 2 = NW) instead of degrees.

### `curvature(elevation, ctype='profc')`

Surface curvature from the Evans 3x3 quadratic fit.

Equivalent to ``GRIDobj/curvature``.

Parameters
----------
ctype : {'profc', 'planc', 'tangc', 'meanc', 'total'}
    Profile, plan, tangential, mean, or total curvature.  ``'profc'``
    (the default) is positive where a slope steepens downhill.

Returns
-------
BaseSpatialGrid
    Units are 1/length, except ``'total'`` which is 1/length^2.

Notes
-----
Cells where the gradient vanishes give 0/0; TopoToolbox sets those to 0,
and so does this.  Cells adjacent to no-data therefore come out as 0
rather than ``NaN`` -- again matching TopoToolbox.

### `cellarea(grid, unit='m')`

Area of each cell of a geographic (lat/lon) grid.

Equivalent to ``GRIDobj/cellarea``.  Cells are treated as spherical
quadrangles on a sphere whose total area is
:data:`EARTH_SURFACE_AREA_KM2`.

Parameters
----------
unit : {'m', 'km'}
    Square metres (the default) or square kilometres.

Returns
-------
BaseSpatialGrid
    Constant along each row, since area depends only on latitude.

### `localtopography(elevation, radius=5000.0, type='range')`

Statistic of elevation within a circular window.

Equivalent to ``GRIDobj/localtopography``.  The window radius is given in
map units and rounded up to a whole number of cells.

Parameters
----------
type : {'range', 'max', 'min', 'mean', 'median', 'std'}

Notes
-----
TopoToolbox's default disk (``strel('disk', r)`` with ``N=8``) is a
periodic-line approximation whose exact shape is an internal MATLAB
detail.  This uses an exact disc -- TopoToolbox's ``N=0`` -- so values
can differ slightly near the window scale.

### `streamorder(flow_direction, stream_mask, order_type='strahler')`

Strahler or Shreve stream order.

Equivalent to ``FLOWobj/streamorder`` / ``STREAMobj/streamorder``.

Parameters
----------
stream_mask : ndarray or BaseSpatialGrid
    Non-zero on the channel network, e.g. ``flowacc(fd)._griddata > 1000``.
order_type : {'strahler', 'shreve'}

Returns
-------
BaseSpatialGrid
    Integer orders; 0 off the network.

### `chitransform(flow_direction, area, mn=0.45, a0=1000000.0, outlets=None, correctcellsize=True, stream_mask=None)`

The chi integral transform.

Equivalent to ``STREAMobj/chitransform``:

.. math::
    \chi = \int_{\mathrm{outlet}}^{x}
           \left(\frac{A_0}{A(x')}\right)^{m/n}\,\mathrm{d}x'

integrated upstream with the trapezoidal rule, so chi is 0 at each outlet
and increases upstream.

Parameters
----------
area : BaseSpatialGrid
    Flow accumulation.  With ``correctcellsize=True`` (the default, as in
    TopoToolbox) it is taken to be in **number of cells** and multiplied
    by ``cellsize**2``; pass ``False`` if it is already in m^2, e.g. if
    it came from :class:`TopoAnalysis.dem.Area`.
mn : float
    The m/n concavity ratio.  TopoToolbox's default is 0.45.
a0 : float
    Reference area in m^2.  TopoToolbox's default is 1e6 (1 km^2).
outlets : sequence of (x, y), optional
    Integration starts here.  Defaults to every terminal cell of the
    network, which is what a STREAMobj built over the whole grid gives.
stream_mask : ndarray or BaseSpatialGrid, optional
    Restrict the integration to the channel network.

Returns
-------
BaseSpatialGrid
    Chi in map units; ``NaN`` away from the integrated network.

### `ksn(flow_direction, elevation, area, theta=0.45, correctcellsize=True, min_gradient=1e-05)`

Normalized channel steepness index, :math:`k_{sn} = S A^{\theta}`.

Equivalent to ``STREAMobj/ksn``, including its unconditional minima
imposition at a gradient of 1e-5 before the slope is measured.

Parameters
----------
area : BaseSpatialGrid
    Flow accumulation in cells (see ``correctcellsize``).
theta : float
    Reference concavity; TopoToolbox's default is 0.45.
min_gradient : float
    Gradient imposed by the carving step.  TopoToolbox hard-codes 1e-5.

Returns
-------
BaseSpatialGrid
    ``ksn`` with a reference area of 1 m^2 -- TopoToolbox's ``ksn`` has
    no ``a0``.  Multiply by ``a0**theta`` to compare against
    ``chiplot``'s ``ks``.

### `mchi(flow_direction, elevation, chi)`

Gradient of elevation in chi space, :math:`M_{\chi}`.

Equivalent to ``STREAMobj/mchi``: a forward difference assigned to the
upstream cell, with no carving and no smoothing.

``M_chi`` scales with the ``a0`` baked into ``chi``: computing chi with
``a0=1`` puts it in the same units as :func:`ksn`, and
``mchi(a0) * a0**mn`` is invariant.

It is *close to* but not equal to :func:`ksn`, despite what TopoToolbox's
own docstring says: ``ksn`` carves the DEM first and evaluates the
integrand at the upstream cell, while ``mchi`` uses the raw elevations
and inherits chi's trapezoidal averaging.

### `getcoordinates(grid)`

``(x, y)`` cell-centre coordinate vectors.

Equivalent to ``GRIDobj/getcoordinates``.  ``x`` runs west to east;
``y`` runs **north to south**, matching the row order of the grid.

### `getextent(grid)`

``(xmin, xmax, ymin, ymax)`` of the outermost **cell centres**.

Equivalent to ``GRIDobj/getextent``.  Note this is half a cell inside
:meth:`TopoAnalysis.dem.BaseSpatialGrid.extent`, which returns the pixel
edges that matplotlib wants.

### `coord2sub(grid, x, y)`

Map coordinates to 0-based ``(row, col)`` subscripts.

Equivalent to ``GRIDobj/coord2sub`` apart from the index base.  Points
outside the grid give ``NaN``.

Notes
-----
Halves are rounded away from zero, as MATLAB does.  Python's built-in
``round`` -- which
:meth:`TopoAnalysis.dem.BaseSpatialGrid._xy_to_rowscols` uses -- rounds
halves to even, so the two disagree for points landing exactly on a cell
boundary.

### `sub2coord(grid, row, col)`

0-based ``(row, col)`` subscripts to cell-centre map coordinates.

Equivalent to ``GRIDobj/sub2coord``.  Subscripts outside the grid are
extrapolated rather than rejected, as in TopoToolbox.

## `TopoAnalysis.kernels` -- backend dispatch

Backend dispatch for the TopoAnalysis grid kernels.

Import the kernels from here rather than from :mod:`TopoAnalysis.fastops` or
the compiled extension directly::

    from TopoAnalysis import kernels
    kernels.priority_flood(z, mode="flat")

If the C++ extension built successfully it is used; otherwise the pure-NumPy
implementations in :mod:`TopoAnalysis.fastops` take over.  Both produce the
same numbers -- the test suite asserts it -- so the choice is purely about
speed.

Set the environment variable ``TOPOANALYSIS_PURE_PYTHON=1`` to force the
NumPy path, which is useful when debugging or comparing the two.

| Function | Summary |
|---|---|
| `accumulate(recv, order, weights, gate=None)` | Downstream accumulation; see :func:`TopoAnalysis.fastops.accumulate`. |
| `backend() -> 'str'` | ``'c++'`` or ``'python'``. |
| `chi(recv, order, area, step_length, outlets, A0, theta, trapezoid=False, max_length=0.0, mask=None)` | ``(chi, distance)``; see :func:`TopoAnalysis.fastops.chi`. |
| `default_seeds(elevations)` | Flat indices of the cells a flood should start from. |
| `downstream_distance(recv, order, step_length)` | See :func:`TopoAnalysis.fastops.downstream_distance`. |
| `drainage_basins(recv, order, outlets=None, valid=None)` | See :func:`TopoAnalysis.fastops.drainage_basins`. |
| `flow_directions(elevations, cellsize_grid=None, cellsize=1.0)` | Steepest-descent D8 codes; see :func:`TopoAnalysis.fastops.flow_directions`. |
| `flow_length(recv, order, step_length)` | ``(length, from_codes)``; see :func:`TopoAnalysis.fastops.flow_length`. |
| `have_extension() -> 'bool'` | True when the compiled kernels are in use. |
| `imposemin(recv, order, step_length, sl, elevations)` | Carve in place; see :func:`TopoAnalysis.fastops.imposemin`. |
| `priority_flood(elevations, closed=None, seeds=None, mode='flat', epsilon=0.0, cellsize=1.0, max_pit_depth=0.0, track_visited=False, reference_algorithm=False)` | Fill depressions in place; see :func:`TopoAnalysis.fastops.priority_flood`. |
| `propagate_along_main_stem(recv, order, main_stem_from_codes, gate, values, mode)` | See :func:`TopoAnalysis.fastops.propagate_along_main_stem`. |
| `receivers(codes)` | Receiver index per cell; see :func:`TopoAnalysis.fastops.receivers`. |
| `stream_order(recv, order, is_stream, kind='strahler')` | See :func:`TopoAnalysis.fastops.stream_order`. |
| `topological_order(recv)` | ``(order, n_cells_on_cycles)``; see :func:`TopoAnalysis.fastops.topological_order`. |
| `upstream_mask(recv, outlets)` | See :func:`TopoAnalysis.fastops.upstream_mask`. |

## `TopoAnalysis.fastops` -- the reference NumPy kernels

Pure-NumPy implementations of the TopoAnalysis grid kernels.

Every function here mirrors, name for name and result for result, one of the
functions in the compiled ``TopoAnalysis._topoanalysis`` extension.  They are
the reference semantics: the C++ code is validated against them in the test
suite, and they are what :mod:`TopoAnalysis.kernels` falls back to when the
extension is unavailable.

The depression-filling routine implements the Priority-Flood algorithms of

    Barnes, R., Lehman, C., Mulla, D. (2014). "Priority-flood: An optimal
    depression-filling and watershed-labeling algorithm for digital elevation
    models." *Computers & Geosciences* 62, 117-127.

These implementations are correct but roughly one to two orders of magnitude
slower than the compiled kernels on real DEMs.  If you are processing
anything larger than a few hundred thousand cells, build the extension::

    pip install -e .

### `default_seeds(elevations)`

Flat indices of the cells a flood should start from.

These are the grid perimeter plus every valid cell touching a no-data
cell -- the same set that TopoToolbox's ``fillsinks`` keeps in its marker
image.

Parameters
----------
elevations : ndarray
    2-D elevation grid; ``NaN`` marks no-data.

Returns
-------
ndarray of int64
    Flat (C-order) indices of the seed cells.

### `priority_flood(elevations, closed=None, seeds=None, mode='flat', epsilon=0.0, cellsize=1.0, max_pit_depth=0.0, track_visited=False, reference_algorithm=False)`

Fill depressions in ``elevations`` in place (Barnes et al. 2014).

Parameters
----------
elevations : ndarray
    2-D float64 grid, modified in place.  ``NaN`` is no-data.
closed : ndarray of uint8, optional
    Cells with a non-zero entry are never processed.  Use this to restrict
    the flood to a mask.
seeds : ndarray of int64, optional
    Flat indices to start from.  Defaults to :func:`default_seeds`.
mode : {'flat', 'epsilon'}
    ``'flat'`` fills depressions to a level surface (Algorithm 2), which
    is what TopoToolbox's ``fillsinks`` does.  ``'epsilon'`` adds a small
    upslope increment so the result has no flats (Algorithm 4).
epsilon : float
    Gradient of the increment used by ``mode='epsilon'``.  Each filled
    cell is raised to ``parent + epsilon * cellsize * step``.  Zero means
    "use the smallest representable increment".
cellsize : float
    Grid spacing, used only to scale ``epsilon``.
max_pit_depth : float
    Depressions deeper than this are left unfilled.  Zero disables the
    check.  Unlike the original TopoAnalysis implementation, the flood
    still traverses such depressions, so terrain beyond them is reached.
track_visited : bool
    Also return a mask of the cells the flood reached.
reference_algorithm : bool
    Use Algorithm 1 (everything through the priority queue) instead of the
    faster Algorithm 2.  The filled surface is identical; this exists for
    testing.

Returns
-------
dict
    ``{'visited': ndarray|None, 'cells_filled': int, 'max_fill_depth': float}``

### `flow_directions(elevations, cellsize_grid=None, cellsize=1.0)`

Steepest-descent D8 flow directions as ArcGIS codes.

The drop to each neighbour is divided by the centre-to-centre distance, so
a diagonal neighbour only wins if its *gradient* is steeper.  Ties go to
the first direction in the order E, SE, S, SW, W, NW, N, NE.

Cells with no lower neighbour -- pits, and no-data cells -- get code 0.

### `receivers(codes)`

Flat index of the cell each cell drains to, or -1 if there is none.

### `topological_order(recv)`

Order the cells donors-first.

Returns ``(order, n_cells_on_cycles)``.  Cells that lie on a cycle -- only
possible in a hand-edited or externally supplied flow-direction grid --
cannot be ordered and are appended at the end in index order.

### `accumulate(recv, order, weights, gate=None)`

Accumulate ``weights`` downstream along the D8 network.

``order`` is accepted for signature compatibility with the compiled
kernel; the result does not depend on which valid topological order is
used, so this implementation derives its own level structure -- which
lets it run vectorised rather than one cell at a time.

### `flow_length(recv, order, step_length)`

Longest upstream flow distance, and the main-stem donor direction.

Returns ``(length, from_codes)``.  ``from_codes[k]`` is the direction code
pointing from ``k`` back at the donor that supplies its longest path.
Ties are settled in favour of the donor with the lower flat index, so the
main stem does not depend on which topological order was used.

### `propagate_along_main_stem(recv, order, main_stem_from_codes, gate, values, mode)`

Carry or sum ``values`` downstream, but only along longest-flow paths.

``mode='carry'`` overwrites the receiver with the donor's value (used by
:class:`~TopoAnalysis.dem.Relief`); ``mode='sum'`` adds it (used by
:class:`~TopoAnalysis.dem.Ksi`).

### `upstream_mask(recv, outlets)`

Mask of every cell draining to any of ``outlets`` (flat indices).

### `chi(recv, order, area, step_length, outlets, A0, theta, trapezoid=False, max_length=0.0, mask=None)`

Integrate the chi coordinate upstream from ``outlets``.

``chi(cell) = chi(receiver) + (A0 / A)^theta * step_length``

With ``trapezoid=True`` the integrand is averaged between the cell and its
receiver, which is the rule TopoToolbox's ``chitransform`` uses.  Cells
outside the outlets' drainage areas stay ``NaN``.

Returns ``(chi, distance_from_outlet)``.

### `drainage_basins(recv, order, outlets=None, valid=None)`

Label drainage basins.

With ``outlets`` given, each cell is labelled with the 1-based position of
the outlet it drains to (0 if it reaches none).  Without, every terminal
cell -- pit or edge outlet -- seeds its own basin, matching TopoToolbox's
``drainagebasins(FD)``.

## `TopoAnalysis.demRecursionTools`

Profile extraction and chi/steepness fitting over recursive basin maps.

These functions consume the nested dictionaries produced by
:meth:`TopoAnalysis.dem.FlowDirectionD8.map_values_to_recursive_list`.

| Function | Summary |
|---|---|
| `area_elevation_for_mainstem_and_tributaries(outlet, flow_direction, elevation, area, theta=0.5, minimum_area=10000000.0)` |  |
| `best_ks_and_theta_with_wrss(elevation, flow_direction_or_length, area, outlet, xo=500)` |  |
| `best_ks_and_theta_with_wrss_list(ld_list, de, xo=500, maxiter=100, maxfun=200)` |  |
| `best_ks_theta(outlet, flow_direction, elevation, area, minimum_area)` |  |
| `best_ks_theta_wrss_for_outlet(outlet, flow_direction, elevation, area, minimum_area=10000000.0)` |  |
| `best_ks_with_r2_list(ld_list, de, theta, xo=500)` |  |
| `best_ks_with_wrss_list(ld_list, de, theta, xo=500)` |  |
| `chi_elevation(ld_list, de, theta, xo=500.0)` |  |
| `extract_chi_elevation_values(ld_list, de, theta, chi_o, elevation, chi, base_elevation, A_mdx=None, xo=500.0)` |  |
| `extract_dA_elevation_values(ld_list)` |  |
| `extract_profile_values(ld_list, xo=500.0, items=())` |  |
| `hi(elevation, flow_direction, dA, outlet)` |  |
| `hi_list(ld_list)` |  |
| `map_chi_profiles(elevation, flow_direction, area, outlet, minimum_area=1000000.0, theta=0.5, start_at=0.0, downstream=True, Ao=1000000.0)` |  |
| `uninformative_SS_list(ld_list, de, xo=500)` |  |

## `TopoAnalysis.plotting`

Longitudinal- and chi-profile plotting helpers.

| Function | Summary |
|---|---|
| `interactive_chi_profiles_and_map_view(prefix, code, plot_code, dem, fd, area, hillshade, minimum_area=10000000.0, theta=0.5, Ao=1000000.0)` |  |
| `plot_chi_profiles(elevation, flow_direction, area, outlet, plot_code, minimum_area=1000000.0, figure=None, theta=0.5, start_at=0.0, downstream=True, Ao=1000000.0)` |  |
| `plot_chi_profiles_with_outlet_code(prefix, code, plot_code, dem, fd, area, minimum_area=10000000.0, theta=0.5)` |  |
| `plot_downstream_profile(elevation, flow_direction, outlet, plot_code, downstream=True, start_at=0.0, mean_pixel_dimension=None, figure=None)` |  |
| `plot_profiles(elevation, flow_direction, area, outlet, plot_code, minimum_area=1000000.0, figure=None)` |  |
| `plot_profiles_with_outlet_code(prefix, code, plot_code, dem, fd, area, minimum_area=10000000.0)` |  |
| `plot_recursive_upstream_profiles(elevation, flow_direction, area, outlet, plot_code, downstream=False, start_at=0.0, figure=None, minimum_area=1000000.0)` |  |

## `TopoAnalysis.datasets`

Utilities for generating synthetic data and loading examples

| Function | Summary |
|---|---|
| `sinusoid_grid(ny, nx, width, amp=1, sig=0, slope_y=None)` | Returns a synthetic landscape with sinusoidal ridges and valleys |
| `triangle_grid(ny, nx, width, amp=1, sig=0, slope_y=None)` | Returns a synthetic landscape with triangular ridges and valleys |

## `TopoAnalysis.utils`

Convenience wrappers for fitting channel steepness at a single outlet.

| Function | Summary |
|---|---|
| `calc_ks_for_outlet(outlet, theta, **kwargs)` | Fit channel steepness to the basin draining to ``outlet``. |

## `TopoAnalysis.analysis`

Quadrat sampling of a raster.

### `Quadrats`



| Method | Summary |
|---|---|
| `load_data(filename, band=1)` | Read a raster band with rasterio. |
| `make_quadrats(dx, dy=None)` |  |
| `map_quadrats(func, **kwargs)` |  |
| `plot(values, ax=None, **kwargs)` |  |
| `quiver(u, v, ax=None, **kwargs)` |  |

## `TopoAnalysis.cli`

Command-line entry point: ``topoanalysis-process``.

Runs the standard pipeline over a DEM and writes the derived grids as
GeoTIFFs beside it::

    topoanalysis-process srtm.tif --outdir results --theta 0.45 --a0 1e6

Use ``topoanalysis-process --info`` to report which compute backend is
active and whether GDAL is available.

| Function | Summary |
|---|---|
| `main(argv=None)` | Run the CLI.  Returns a process exit status. |

## `TopoAnalysis.MovingWindow`

Moving-window filters over a grid.

A :class:`MovingWindow` applies a reducing function to the cells inside a
window centred on every cell.  Subclasses choose the window shape; a
concrete class supplies the reduction as :attr:`~MovingWindow.function`::

    class WindowMean(CircularMovingWindow):
        function = staticmethod(np.mean)

    grid.apply_moving_window(WindowMean(window_radius=200.0))

Cells outside the grid, and no-data cells, are dropped from the window
rather than filled, so the function always sees real data.

For the common reductions (min, max, mean, median, standard deviation)
:func:`TopoAnalysis.topotoolbox.localtopography` is far faster, because it
uses SciPy's separable filters instead of a Python loop.

### `CircularMovingWindow`

A disc-shaped window of radius ``round(window_radius / dx)`` cells.

### `MovingWindow`

Abstract base for moving-window filters.

Parameters
----------
window_radius : float
    Radius of the window in map units.  ``window_dimension`` is accepted
    as a synonym for backward compatibility.

| Method | Summary |
|---|---|
| `apply_moving_window(grid, dx, dtype)` | Run the window over every cell of ``grid``. |

### `RectangularMovingWindow`

A square window of side ``2 * round(window_radius / dx) + 1`` cells.

## `TopoAnalysis.error`

### `Error`

Base class for exceptions in this module.

### `InputError`

Exception raised for errors in the input.

Attributes:
    expr -- input expression in which the error occurred
    msg  -- explanation of the error

### `TransitionError`

Raised when an operation attempts a state transition that's not
allowed.

Attributes:
    prev -- state at beginning of transition
    next -- attempted new state
    msg  -- explanation of why the specific transition is not allowed

