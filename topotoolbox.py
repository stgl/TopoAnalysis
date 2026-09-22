"""TopoToolbox-equivalent terrain analysis.

Every function here reproduces the numerical result of the identically named
function in TopoToolbox 2 (Schwanghart & Scherler 2014), the MATLAB package,
so that a result computed in TopoAnalysis can be checked against one computed
in MATLAB.

The rest of TopoAnalysis keeps its own idioms -- drainage area in square
metres, ESRI-style hillshade, chi integrated with a left-endpoint rule.  This
module exists for the cases where you need the *same number*, not just the
same quantity.  :doc:`docs/topotoolbox_parity` tabulates the correspondence
and lists the places where exact agreement is not achievable.

Example
-------
::

    from TopoAnalysis import dem, topotoolbox as tt

    z  = dem.Elevation(gdal_filename='srtm.tif')
    zf = tt.fillsinks(z)                     # == GRIDobj/fillsinks
    fd = tt.flowobj(z, preprocess='fill')    # == FLOWobj(DEM,'preprocess','fill')
    a  = tt.flowacc(fd)                      # == flowacc(FD), in CELLS
    c  = tt.chitransform(fd, a, mn=0.45)     # == chitransform(S,A)

Conventions
-----------
* Grids are TopoAnalysis grid objects; results are returned as grid objects
  of the same shape and georeferencing.
* Drainage area is in **number of cells**, as TopoToolbox reports it -- not
  in square metres, which is what :class:`TopoAnalysis.dem.Area` returns.
* Row 0 is the north edge, matching both packages.

References
----------
Schwanghart, W. and Scherler, D. (2014). "TopoToolbox 2 -- MATLAB-based
software for topographic analysis and modeling in Earth surface sciences."
*Earth Surface Dynamics* 2, 1-7.
"""

from __future__ import annotations

import numpy as np

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import dem as _dem
    from . import error as Error
    from . import kernels
except ImportError:  # pragma: no cover
    import dem as _dem
    import error as Error
    import kernels

__all__ = [
    "fillsinks",
    "identifyflats",
    "imposemin",
    "flowobj",
    "flowacc",
    "flowdistance",
    "drainagebasins",
    "gradient8",
    "arcslope",
    "hillshade",
    "aspect",
    "curvature",
    "cellarea",
    "localtopography",
    "streamorder",
    "chitransform",
    "ksn",
    "mchi",
    "getcoordinates",
    "getextent",
    "coord2sub",
    "sub2coord",
    "STREAM_UNITS",
]

#: Unit conversions shared by ``gradient8`` and ``arcslope``.
_UNIT_FUNCS = {
    'tangent': lambda g: g,
    'degree': lambda g: np.degrees(np.arctan(g)),
    'radian': lambda g: np.arctan(g),
    'sine': lambda g: np.sin(np.arctan(g)),
    'percent': lambda g: g * 100.0,
}

STREAM_UNITS = tuple(sorted(_UNIT_FUNCS))


def _like(source, values, cls=None):
    """A new grid with ``source``'s georeferencing and ``values`` as data."""
    target = (cls or _dem.BaseSpatialGrid)()
    target._copy_info_from_grid(source, set_zeros=True)
    target._griddata = values
    return target


def _convert_unit(values, unit):
    try:
        return _UNIT_FUNCS[unit](values)
    except KeyError:
        raise Error.InputError(
            'unit', "unit must be one of {0}; got {1!r}".format(STREAM_UNITS, unit))


# ---------------------------------------------------------------------------
# Depression handling
# ---------------------------------------------------------------------------


def fillsinks(elevation, maxdepth=None):
    """Fill topographic depressions to a flat surface.

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
    """
    from scipy.ndimage import generate_binary_structure, label

    original = np.array(elevation._griddata, dtype=np.float64, copy=True)

    filled = _dem.Elevation()
    filled._copy_info_from_grid(elevation)
    # An explicit copy: np.ascontiguousarray returns the *same* array when the
    # input is already float64 and contiguous, and the flood works in place.
    filled._griddata = original.copy()

    if maxdepth is None:
        kernels.priority_flood(filled._griddata, mode='flat')
        return filled

    # MATLAB's depth-limited branch is all-or-nothing per *connected
    # depression*, not per cell: whenever a depression's deepest point
    # exceeds maxdepth, that cell is added to the marker set, which restores
    # the whole depression to its original elevations. Unblocking it can
    # split the depression into shallower ones, which then do get filled, so
    # the process repeats until nothing exceeds the limit.
    maxdepth = float(maxdepth)
    structure = generate_binary_structure(2, 2)
    seeds = list(kernels.default_seeds(original).tolist())

    for _ in range(original.size + 1):
        surface = original.copy()
        kernels.priority_flood(surface, seeds=np.array(seeds, dtype=np.int64),
                               mode='flat')
        depth = surface - original
        labels, count = label(depth > 0, structure=structure)
        if count == 0:
            break

        added = []
        for component in range(1, count + 1):
            inside = labels == component
            if depth[inside].max() <= maxdepth:
                continue
            # Seed the deepest cell of the component.
            flat_index = int(np.flatnonzero((labels == component).reshape(-1))[
                int(np.argmax(depth[inside]))])
            added.append(flat_index)
        if not added:
            break
        seeds.extend(added)

    filled._griddata = surface
    return filled


def identifyflats(elevation, output='both'):
    """Locate flat areas, their sills and the closed basins.

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
    """
    from scipy.ndimage import binary_dilation, generate_binary_structure, label

    z = np.asarray(elevation._griddata, dtype=np.float64)
    nodata = np.isnan(z)
    work = np.where(nodata, -np.inf, z)
    structure = generate_binary_structure(2, 2)      # ones(3)

    # imerode with a 3x3 of ones is the local minimum over the 3x3 INCLUDING
    # the centre, padded with +inf; `== dem` therefore means "no strictly
    # lower neighbour".  There is deliberately no requirement of an *equal*
    # neighbour, so a cell lower than all eight of its neighbours is a flat.
    local_min = work.copy()
    for d in range(8):
        di, dj = int(kernels.D8_DI[d]), int(kernels.D8_DJ[d])
        local_min = np.minimum(local_min, _shift(work, di, dj, np.inf))
    flats = (local_min == work)
    if nodata.any():
        flats &= ~nodata

    # Regional minima, before the rim is cleared: connected components of
    # `flats` are equal-valued plateaus with no lower neighbour outside them.
    regional_minima = flats.copy()

    flats[0, :] = False
    flats[-1, :] = False
    flats[:, 0] = False
    flats[:, -1] = False
    if nodata.any():
        flats[binary_dilation(nodata, structure=structure)] = False

    # Sills: the 3x3 maximum taken over the flat cells only, compared with
    # the cell's own elevation, excluding the flats themselves.
    marker = np.full(z.shape, -np.inf)
    marker[flats] = work[flats]
    dilated = marker.copy()
    for d in range(8):
        di, dj = int(kernels.D8_DI[d]), int(kernels.D8_DJ[d])
        dilated = np.maximum(dilated, _shift(marker, di, dj, -np.inf))
    sills = (dilated == work) & ~flats
    if nodata.any():
        sills[nodata] = False

    if output in ('closedbasins', 'all'):
        # imclearborder: drop the components that touch the grid edge.
        labels, count = label(regional_minima, structure=structure)
        on_border = set(labels[0, :].tolist()) | set(labels[-1, :].tolist()) \
            | set(labels[:, 0].tolist()) | set(labels[:, -1].tolist())
        on_border.discard(0)
        closed = regional_minima & ~np.isin(labels, list(on_border))
        if nodata.any():
            closed |= nodata
            labels, _ = label(closed, structure=structure)
            on_border = set(labels[0, :].tolist()) | set(labels[-1, :].tolist()) \
                | set(labels[:, 0].tolist()) | set(labels[:, -1].tolist())
            on_border.discard(0)
            closed &= ~np.isin(labels, list(on_border))
            closed[nodata] = False

    flats_grid = _like(elevation, flats)
    sills_grid = _like(elevation, sills)
    if output == 'flats':
        return flats_grid
    if output == 'sills':
        return sills_grid
    if output == 'closedbasins':
        return _like(elevation, closed)
    if output == 'all':
        return flats_grid, sills_grid, _like(elevation, closed)
    return flats_grid, sills_grid


def _shift(grid, di, dj, fill):
    """``grid`` translated so element (i, j) holds ``grid[i+di, j+dj]``."""
    out = np.full_like(grid, fill)
    ny, nx = grid.shape
    src_i = slice(max(di, 0), ny + min(di, 0))
    src_j = slice(max(dj, 0), nx + min(dj, 0))
    dst_i = slice(max(-di, 0), ny + min(-di, 0))
    dst_j = slice(max(-dj, 0), nx + min(-dj, 0))
    out[dst_i, dst_j] = grid[src_i, src_j]
    return out


def imposemin(flow_direction, elevation, sl=0.0):
    """Carve elevations so every flow path descends at gradient >= ``sl``.

    Equivalent to ``FLOWobj/imposemin(FD, DEM, sl)``.  Cells are only ever
    lowered, so channels are cut through obstructions instead of being
    drowned by them.

    Returns
    -------
    Elevation
        A new grid.
    """
    carved = _dem.Elevation()
    carved._copy_info_from_grid(elevation)
    carved._griddata = np.array(elevation._griddata, dtype=np.float64, copy=True)
    kernels.imposemin(
        flow_direction.receivers,
        flow_direction.topological_order,
        np.ascontiguousarray(flow_direction.step_lengths(), dtype=np.float64),
        float(sl),
        carved._griddata)
    return carved


# ---------------------------------------------------------------------------
# Flow
# ---------------------------------------------------------------------------


def flowobj(elevation, preprocess='fill'):
    """Build a D8 flow-direction grid the way ``FLOWobj`` does.

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
    """
    if preprocess not in ('fill', 'none'):
        raise Error.InputError(
            'preprocess', "preprocess must be 'fill' or 'none'; 'carve' is not "
                          "implemented (see the module docstring)")

    if preprocess == 'fill':
        # A small gradient across filled flats, so every cell has a defined
        # downslope neighbour.  A pure flat fill would leave plateaus with
        # no direction at all.
        surface = _dem.FilledElevation(elevation=elevation)
    else:
        surface = elevation

    return _dem.FlowDirectionD8(flooded_dem=surface)


def flowacc(flow_direction, weights=None):
    """Flow accumulation in **number of cells**.

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
    """
    shape = flow_direction._griddata.shape
    if weights is None:
        w = np.ones(shape, dtype=np.float64)
    else:
        w = np.ascontiguousarray(
            getattr(weights, '_griddata', weights), dtype=np.float64)

    accumulated = kernels.accumulate(
        flow_direction.receivers, flow_direction.topological_order, w)
    return _like(flow_direction, np.asarray(accumulated))


def flowdistance(flow_direction, direction='upstream'):
    """Along-flow distance in map units.

    Equivalent to ``FLOWobj/flowdistance(FD, direction)`` with no seed
    points.

    Parameters
    ----------
    direction : {'upstream', 'downstream', 'maxdownstream'}
        ``'upstream'`` gives the distance from each cell down to its outlet.
        ``'downstream'`` (a synonym for ``'maxdownstream'``) gives the length
        of the longest flow path reaching each cell -- the quantity
        :class:`TopoAnalysis.dem.FlowLength` computes.
    """
    step = np.ascontiguousarray(flow_direction.step_lengths(), dtype=np.float64)
    receivers = flow_direction.receivers
    order = flow_direction.topological_order

    if direction == 'upstream':
        distance = kernels.downstream_distance(receivers, order, step)
    elif direction in ('downstream', 'maxdownstream'):
        distance, _ = kernels.flow_length(receivers, order, step)
    else:
        raise Error.InputError(
            'direction',
            "direction must be 'upstream', 'downstream' or 'maxdownstream'")

    return _like(flow_direction, np.asarray(distance))


def drainagebasins(flow_direction, outlets=None, valid=None):
    """Label drainage basins.

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
    """
    outlet_indices = None
    if outlets is not None:
        flat = []
        for (row, col) in flow_direction._xy_to_rowscols(outlets):
            if row is not None:
                flat.append(row * flow_direction._georef_info.nx + col)
        outlet_indices = np.array(flat, dtype=np.int64)

    valid_mask = None
    if valid is not None:
        valid_mask = np.ascontiguousarray(
            (np.asarray(getattr(valid, '_griddata', valid)) != 0).astype(np.uint8))

    labels = kernels.drainage_basins(
        flow_direction.receivers, flow_direction.topological_order,
        outlet_indices, valid_mask)
    return _like(flow_direction, np.asarray(labels))


def streamorder(flow_direction, stream_mask, order_type='strahler'):
    """Strahler or Shreve stream order.

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
    """
    if order_type not in ('strahler', 'shreve'):
        raise Error.InputError('order_type', "order_type must be 'strahler' or 'shreve'")

    mask = np.ascontiguousarray(
        (np.asarray(getattr(stream_mask, '_griddata', stream_mask)) != 0).astype(np.uint8))
    orders = kernels.stream_order(
        flow_direction.receivers, flow_direction.topological_order, mask, order_type)
    return _like(flow_direction, np.asarray(orders))


# ---------------------------------------------------------------------------
# Terrain attributes
# ---------------------------------------------------------------------------


def gradient8(elevation, unit='tangent'):
    """Steepest downward gradient over the 8 neighbours.

    Equivalent to ``GRIDobj/gradient8``.  The cardinal and diagonal drops are
    divided by ``cellsize`` and ``sqrt(2)*cellsize`` respectively, and the
    result is clamped at 0, so pits and flats give exactly 0 rather than a
    negative gradient.

    Parameters
    ----------
    unit : {'tangent', 'degree', 'radian', 'sine', 'percent'}
    """
    z = np.asarray(elevation._griddata, dtype=np.float64)
    cs = elevation._georef_info.dx
    nodata = np.isnan(z)
    work = np.where(nodata, np.inf, z)

    gradient = np.zeros(z.shape, dtype=np.float64)
    for d in range(8):
        di, dj = int(kernels.D8_DI[d]), int(kernels.D8_DJ[d])
        # +inf outside the grid, so an off-grid neighbour never wins the
        # minimum (TopoToolbox erodes with an infinite pad for the same
        # reason).
        neighbour = _shift(work, di, dj, np.inf)
        distance = cs * kernels.D8_DIST[d]
        gradient = np.maximum(gradient, (work - neighbour) / distance)

    gradient[nodata] = np.nan
    return _like(elevation, _convert_unit(gradient, unit))


def arcslope(elevation, unit='tangent'):
    """ArcGIS/Horn 3x3 slope magnitude.

    Equivalent to ``GRIDobj/arcslope``.  Uses the 1-2-1 weighted 3x3 stencil
    with an ``8 * cellsize`` denominator -- there is no ``sqrt(2)`` anywhere;
    the diagonals enter through the weights.

    ``NaN`` holes are filled from their nearest valid neighbour before
    differencing, then restored, which is how TopoToolbox handles them here
    (and is the opposite of what :func:`gradient8` does).
    """
    from scipy.ndimage import distance_transform_edt
    from scipy.signal import convolve2d

    z = np.asarray(elevation._griddata, dtype=np.float64)
    cs = elevation._georef_info.dx
    nodata = np.isnan(z)

    padded = np.pad(z, 1, mode='edge')
    padded_nodata = np.pad(nodata, 1, mode='edge')
    if padded_nodata.any():
        _, indices = distance_transform_edt(padded_nodata, return_indices=True)
        padded = padded[tuple(indices)]

    kernel_y = np.array([[-1.0, -2.0, -1.0], [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]])
    kernel_x = kernel_y.T
    gx = convolve2d(padded, kernel_x, mode='valid') / (8.0 * cs)
    gy = convolve2d(padded, kernel_y, mode='valid') / (8.0 * cs)

    slope = np.hypot(gx, gy)
    slope[nodata] = np.nan
    return _like(elevation, _convert_unit(slope, unit))


def _surfnorm(z):
    """MATLAB ``surfnorm(Z)`` for the single-argument case.

    Returns ``(Nx, Ny, Nz)``, the unit surface normal with derivatives taken
    per cell index: ``Nx`` positive where elevation decreases eastward,
    ``Ny`` positive where it decreases southward (row-increasing), ``Nz``
    always positive.

    The edge treatment is the second-order one-sided difference, which is
    what MATLAB's quadratic ghost-cell extrapolation reduces to.
    """
    z = np.asarray(z, dtype=np.float64)
    if min(z.shape) < 3:
        raise Error.InputError('surfnorm', 'the grid must be at least 3x3')

    dz_drow, dz_dcol = np.gradient(z, edge_order=2)
    nx = -dz_dcol
    ny = -dz_drow
    nz = np.ones_like(z)

    magnitude = np.sqrt(nx * nx + ny * ny + nz * nz)
    magnitude[magnitude == 0] = np.finfo(np.float64).eps
    return nx / magnitude, ny / magnitude, nz / magnitude


def hillshade(elevation, azimuth=315.0, altitude=60.0, exaggerate=1.0):
    """Shaded relief as the cosine of the illumination incidence angle.

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
    """
    cs = elevation._georef_info.dx
    z = np.asarray(elevation._griddata, dtype=np.float64)

    # The compass bearing is converted by subtracting 90, which is correct
    # here because surfnorm's +y axis points along increasing row (south).
    azimuth_rad = np.radians(azimuth - 90.0)
    altitude_rad = np.radians(altitude)
    cos_elevation = np.cos(altitude_rad)
    sx = cos_elevation * np.cos(azimuth_rad)
    sy = cos_elevation * np.sin(azimuth_rad)
    sz = np.sin(altitude_rad)

    nx, ny, nz = _surfnorm(z / cs * exaggerate)
    return _like(elevation, nx * sx + ny * sy + nz * sz)


def aspect(elevation, classify=False):
    """Down-slope direction in degrees clockwise from north.

    Equivalent to ``GRIDobj/aspect``.  0 = north, 90 = east, 180 = south,
    270 = west.  Perfectly flat cells come out as 90, exactly as they do in
    TopoToolbox.

    Parameters
    ----------
    classify : bool
        Return the 8-class Gomez-Plaza wetness ranking (1 = N, 3 = NE,
        5 = E, 7 = SE, 8 = S, 6 = SW, 4 = W, 2 = NW) instead of degrees.
    """
    nx, ny, _ = _surfnorm(np.asarray(elevation._griddata, dtype=np.float64))
    degrees = np.mod(90.0 + np.degrees(np.arctan2(ny, nx)), 360.0)
    # atan2(0, 0) is 0 in MATLAB, so a flat cell comes out as 90 (east).
    # Reproducing that through NumPy depends on the sign of zero, which
    # np.gradient does not preserve, so set it explicitly.
    degrees[(nx == 0.0) & (ny == 0.0)] = 90.0

    if not classify:
        return _like(elevation, degrees)

    edges = np.arange(0, 361, 45)
    classes = np.array([1, 3, 5, 7, 8, 6, 4, 2], dtype=np.uint8)
    bins = np.digitize(degrees, edges, right=False) - 1
    out = np.zeros(degrees.shape, dtype=np.uint8)
    inside = (bins >= 0) & (bins < classes.size) & ~np.isnan(degrees)
    out[inside] = classes[bins[inside]]
    return _like(elevation, out)


def curvature(elevation, ctype='profc'):
    """Surface curvature from the Evans 3x3 quadratic fit.

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
    """
    from scipy.signal import convolve2d

    valid_types = ('profc', 'planc', 'tangc', 'meanc', 'total')
    if ctype not in valid_types:
        raise Error.InputError('ctype', 'ctype must be one of {0}'.format(valid_types))

    z = np.asarray(elevation._griddata, dtype=np.float64)
    cs = elevation._georef_info.dx
    padded = np.pad(z, 1, mode='symmetric')

    def conv(kernel):
        return convolve2d(padded, kernel, mode='valid')

    fx = conv(np.array([[-1.0, 0, 1.0]] * 3) / (6 * cs))
    fy = conv(np.array([[1.0, 1.0, 1.0], [0, 0, 0], [-1.0, -1.0, -1.0]]) / (6 * cs))
    kxx = np.array([[1.0, -2.0, 1.0]] * 3) / (3 * cs ** 2)
    fxx = conv(kxx)
    fyy = conv(kxx.T)
    fxy = conv(np.array([[-1.0, 0, 1.0], [0, 0, 0], [1.0, 0, -1.0]]) / (4 * cs ** 2))

    p2 = fx ** 2 + fy ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        if ctype == 'profc':
            curv = -(fx ** 2 * fxx + 2 * fx * fy * fxy + fy ** 2 * fyy) / (
                p2 * (1 + p2) ** 1.5)
        elif ctype == 'tangc':
            curv = -(fy ** 2 * fxx - 2 * fx * fy * fxy + fx ** 2 * fyy) / (
                p2 * (1 + p2) ** 0.5)
        elif ctype == 'planc':
            curv = -(fy ** 2 * fxx - 2 * fx * fy * fxy + fx ** 2 * fyy) / (p2 ** 1.5)
        elif ctype == 'meanc':
            curv = -((1 + fy ** 2) * fxx - 2 * fxy * fx * fy + (1 + fx ** 2) * fyy) / (
                2 * (p2 + 1) ** 1.5)
        else:
            curv = fxx ** 2 + 2 * fxy ** 2 + fyy ** 2

    # Order matters: zero the 0/0 cells first, then restore NaN only where
    # the DEM itself was no-data.
    curv[~np.isfinite(curv)] = 0.0
    curv[np.isnan(z)] = np.nan
    return _like(elevation, curv)


#: Earth's total surface area in km^2, the constant TopoToolbox's ``cellarea``
#: uses.  It implies a spherical radius of 6 371 047 m.
EARTH_SURFACE_AREA_KM2 = 510072000.0


def cellarea(grid, unit='m'):
    """Area of each cell of a geographic (lat/lon) grid.

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
    """
    unit = unit.lower()
    if unit not in ('m', 'km'):
        raise Error.InputError('unit', "unit must be 'm' or 'km'")
    scale = 1e6 if unit == 'm' else 1.0

    info = grid._georef_info
    dlon = abs(info.geoTransform[1])
    dlat = abs(info.geoTransform[5])

    # Row 0 is the north edge, so latitudes descend with row index.
    lat_center = info.yllcenter + (info.ny - 1 - np.arange(info.ny)) * dlat
    lat_north = np.radians(lat_center + dlat / 2.0)
    lat_south = np.radians(lat_center - dlat / 2.0)

    fraction = (np.radians(dlon) * (np.sin(lat_north) - np.sin(lat_south))) / (4 * np.pi)
    per_row = EARTH_SURFACE_AREA_KM2 * fraction * scale
    return _like(grid, np.repeat(per_row[:, None], info.nx, axis=1))


def localtopography(elevation, radius=5000.0, type='range'):
    """Statistic of elevation within a circular window.

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
    """
    from scipy.ndimage import (convolve, distance_transform_edt, generic_filter,
                               maximum_filter, median_filter, minimum_filter)

    z = np.asarray(elevation._griddata, dtype=np.float64)
    cs = elevation._georef_info.dx
    radius_px = int(np.ceil(radius / cs))
    yy, xx = np.ogrid[-radius_px:radius_px + 1, -radius_px:radius_px + 1]
    footprint = (xx * xx + yy * yy) <= radius_px * radius_px

    nodata = np.isnan(z)
    if nodata.any() and type in ('mean', 'median', 'std'):
        _, indices = distance_transform_edt(nodata, return_indices=True)
        work = z[tuple(indices)]
    else:
        work = z

    if type == 'max':
        result = maximum_filter(np.where(nodata, -np.inf, work), footprint=footprint)
    elif type == 'min':
        result = minimum_filter(np.where(nodata, np.inf, work), footprint=footprint)
    elif type == 'range':
        result = (maximum_filter(np.where(nodata, -np.inf, work), footprint=footprint)
                  - minimum_filter(np.where(nodata, np.inf, work), footprint=footprint))
    elif type in ('mean', 'average'):
        # The disc, not the bounding square: uniform_filter(size=2r+1) would
        # include the four corners, which lie outside the window.
        kernel = footprint.astype(np.float64)
        kernel /= kernel.sum()
        result = convolve(work, kernel, mode='reflect')
    elif type == 'median':
        result = median_filter(work, footprint=footprint, mode='reflect')
    elif type == 'std':
        # ddof=1: MATLAB's stdfilt is the sample standard deviation.
        result = generic_filter(work, lambda v: np.std(v, ddof=1),
                                footprint=footprint, mode='reflect')
    else:
        raise Error.InputError(
            'type', "type must be one of 'range', 'max', 'min', 'mean', 'median', 'std'")

    result = np.asarray(result, dtype=np.float64)
    result[nodata] = np.nan
    return _like(elevation, result)


# ---------------------------------------------------------------------------
# Chi and channel steepness
# ---------------------------------------------------------------------------


def chitransform(flow_direction, area, mn=0.45, a0=1e6, outlets=None,
                 correctcellsize=True, stream_mask=None):
    r"""The chi integral transform.

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
    """
    a = np.asarray(area._griddata, dtype=np.float64)
    if correctcellsize:
        # TopoToolbox multiplies by S.cellsize^2 unconditionally; feeding it
        # an area that is already in m^2 squares the cell size twice.
        a = a * (flow_direction._georef_info.dx ** 2)

    receivers = flow_direction.receivers
    order = flow_direction.topological_order
    step = np.ascontiguousarray(flow_direction.step_lengths(), dtype=np.float64)

    if outlets is None:
        outlet_indices = np.flatnonzero(receivers.reshape(-1) < 0).astype(np.int64)
    else:
        flat = []
        for (row, col) in flow_direction._xy_to_rowscols(outlets):
            if row is not None:
                flat.append(row * flow_direction._georef_info.nx + col)
        outlet_indices = np.array(flat, dtype=np.int64)

    mask = None
    if stream_mask is not None:
        mask = np.ascontiguousarray(
            (np.asarray(getattr(stream_mask, '_griddata', stream_mask)) != 0).astype(np.uint8))

    chi_values, _ = kernels.chi(
        receivers, order, np.ascontiguousarray(a), step, outlet_indices,
        float(a0), float(mn), True, 0.0, mask)
    return _like(flow_direction, np.asarray(chi_values))


def ksn(flow_direction, elevation, area, theta=0.45, correctcellsize=True,
        min_gradient=1e-5):
    r"""Normalized channel steepness index, :math:`k_{sn} = S A^{\theta}`.

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
    """
    a = np.asarray(area._griddata, dtype=np.float64)
    if correctcellsize:
        a = a * (flow_direction._georef_info.dx ** 2)

    carved = imposemin(flow_direction, elevation, sl=min_gradient)

    receivers = flow_direction.receivers.reshape(-1)
    z = np.asarray(carved._griddata, dtype=np.float64).reshape(-1)
    step = flow_direction.step_lengths().reshape(-1)

    gradient = np.zeros(z.shape, dtype=np.float64)
    donors = np.flatnonzero(receivers >= 0)
    gradient[donors] = (z[donors] - z[receivers[donors]]) / step[donors]

    shape = flow_direction._griddata.shape
    with np.errstate(divide='ignore', invalid='ignore'):
        values = gradient.reshape(shape) / (a ** (-theta))
    return _like(flow_direction, values)


def mchi(flow_direction, elevation, chi):
    r"""Gradient of elevation in chi space, :math:`M_{\chi}`.

    Equivalent to ``STREAMobj/mchi``: a forward difference assigned to the
    upstream cell, with no carving and no smoothing.

    ``M_chi`` scales with the ``a0`` baked into ``chi``: computing chi with
    ``a0=1`` puts it in the same units as :func:`ksn`, and
    ``mchi(a0) * a0**mn`` is invariant.

    It is *close to* but not equal to :func:`ksn`, despite what TopoToolbox's
    own docstring says: ``ksn`` carves the DEM first and evaluates the
    integrand at the upstream cell, while ``mchi`` uses the raw elevations
    and inherits chi's trapezoidal averaging.
    """
    receivers = flow_direction.receivers.reshape(-1)
    z = np.asarray(elevation._griddata, dtype=np.float64).reshape(-1)
    c = np.asarray(getattr(chi, '_griddata', chi), dtype=np.float64).reshape(-1)

    values = np.zeros(z.shape, dtype=np.float64)
    donors = np.flatnonzero(receivers >= 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        values[donors] = ((z[donors] - z[receivers[donors]])
                          / (c[donors] - c[receivers[donors]]))
    return _like(flow_direction, values.reshape(flow_direction._griddata.shape))


# ---------------------------------------------------------------------------
# Referencing
# ---------------------------------------------------------------------------


def getcoordinates(grid):
    """``(x, y)`` cell-centre coordinate vectors.

    Equivalent to ``GRIDobj/getcoordinates``.  ``x`` runs west to east;
    ``y`` runs **north to south**, matching the row order of the grid.
    """
    info = grid._georef_info
    x = info.xllcenter + np.arange(info.nx) * info.dx
    y = info.yllcenter + (info.ny - 1 - np.arange(info.ny)) * info.dx
    return x, y


def getextent(grid):
    """``(xmin, xmax, ymin, ymax)`` of the outermost **cell centres**.

    Equivalent to ``GRIDobj/getextent``.  Note this is half a cell inside
    :meth:`TopoAnalysis.dem.BaseSpatialGrid.extent`, which returns the pixel
    edges that matplotlib wants.
    """
    x, y = getcoordinates(grid)
    return (float(x.min()), float(x.max()), float(y.min()), float(y.max()))


def _round_half_away(values):
    """MATLAB's ``round``: halves go away from zero, not to even."""
    values = np.asarray(values, dtype=np.float64)
    return np.where(values >= 0, np.floor(values + 0.5), np.ceil(values - 0.5))


def coord2sub(grid, x, y):
    """Map coordinates to 0-based ``(row, col)`` subscripts.

    Equivalent to ``GRIDobj/coord2sub`` apart from the index base.  Points
    outside the grid give ``NaN``.

    Notes
    -----
    Halves are rounded away from zero, as MATLAB does.  Python's built-in
    ``round`` -- which
    :meth:`TopoAnalysis.dem.BaseSpatialGrid._xy_to_rowscols` uses -- rounds
    halves to even, so the two disagree for points landing exactly on a cell
    boundary.
    """
    info = grid._georef_info
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    y = np.atleast_1d(np.asarray(y, dtype=np.float64))

    north = info.yllcenter + (info.ny - 1) * info.dx
    col = _round_half_away((x - info.xllcenter) / info.dx)
    row = _round_half_away((north - y) / info.dx)

    outside = (col < 0) | (col > info.nx - 1) | (row < 0) | (row > info.ny - 1) \
        | np.isnan(x) | np.isnan(y)
    row = np.where(outside, np.nan, row)
    col = np.where(outside, np.nan, col)
    return row, col


def sub2coord(grid, row, col):
    """0-based ``(row, col)`` subscripts to cell-centre map coordinates.

    Equivalent to ``GRIDobj/sub2coord``.  Subscripts outside the grid are
    extrapolated rather than rejected, as in TopoToolbox.
    """
    info = grid._georef_info
    row = np.asarray(row, dtype=np.float64)
    col = np.asarray(col, dtype=np.float64)
    x = info.xllcenter + col * info.dx
    y = info.yllcenter + (info.ny - 1 - row) * info.dx
    return x, y
