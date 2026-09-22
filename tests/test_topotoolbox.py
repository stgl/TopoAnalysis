"""Tests for the TopoToolbox-equivalent functions.

Without MATLAB on hand, parity is established three ways:

* against closed-form answers on surfaces whose value is known exactly
  (planes, quadratics, spherical cells);
* against the documented definitions transcribed from the MATLAB source
  (stencils, constants, conventions); and
* against TopoAnalysis' own results where the two are meant to agree, with
  the documented conversion applied.
"""

from __future__ import annotations

import numpy as np
import pytest

from TopoAnalysis import dem, error, topotoolbox as tt


# ---------------------------------------------------------------------------
# Filling and flats
# ---------------------------------------------------------------------------


def test_fillsinks_is_a_flat_fill(single_pit):
    filled = tt.fillsinks(single_pit)
    # TopoToolbox's fillsinks adds no epsilon: the depression becomes an
    # exactly flat plateau at the spill elevation.
    assert np.all(filled._griddata == 10.0)
    assert single_pit._griddata[2, 2] == 1.0, 'the input must not be modified'


def test_fillsinks_respects_maxdepth(single_pit):
    filled = tt.fillsinks(single_pit, maxdepth=2.0)
    assert filled._griddata[2, 2] == pytest.approx(1.0)


def test_fillsinks_treats_nodata_as_an_outlet():
    z = np.array([
        [10, 10, 10, 10, 10],
        [10, 5, 5, 5, 10],
        [np.nan, 5, 1, 5, 10],
        [10, 5, 5, 5, 10],
        [10, 10, 10, 10, 10]], dtype=float)
    filled = tt.fillsinks(dem.Elevation(dx=1.0, grid=z))
    assert filled._griddata[2, 1] == pytest.approx(5.0)
    assert np.isnan(filled._griddata[2, 0])


def test_identifyflats_finds_the_plateau_and_its_sills(plateau):
    flats, sills = tt.identifyflats(plateau)
    # Every plateau cell has an equal neighbour and no lower one ...
    assert flats._griddata[1, 1]
    assert flats._griddata[1, 2]
    # ... except those beside the 2 m outlet, which do have a lower one.
    assert not flats._griddata[3, 2]
    assert sills._griddata[3, 2]
    assert not sills._griddata[1, 1]


def test_identifyflats_single_output(plateau):
    only_flats = tt.identifyflats(plateau, output='flats')
    both = tt.identifyflats(plateau)
    assert np.array_equal(only_flats._griddata, both[0]._griddata)


def test_imposemin_only_lowers(v_valley):
    fd = tt.flowobj(v_valley)
    carved = tt.imposemin(fd, v_valley, sl=1e-3)
    assert np.all(carved._griddata <= v_valley._griddata + 1e-12)


def test_imposemin_enforces_the_gradient(v_valley):
    fd = tt.flowobj(v_valley)
    sl = 1e-3
    carved = tt.imposemin(fd, v_valley, sl=sl)

    z = carved._griddata.reshape(-1)
    step = fd.step_lengths().reshape(-1)
    for k, r in enumerate(fd.receivers.reshape(-1)):
        if r >= 0:
            assert z[k] - z[r] >= sl * step[k] - 1e-9


# ---------------------------------------------------------------------------
# Flow
# ---------------------------------------------------------------------------


def test_flowacc_counts_cells_not_area(east_plane):
    fd = tt.flowobj(east_plane)
    counts = tt.flowacc(fd)
    # TopoToolbox's flowacc has unit weight per cell and counts the cell
    # itself, so the minimum is 1 -- not dx**2.
    assert counts._griddata.min() == 1.0
    assert counts._griddata[2].tolist() == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_flowacc_relates_to_Area_by_the_cell_area(rough_dome):
    fd = tt.flowobj(rough_dome)
    counts = tt.flowacc(fd)
    area = dem.Area(flow_direction=fd)
    cell = rough_dome._georef_info.dx ** 2
    assert np.allclose(counts._griddata * cell, area._griddata)


def test_flowacc_accepts_weights(east_plane):
    fd = tt.flowobj(east_plane)
    weights = np.full(east_plane._griddata.shape, 2.0)
    assert tt.flowacc(fd, weights)._griddata[2, 4] == pytest.approx(10.0)


def test_flowdistance_upstream_is_zero_at_the_outlet(east_plane):
    fd = tt.flowobj(east_plane)
    distance = tt.flowdistance(fd, 'upstream')
    # Distance to the outlet: 4 cells away at the western edge, 0 at the east.
    assert distance._griddata[2].tolist() == [4.0, 3.0, 2.0, 1.0, 0.0]


def test_flowdistance_downstream_is_the_flow_length(east_plane):
    fd = tt.flowobj(east_plane)
    tt_length = tt.flowdistance(fd, 'downstream')
    ta_length = dem.FlowLength(flow_direction=fd)
    assert np.allclose(tt_length._griddata, ta_length._griddata)


def test_flowdistance_rejects_an_unknown_direction(east_plane):
    fd = tt.flowobj(east_plane)
    with pytest.raises(error.InputError):
        tt.flowdistance(fd, 'sideways')


def test_drainagebasins_partition_matches_the_flow_network(rough_dome):
    fd = tt.flowobj(rough_dome)
    labels = tt.drainagebasins(fd)._griddata

    flat = labels.reshape(-1)
    for k, r in enumerate(fd.receivers.reshape(-1)):
        if r >= 0:
            assert flat[k] == flat[r]
    assert labels.min() >= 1


def test_drainagebasins_from_outlets(rough_dome):
    fd = tt.flowobj(rough_dome)
    outlets = fd._rowscols_to_xy(((20, 20), (10, 30)))
    labels = tt.drainagebasins(fd, outlets=outlets)._griddata

    assert set(np.unique(labels).tolist()) <= {0, 1, 2}
    assert labels[20, 20] == 1
    assert labels[10, 30] == 2


def test_streamorder_of_a_confluence():
    z = np.array([[3.0, 9.0, 3.0],
                  [9.0, 2.0, 9.0],
                  [9.0, 1.0, 9.0]])
    grid = dem.Elevation(dx=1.0, grid=z)
    fd = dem.FlowDirectionD8(elevation=grid)
    streams = np.array([[1, 0, 1], [0, 1, 0], [0, 1, 0]], dtype=np.uint8)

    strahler = tt.streamorder(fd, streams)._griddata
    assert strahler[0, 0] == 1 and strahler[1, 1] == 2
    shreve = tt.streamorder(fd, streams, 'shreve')._griddata
    assert shreve[1, 1] == 2


def test_streamorder_rejects_an_unknown_type(east_plane):
    fd = tt.flowobj(east_plane)
    with pytest.raises(error.InputError):
        tt.streamorder(fd, np.ones((5, 5)), 'horton')


def test_flowobj_rejects_carve(east_plane):
    """'carve' is not implemented and must say so rather than guess."""
    with pytest.raises(error.InputError):
        tt.flowobj(east_plane, preprocess='carve')


# ---------------------------------------------------------------------------
# Terrain attributes against closed-form answers
# ---------------------------------------------------------------------------


def test_gradient8_on_a_cardinal_plane(east_plane):
    assert tt.gradient8(east_plane)._griddata[2, 2] == pytest.approx(1.0)


def test_gradient8_on_a_diagonal_plane():
    """A 2 m drop over a sqrt(2) m step is a gradient of sqrt(2)."""
    z = np.add.outer(np.arange(5, 0, -1, dtype=float),
                     np.arange(5, 0, -1, dtype=float))
    grid = dem.Elevation(dx=1.0, grid=z)
    assert tt.gradient8(grid)._griddata[2, 2] == pytest.approx(np.sqrt(2.0))


def test_gradient8_is_never_negative(rough_dome):
    """The structuring elements include the centre, so pits give 0, not < 0."""
    assert tt.gradient8(rough_dome)._griddata.min() >= 0.0


def test_gradient8_is_zero_on_a_flat():
    flat = dem.Elevation(dx=1.0, grid=np.full((5, 5), 7.0))
    assert np.all(tt.gradient8(flat)._griddata == 0.0)


@pytest.mark.parametrize('unit,expected', [
    ('tangent', 1.0),
    ('degree', 45.0),
    ('radian', np.pi / 4),
    ('sine', np.sqrt(0.5)),
    ('percent', 100.0),
])
def test_gradient8_units(east_plane, unit, expected):
    assert tt.gradient8(east_plane, unit)._griddata[2, 2] == pytest.approx(expected)


def test_gradient8_rejects_an_unknown_unit(east_plane):
    with pytest.raises(error.InputError):
        tt.gradient8(east_plane, 'furlongs')


def test_arcslope_on_a_plane(east_plane):
    """The Horn 1-2-1 stencil reproduces a plane exactly."""
    assert tt.arcslope(east_plane)._griddata[2, 2] == pytest.approx(1.0)


def test_arcslope_has_no_sqrt2_factor():
    """The 8*cellsize denominator carries all the geometry."""
    z = np.add.outer(np.arange(5, 0, -1, dtype=float),
                     np.arange(5, 0, -1, dtype=float))
    grid = dem.Elevation(dx=1.0, grid=z)
    # dz/dx = dz/dy = 1, so the Horn magnitude is sqrt(2).
    assert tt.arcslope(grid)._griddata[2, 2] == pytest.approx(np.sqrt(2.0))


@pytest.mark.parametrize('descent,expected', [
    ('north', 0.0),
    ('east', 90.0),
    ('south', 180.0),
    ('west', 270.0),
])
def test_aspect_points_down_slope(descent, expected):
    """Aspect is the compass bearing the surface descends towards."""
    n = 7
    y, x = np.mgrid[0:n, 0:n].astype(float)
    surface = {
        'north': y,        # elevation grows with row, i.e. downhill to north
        'south': -y,
        'east': -x,
        'west': x,
    }[descent]
    grid = dem.Elevation(dx=1.0, grid=surface)
    assert tt.aspect(grid)._griddata[3, 3] == pytest.approx(expected)


def test_aspect_of_a_flat_cell_is_ninety():
    """atan2(0, 0) is 0, so TopoToolbox reports flats as east-facing."""
    flat = dem.Elevation(dx=1.0, grid=np.zeros((5, 5)))
    assert tt.aspect(flat)._griddata[2, 2] == pytest.approx(90.0)


def test_aspect_classification():
    n = 7
    y, x = np.mgrid[0:n, 0:n].astype(float)
    grid = dem.Elevation(dx=1.0, grid=-x)       # descends east, aspect 90
    classes = tt.aspect(grid, classify=True)._griddata
    assert classes.dtype == np.uint8
    assert classes[3, 3] == 5                   # [90, 135) -> class 5 (E)


def test_aspect_is_independent_of_cell_size():
    n = 7
    y, x = np.mgrid[0:n, 0:n].astype(float)
    fine = dem.Elevation(dx=1.0, grid=-x)
    coarse = dem.Elevation(dx=30.0, grid=-x)
    assert np.allclose(tt.aspect(fine)._griddata, tt.aspect(coarse)._griddata)


def test_curvature_of_a_plane_is_zero(east_plane):
    for ctype in ('profc', 'planc', 'tangc', 'meanc', 'total'):
        curv = tt.curvature(east_plane, ctype)._griddata
        assert np.allclose(curv[1:-1, 1:-1], 0.0, atol=1e-12), ctype


def test_total_curvature_of_a_quadratic():
    """'total' is fxx^2 + 2 fxy^2 + fyy^2, with no minus and no denominator."""
    n = 15
    y, x = np.mgrid[0:n, 0:n].astype(float)
    a = 0.02
    grid = dem.Elevation(dx=1.0, grid=a * x ** 2)
    curv = tt.curvature(grid, 'total')._griddata
    # fxx = 2a, fyy = fxy = 0.
    assert np.allclose(curv[2:-2, 2:-2], (2 * a) ** 2)


def test_curvature_rejects_an_unknown_type(east_plane):
    with pytest.raises(error.InputError):
        tt.curvature(east_plane, 'gaussian')


def test_curvature_zeroes_the_zero_gradient_cells():
    """0/0 becomes 0, not NaN -- matching TopoToolbox."""
    flat = dem.Elevation(dx=1.0, grid=np.zeros((7, 7)))
    assert np.all(tt.curvature(flat, 'profc')._griddata == 0.0)


def test_hillshade_is_a_cosine_not_a_byte(east_plane):
    shade = tt.hillshade(east_plane, azimuth=90, altitude=45)._griddata
    # TopoToolbox scales to 0-255 only when plotting.
    assert shade.dtype == np.float64
    assert -1.0 - 1e-9 <= shade.min() and shade.max() <= 1.0 + 1e-9


def test_hillshade_is_brightest_facing_the_light(east_plane):
    lit = tt.hillshade(east_plane, azimuth=90, altitude=45)._griddata[2, 2]
    shadowed = tt.hillshade(east_plane, azimuth=270, altitude=45)._griddata[2, 2]
    assert lit == pytest.approx(1.0)
    assert shadowed == pytest.approx(0.0, abs=1e-9)


def test_hillshade_defaults_match_topotoolbox(rough_dome):
    """azimuth 315, altitude 60 -- the code's defaults, not the docstring's."""
    default = tt.hillshade(rough_dome)._griddata
    explicit = tt.hillshade(rough_dome, azimuth=315, altitude=60)._griddata
    assert np.array_equal(default, explicit)


def test_hillshade_of_a_flat_surface_is_the_sine_of_the_altitude():
    flat = dem.Elevation(dx=1.0, grid=np.zeros((7, 7)))
    for altitude in (30.0, 45.0, 60.0):
        shade = tt.hillshade(flat, altitude=altitude)._griddata[3, 3]
        assert shade == pytest.approx(np.sin(np.radians(altitude)))


def test_hillshade_agrees_in_sign_with_the_esri_version(rough_dome):
    """The two formulations must at least light the same slopes."""
    tt_shade = tt.hillshade(rough_dome, azimuth=315, altitude=45)._griddata
    esri = dem.Hillshade(elevation=rough_dome, azimuth=315,
                         inclination=45)._griddata.astype(float)
    interior = (slice(2, -2), slice(2, -2))
    correlation = np.corrcoef(tt_shade[interior].ravel(), esri[interior].ravel())[0, 1]
    assert correlation > 0.97


def test_localtopography_range_matches_local_relief(rough_dome):
    tt_relief = tt.localtopography(rough_dome, radius=60.0, type='range')._griddata
    ta_relief = dem.LocalRelief(elevation=rough_dome, pixel_radius=2)._griddata
    assert np.allclose(tt_relief, ta_relief)


# ---------------------------------------------------------------------------
# Cell area
# ---------------------------------------------------------------------------


def test_cellarea_uses_the_topotoolbox_sphere(geographic_grid):
    areas = tt.cellarea(geographic_grid)._griddata
    radius_km = np.sqrt(tt.EARTH_SURFACE_AREA_KM2 / (4 * np.pi))
    dlon = dlat = 0.01
    lat = geographic_grid._georef_info.yllcenter + 9 * dlat
    expected = (radius_km * 1000.0) ** 2 * np.radians(dlon) * (
        np.sin(np.radians(lat + dlat / 2)) - np.sin(np.radians(lat - dlat / 2)))
    assert areas[0, 0] == pytest.approx(expected, rel=1e-9)


def test_cellarea_km_units(geographic_grid):
    metres = tt.cellarea(geographic_grid, 'm')._griddata
    kilometres = tt.cellarea(geographic_grid, 'km')._griddata
    assert np.allclose(metres, kilometres * 1e6)


def test_cellarea_matches_topoanalysis_to_the_radius_difference(geographic_grid):
    """The two differ only through their Earth radii, by ~1.5e-5."""
    tt_area = tt.cellarea(geographic_grid)._griddata
    ta_area = geographic_grid._area_per_pixel()
    ratio = (np.sqrt(tt.EARTH_SURFACE_AREA_KM2 / (4 * np.pi)) * 1000.0
             / geographic_grid.earth_radius) ** 2
    # TopoAnalysis evaluates the longitude difference per cell, so its
    # areas vary in the last few bits along a row; TopoToolbox's do not.
    assert np.allclose(tt_area, ta_area * ratio, rtol=1e-9)


def test_cellarea_rejects_an_unknown_unit(geographic_grid):
    with pytest.raises(error.InputError):
        tt.cellarea(geographic_grid, 'miles')


# ---------------------------------------------------------------------------
# Chi and steepness
# ---------------------------------------------------------------------------


def test_chitransform_is_zero_at_the_outlet(v_valley):
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    chi = tt.chitransform(fd, counts, mn=0.45, a0=1e6)._griddata

    outlets = fd.receivers < 0
    assert np.allclose(chi[outlets], 0.0)
    assert np.nanmax(chi) > 0.0


def test_chitransform_correctcellsize(v_valley):
    """With correctcellsize the area is taken to be in cells, not m^2."""
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    area = dem.Area(flow_direction=fd)      # already m^2

    corrected = tt.chitransform(fd, counts, correctcellsize=True)._griddata
    direct = tt.chitransform(fd, area, correctcellsize=False)._griddata
    assert np.allclose(corrected, direct, equal_nan=True)


def test_chitransform_scales_with_a0(v_valley):
    """chi is proportional to a0**mn."""
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    mn = 0.45
    one = tt.chitransform(fd, counts, mn=mn, a0=1.0)._griddata
    million = tt.chitransform(fd, counts, mn=mn, a0=1e6)._griddata
    valid = np.isfinite(one) & (one > 0)
    assert np.allclose(million[valid] / one[valid], 1e6 ** mn)


def test_chitransform_uses_the_trapezoid_rule(v_valley):
    """It must agree with Chi(trapezoid=True), not with the default rule."""
    fd = tt.flowobj(v_valley)
    area = dem.Area(flow_direction=fd)
    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlets = area._rowscols_to_xy((outlet_rc,))

    tt_chi = tt.chitransform(fd, area, mn=0.45, a0=1e6, outlets=outlets,
                             correctcellsize=False)._griddata
    ta_chi = dem.Chi(area=area, flow_direction=fd, theta=0.45, Ao=1e6,
                     outlets=outlets, trapezoid=True)._griddata

    finite = np.isfinite(tt_chi)
    assert np.allclose(tt_chi[finite], ta_chi[finite])


def test_ksn_is_slope_times_area_to_theta(v_valley):
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    theta = 0.45
    k = tt.ksn(fd, v_valley, counts, theta=theta)._griddata

    carved = tt.imposemin(fd, v_valley, sl=1e-5)
    z = carved._griddata.reshape(-1)
    step = fd.step_lengths().reshape(-1)
    receivers = fd.receivers.reshape(-1)
    area = (counts._griddata * fd._georef_info.dx ** 2).reshape(-1)

    donors = np.flatnonzero(receivers >= 0)
    gradient = (z[donors] - z[receivers[donors]]) / step[donors]
    assert np.allclose(k.reshape(-1)[donors], gradient * area[donors] ** theta)


def test_mchi_is_the_elevation_gradient_in_chi_space(v_valley):
    """M_chi is dz/dchi, a plain forward difference with no carving."""
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    chi = tt.chitransform(fd, counts, mn=0.45, a0=1.0)
    mc = tt.mchi(fd, v_valley, chi)._griddata.reshape(-1)

    z = np.asarray(v_valley._griddata, dtype=float).reshape(-1)
    c = chi._griddata.reshape(-1)
    receivers = fd.receivers.reshape(-1)
    donors = np.flatnonzero(receivers >= 0)
    expected = (z[donors] - z[receivers[donors]]) / (c[donors] - c[receivers[donors]])

    finite = np.isfinite(expected)
    assert np.allclose(mc[donors][finite], expected[finite])


def test_mchi_scales_inversely_with_a0(v_valley):
    """M_chi depends on the a0 baked into chi: mchi(a0) * a0**mn is invariant."""
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    mn = 0.45

    one = tt.mchi(fd, v_valley, tt.chitransform(fd, counts, mn=mn, a0=1.0))._griddata
    million = tt.mchi(fd, v_valley, tt.chitransform(fd, counts, mn=mn, a0=1e6))._griddata

    finite = np.isfinite(one) & np.isfinite(million) & (np.abs(one) > 1e-12)
    assert np.allclose(million[finite] * (1e6 ** mn), one[finite], rtol=1e-9)


def test_mchi_and_ksn_are_close_but_not_identical(v_valley):
    """TopoToolbox's docstring claims they are the same; the code differs.

    ``ksn`` carves the DEM first and takes the integrand at the upstream
    cell, while ``mchi`` uses the raw DEM and chi's trapezoidal averaging.
    """
    fd = tt.flowobj(v_valley)
    counts = tt.flowacc(fd)
    mn = 0.45

    carved = tt.imposemin(fd, v_valley, sl=1e-5)
    mc = tt.mchi(fd, carved, tt.chitransform(fd, counts, mn=mn, a0=1.0))._griddata
    k = tt.ksn(fd, v_valley, counts, theta=mn)._griddata

    receivers = fd.receivers.reshape(-1)
    donors = np.flatnonzero(receivers >= 0)
    a = mc.reshape(-1)[donors]
    b = k.reshape(-1)[donors]
    finite = np.isfinite(a) & np.isfinite(b) & (b > 0)

    ratio = a[finite] / b[finite]
    assert np.all(ratio > 0.5) and np.all(ratio < 2.0)
    assert not np.allclose(a[finite], b[finite])


# ---------------------------------------------------------------------------
# Referencing
# ---------------------------------------------------------------------------


def test_getcoordinates_runs_north_to_south():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((4, 6)))
    x, y = tt.getcoordinates(grid)
    assert x.size == 6 and y.size == 4
    assert y[0] > y[-1], 'row 0 is the north edge'
    assert x[1] - x[0] == pytest.approx(10.0)


def test_getextent_is_cell_centres_not_edges():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((4, 6)))
    centres = tt.getextent(grid)
    edges = grid.extent()
    assert centres[0] == pytest.approx(edges[0] + 5.0)
    assert centres[1] == pytest.approx(edges[1] - 5.0)


def test_coord2sub_and_sub2coord_round_trip():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((7, 5)))
    for row in range(7):
        for col in range(5):
            x, y = tt.sub2coord(grid, row, col)
            r, c = tt.coord2sub(grid, x, y)
            assert (int(r[0]), int(c[0])) == (row, col)


def test_coord2sub_returns_nan_outside_the_grid():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((3, 3)))
    r, c = tt.coord2sub(grid, [1000.0], [1000.0])
    assert np.isnan(r[0]) and np.isnan(c[0])


def test_coord2sub_rounds_halves_away_from_zero():
    """MATLAB's round, not Python's round-half-to-even."""
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((5, 5)))
    info = grid._georef_info
    # Exactly half a cell east of column 0.
    r, c = tt.coord2sub(grid, [info.xllcenter + 0.5], [info.yllcenter])
    assert int(c[0]) == 1


# ---------------------------------------------------------------------------
# Checks against literal ports of the MATLAB source
#
# These transcribe the .m files directly rather than restating the intent, so
# they catch a drift in the Python that a behavioural test would not.
# ---------------------------------------------------------------------------


def matlab_surfnorm(z):
    """A literal port of MATLAB's ``surfnorm(Z)`` (graph3d/surfnorm.m).

    The boundary is expanded by the quadratic ghost cell ``3*a1 - 3*a2 + a3``
    -- *not* the linear ``2*a1 - a2`` -- which makes the edge derivative the
    second-order one-sided difference.
    """
    m, n = z.shape
    x, y = np.meshgrid(np.arange(1.0, n + 1), np.arange(1.0, m + 1))

    def expand(a):
        a = np.vstack([3 * a[0:1] - 3 * a[1:2] + a[2:3], a,
                       3 * a[m - 1:m] - 3 * a[m - 2:m - 1] + a[m - 3:m - 2]])
        return np.hstack([3 * a[:, 0:1] - 3 * a[:, 1:2] + a[:, 2:3], a,
                          3 * a[:, n - 1:n] - 3 * a[:, n - 2:n - 1] + a[:, n - 3:n - 2]])

    xx, yy, zz = expand(x), expand(y), expand(z)

    def filter_row(a):          # filter2([1 0 -1]/2, a) -- correlation
        out = np.zeros_like(a)
        out[:, 1:-1] = 0.5 * (a[:, :-2] - a[:, 2:])
        return out

    def filter_col(a):          # filter2([-1;0;1]/2, a)
        out = np.zeros_like(a)
        out[1:-1, :] = 0.5 * (a[2:, :] - a[:-2, :])
        return out

    rows, cols = slice(1, m + 1), slice(1, n + 1)
    ax, ay, az = (filter_row(v)[rows, cols] for v in (xx, yy, zz))
    bx, by, bz = (filter_col(v)[rows, cols] for v in (xx, yy, zz))

    nx = -(ay * bz - az * by)
    ny = -(az * bx - ax * bz)
    nz = -(ax * by - ay * bx)
    magnitude = np.sqrt(nx * nx + ny * ny + nz * nz)
    magnitude[magnitude == 0] = np.finfo(float).eps
    return nx / magnitude, ny / magnitude, nz / magnitude


@pytest.mark.parity
@pytest.mark.parametrize('shape', [(7, 8), (12, 12), (5, 20)])
def test_surfnorm_matches_matlab_including_the_border(shape):
    from TopoAnalysis.topotoolbox import _surfnorm

    z = np.random.default_rng(7).random(shape) * 100
    for ours, theirs in zip(_surfnorm(z), matlab_surfnorm(z)):
        assert np.allclose(ours, theirs, atol=1e-12)


def matlab_identifyflats(z):
    """A literal port of ``@GRIDobj/identifyflats.m``."""
    from scipy.ndimage import binary_dilation, grey_dilation, grey_erosion

    log_nans = np.isnan(z)
    dem_ = np.where(log_nans, -np.inf, z)

    flats = grey_erosion(dem_, size=(3, 3), mode='constant', cval=np.inf) == dem_
    if log_nans.any():
        flats &= ~log_nans
    flats[:, [0, -1]] = False
    flats[[0, -1], :] = False
    if log_nans.any():
        flats[binary_dilation(log_nans, structure=np.ones((3, 3), bool))] = False

    marker = np.full(dem_.shape, -np.inf)
    marker[flats] = dem_[flats]
    sills = (grey_dilation(marker, size=(3, 3), mode='constant',
                           cval=-np.inf) == dem_) & ~flats
    if log_nans.any():
        sills[log_nans] = False
    return flats, sills


@pytest.mark.parity
@pytest.mark.parametrize('seed', [0, 1, 2, 3])
def test_identifyflats_matches_matlab(seed):
    """Including the cleared rim, the NaN halo, and isolated pits."""
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(seed)
    n = int(rng.integers(20, 45))
    # Rounded, so plateaus and ties actually occur.
    z = np.round(gaussian_filter(rng.random((n, n)) * 50, 2.0), 1)
    if seed % 2:
        z[rng.random(z.shape) < 0.05] = np.nan

    filled = tt.fillsinks(dem.Elevation(dx=1.0, grid=z))
    expected_flats, expected_sills = matlab_identifyflats(filled._griddata)
    flats, sills = tt.identifyflats(filled)

    assert np.array_equal(flats._griddata, expected_flats)
    assert np.array_equal(sills._griddata, expected_sills)


@pytest.mark.parity
def test_identifyflats_clears_the_rim():
    """A plateau reaching the grid edge is flat only in the interior."""
    z = np.full((6, 6), 3.0)
    z[5, 5] = 2.0                       # the plateau's only outlet, at a corner
    flats = tt.identifyflats(dem.Elevation(dx=1.0, grid=z), output='flats')._griddata

    # The rim is cleared unconditionally ...
    assert not flats[0, :].any() and not flats[-1, :].any()
    assert not flats[:, 0].any() and not flats[:, -1].any()
    # ... and interior cells away from the outlet are flat.
    assert flats[1, 1] and flats[2, 3]
    # (4, 4) touches the 2, so it has a lower neighbour and is not flat.
    assert not flats[4, 4]


@pytest.mark.parity
def test_identifyflats_clears_the_nodata_halo():
    """MATLAB drops the 3x3 dilation of the no-data mask from the flats."""
    z = np.full((7, 7), 3.0)
    z[3, 3] = np.nan
    flats = tt.identifyflats(dem.Elevation(dx=1.0, grid=z), output='flats')._griddata

    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            assert not flats[3 + di, 3 + dj], 'the NaN halo must be cleared'
    assert flats[1, 1], 'the rest of the plateau is still flat'


@pytest.mark.parity
def test_identifyflats_counts_an_isolated_pit_as_a_flat():
    """MATLAB's test is 'no strictly lower neighbour', with no equal-neighbour rule."""
    z = np.full((5, 5), 99.0)
    z[1:4, 1:4] = 5.0
    z[2, 2] = 1.0                      # lower than all eight neighbours
    flats = tt.identifyflats(dem.Elevation(dx=1.0, grid=z), output='flats')._griddata
    assert flats[2, 2]


@pytest.mark.parity
def test_closed_basins_exclude_the_grid_edge():
    z = np.full((7, 7), 50.0)
    z[3, 3] = 10.0                      # an interior closed basin
    z[0, 3] = 1.0                       # touches the edge, so not closed
    closed = tt.identifyflats(dem.Elevation(dx=1.0, grid=z),
                              output='closedbasins')._griddata
    assert closed[3, 3]
    assert not closed[0, 3]


@pytest.mark.parity
def test_fillsinks_maxdepth_is_all_or_nothing_per_depression():
    """MATLAB restores the whole depression, not just its deepest cells."""
    z = np.full((7, 7), 10.0)
    z[1:6, 1:6] = 9.0
    z[2:5, 2:5] = 6.0
    z[3, 3] = 4.0                       # one depression, spill 10, depth 6
    grid = dem.Elevation(dx=1.0, grid=z)

    # Below the depression's depth: nothing is touched at all.
    for maxdepth in (1.0, 3.0, 5.0):
        assert np.array_equal(tt.fillsinks(grid, maxdepth=maxdepth)._griddata, z), \
            'maxdepth={0} must leave the whole depression alone'.format(maxdepth)

    # At or above it: the depression fills completely.
    for maxdepth in (6.0, 7.0):
        filled = tt.fillsinks(grid, maxdepth=maxdepth)._griddata
        assert np.all(filled[1:6, 1:6] == 10.0)


@pytest.mark.parity
def test_localtopography_mean_uses_the_disc_not_its_bounding_box():
    data = np.zeros((15, 15))
    data[4, 4] = 100.0                  # 4.24 cells from the centre
    grid = dem.Elevation(dx=1.0, grid=data)
    # A radius-3 disc excludes it; the 7x7 bounding box does not.
    assert tt.localtopography(grid, radius=3.0, type='mean')._griddata[7, 7] == \
        pytest.approx(0.0)
    assert tt.localtopography(grid, radius=6.0, type='mean')._griddata[7, 7] > 0.0


@pytest.mark.parity
def test_localtopography_std_is_the_sample_standard_deviation():
    """MATLAB's stdfilt normalises by N-1."""
    values = np.arange(25.0).reshape(5, 5)
    grid = dem.Elevation(dx=1.0, grid=values)
    window = [7.0, 11.0, 12.0, 13.0, 17.0]      # the radius-1 disc at (2, 2)
    assert tt.localtopography(grid, radius=1.0, type='std')._griddata[2, 2] == \
        pytest.approx(np.std(window, ddof=1))
