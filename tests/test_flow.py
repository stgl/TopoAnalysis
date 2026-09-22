"""Tests for the flow-routing grid classes."""

from __future__ import annotations

import numpy as np
import pytest

from TopoAnalysis import dem


# ---------------------------------------------------------------------------
# Filling
# ---------------------------------------------------------------------------


def test_filled_elevation_removes_the_pit(single_pit):
    filled = dem.FilledElevation(elevation=single_pit)
    assert np.all(filled._griddata >= single_pit._griddata - 1e-12)
    # With the default aggradation slope the fill is just above the rim.
    assert filled._griddata[2, 2] > 10.0
    assert filled._griddata[2, 2] == pytest.approx(10.0, abs=1e-6)
    assert filled.fill_report['cells_filled'] == 9


def test_flat_fill_option_matches_the_rim_exactly(single_pit):
    filled = dem.FilledElevation(elevation=single_pit, aggradation_slope=0.0)
    assert np.all(filled._griddata == 10.0)


def test_fill_restricted_to_a_mask(single_pit):
    """Only masked cells take part; the rest keep their exact values."""
    mask = dem.Mask(dx=1.0, grid=np.zeros((5, 5), dtype=np.uint8))
    mask._griddata[1:4, 1:4] = 1
    masked = dem.FilledElevation(elevation=single_pit, mask=mask)
    unmasked = dem.FilledElevation(elevation=single_pit)

    outside = mask._griddata == 0
    assert np.array_equal(masked._griddata[outside], single_pit._griddata[outside])
    # And the mask must actually change the answer.  Confined to the 3x3
    # interior, the flood spills over that block's own rim of 4 rather than
    # over the grid's 10 m rim.
    assert not np.allclose(masked._griddata, unmasked._griddata)
    assert masked._griddata[2, 2] == pytest.approx(4.0, abs=1e-6)
    assert unmasked._griddata[2, 2] == pytest.approx(10.0, abs=1e-6)


def test_fill_from_explicit_outlets(single_pit):
    outlet = single_pit._rowscols_to_xy(((2, 0),))
    filled = dem.FilledElevation(elevation=single_pit, outlets=outlet)
    assert np.all(filled._griddata >= single_pit._griddata - 1e-12)


def test_clip_to_fill_marks_unreached_cells():
    """Cells the flood cannot reach become no-data."""
    # Two basins separated by a wall of NaN: flooding from the left outlet
    # cannot reach the right half.
    z = np.array([
        [5., 4., 3., np.nan, 3., 4., 5.],
        [5., 4., 3., np.nan, 3., 4., 5.],
        [5., 4., 3., np.nan, 3., 4., 5.]])
    grid = dem.Elevation(dx=1.0, grid=z)
    outlet = grid._rowscols_to_xy(((1, 2),))

    filled = dem.FilledElevation(elevation=grid, outlets=outlet, clip_to_fill=True)
    assert np.isfinite(filled._griddata[:, :3]).all(), 'the left basin is reached'
    assert np.isnan(filled._griddata[:, 4:]).all(), 'the right basin is not'


def test_clip_to_fill_excises_depressions_left_by_a_depth_limit():
    """A depression too deep to fill is marked no-data, as it always was."""
    ny, nx = 11, 21
    z = np.tile(np.linspace(100.0, 80.0, nx), (ny, 1))
    z[4:7, 9:12] = 50.0                      # a 40 m pit
    grid = dem.Elevation(dx=10.0, grid=z)

    clipped = dem.FilledElevation(elevation=grid, maximum_pit_depth=5.0,
                                  clip_to_fill=True)
    assert int(np.isnan(clipped._griddata).sum()) == 9
    assert np.isnan(clipped._griddata[4:7, 9:12]).all()

    # Without the limit the pit is filled and nothing is clipped.
    filled = dem.FilledElevation(elevation=grid, clip_to_fill=True)
    assert np.isfinite(filled._griddata).all()
    assert filled._griddata[5, 10] > 80.0


def test_findDEMedge_matches_the_flood_seeds(rough_dome):
    rows, cols = rough_dome.findDEMedge()
    n = rough_dome._georef_info.nx
    assert set((rows * n + cols).tolist()) == set(
        dem.kernels.default_seeds(rough_dome._griddata).tolist())


def test_randomize_subbasins_with_mask_runs(single_pit):
    """The method used to call a name that Python had mangled away."""
    filled = dem.FilledElevation(elevation=single_pit)
    mask = dem.Mask(dx=1.0, grid=np.ones((5, 5), dtype=np.uint8))
    outlets = single_pit._rowscols_to_xy(((0, 0),))
    filled.randomize_subbasins_with_mask(mask, outlets)
    assert filled._griddata.shape == (5, 5)


# ---------------------------------------------------------------------------
# Flow direction
# ---------------------------------------------------------------------------


def test_flow_direction_on_a_plane(east_plane):
    filled = dem.FilledElevation(elevation=east_plane)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    # Everything flows east except the last column, which leaves the grid.
    assert np.all(d8._griddata[:, :-1] == 1)
    assert np.all(d8._griddata[:, -1] == 0)


def test_flow_direction_covers_the_last_row_and_column(rough_dome):
    """The old windowed implementation left the last row and column unrouted.

    On a dome every cell drains outwards, so the only cells without a
    direction are those on the perimeter whose steepest neighbour is off the
    grid -- never a whole row or column.
    """
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)

    assert (d8._griddata[1:-1, 1:-1] == 0).sum() == 0
    # The last row and the last column are what the old code skipped
    # entirely; each must now be routed except where flow leaves the grid.
    for label, line in (('last row', d8._griddata[-1, :]),
                        ('last column', d8._griddata[:, -1]),
                        ('first row', d8._griddata[0, :]),
                        ('first column', d8._griddata[:, 0])):
        assert (line != 0).any(), '{0} is entirely unrouted'.format(label)


def test_flow_direction_from_an_unfilled_elevation(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    assert np.all(d8._griddata[:, :-1] == 1)


def test_get_flow_to_cell_round_trips(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    receivers = d8.receivers
    nx = d8._georef_info.nx

    for i in (5, 17, 30):
        for j in (4, 22, 38):
            row, col, good = d8.get_flow_to_cell(i, j)
            if good:
                assert receivers[i, j] == row * nx + col
            else:
                assert receivers[i, j] == -1


def test_upstream_cells_are_the_inverse_of_downstream(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)

    for (i, j) in ((10, 10), (20, 25), (33, 7)):
        for (ui, uj) in d8.get_upstream_cell_indexes(i, j):
            assert d8.get_flow_to_cell(ui, uj)[:2] == (i, j)


def test_upstream_cell_set_is_closed(rough_dome):
    """Everything upstream of a cell must itself drain to that cell."""
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)

    upstream = set(d8.get_indexes_of_upstream_cells(20, 20))
    assert (20, 20) in upstream
    for (i, j) in upstream:
        if (i, j) == (20, 20):
            continue
        row, col, good = d8.get_flow_to_cell(i, j)
        assert good and (row, col) in upstream


def test_get_indexes_of_upstream_cells_can_be_iterated_twice(rough_dome):
    """The result used to be a single-use zip object."""
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    indexes = d8.get_indexes_of_upstream_cells(20, 20)
    assert len(list(indexes)) == len(list(indexes)) > 0


def test_deep_basins_do_not_need_recursion():
    """A 4000-cell snake would have overflowed the recursive traversal."""
    n = 4000
    z = np.arange(n, 0, -1, dtype=float).reshape(1, n)
    grid = dem.Elevation(dx=1.0, grid=z)
    d8 = dem.FlowDirectionD8(elevation=grid)
    assert len(d8.get_indexes_of_upstream_cells(0, n - 1)) == n


def test_pixel_scale_uses_the_true_diagonal():
    codes = np.array([[1, 2], [4, 32]], dtype=np.uint8)
    d8 = dem.FlowDirectionD8(dx=1.0, grid=codes)
    scale = d8.pixel_scale(np.float64)
    assert scale[0, 0] == 1.0 and scale[1, 0] == 1.0
    assert scale[0, 1] == pytest.approx(np.sqrt(2.0))
    assert scale[1, 1] == pytest.approx(np.sqrt(2.0))


def test_rivertools_conversion_rotates_the_codes():
    # RiverTools starts its numbering one position anticlockwise of ArcGIS.
    codes = np.array([[1, 2, 4, 8], [16, 32, 64, 128]], dtype=np.uint8)
    d8 = dem.FlowDirectionD8(dx=1.0, grid=codes.copy())
    d8.convert_rivertools_directions_to_arc()
    assert d8._griddata.tolist() == [[128, 1, 2, 4], [8, 16, 32, 64]]


def test_search_down_flow_direction_follows_the_plane(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    start = east_plane._rowscols_to_xy(((2, 0),))[0]
    cells, lengths = d8.search_down_flow_direction_with_length(start)
    assert [c[1] for c in cells] == [0, 1, 2, 3, 4]
    assert list(lengths) == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0])


def test_divides_have_no_upstream_neighbours(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    divides = d8.divides()

    rows, cols = np.where(divides._griddata == 1)
    for (i, j) in list(zip(rows.tolist(), cols.tolist()))[:50]:
        assert d8.get_upstream_cell_indexes(i, j) == []


def test_bounds_of_basin_for_outlet(rough_dome):
    """Comparing a float with None used to raise on the first iteration."""
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    ((x, y),) = d8._rowscols_to_xy(((20, 20),))

    (xr, yr) = d8.bounds_of_basin_for_outlet((y, x))
    assert xr[0] <= x <= xr[1]
    assert yr[0] <= y <= yr[1]


def test_divides_for_outlets_returns_reusable_sequences(v_valley):
    """Two adjoining basins must share a non-empty divide, listed twice over."""
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)

    # The valley's outlet, and a point part-way up it: the second basin is
    # nested inside the first, so the two adjoin along its whole rim.
    a = d8._rowscols_to_xy(((8, 3),))[0]
    b = d8._rowscols_to_xy(((0, 3),))[0]

    first, second = d8.divides_for_outlets(a, b)
    assert len(first) > 0 and len(second) > 0, 'the basins must actually meet'
    # Both must survive a second pass; they used to be single-use zips.
    assert len(list(first)) == len(list(first)) == len(first)
    assert len(list(second)) == len(list(second)) == len(second)
    # Every reported cell belongs to the basin it was reported for.
    in_a = d8.basin_mask((a,))
    in_b = d8.basin_mask((b,))
    assert all(in_a[i, j] for (i, j) in first)
    assert all(in_b[i, j] for (i, j) in second)


# ---------------------------------------------------------------------------
# Accumulation and length
# ---------------------------------------------------------------------------


def test_area_conserves_the_total(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)

    cell_area = rough_dome._georef_info.dx ** 2
    assert area._griddata.min() == pytest.approx(cell_area)

    outlets = d8.receivers < 0
    assert area._griddata[outlets].sum() == pytest.approx(
        cell_area * rough_dome._griddata.size)


def test_area_on_a_plane_grows_downstream(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    area = dem.Area(flow_direction=d8)
    # Each row is an independent flow line: 1, 2, 3, 4, 5 cells.
    assert area._griddata[2].tolist() == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_area_ignores_a_loaded_grids_missing_sort_order(rough_dome, tmp_path):
    """Accumulation must not depend on an elevation sort being cached."""
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    expected = dem.Area(flow_direction=d8)._griddata

    # A flow-direction grid built straight from codes has no sort cache, so
    # sort() would have ordered cells by their direction *code*.
    fresh = dem.FlowDirectionD8(dx=rough_dome._georef_info.dx,
                                grid=d8._griddata.copy())
    assert np.allclose(dem.Area(flow_direction=fresh)._griddata, expected)


def test_area_gate_restricts_accumulation(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)

    gate = dem.Mask(dx=rough_dome._georef_info.dx,
                    grid=np.zeros(rough_dome._griddata.shape, dtype=np.uint8))
    gate._griddata[10:30, 10:30] = 1

    gated = dem.Area(flow_direction=d8, evaluate_at=gate)
    plain = dem.Area(flow_direction=d8)

    # A mask that used to be ignored entirely must now reduce the totals ...
    assert gated._griddata.max() < plain._griddata.max()
    # ... and no cell can accumulate more than the gate could deliver.
    cell = rough_dome._georef_info.dx ** 2
    assert gated._griddata.max() <= cell * int(gate._griddata.sum())
    # Cells outside the gate contribute nothing of their own, so any value
    # they hold came from inside it.
    assert np.all(gated._griddata <= plain._griddata + 1e-9)


def test_areas_between_finds_channel_heads(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)

    cell = rough_dome._georef_info.dx ** 2
    heads = area.areas_between(d8, 5 * cell, 20 * cell)
    assert len(heads) > 0
    for (x, y) in heads:
        ((i, j),) = area._xy_to_rowscols(((x, y),))
        assert 5 * cell <= area._griddata[i, j] <= 20 * cell


def test_flow_length_on_a_plane(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    length = dem.FlowLength(flow_direction=d8)
    assert length._griddata[2].tolist() == [0.0, 1.0, 2.0, 3.0, 4.0]


def test_flow_length_main_stem_uses_the_shared_code_convention(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    length = dem.FlowLength(flow_direction=d8)
    # Cell (2,2)'s longest path comes from the west, so the code stored
    # there must be 16 (W) in the same convention FlowDirectionD8 uses.
    assert length.main_stem_directions[2, 2] == 16


def test_flow_length_path_walks_upstream(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    length = dem.FlowLength(flow_direction=d8)
    outlet = east_plane._rowscols_to_xy(((2, 4),))[0]
    path = length.indexes_along_flow_path_from_outlet(outlet)
    assert [c[1] for c in path] == [4, 3, 2, 1, 0]


def test_is_along_flow_length_agrees_with_the_stored_codes(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    length = dem.FlowLength(flow_direction=d8)

    hits = 0
    for i in range(5, 35, 7):
        for j in range(5, 35, 7):
            row, col, good = d8.get_flow_to_cell(i, j)
            if good and length.is_along_flow_length((i, j), (row, col)):
                hits += 1
    assert hits > 0


@pytest.mark.gdal
def test_flow_length_round_trip(rough_dome, tmp_path):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    length = dem.FlowLength(flow_direction=d8)

    path = str(tmp_path / 'length.tif')
    length.save(path)
    back = dem.FlowLength.load(path)

    assert np.allclose(back._griddata, length._griddata)
    assert np.array_equal(back.main_stem_directions, length.main_stem_directions)


# ---------------------------------------------------------------------------
# Chi, relief and steepness
# ---------------------------------------------------------------------------


def test_chi_is_zero_at_the_outlet(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)

    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlet = area._rowscols_to_xy((outlet_rc,))
    chi = dem.Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e4, outlets=outlet)

    assert chi._griddata[outlet_rc] == 0.0
    assert chi._griddata.max() > 0.0


def test_chi_trapezoid_option(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlet = area._rowscols_to_xy((outlet_rc,))

    left = dem.Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e4, outlets=outlet)
    trap = dem.Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e4, outlets=outlet,
                   trapezoid=True)
    assert not np.allclose(left._griddata, trap._griddata)
    assert trap._griddata.max() <= left._griddata.max() + 1e-9


def test_ksi_uses_the_ratio_of_areas(v_valley):
    """chi must integrate (Ao/A)^theta, not (Ao/(A-Ao))^theta."""
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    length = dem.FlowLength(flow_direction=d8)

    Ao, theta = 200.0, 0.45
    ksi = dem.Ksi(area=area, flow_direction=d8, flow_length=length,
                  theta=theta, Ao=Ao)

    # Pin the integrand exactly.  At a cell where the channel network
    # *starts* -- nothing above it on the main stem is on the network -- chi
    # is one step of (Ao/A)^theta * step_length and nothing more.  The old
    # form, (Ao/(A - Ao))^theta, is a different number and diverges as A
    # approaches Ao.
    step = d8.step_lengths()
    a = area._griddata
    on_network = a > Ao

    sources = on_network.copy()
    codes = length.main_stem_directions
    ny, nx = a.shape
    for i in range(ny):
        for j in range(nx):
            code = codes[i, j]
            if not code:
                continue
            for d in range(8):
                if code == dem.kernels.D8_CODES[d]:
                    ui = i + int(dem.kernels.D8_DI[d])
                    uj = j + int(dem.kernels.D8_DJ[d])
                    if 0 <= ui < ny and 0 <= uj < nx and on_network[ui, uj]:
                        sources[i, j] = False
                    break
    assert sources.any(), 'the test needs at least one network source'

    expected = (Ao / a[sources]) ** theta * step[sources]
    assert np.allclose(ksi._griddata[sources], expected)

    buggy = (Ao / (a[sources] - Ao)) ** theta * step[sources]
    assert not np.allclose(ksi._griddata[sources], buggy)
    assert np.isfinite(ksi._griddata[on_network]).all()


def test_relief_measures_the_drop_to_the_path_head(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    length = dem.FlowLength(flow_direction=d8)

    relief = dem.Relief(flow_direction=d8, elevation=v_valley, flow_length=length)
    assert relief._griddata.min() >= -1e-9
    # The outlet of the longest path sees the largest relief.
    assert relief._griddata.max() > 0.0


def test_scaled_relief_accepts_a_flooded_dem(v_valley):
    """The flooded_dem branch called a method that did not exist."""
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    length = dem.FlowLength(flow_direction=d8)

    scaled = dem.ScaledRelief(flow_direction=d8, flooded_dem=filled,
                              elevation=v_valley, flow_length=length,
                              Ao=1e4, theta=0.45)
    assert np.isfinite(scaled._griddata).all()


def test_chi_scaled_relief(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlet = area._rowscols_to_xy((outlet_rc,))

    relief = dem.ChiScaledRelief(elevation=v_valley, flow_direction=d8,
                                 theta=0.45, Ao=1e4, outlets=outlet)
    assert relief._griddata[outlet_rc] == pytest.approx(0.0)
    assert relief._griddata.max() > 0.0


def test_cross_divide_dchi_runs(rough_dome):
    """zip() objects are not subscriptable; this used to raise immediately."""
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    outlets = area.areas_greater_than(0.5 * area._griddata.max())
    chi = dem.Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e4,
                  outlets=outlets[:1])

    result = dem.CrossDivideDChi(chi=chi, flow_direction=d8)
    assert result._griddata.shape == rough_dome._griddata.shape
    normalized = dem.NormalizedCrossDivideDChi(chi=chi, flow_direction=d8)
    assert np.isfinite(normalized._griddata).all()


# ---------------------------------------------------------------------------
# Masks and basins
# ---------------------------------------------------------------------------


def test_mask_matches_the_upstream_set(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    outlet = d8._rowscols_to_xy(((25, 25),))

    mask = dem.Mask(flow_direction=d8, outlets=outlet)
    expected = set(d8.get_indexes_of_upstream_cells(25, 25))
    assert int(mask._griddata.sum()) == len(expected)
    for (i, j) in list(expected)[:100]:
        assert mask._griddata[i, j] == 1


def test_channel_slope_is_positive_downhill(east_plane):
    d8 = dem.FlowDirectionD8(elevation=east_plane)
    slope = dem.ChannelSlope(flow_direction=d8, elevation=east_plane)
    assert np.all(slope._griddata[:, :-1] == pytest.approx(1.0))
    assert np.all(slope._griddata[:, -1] == 0.0)


def test_discrete_flow_accumulation_walks_downhill(east_plane):
    # The walk needs a full 8-neighbourhood at every step, so it must start
    # away from the grid edge.
    outlet = east_plane._rowscols_to_xy(((2, 1),))
    dfa = dem.DiscreteFlowAccumulation(elevation=east_plane, outlets=outlet)
    # The walk descends east, ties going to the first neighbour in the
    # search order, so the path is diagonal; all that matters is that it
    # advanced and that area accumulated along it.
    visited = dfa._griddata > 0
    assert visited.sum() > 1
    assert dfa._griddata.max() > dfa._griddata[2, 1]


def test_restored_elevation_converges(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlet = area._rowscols_to_xy((outlet_rc,))

    restored = dem.RestoredElevation(flow_direction=d8, elevation=v_valley,
                                     area=area, theta=0.45, ks=20.0,
                                     outlets=outlet, iterations=2, verbose=False)
    assert len(restored.convergence) == 2
    assert np.isfinite(restored._griddata).all()


def test_deflection_of_a_uniform_load():
    """A flat load must deflect uniformly, at the Airy value."""
    load = dem.Elevation(dx=1000.0, grid=np.full((32, 32), 1000.0))
    rho_c, rho_m, g = 2700.0, 3300.0, 9.81

    w = dem.Deflection(elevation=load, D=1e23, rho_m=rho_m, rho_c=rho_c, g=g)
    airy = -1000.0 * rho_c / (rho_m - rho_c)
    assert np.allclose(w._griddata, airy, rtol=1e-6)
