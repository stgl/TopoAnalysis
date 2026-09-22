"""Tests for the supporting modules: windows, datasets, profiles, the CLI."""

from __future__ import annotations

import os

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')  # no display in CI

from TopoAnalysis import (MovingWindow, analysis, datasets, dem,
                          demRecursionTools, error, plotting, utils)
from TopoAnalysis import cli


# ---------------------------------------------------------------------------
# MovingWindow
# ---------------------------------------------------------------------------


class WindowMean(MovingWindow.CircularMovingWindow):
    function = staticmethod(np.mean)


class WindowRange(MovingWindow.RectangularMovingWindow):
    function = staticmethod(lambda values: values.max() - values.min())


def test_moving_window_needs_a_radius():
    with pytest.raises(error.InputError):
        WindowMean()


def test_abstract_windows_cannot_be_instantiated():
    with pytest.raises(error.Error):
        MovingWindow.MovingWindow(window_radius=1.0)
    with pytest.raises(error.Error):
        MovingWindow.RectangularMovingWindow(window_radius=1.0)
    with pytest.raises(error.Error):
        MovingWindow.CircularMovingWindow(window_radius=1.0)


def test_circular_window_is_a_disc_not_a_band():
    """The kernel used to square the row offset twice, giving a column."""
    window = WindowMean(window_radius=2.0)
    rows, cols = window._build_search_kernel(1.0)
    offsets = set(zip(rows.tolist(), cols.tolist()))
    assert (0, 2) in offsets and (2, 0) in offsets
    assert (2, 2) not in offsets, 'the corner is outside a radius-2 disc'


def test_moving_window_mean_of_a_constant_field():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.full((9, 9), 4.0))
    # The method used to call a mangled name that did not exist, and to pass
    # the grid where the cell size was expected.
    result = grid.apply_moving_window(WindowMean(window_radius=2.0))
    assert np.allclose(result._griddata, 4.0)


def test_moving_window_ignores_nodata():
    data = np.ones((7, 7))
    data[3, 3] = np.nan
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)
    result = grid.apply_moving_window(WindowMean(window_radius=1.0))
    # The NaN is dropped from the window rather than poisoning it.
    assert result._griddata[3, 3] == pytest.approx(1.0)


def test_rectangular_window_range():
    data = np.tile(np.arange(9, dtype=float), (9, 1))
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)
    result = grid.apply_moving_window(WindowRange(window_radius=1.0))
    assert result._griddata[4, 4] == pytest.approx(2.0)


def test_window_dimension_is_still_accepted():
    """The original keyword name must keep working."""
    window = WindowMean(window_dimension=3.0)
    assert window.window_radius == 3.0


# ---------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------


def test_triangle_grid_has_the_requested_shape_and_relief():
    grid = datasets.triangle_grid(30, 40, width=10, amp=5.0)
    assert grid._griddata.shape == (30, 40)
    assert grid._griddata.min() == pytest.approx(0.0)
    assert grid._griddata.max() == pytest.approx(10.0, rel=0.1)


def test_sinusoid_grid_spans_the_requested_amplitude():
    grid = datasets.sinusoid_grid(20, 40, width=10, amp=3.0)
    row = grid._griddata[0]
    assert np.allclose(grid._griddata, row, atol=1e-12), 'every row is the same'
    assert row.min() == pytest.approx(0.0)
    assert row.max() == pytest.approx(6.0, rel=0.05)
    # Four ridges across 40 columns at a spacing of 10.
    peaks = np.flatnonzero((row[1:-1] > row[:-2]) & (row[1:-1] > row[2:])) + 1
    assert len(peaks) == 4


def test_slope_adds_a_ramp_rather_than_scaling_the_relief():
    """`triangle *= tilt` zeroed row 0 and amplified the rest."""
    flat = datasets.triangle_grid(20, 20, width=5, amp=1.0)
    sloped = datasets.triangle_grid(20, 20, width=5, amp=1.0, slope_y=2.0)
    assert not np.allclose(sloped._griddata[0], 0.0)
    # The added ramp is the only difference between the two.
    difference = sloped._griddata - flat._griddata
    assert np.allclose(difference, difference[:, [0]])
    assert difference[1, 0] - difference[0, 0] == pytest.approx(
        2.0 * 20.0 / 19.0, rel=1e-6)


def test_synthetic_grids_are_georeferenced_enough_to_route():
    grid = datasets.triangle_grid(24, 24, width=6, amp=4.0, sig=0.05)
    filled = dem.FilledElevation(elevation=grid)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    assert dem.Area(flow_direction=d8)._griddata.max() > 1.0
    # A grid built in memory used to have no geotransform, so save() failed.
    assert grid._georef_info.geoTransform != 0


# ---------------------------------------------------------------------------
# Recursive profile tools
# ---------------------------------------------------------------------------


@pytest.fixture
def profile_network(v_valley):
    filled = dem.FilledElevation(elevation=v_valley)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    outlet_rc = np.unravel_index(int(np.argmax(area._griddata)), area._griddata.shape)
    outlet = area._rowscols_to_xy((outlet_rc,))[0]
    return v_valley, d8, area, outlet


def test_map_values_to_recursive_list_structure(profile_network):
    elevation, d8, area, outlet = profile_network
    tree = d8.map_values_to_recursive_list(outlet, elevation=elevation, area=area)

    assert 'index' in tree and 'distance_scale' in tree
    assert tree['distance_scale'] == 1.0, 'the root has no parent step'
    assert tree['area'] == pytest.approx(area._griddata.max())


def test_distance_scale_describes_the_step_from_the_parent(profile_network):
    """It used to be overwritten on the parent by each child in turn."""
    elevation, d8, area, outlet = profile_network
    tree = d8.map_values_to_recursive_list(outlet, elevation=elevation, area=area)

    stack = [tree]
    checked = 0
    while stack:
        node = stack.pop()
        parent_index = node['index']
        for child in node.get('next', []):
            di = child['index'][0] - parent_index[0]
            dj = child['index'][1] - parent_index[1]
            expected = np.sqrt(2.0) if (di and dj) else 1.0
            assert child['distance_scale'] == pytest.approx(expected)
            checked += 1
            stack.append(child)
    assert checked > 0


def test_recursive_map_covers_the_whole_basin(profile_network):
    elevation, d8, area, outlet = profile_network
    tree = d8.map_values_to_recursive_list(outlet, elevation=elevation, area=area)

    seen = set()
    stack = [tree]
    while stack:
        node = stack.pop()
        seen.add(tuple(node['index']))
        stack.extend(node.get('next', []))

    ((row, col),) = d8._xy_to_rowscols((outlet,))
    assert seen == set(d8.get_indexes_of_upstream_cells(row, col))


def test_chi_elevation_and_ks_fit(profile_network):
    elevation, d8, area, outlet = profile_network
    de = area._mean_pixel_dimension()
    tree = d8.map_values_to_recursive_list(outlet, elevation=elevation, area=area)

    e, c = demRecursionTools.chi_elevation(tree, de, np.array([0.45]), xo=10.0)
    assert e.size == c.size > 0
    assert c.min() == pytest.approx(0.0)

    ks, r2 = demRecursionTools.best_ks_with_r2_list(tree, de, np.array([0.45]), xo=10.0)
    assert np.isfinite(ks).all()


def test_hypsometric_integral_is_between_zero_and_one(profile_network):
    elevation, d8, area, outlet = profile_network
    dA = dem.BaseSpatialGrid(dx=area._georef_info.dx,
                             grid=area._area_per_pixel())
    hi = demRecursionTools.hi(elevation, d8, dA, outlet)
    assert 0.0 <= hi <= 1.0


def test_calc_ks_for_outlet_without_an_xo(profile_network):
    """kwargs.pop('xo') used to raise KeyError when xo was left at default."""
    elevation, d8, area, outlet = profile_network
    ks, r2 = utils.calc_ks_for_outlet(outlet, 0.45, xo=10.0, flow_direction=d8,
                                      elevation=elevation, area=area)
    assert np.isfinite(ks).all()


def test_map_chi_profiles(profile_network):
    elevation, d8, area, outlet = profile_network
    chi_map = demRecursionTools.map_chi_profiles(
        elevation, d8, area, outlet, minimum_area=100.0, theta=0.45)
    assert len(chi_map) > 0
    for (chi, relief) in chi_map.values():
        assert np.isfinite(chi)


# ---------------------------------------------------------------------------
# Plotting (smoke tests -- Agg backend, no display)
# ---------------------------------------------------------------------------


def test_plot_downstream_profile(profile_network):
    elevation, d8, area, outlet = profile_network
    rows, cols, length, profile = plotting.plot_downstream_profile(
        elevation, d8, outlet, 'k-')
    assert rows.size == cols.size == length.size == profile.size > 0


def test_plot_two_grids(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    slope = dem.MaxSlope(elevation=rough_dome)

    # Passing xlabel used to leave the local name unbound and raise.
    axes = dem.plot(area, slope, xlabel='area', ylabel='slope',
                    interactive=False, decimation_factor=7)
    assert axes.get_xlabel() == 'area'


def test_grid_plot_accepts_an_axes(rough_dome):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    returned = rough_dome.plot(ax=ax, interactive=False, colorbar=False)
    assert returned is fig
    plt.close(fig)


# ---------------------------------------------------------------------------
# Quadrats
# ---------------------------------------------------------------------------


def test_quadrats_partition_the_data():
    data = np.arange(100, dtype=float).reshape(10, 10)
    quadrats = analysis.Quadrats(data=data, dx=5, dy=5)
    assert len(quadrats.quadrats) == 4
    assert all(q.shape == (5, 5) for q in quadrats.quadrats)
    assert quadrats.map_quadrats(np.sum) == [
        float(data[r:r + 5, c:c + 5].sum()) for r in (0, 5) for c in (0, 5)]


def test_quadrats_without_rasterio_message():
    quadrats = analysis.Quadrats(data=np.zeros((4, 4)))
    try:
        import rasterio  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match='rasterio'):
            quadrats.load_data('nowhere.tif')


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_cli_info(capsys):
    assert cli.main(['--info']) == 0
    out = capsys.readouterr().out
    assert 'compute backend' in out
    assert 'Barnes' in out


def test_cli_rejects_an_unknown_product():
    with pytest.raises(SystemExit):
        cli.main(['dem.tif', '--products', 'not-a-product'])


def test_cli_requires_a_dem():
    with pytest.raises(SystemExit):
        cli.main([])


@pytest.mark.gdal
def test_cli_end_to_end(tmp_path, rough_dome):
    source = str(tmp_path / 'dem.tif')
    rough_dome.save(source)

    outdir = str(tmp_path / 'out')
    assert cli.main([source, '-o', outdir, '--products',
                     'filled,flowdir,area,length,slope,hillshade',
                     '--min-area', '9000']) == 0

    produced = sorted(os.listdir(outdir))
    for name in ('dem_filled.tif', 'dem_flowdir.tif', 'dem_area.tif',
                 'dem_length.tif', 'dem_slope.tif', 'dem_hillshade.tif'):
        assert name in produced

    back = dem.Area.load(os.path.join(outdir, 'dem_area.tif'))
    assert back._griddata.max() > 0.0
