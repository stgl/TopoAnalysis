"""Tests for :class:`TopoAnalysis.dem.BaseSpatialGrid` and its geometry."""

from __future__ import annotations

import os

import numpy as np
import pytest

from TopoAnalysis import dem, error


# ---------------------------------------------------------------------------
# Construction and dispatch
# ---------------------------------------------------------------------------


def test_grid_from_array_sets_a_usable_geotransform():
    data = np.arange(12, dtype=float).reshape(3, 4)
    grid = dem.BaseSpatialGrid(dx=2.0, grid=data)

    assert grid._georef_info.nx == 4
    assert grid._georef_info.ny == 3
    assert grid._georef_info.dx == 2.0
    # A grid built in memory must still be georeferenced well enough to save.
    gt = grid._georef_info.geoTransform
    assert gt[1] == 2.0 and gt[5] == -2.0
    assert gt[0] == pytest.approx(grid._georef_info.xllcenter - 1.0)
    assert gt[3] == pytest.approx(grid._georef_info.yllcenter + 2.5 * 2.0)


def test_unsatisfiable_keywords_raise_a_helpful_error():
    with pytest.raises(error.InputError) as excinfo:
        dem.BaseSpatialGrid(not_a_real_keyword=1)
    # The message must list what the class does accept.
    assert 'gdal_filename' in excinfo.value.msg


def test_random_grid_is_bounded_and_seedable():
    """Values lie in [0, 1); seeding NumPy makes the grid reproducible."""
    grid = dem.BaseSpatialGrid(nx=6, ny=4, dx=1.0)
    assert grid._griddata.shape == (4, 6)
    assert grid._griddata.min() >= 0.0 and grid._griddata.max() < 1.0
    # Two grids drawn in a row differ -- it is a random grid, after all ...
    assert not np.array_equal(grid._griddata,
                              dem.BaseSpatialGrid(nx=6, ny=4, dx=1.0)._griddata)
    # ... but reseeding NumPy's global generator reproduces one exactly.
    np.random.seed(1234)
    first = dem.BaseSpatialGrid(nx=6, ny=4, dx=1.0)._griddata.copy()
    np.random.seed(1234)
    second = dem.BaseSpatialGrid(nx=6, ny=4, dx=1.0)._griddata
    assert np.array_equal(first, second)


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------


def test_xy_and_rowcol_round_trip():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((7, 5)))
    for row in range(7):
        for col in range(5):
            ((x, y),) = grid._rowscols_to_xy(((row, col),))
            assert grid._xy_to_rowscols(((x, y),)) == ((row, col),)


def test_row_zero_is_the_north_edge():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((3, 3)))
    ((_, y_top),) = grid._rowscols_to_xy(((0, 0),))
    ((_, y_bottom),) = grid._rowscols_to_xy(((2, 0),))
    assert y_top > y_bottom


def test_points_outside_the_grid_are_rejected():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((3, 3)))
    info = grid._georef_info
    # One whole cell beyond the eastern edge -- this used to come back as a
    # valid subscript of nx, which then raised IndexError downstream.
    outside = (info.xllcenter + info.nx * info.dx, info.yllcenter)
    assert grid._xy_to_rowscols((outside,)) == ((None, None),)


def test_getitem_returns_none_outside_the_grid():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.ones((3, 3)))
    assert grid[1, 1] == 1.0
    assert grid[-1, 0] is None
    assert grid[3, 0] is None
    assert grid[0, 3] is None


def test_setitem_ignores_out_of_range_writes():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((3, 3)))
    grid[3, 0] = 5.0      # one past the last row
    grid[0, 3] = 5.0      # one past the last column
    grid[-1, 0] = 5.0
    assert not grid._griddata.any()


def test_coordinate_vectors_have_one_entry_per_cell():
    grid = dem.BaseSpatialGrid(dx=2.0, grid=np.zeros((4, 6)))
    x, y = grid.coordinate_vectors()
    # np.arange(start, start + (n-1)*dx, dx) yields n-1 values, which is what
    # get_XY_matricies used to return.
    assert x.size == 6 and y.size == 4
    X, Y = grid.get_XY_matricies()
    assert X.shape == (4, 6) and Y.shape == (4, 6)
    assert grid.get_XY_matrices()[0].shape == (4, 6)


def test_extent_covers_the_pixel_edges():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((3, 4)))
    xmin, xmax, ymin, ymax = grid.extent()
    assert xmax - xmin == pytest.approx(4 * 10.0)
    assert ymax - ymin == pytest.approx(3 * 10.0)


# ---------------------------------------------------------------------------
# Clipping, tiling, resampling
# ---------------------------------------------------------------------------


def test_clip_to_extent_keeps_the_requested_cells():
    data = np.arange(100, dtype=float).reshape(10, 10)
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    ((x0, y0),) = grid._rowscols_to_xy(((7, 2),))   # south-west corner
    ((x1, y1),) = grid._rowscols_to_xy(((3, 6),))   # north-east corner
    clipped = grid.clip_to_extent((x0, x1, y0, y1))

    assert clipped._griddata.shape == (5, 5)
    assert np.array_equal(clipped._griddata, data[3:8, 2:7])
    # The clip must be georeferenced where it was cut from.
    assert clipped._rowscols_to_xy(((0, 0),))[0] == pytest.approx((x0, y1))


def test_clip_rejects_an_extent_that_misses_the_grid():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((5, 5)))
    with pytest.raises(error.InputError):
        grid.clip_to_extent((1000.0, 1010.0, 1000.0, 1010.0))


def test_extent_of_data_ignores_nodata_margins():
    data = np.full((8, 10), np.nan)
    data[2:5, 3:7] = 1.0
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    xmin, xmax, ymin, ymax = grid.extent_of_data()
    ((x_left, y_bottom),) = grid._rowscols_to_xy(((4, 3),))
    ((x_right, y_top),) = grid._rowscols_to_xy(((2, 6),))
    assert (xmin, xmax, ymin, ymax) == pytest.approx(
        (x_left, x_right, y_bottom, y_top))


def test_tile_then_mosaic_reproduces_the_grid():
    rng = np.random.default_rng(4)
    data = rng.random((37, 41))
    grid = dem.BaseSpatialGrid(dx=5.0, grid=data)

    padding = 4
    tiles = grid.tile(tile_xdim=16, tile_ydim=16,
                      tile_xpadding=padding, tile_ypadding=padding)
    unpadded = [t.remove_padding(xpadding=padding, ypadding=padding) for t in tiles]
    rebuilt = dem.BaseSpatialGrid.mosaic(unpadded)

    assert rebuilt._griddata.shape == data.shape
    assert np.allclose(rebuilt._griddata, data)
    assert rebuilt._georef_info.xllcenter == pytest.approx(grid._georef_info.xllcenter)
    assert rebuilt._georef_info.yllcenter == pytest.approx(grid._georef_info.yllcenter)


def test_resample_preserves_a_linear_ramp():
    """A plane must resample to the same plane at any spacing."""
    nx = ny = 21
    x = np.arange(nx) * 4.0
    data = np.tile(x, (ny, 1))
    grid = dem.BaseSpatialGrid(dx=4.0, grid=data)

    coarse = grid.resample(8.0, interpolation='linear')
    assert coarse._georef_info.dx == 8.0
    assert coarse._griddata.shape[1] == pytest.approx(nx // 2, abs=1)

    # The value at each new cell centre must equal the plane's value there.
    xs, _ = coarse.coordinate_vectors()
    expected = np.tile(xs - grid._georef_info.xllcenter, (coarse._georef_info.ny, 1))
    assert np.allclose(coarse._griddata, expected, atol=1e-8)


def test_resample_rejects_an_unknown_interpolation():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((6, 6)))
    with pytest.raises(error.InputError):
        grid.resample(2.0, interpolation='bilinear-ish')


# ---------------------------------------------------------------------------
# Sorting and value lookup
# ---------------------------------------------------------------------------


def test_sort_orders_cells_by_value():
    data = np.array([[3.0, 1.0], [4.0, 2.0]])
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    ascending = grid.sort(reverse=False)
    assert data.reshape(-1)[ascending].tolist() == [1.0, 2.0, 3.0, 4.0]
    descending = grid.sort(reverse=True)
    assert data.reshape(-1)[descending].tolist() == [4.0, 3.0, 2.0, 1.0]


def test_sort_with_a_mask_keeps_only_masked_cells():
    data = np.array([[3.0, 1.0], [4.0, 2.0]])
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)
    mask = dem.BaseSpatialGrid(dx=1.0, grid=np.array([[1.0, 0.0], [1.0, 0.0]]))

    indices = grid.sort(reverse=False, force=True, mask=mask)
    # Only the left column (values 3 and 4) survives, in ascending order.
    assert data.reshape(-1)[indices].tolist() == [3.0, 4.0]


def test_find_nearest_cell_searches_the_full_radius():
    data = np.zeros((9, 9))
    data[4, 8] = 10.0            # four cells east of the centre
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    # range(i-r, i+r) stopped one short of the radius, so the far cell was
    # invisible at exactly r = 4.
    assert grid.find_nearest_cell_with_greatest_value((4, 4), pixel_radius=4) == (4, 8)


def test_find_nearest_cell_does_not_wrap_around_the_grid():
    data = np.zeros((5, 5))
    data[0, 0] = 1.0
    data[4, 4] = 99.0            # would be reached by a negative index
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    assert grid.find_nearest_cell_with_greatest_value((0, 0), pixel_radius=2) == (0, 0)


def test_find_nearest_cell_returns_none_when_there_is_no_data():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.full((5, 5), np.nan))
    assert grid.find_nearest_cell_with_greatest_value((2, 2), 2) == (None, None)


def test_set_value_at_rowscols_handles_an_empty_list():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((3, 3)))
    grid.set_value_at_rowscols(1.0, [])
    assert not grid._griddata.any()
    grid.set_value_at_rowscols(1.0, [(0, 0), (2, 2)])
    assert grid._griddata[0, 0] == 1.0 and grid._griddata[2, 2] == 1.0


def test_value_grid_set_value_at_indexes():
    grid = dem.ValueGrid(dx=1.0, grid=np.zeros((3, 3)))
    # zip() is a lazy iterator in Python 3; using one as a subscript raised.
    grid.set_value_at_indexes([(0, 1), (2, 2)], 7.0)
    assert grid._griddata[0, 1] == 7.0 and grid._griddata[2, 2] == 7.0


# ---------------------------------------------------------------------------
# Smoothing and derivatives
# ---------------------------------------------------------------------------


def test_average_over_distance_preserves_a_constant_field():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.full((21, 21), 5.0))
    averaged = grid.average_over_distance(3.0)
    assert averaged.shape == (21, 21)
    assert np.allclose(averaged, 5.0)


@pytest.mark.parametrize('shape', [(21, 21), (20, 20), (20, 25)])
def test_average_over_distance_is_centred(shape):
    """A point source must smear symmetrically about its own cell."""
    data = np.zeros(shape)
    centre = (shape[0] // 2, shape[1] // 2)
    data[centre] = 1.0
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    averaged = grid.average_over_distance(2.0)

    # The answer is the normalised disc, centred on the source cell.  An
    # off-by-half-a-cell kernel would place it one cell away.
    ny, nx = shape
    di = np.abs(np.arange(ny) - centre[0])
    dj = np.abs(np.arange(nx) - centre[1])
    di = np.minimum(di, ny - di)
    dj = np.minimum(dj, nx - dj)
    inside = np.hypot(di[:, None], dj[None, :]) <= 2.0
    expected = inside / inside.sum()

    assert np.allclose(averaged, expected, atol=1e-12)
    assert averaged.sum() == pytest.approx(1.0)


def test_average_over_distance_rejects_a_radius_that_selects_nothing():
    grid = dem.BaseSpatialGrid(dx=10.0, grid=np.zeros((5, 5)))
    with pytest.raises(error.InputError):
        grid.average_over_distance(-1.0)


def test_gradient_over_length_scale_on_a_plane():
    """A plane of slope 0.5 must give 0.5 whatever the measuring length."""
    nx = ny = 21
    data = np.tile(np.arange(nx) * 0.5, (ny, 1))
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    for length in (1.0, 3.0, 5.0):
        sx, sy = grid.calculate_gradient_over_length_scale(length)
        n = int(np.ceil(length))
        assert np.allclose(sx[n:-n, n:-n], 0.5)
        assert np.allclose(sy[n:-n, n:-n], 0.0, atol=1e-12)
        assert np.isnan(sx[0, 0])


def test_laplacian_over_length_scale_on_a_quadratic():
    """d2/dx2 of x^2 is 2, at any measuring length."""
    nx = ny = 25
    x = np.arange(nx, dtype=float)
    data = np.tile(x ** 2, (ny, 1))
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    for length in (1.0, 2.0, 4.0):
        curvature = grid.calculate_laplacian_over_length_scale(length)
        n = int(np.ceil(length))
        assert np.allclose(curvature[n:-n, n:-n], 2.0)


def test_length_scale_larger_than_the_grid_is_rejected():
    grid = dem.BaseSpatialGrid(dx=1.0, grid=np.zeros((5, 5)))
    with pytest.raises(error.InputError):
        grid.calculate_gradient_over_length_scale(10.0)


def test_principal_curvatures_of_a_sphere_cap():
    """Both principal curvatures of a paraboloid of revolution are equal."""
    n = 41
    y, x = np.mgrid[-20:21, -20:21].astype(float)
    data = -0.001 * (x ** 2 + y ** 2)
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    k1, k2 = grid.principal_curvatures()
    centre = (slice(15, 26), slice(15, 26))
    assert np.allclose(k1._griddata[centre], k2._griddata[centre], atol=2e-4)
    assert np.all(np.isfinite(k1._griddata[centre]))


def test_assign_bcs_replicates_every_corner():
    grid = dem.Elevation(dx=1.0, grid=np.arange(9, dtype=float).reshape(3, 3))
    padded = grid.assignBCs(grid._griddata, 3, 3)
    assert padded[0, 0] == grid._griddata[0, 0]
    assert padded[0, -1] == grid._griddata[0, -1]
    assert padded[-1, 0] == grid._griddata[-1, 0]
    # The south-east corner used to copy the south-west one.
    assert padded[-1, -1] == grid._griddata[-1, -1]


# ---------------------------------------------------------------------------
# File round-trips
# ---------------------------------------------------------------------------


@pytest.mark.gdal
def test_geotiff_round_trip(tmp_path):
    rng = np.random.default_rng(9)
    data = rng.random((11, 13)) * 100.0
    data[3, 4] = np.nan
    grid = dem.Elevation(dx=25.0, grid=data)

    path = str(tmp_path / 'grid.tif')
    grid.save(path)
    back = dem.Elevation.load(path)

    assert np.array_equal(back._griddata, data, equal_nan=True)
    assert back._georef_info.nx == 13 and back._georef_info.ny == 11
    assert back._georef_info.dx == pytest.approx(25.0)
    assert back._georef_info.xllcenter == pytest.approx(grid._georef_info.xllcenter)
    assert back._georef_info.yllcenter == pytest.approx(grid._georef_info.yllcenter)


@pytest.mark.gdal
def test_ascii_grid_round_trip_maps_nodata(tmp_path):
    data = np.array([[1.0, 2.0], [np.nan, 4.0]])
    grid = dem.Elevation(dx=1.0, grid=data)

    path = str(tmp_path / 'grid.asc')
    grid.write_to_ai(path)
    text = open(path).read()
    assert 'NODATA_value -9999' in text
    assert 'nan' not in text.lower()

    back = dem.Elevation(ai_ascii_filename=path, EPSGprojectionCode=32610)
    # The no-data flag must come back as NaN, not as an elevation of -9999.
    assert np.isnan(back._griddata[1, 0])
    assert back._griddata[0, 0] == pytest.approx(1.0)


@pytest.mark.gdal
def test_save_flushes_to_disk(tmp_path):
    grid = dem.Elevation(dx=1.0, grid=np.ones((4, 4)))
    path = str(tmp_path / 'flushed.tif')
    grid.save(path)
    assert os.path.getsize(path) > 0
