"""Tests for the terrain-attribute grids, against analytic answers."""

from __future__ import annotations

import numpy as np
import pytest

from TopoAnalysis import dem


# ---------------------------------------------------------------------------
# Slope and hillshade
# ---------------------------------------------------------------------------


def test_max_slope_on_a_plane(east_plane):
    slope = dem.MaxSlope(elevation=east_plane)
    # Centred differences reproduce a plane exactly in the interior.
    assert np.allclose(slope._griddata[1:-1, 1:-1], 1.0)


def test_max_slope_scales_with_cell_size():
    data = np.tile(np.arange(7, dtype=float), (7, 1))
    fine = dem.MaxSlope(elevation=dem.Elevation(dx=1.0, grid=data))
    coarse = dem.MaxSlope(elevation=dem.Elevation(dx=2.0, grid=data))
    assert np.allclose(fine._griddata[1:-1, 1:-1],
                       2.0 * coarse._griddata[1:-1, 1:-1])


def test_hillshade_lights_the_slope_that_faces_the_sun(east_plane):
    """An east-facing slope is lit from the east and shadowed from the west."""
    lit = dem.Hillshade(elevation=east_plane, azimuth=90, inclination=45)
    shadowed = dem.Hillshade(elevation=east_plane, azimuth=270, inclination=45)
    assert lit._griddata[2, 2] == 255
    assert shadowed._griddata[2, 2] == 0


def test_hillshade_of_a_south_facing_slope(south_plane):
    lit = dem.Hillshade(elevation=south_plane, azimuth=180, inclination=45)
    shadowed = dem.Hillshade(elevation=south_plane, azimuth=0, inclination=45)
    assert lit._griddata[2, 2] > shadowed._griddata[2, 2]


def test_hillshade_is_a_byte_grid_without_wraparound(rough_dome):
    hs = dem.Hillshade(elevation=rough_dome, azimuth=315, inclination=45)
    assert hs._griddata.dtype == np.uint8
    # Negative cosines used to wrap round the byte into bright speckle.
    assert hs._griddata.min() >= 0 and hs._griddata.max() <= 255


def test_hillshade_of_a_flat_surface_is_the_sine_of_the_altitude():
    flat = dem.Elevation(dx=1.0, grid=np.zeros((7, 7)))
    hs = dem.Hillshade(elevation=flat, azimuth=315, inclination=30)
    assert hs._griddata[3, 3] == pytest.approx(255 * np.cos(np.radians(60)), abs=1)


# ---------------------------------------------------------------------------
# Curvature
# ---------------------------------------------------------------------------


def test_laplacian_of_a_quadratic_bowl():
    """The Laplacian of a(x^2 + y^2) is 4a."""
    n = 21
    y, x = np.mgrid[0:n, 0:n].astype(float)
    a = 0.25
    grid = dem.Elevation(dx=1.0, grid=a * (x ** 2 + y ** 2))

    curvature = dem.Laplacian(elevation=grid)
    assert np.allclose(curvature._griddata[2:-2, 2:-2], 4 * a)


def test_laplacian_scales_with_cell_size():
    n = 15
    y, x = np.mgrid[0:n, 0:n].astype(float)
    data = x ** 2
    fine = dem.Laplacian(elevation=dem.Elevation(dx=1.0, grid=data))
    coarse = dem.Laplacian(elevation=dem.Elevation(dx=2.0, grid=data))
    assert np.allclose(fine._griddata[2:-2, 2:-2],
                       4.0 * coarse._griddata[2:-2, 2:-2])


def test_local_relief_window_is_centred():
    """A single peak must raise the relief symmetrically about itself."""
    data = np.zeros((21, 21))
    data[10, 10] = 100.0
    grid = dem.Elevation(dx=1.0, grid=data)

    relief = dem.LocalRelief(elevation=grid, pixel_radius=3)
    # np.ogrid[-r:r] gave a window one cell short on the far side, so the
    # footprint was off-centre.
    assert relief._griddata[10, 7] == pytest.approx(relief._griddata[10, 13])
    assert relief._griddata[7, 10] == pytest.approx(relief._griddata[13, 10])
    assert relief._griddata[10, 10] == pytest.approx(100.0)


def test_local_relief_equals_the_range_over_the_disc(rough_dome):
    """The value is exactly the range over the disc, not over its bounding box."""
    radius = 2
    relief = dem.LocalRelief(elevation=rough_dome, pixel_radius=radius)

    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    disc = (xx * xx + yy * yy) <= radius * radius
    for (i, j) in ((10, 10), (20, 15), (31, 27)):
        window = rough_dome._griddata[i - radius:i + radius + 1,
                                      j - radius:j + radius + 1][disc]
        assert relief._griddata[i, j] == pytest.approx(window.max() - window.min())

    # The bounding box is strictly larger, so it is a different answer.
    box = rough_dome._griddata[8:13, 8:13]
    assert relief._griddata[10, 10] <= box.max() - box.min() + 1e-9
    assert relief._griddata.min() >= 0.0


def test_log_area_is_a_float_grid(rough_dome):
    filled = dem.FilledElevation(elevation=rough_dome)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)

    log_area = dem.LogArea(area=area)
    # Declaring this uint8 truncated log10(area) to a handful of integers.
    assert log_area._griddata.dtype == np.float64
    finite = np.isfinite(log_area._griddata)
    assert np.allclose(10 ** log_area._griddata[finite], area._griddata[finite])


def test_gradient_grid_components(east_plane):
    gradient = dem.Gradient(elevation=east_plane)
    assert np.allclose(gradient._gx[1:-1, 1:-1], -1.0)
    assert np.allclose(gradient._gy[1:-1, 1:-1], 0.0)


# ---------------------------------------------------------------------------
# Geographic grids
# ---------------------------------------------------------------------------


def test_geographic_cell_area_shrinks_polewards(geographic_grid):
    areas = geographic_grid._area_per_pixel()
    # Row 0 is the north edge, so its cells are the smallest.
    assert areas[0, 0] < areas[-1, 0]
    assert np.allclose(areas[:, 0], areas[:, -1])


def test_geographic_cell_area_matches_the_spherical_formula(geographic_grid):
    areas = geographic_grid._area_per_pixel()
    re = geographic_grid.earth_radius
    dlon = dlat = 0.01
    lat = geographic_grid._georef_info.yllcenter + 9 * dlat   # northernmost row
    expected = (re ** 2) * np.radians(dlon) * (
        np.sin(np.radians(lat + dlat / 2)) - np.sin(np.radians(lat - dlat / 2)))
    assert areas[0, 0] == pytest.approx(expected, rel=1e-9)


def test_geographic_mean_pixel_dimension_is_the_area_root(geographic_grid):
    assert np.allclose(geographic_grid._mean_pixel_dimension() ** 2,
                       geographic_grid._area_per_pixel())


def test_geographic_area_uses_spherical_cells():
    """Drainage area on a lat/lon grid must use the true cell areas."""
    n = 12
    z = np.tile(np.arange(n, 0, -1, dtype=float), (n, 1))
    grid = dem.GeographicElevation(dx=0.01, grid=z)
    grid._georef_info.geoTransform = (-120.0, 0.01, 0.0, 40.12, 0.0, -0.01)
    grid._georef_info.xllcenter = -119.995
    grid._georef_info.yllcenter = 40.005

    d8 = dem.FlowDirectionD8(flooded_dem=dem.FilledElevation(elevation=grid))
    area = dem.GeographicArea(flow_direction=d8)

    per_pixel = grid._area_per_pixel()
    assert area._griddata[5, 0] == pytest.approx(per_pixel[5, 0])
    # Cells accumulate the true, latitude-dependent areas.
    assert area._griddata[5, 4] == pytest.approx(per_pixel[5, 0:5].sum())
