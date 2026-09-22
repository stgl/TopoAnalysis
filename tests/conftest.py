"""Shared fixtures: small synthetic landscapes with known answers.

This file also makes ``import TopoAnalysis`` resolve when the suite is run
from a checkout rather than from an installed copy.  It lives in ``tests/``
rather than beside ``__init__.py`` deliberately: a ``conftest.py`` next to a
package's ``__init__.py`` makes pytest import that package by the checkout
directory's *name*, which fails whenever the directory is called anything
other than ``TopoAnalysis`` -- a GitHub zip unpacks to ``TopoAnalysis-main``.
"""

from __future__ import annotations

import importlib.util
import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_PACKAGE_DIR = os.path.dirname(_TESTS_DIR)
_PARENT = os.path.dirname(_PACKAGE_DIR)


def _register_package():
    """Make ``import TopoAnalysis`` resolve to this checkout."""
    try:
        import TopoAnalysis  # noqa: F401
        return
    except ImportError:
        pass

    if os.path.basename(_PACKAGE_DIR) == "TopoAnalysis":
        if _PARENT not in sys.path:
            sys.path.insert(0, _PARENT)
        return

    # The checkout has some other name; load it under the one the package
    # expects rather than asking the user to rename their directory.
    spec = importlib.util.spec_from_file_location(
        "TopoAnalysis", os.path.join(_PACKAGE_DIR, "__init__.py"),
        submodule_search_locations=[_PACKAGE_DIR])
    module = importlib.util.module_from_spec(spec)
    sys.modules["TopoAnalysis"] = module
    spec.loader.exec_module(module)


_register_package()

import numpy as np                                            # noqa: E402
import pytest                                                 # noqa: E402

from TopoAnalysis import dem, kernels                          # noqa: E402


# ---------------------------------------------------------------------------
# Surfaces whose correct answers can be written down
# ---------------------------------------------------------------------------


@pytest.fixture
def east_plane():
    """A 5x5 plane dropping 1 m per 1 m eastward.

    Everything flows east, gradient is exactly 1, aspect is exactly 90.
    """
    return dem.Elevation(dx=1.0, grid=np.tile(np.arange(5, 0, -1, dtype=float), (5, 1)))


@pytest.fixture
def south_plane():
    """A 5x5 plane dropping 1 m per 1 m southward (row-increasing)."""
    return dem.Elevation(dx=1.0, grid=np.tile(
        np.arange(5, 0, -1, dtype=float).reshape(-1, 1), (1, 5)))


@pytest.fixture
def single_pit():
    """A 5x5 bowl: a rim at 10, a ring at 4-5 and a 1 m pit in the middle.

    Filling raises every interior cell to the 10 m rim.
    """
    return dem.Elevation(dx=1.0, grid=np.array([
        [10, 10, 10, 10, 10],
        [10, 5, 4, 5, 10],
        [10, 4, 1, 4, 10],
        [10, 5, 4, 5, 10],
        [10, 10, 10, 10, 10]], dtype=float))


@pytest.fixture
def plateau():
    """A 5x5 grid with a flat plateau at 3 and one outlet cell at 2."""
    return dem.Elevation(dx=1.0, grid=np.array([
        [5, 5, 5, 5, 5],
        [5, 3, 3, 3, 5],
        [5, 3, 3, 3, 5],
        [5, 3, 3, 2, 5],
        [5, 5, 5, 5, 5]], dtype=float))


@pytest.fixture
def v_valley():
    """A 7-column V-shaped valley draining south along the middle column.

    The channel is column 3; every row drains inward then down.
    """
    ny, nx = 9, 7
    col = np.abs(np.arange(nx) - 3.0)
    z = np.tile(col, (ny, 1)) * 2.0
    z += np.arange(ny, dtype=float).reshape(-1, 1)[::-1] * 0.5
    return dem.Elevation(dx=10.0, grid=z)


@pytest.fixture
def rough_dome():
    """A 40x40 noisy dome with real depressions.

    The noise amplitude is chosen so that the surface genuinely has closed
    depressions -- a smoother dome makes every filling test vacuous, because
    there is nothing to fill.  ``test_dome_fixture_has_depressions`` pins
    that down.
    """
    rng = np.random.default_rng(20240921)
    n = 40
    y, x = np.mgrid[0:n, 0:n] / n
    z = (300.0 * (1.0 - ((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.5)
         + 12.0 * np.sin(5 * np.pi * x) * np.sin(5 * np.pi * y)
         + rng.normal(0.0, 6.0, (n, n)))
    # A handful of deliberate pits, so the depressions are not an accident of
    # the random seed.
    for (i, j) in ((9, 12), (22, 8), (30, 27), (15, 31)):
        z[i, j] -= 40.0
    return dem.Elevation(dx=30.0, grid=z)


@pytest.fixture
def dome_network(rough_dome):
    """``(elevation, filled, flow_direction, area, flow_length)`` for the dome."""
    filled = dem.FilledElevation(elevation=rough_dome)
    flow_direction = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=flow_direction)
    flow_length = dem.FlowLength(flow_direction=flow_direction)
    return rough_dome, filled, flow_direction, area, flow_length


@pytest.fixture
def geographic_grid():
    """A 10x10 lat/lon grid near 40 degrees north, 0.01 degree cells."""
    grid = dem.GeographicElevation(dx=0.01, grid=np.zeros((10, 10)))
    grid._georef_info.geoTransform = (-120.0, 0.01, 0.0, 40.1, 0.0, -0.01)
    grid._georef_info.xllcenter = -119.995
    grid._georef_info.yllcenter = 40.005
    return grid


# ---------------------------------------------------------------------------
# Environment markers
# ---------------------------------------------------------------------------


def pytest_runtest_setup(item):
    """Skip GDAL-dependent tests when the bindings are missing."""
    if any(mark.name == 'gdal' for mark in item.iter_markers()):
        if not dem.gdal_is_available():
            pytest.skip('the GDAL Python bindings are not installed')


@pytest.fixture
def backend_name():
    """Which kernel backend the test run is exercising."""
    return kernels.backend()
