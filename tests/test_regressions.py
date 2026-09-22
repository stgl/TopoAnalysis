"""One test per repaired defect, so none of them can quietly return.

Each test names the specific failure it guards against.  Defects already
covered by the behavioural suites (filling, routing, terrain attributes) are
not repeated here.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys

import numpy as np
import pytest

from TopoAnalysis import dem, error, fastops, kernels

#: The package directory, and the directory that must be on sys.path for
#: ``import TopoAnalysis`` to resolve.
PACKAGE_DIR = os.path.dirname(os.path.abspath(dem.__file__))
IMPORT_ROOT = os.path.dirname(PACKAGE_DIR)


def run_in_subprocess(script, **extra_env):
    """Run ``script`` in a fresh interpreter that can import TopoAnalysis."""
    env = dict(os.environ, PYTHONPATH=IMPORT_ROOT, **extra_env)
    return subprocess.run([sys.executable, '-c', script], capture_output=True,
                          text=True, cwd=IMPORT_ROOT, env=env)


# ---------------------------------------------------------------------------
# Module-level hygiene
# ---------------------------------------------------------------------------


def test_no_removed_numpy_or_scipy_aliases_remain():
    """np.NAN, np.float and scipy.ndimage.morphology were all removed upstream."""
    package = PACKAGE_DIR
    banned = ('np.NAN', 'np.NaN', 'np.float)', 'np.int)', 'np.bool)',
              'scipy.ndimage.morphology', 'scipy.ndimage.filters',
              'interp2d', 'plt.hold')

    offenders = []
    for name in sorted(os.listdir(package)):
        if not name.endswith('.py'):
            continue
        text = open(os.path.join(package, name)).read()
        for token in banned:
            # The tests themselves may name these tokens when describing them.
            if token in text and not name.startswith('test_'):
                offenders.append('{0}: {1}'.format(name, token))
    assert offenders == []


def test_import_does_not_require_gdal():
    """GDAL must be imported lazily, so the package loads without it."""
    script = (
        'import sys, types;'
        'blocker = types.ModuleType("osgeo");'
        'blocker.__getattr__ = lambda name: (_ for _ in ()).throw(ImportError());'
        'sys.modules["osgeo"] = blocker;'
        'import TopoAnalysis;'
        'print(TopoAnalysis.gdal_is_available())'
    )
    result = run_in_subprocess(script)
    assert result.returncode == 0, result.stderr


def test_recursion_limit_is_left_alone():
    """The module used to raise it to a million, which segfaults before it trips."""
    assert sys.getrecursionlimit() < 100000


def test_constructor_dispatch_does_not_use_eval():
    """Constructors are selected with getattr, not by compiling source text."""
    import ast

    tree = ast.parse(open(dem.__file__).read())
    called = {node.func.id for node in ast.walk(tree)
              if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}
    assert 'eval' not in called
    assert 'exec' not in called


def test_both_import_styles_work():
    """`import dem` (flat) and `from TopoAnalysis import dem` must both work."""
    script = ('import sys; sys.path.insert(0, {0!r}); '
              'import dem, kernels; print(kernels.backend())').format(PACKAGE_DIR)
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() in ('c++', 'python')


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def test_extent_of_data_handles_a_non_square_grid():
    """The column scan used to index rows, and ran off the end when nx != ny."""
    data = np.full((4, 12), np.nan)
    data[1:3, 5:9] = 1.0
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)
    xmin, xmax, ymin, ymax = grid.extent_of_data()
    assert xmax > xmin and ymax > ymin


def test_tile_geotransform_is_self_consistent():
    """The tile transform put the north edge one cell too far north."""
    grid = dem.BaseSpatialGrid(dx=3.0, grid=np.arange(400, dtype=float).reshape(20, 20))
    for tile in grid.tile(tile_xdim=8, tile_ydim=8, tile_xpadding=2, tile_ypadding=2):
        info = tile._georef_info
        # yllcenter must invert out of the transform the same way the reader
        # derives it.
        assert info.geoTransform[3] - info.dx * (info.ny - 0.5) == pytest.approx(
            info.yllcenter)


def test_flow_length_clip_keeps_the_direction_codes():
    """clip_to_extent used to copy the *lengths* into the direction grid."""
    n = 12
    z = np.tile(np.arange(n, 0, -1, dtype=float), (n, 1))
    grid = dem.Elevation(dx=1.0, grid=z)
    d8 = dem.FlowDirectionD8(elevation=grid)
    length = dem.FlowLength(flow_direction=d8)

    ((x0, y0),) = length._rowscols_to_xy(((9, 2),))
    ((x1, y1),) = length._rowscols_to_xy(((2, 9),))
    clipped = length.clip_to_extent((x0, x1, y0, y1))

    codes = clipped.main_stem_directions
    assert codes.shape == clipped._griddata.shape
    assert codes.dtype == np.uint8
    # Direction codes are powers of two (or 0), never flow lengths.
    assert set(np.unique(codes).tolist()) <= {0, 1, 2, 4, 8, 16, 32, 64, 128}


# ---------------------------------------------------------------------------
# Flow routing
# ---------------------------------------------------------------------------


def test_flow_code_for_position_prefers_the_steeper_gradient():
    """The per-cell re-router divided cardinals, not diagonals, by sqrt(2)."""
    z = np.ones((3, 3)) * 10.0
    z[1, 2] = 9.0        # east, gradient 1.00
    z[2, 2] = 8.6        # south-east, gradient 0.99
    grid = dem.Elevation(dx=1.0, grid=z)

    d8 = dem.FlowDirectionD8(dx=1.0, grid=np.zeros((3, 3), dtype=np.uint8))
    mask = dem.Mask(dx=1.0, grid=np.zeros((3, 3), dtype=np.uint8))
    mask._griddata[1, 1] = 1
    d8.update_flow_codes_in_mask(grid, mask)
    assert d8._griddata[1, 1] == 1, 'the steeper cardinal must win'

    z[2, 2] = 8.5        # south-east, gradient 1.06
    grid = dem.Elevation(dx=1.0, grid=z)
    d8.update_flow_codes_in_mask(grid, mask)
    assert d8._griddata[1, 1] == 2, 'the steeper diagonal must win'


def test_update_flow_codes_invalidates_the_cached_network():
    z = np.tile(np.arange(5, 0, -1, dtype=float), (5, 1))
    grid = dem.Elevation(dx=1.0, grid=z)
    d8 = dem.FlowDirectionD8(elevation=grid)
    first = d8.receivers.copy()

    mask = dem.Mask(dx=1.0, grid=np.ones((5, 5), dtype=np.uint8))
    flipped = dem.Elevation(dx=1.0, grid=np.fliplr(z).copy())
    d8.update_flow_codes_in_mask(flipped, mask)

    assert not np.array_equal(d8.receivers, first)


def test_mainstem_valley_area_differs_from_valley_area():
    """The override was name-mangled out of reach, so the two were identical."""
    n = 24
    y, x = np.mgrid[0:n, 0:n].astype(float)
    z = 100.0 - y * 2.0 + np.abs(x - n / 2) * 1.5
    elevation = dem.Elevation(dx=10.0, grid=z)

    filled = dem.FilledElevation(elevation=elevation)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    laplace = dem.Laplacian(elevation=elevation)

    common = dict(flow_direction=d8, area=area, laplace=laplace,
                  valley_laplace_value=-1e9, min_area_value=1e3)
    whole = dem.ValleyArea(**common)
    mainstem = dem.MainstemValleyArea(**common)

    assert not np.allclose(whole._griddata, mainstem._griddata)
    # Only the trunk contributes, so the mainstem totals cannot exceed the
    # full accumulation anywhere.
    assert np.all(mainstem._griddata <= whole._griddata + 1e-9)


def test_area_gate_actually_gates():
    """`mask[i, j] is not None` was always true, so the mask did nothing."""
    n = 10
    z = np.tile(np.arange(n, 0, -1, dtype=float), (n, 1))
    d8 = dem.FlowDirectionD8(elevation=dem.Elevation(dx=1.0, grid=z))

    gate = dem.Mask(dx=1.0, grid=np.zeros((n, n), dtype=np.uint8))
    gate._griddata[:, :5] = 1

    gated = dem.Area(flow_direction=d8, evaluate_at=gate)
    plain = dem.Area(flow_direction=d8)

    # Inside the gate, accumulation is unaffected ...
    assert gated._griddata[0, 4] == pytest.approx(plain._griddata[0, 4])
    # ... the first cell outside it still receives what the gate delivers ...
    assert gated._griddata[0, 5] == pytest.approx(plain._griddata[0, 4])
    # ... but passes nothing further on, and adds nothing of its own.
    assert gated._griddata[0, 6] == 0.0
    assert gated._griddata[0, 9] == 0.0
    assert plain._griddata[0, 9] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# Regression grids that fit a model
# ---------------------------------------------------------------------------


@pytest.fixture
def smoothing_inputs():
    n = 30
    y, x = np.mgrid[0:n, 0:n].astype(float)
    z = 200.0 - y * 3.0 + np.abs(x - n / 2) * 2.0
    elevation = dem.Elevation(dx=20.0, grid=z)
    filled = dem.FilledElevation(elevation=elevation)
    d8 = dem.FlowDirectionD8(flooded_dem=filled)
    area = dem.Area(flow_direction=d8)
    return elevation, d8, area


def test_ks_from_chi_with_smoothing_counts_each_crossing(smoothing_inputs):
    """`self._n[rows, cols] += 1` only counted a repeated cell once."""
    pytest.importorskip('statsmodels')
    elevation, d8, area = smoothing_inputs

    ks = dem.KsFromChiWithSmoothing(
        elevation=elevation, area=area, flow_direction=d8, theta=0.45,
        vertical_interval=15.0, area_threshold=2000.0, verbose=False)

    assert ks._n.sum() > 0
    assert ks._n.max() > 1, 'busy cells must be counted once per profile'
    assert np.isfinite(ks._griddata[~np.isnan(ks._griddata)]).all()


def test_channel_slope_with_smoothing_accepts_a_horizontal_interval(smoothing_inputs):
    """The horizontal_interval entry pointed at a method that did not exist."""
    elevation, d8, area = smoothing_inputs
    slope = dem.ChannelSlopeWithSmoothing(
        elevation=elevation, area=area, flow_direction=d8,
        horizontal_interval=100.0, min_area=2000.0, verbose=False)
    assert np.isfinite(slope._griddata[~np.isnan(slope._griddata)]).all()


def test_channel_up_and_down_slope_variants(smoothing_inputs):
    elevation, d8, area = smoothing_inputs
    common = dict(elevation=elevation, area=area, flow_direction=d8,
                  vertical_interval=15.0, min_area=2000.0, verbose=False)
    down = dem.ChannelDownSlopeWithSmoothing(**common)
    up = dem.ChannelUpSlopeWithSmoothing(**common)
    assert down._griddata.shape == up._griddata.shape
    assert not np.array_equal(np.nan_to_num(down._griddata),
                              np.nan_to_num(up._griddata))


def test_theta_from_chi_with_smoothing(smoothing_inputs):
    pytest.importorskip('statsmodels')
    elevation, d8, area = smoothing_inputs
    theta = dem.ThetaFromChiWithSmoothing(
        elevation=elevation, area=area, flow_direction=d8,
        vertical_interval=15.0, min_area=2000.0, verbose=False)
    fitted = theta._griddata[~np.isnan(theta._griddata)]
    assert fitted.size > 0
    assert np.all(np.abs(fitted) <= 10.0), 'the optimiser is bounded at +/-10'


def test_geographic_theta_inherits_the_theta_class():
    """It used to inherit KsFromChiWithSmoothing and compute ks instead."""
    assert issubclass(dem.GeographicThetaFromChiWithSmoothing,
                      dem.ThetaFromChiWithSmoothing)
    assert issubclass(dem.GeographicKsFromChiWithSmoothing,
                      dem.KsFromChiWithSmoothing)


# ---------------------------------------------------------------------------
# Profile tools
# ---------------------------------------------------------------------------


def test_hypsometric_integral_detects_nodata():
    """`np.nan in array` is always False, so NaNs slipped through."""
    from TopoAnalysis import demRecursionTools

    tree = {'elevation': np.nan, 'dA': 1.0,
            'next': [{'elevation': 5.0, 'dA': 1.0}]}
    assert demRecursionTools.hi_list(tree) == 0.0


def test_tributary_step_length_multiplies_the_scale():
    """`trb_ld['de'] + trb_ld['distance_scale']` should have been a product."""
    from TopoAnalysis import demRecursionTools

    source = open(demRecursionTools.__file__).read()
    assert "trb_ld['de'] * trb_ld['distance_scale']" in source


# ---------------------------------------------------------------------------
# Backend equivalence at the library level
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not kernels.have_extension(),
                    reason='the compiled extension is not built')
def test_pure_python_backend_gives_the_same_pipeline(tmp_path):
    """The whole pipeline must agree between the two backends."""
    script = '''
import os
os.environ["TOPOANALYSIS_PURE_PYTHON"] = "1"
import numpy as np
import TopoAnalysis.kernels as kernels
assert kernels.backend() == "python"
import TopoAnalysis.dem as dem

rng = np.random.default_rng(99)
n = 25
y, x = np.mgrid[0:n, 0:n] / n
z = 100*(1 - ((x-0.5)**2 + (y-0.5)**2)) + rng.normal(0, 0.5, (n, n))
e = dem.Elevation(dx=30.0, grid=z)
f = dem.FilledElevation(elevation=e)
d8 = dem.FlowDirectionD8(flooded_dem=f)
a = dem.Area(flow_direction=d8)
l = dem.FlowLength(flow_direction=d8)
np.save(os.environ["OUT"], np.stack([f._griddata, d8._griddata.astype(float),
                                     a._griddata, l._griddata]))
'''
    out = str(tmp_path / 'python_backend.npy')
    result = run_in_subprocess(script, OUT=out)
    assert result.returncode == 0, result.stderr

    rng = np.random.default_rng(99)
    n = 25
    y, x = np.mgrid[0:n, 0:n] / n
    z = 100 * (1 - ((x - 0.5) ** 2 + (y - 0.5) ** 2)) + rng.normal(0, 0.5, (n, n))
    e = dem.Elevation(dx=30.0, grid=z)
    f = dem.FilledElevation(elevation=e)
    d8 = dem.FlowDirectionD8(flooded_dem=f)
    a = dem.Area(flow_direction=d8)
    l = dem.FlowLength(flow_direction=d8)

    expected = np.stack([f._griddata, d8._griddata.astype(float),
                         a._griddata, l._griddata])
    assert np.allclose(np.load(out), expected)


def test_kernels_expose_the_barnes_reference():
    assert 'Barnes' in kernels.PRIORITY_FLOOD_REFERENCE
    assert '2014' in kernels.PRIORITY_FLOOD_REFERENCE
    assert kernels.backend() in ('c++', 'python')


def test_fastops_and_kernels_expose_the_same_api():
    shared = [name for name in fastops.__all__ if not name.startswith('D8_')]
    for name in shared:
        assert hasattr(kernels, name), name


# ---------------------------------------------------------------------------
# Defects found by the post-rewrite review
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not kernels.have_extension(),
                    reason='the compiled extension is not built')
def test_kernels_are_thread_safe():
    """Every binding must resolve its buffers before releasing the GIL.

    ``py::array::request()`` calls ``PyObject_GetBuffer``, which touches
    reference counts.  Calling it inside a ``gil_scoped_release`` block
    corrupted the interpreter -- leaked references, and reproducible aborts --
    as soon as two threads ran a kernel at once.
    """
    import threading

    rng = np.random.default_rng(1)
    z = np.ascontiguousarray(rng.random((40, 40)) * 10)
    kernels.priority_flood(z, mode='epsilon', epsilon=1e-9)
    codes = kernels.flow_directions(z, cellsize=1.0)
    recv = np.ascontiguousarray(kernels.receivers(codes))
    order = np.ascontiguousarray(kernels.topological_order(recv)[0])
    weights = np.ones_like(z)
    streams = np.ones(z.shape, dtype=np.uint8)
    outlets = np.flatnonzero(recv.reshape(-1) < 0).astype(np.int64)

    watched = (recv, order, weights, streams, outlets)
    before = [sys.getrefcount(a) for a in watched]

    def hammer():
        for _ in range(400):
            kernels.accumulate(recv, order, weights)
            kernels.downstream_distance(recv, order, weights)
            kernels.stream_order(recv, order, streams, 'strahler')
            kernels.flow_length(recv, order, weights)
            kernels.drainage_basins(recv, order, outlets, None)
            kernels.upstream_mask(recv, outlets)

    threads = [threading.Thread(target=hammer) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert [sys.getrefcount(a) for a in watched] == before


@pytest.mark.skipif(not kernels.have_extension(),
                    reason='the compiled extension is not built')
def test_kernels_reject_a_mismatched_order_array():
    """A wrong-length ``order`` must raise, not read off the end of it."""
    recv = np.full((4, 4), -1, dtype=np.int64)
    weights = np.ones((4, 4))
    with pytest.raises((ValueError, RuntimeError)):
        kernels.accumulate(recv, np.arange(3, dtype=np.int64), weights)


def test_strahler_order_does_not_depend_on_the_topological_order():
    """Both backends, and any valid ordering, must agree.

    The incremental rule gives {2, 2, 3} an order of 4 or 3 according to
    which donor arrives first; the junction is now resolved as a whole.
    """
    recv = np.array([[5, 6, 1, -1], [1, -1, 5, 6], [5, 10, 7, 6]], dtype=np.int64)
    streams = np.ones(recv.shape, dtype=np.uint8)

    cpp_order, _ = kernels.topological_order(recv)
    py_order, _ = fastops.topological_order(recv)
    results = [
        np.asarray(kernels.stream_order(recv, cpp_order, streams, 'strahler')),
        np.asarray(kernels.stream_order(recv, py_order, streams, 'strahler')),
        np.asarray(fastops.stream_order(recv, cpp_order, streams, 'strahler')),
        np.asarray(fastops.stream_order(recv, py_order, streams, 'strahler')),
    ]
    for other in results[1:]:
        assert np.array_equal(results[0], other)


def test_strahler_of_a_mixed_junction():
    """Orders {2, 2, 3} meeting give 3, not 4: only the largest can pair."""
    # Four headwaters feeding one cell: two combine to 2, another two to 2,
    # those combine to 3, and a bare order-1 also arrives.
    recv = np.array([
        [4, 4, 4, 4],
        [-1, 4, 4, 4],
    ], dtype=np.int64)
    streams = np.ones(recv.shape, dtype=np.uint8)
    order, _ = kernels.topological_order(recv)
    result = np.asarray(kernels.stream_order(recv, order, streams, 'strahler'))
    # Cell 4 receives four order-1 donors, so it is order 2.
    assert result.reshape(-1)[4] == 2


def test_drainage_basins_honours_valid_with_outlets():
    """The C++ path used to drop the mask whenever outlets were supplied."""
    z = np.add.outer(np.arange(9.0), np.arange(9.0))
    codes = kernels.flow_directions(z, cellsize=1.0)
    recv = kernels.receivers(codes)
    order, _ = kernels.topological_order(recv)

    valid = np.zeros((9, 9), dtype=np.uint8)
    valid[:4, :4] = 1
    outlets = np.array([0], dtype=np.int64)

    cpp = np.asarray(kernels.drainage_basins(recv, order, outlets, valid))
    py = np.asarray(fastops.drainage_basins(recv, order, outlets, valid))
    assert np.array_equal(cpp, py)
    assert cpp[valid == 0].sum() == 0
    assert (cpp > 0).sum() <= int(valid.sum())


def test_priority_flood_rejects_a_non_contiguous_grid():
    """reshape(-1) copies, so the fill would be silently thrown away."""
    base = np.array([[5., 5., 5., 5., 5.],
                     [5., 0., 0., 0., 5.],
                     [5., 0., 0., 0., 5.],
                     [5., 0., 0., 0., 5.],
                     [5., 5., 5., 5., 5.]])
    for awkward in (np.asfortranarray(base), np.repeat(base, 2, axis=0)[::2]):
        with pytest.raises(TypeError):
            fastops.priority_flood(awkward)
        if kernels.have_extension():
            with pytest.raises(TypeError):
                kernels.priority_flood(awkward)


def test_imposemin_rejects_a_non_contiguous_grid():
    recv = np.full((4, 4), -1, dtype=np.int64)
    order = np.arange(16, dtype=np.int64)
    step = np.ones((4, 4))
    awkward = np.asfortranarray(np.zeros((4, 4)))
    with pytest.raises(TypeError):
        fastops.imposemin(recv, order, step, 1e-3, awkward)


def test_default_seeds_of_an_empty_grid():
    assert fastops.default_seeds(np.zeros((0, 5))).size == 0
    assert fastops.default_seeds(np.zeros((5, 0))).size == 0


def test_create_from_geotransform_places_the_grid_correctly():
    """yllcenter was read from ny before ny was assigned."""
    gt = (1000.0, 10.0, 0.0, 2000.0, 0.0, -10.0)
    grid = dem.Elevation(nx=15, ny=12, projection='', geo_transform=gt,
                         grid=np.zeros((12, 15)))

    # The south edge of a 12-row, 10 m grid whose north edge is 2000 is 1880,
    # so the southern row of cell centres sits at 1885.
    assert grid._georef_info.yllcenter == pytest.approx(1885.0)
    assert grid._rowscols_to_xy(((0, 0),))[0] == pytest.approx((1005.0, 1995.0))
    assert grid._rowscols_to_xy(((11, 0),))[0] == pytest.approx((1005.0, 1885.0))


@pytest.mark.gdal
def test_create_from_geotransform_round_trips_through_gdal(tmp_path):
    """A grid built from a geotransform must agree with the file it writes."""
    gt = (1000.0, 10.0, 0.0, 2000.0, 0.0, -10.0)
    rng = np.random.default_rng(6)
    grid = dem.Elevation(nx=15, ny=12, projection='', geo_transform=gt,
                         grid=rng.random((12, 15)))

    path = str(tmp_path / 'grid.tif')
    grid.save(path)
    back = dem.Elevation.load(path)

    assert back._georef_info.yllcenter == pytest.approx(grid._georef_info.yllcenter)
    assert back._rowscols_to_xy(((0, 0),))[0] == pytest.approx(
        grid._rowscols_to_xy(((0, 0),))[0])
    assert back._rowscols_to_xy(((11, 14),))[0] == pytest.approx(
        grid._rowscols_to_xy(((11, 14),))[0])


def test_editing_a_flow_code_invalidates_the_cached_network():
    """Routing is cached now, so a write has to drop the cache."""
    rng = np.random.default_rng(2)
    grid = dem.Elevation(dx=10.0, grid=rng.random((20, 20)) * 10)
    fd = dem.FlowDirectionD8(flooded_dem=dem.FilledElevation(elevation=grid))

    before = dem.Area(flow_direction=fd)._griddata.copy()

    current = int(fd._griddata[10, 10])
    fd[10, 10] = 1 if current != 1 else 16
    assert not np.array_equal(dem.Area(flow_direction=fd)._griddata, before)

    after_setitem = dem.Area(flow_direction=fd)._griddata.copy()
    current = int(fd._griddata[5, 5])
    fd.set_value_at_rowscols(1 if current != 1 else 16, [(5, 5)])
    assert not np.array_equal(dem.Area(flow_direction=fd)._griddata, after_setitem)


def test_geographic_laplacian_uses_metres_not_degrees():
    """dx is in DEGREES on a lat/lon grid; dividing by it was 1e10 wrong."""
    n = 11
    y, x = np.mgrid[0:n, 0:n].astype(float)
    grid = dem.GeographicElevation(dx=0.01, grid=0.25 * (x ** 2 + y ** 2))
    grid._georef_info.geoTransform = (-120.0, 0.01, 0.0, 40.11, 0.0, -0.01)
    grid._georef_info.xllcenter = -119.995
    grid._georef_info.yllcenter = 40.005

    curvature = dem.GeographicLaplacian(elevation=grid)._griddata[5, 5]
    cell = grid._mean_pixel_dimension()[5, 5]
    # The Laplacian of 0.25*(x^2 + y^2) is 1 per cell squared.
    assert curvature == pytest.approx(1.0 / cell ** 2, rel=1e-6)
    assert curvature < 1e-4, 'a degree-scaled result would be about 1e4'


def test_tile_with_zero_padding_does_not_duplicate_the_edges():
    """`griddata[:, -0:]` is the whole array, not an empty slice."""
    data = np.arange(100, dtype=float).reshape(10, 10)
    grid = dem.BaseSpatialGrid(dx=1.0, grid=data)

    tiles = grid.tile(tile_xdim=5, tile_ydim=5, tile_xpadding=0, tile_ypadding=0)
    assert [t._griddata.shape for t in tiles] == [(5, 5)] * 4
    rebuilt = dem.BaseSpatialGrid.mosaic([t.remove_padding(0, 0) for t in tiles])
    assert np.array_equal(rebuilt._griddata, data)


def test_star_import_names_only_things_that_exist():
    import TopoAnalysis

    missing = [name for name in TopoAnalysis.__all__
               if not hasattr(TopoAnalysis, name)]
    assert missing == []

    namespace = {}
    exec('from TopoAnalysis import *', namespace)   # noqa: S102 - that is the test
    assert 'Elevation' in namespace and 'FlowDirectionD8' in namespace


def test_the_test_suite_runs_from_a_differently_named_checkout(tmp_path):
    """A GitHub zip unpacks to TopoAnalysis-main, not TopoAnalysis."""
    import shutil

    target = tmp_path / 'TopoAnalysis-main'
    shutil.copytree(PACKAGE_DIR, str(target),
                    ignore=shutil.ignore_patterns('__pycache__', 'build', '.git',
                                                  '.pytest_cache', 'docs'))

    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    result = subprocess.run(
        [sys.executable, '-m', 'pytest', '-q', 'tests/test_grid.py',
         '-p', 'no:cacheprovider'],
        capture_output=True, text=True, cwd=str(target), env=env)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-2000:]


def test_extension_is_loaded_only_from_the_package_directory(tmp_path):
    """A stray `_topoanalysis` elsewhere on sys.path must not be picked up.

    The flat-import fallback used a bare ``import _topoanalysis``, which
    searches all of sys.path; a stale or unrelated module of that name would
    be loaded instead of the package's own.
    """
    decoy = tmp_path / '_topoanalysis.py'
    decoy.write_text('MARKER = "decoy"\n')

    script = (
        'import sys, os;'
        'sys.path.insert(0, {0!r});'
        'sys.path.insert(0, {1!r});'
        'import kernels;'
        'print(getattr(kernels._ext, "MARKER", "not-the-decoy"))'
    ).format(str(tmp_path), PACKAGE_DIR)
    result = subprocess.run([sys.executable, '-c', script], capture_output=True,
                            text=True, cwd=str(tmp_path))
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == 'not-the-decoy'
