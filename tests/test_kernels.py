"""Tests for the compute kernels.

Two things are checked here:

1. that the compiled C++ kernels and the pure-NumPy ones agree exactly, so
   the two backends are interchangeable; and
2. that the Priority-Flood implementation has the properties Barnes et al.
   (2014) prove it has -- the filled surface is minimal, monotone, and
   independent of which of the paper's algorithms produced it.
"""

from __future__ import annotations

import numpy as np
import pytest

from TopoAnalysis import fastops, kernels


FILL_CASES = [
    ('flat', 0.0, 0.0),
    ('epsilon', 0.0, 0.0),
    ('epsilon', 1e-3, 0.0),
    ('flat', 0.0, 3.0),
    ('epsilon', 1e-4, 5.0),
]


def random_surface(seed, ny=17, nx=23, nodata_fraction=0.0):
    rng = np.random.default_rng(seed)
    z = rng.random((ny, nx)) * 20.0
    if nodata_fraction:
        z[rng.random(z.shape) < nodata_fraction] = np.nan
    return z


def build_network(z, cellsize=1.0):
    """``(filled, codes, receivers, order)`` for a surface."""
    filled = np.ascontiguousarray(z, dtype=float)
    kernels.priority_flood(filled, mode='epsilon', epsilon=1e-9, cellsize=cellsize)
    codes = kernels.flow_directions(filled, cellsize=cellsize)
    receivers = kernels.receivers(codes)
    order, cycles = kernels.topological_order(receivers)
    assert cycles == 0
    return filled, codes, np.asarray(receivers), np.asarray(order)


# ---------------------------------------------------------------------------
# Backend agreement
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('seed', [0, 1, 2, 3])
@pytest.mark.parametrize('mode,epsilon,max_depth', FILL_CASES)
def test_fill_backends_agree(seed, mode, epsilon, max_depth):
    z = random_surface(seed, nodata_fraction=0.08 if seed % 2 else 0.0)

    a, b = z.copy(), z.copy()
    ra = kernels.priority_flood(a, mode=mode, epsilon=epsilon,
                                max_pit_depth=max_depth, track_visited=True)
    rb = fastops.priority_flood(b, mode=mode, epsilon=epsilon,
                                max_pit_depth=max_depth, track_visited=True)

    assert np.array_equal(a, b, equal_nan=True)
    assert np.array_equal(np.asarray(ra['visited']), np.asarray(rb['visited']))
    assert ra['cells_filled'] == rb['cells_filled']
    assert ra['max_fill_depth'] == pytest.approx(rb['max_fill_depth'])


@pytest.mark.parametrize('seed', [4, 5])
def test_routing_backends_agree(seed):
    z = random_surface(seed)
    filled, codes, receivers, order = build_network(z)

    assert np.array_equal(codes, fastops.flow_directions(filled, cellsize=1.0))
    assert np.array_equal(receivers, fastops.receivers(codes))

    weights = np.ones_like(filled)
    gate = (filled > np.median(filled)).astype(np.uint8)
    py_order, _ = fastops.topological_order(receivers)

    assert np.allclose(kernels.accumulate(receivers, order, weights),
                       fastops.accumulate(receivers, py_order, weights))
    assert np.allclose(kernels.accumulate(receivers, order, weights, gate),
                       fastops.accumulate(receivers, py_order, weights, gate))

    length, codes_c = kernels.flow_length(receivers, order, weights)
    length_py, codes_py = fastops.flow_length(receivers, py_order, weights)
    assert np.allclose(length, length_py)
    assert np.array_equal(codes_c, codes_py)

    outlets = np.flatnonzero(receivers.reshape(-1) < 0).astype(np.int64)
    assert np.array_equal(kernels.upstream_mask(receivers, outlets),
                          fastops.upstream_mask(receivers, outlets))
    assert np.array_equal(kernels.drainage_basins(receivers, order),
                          fastops.drainage_basins(receivers, py_order))

    accumulated = np.asarray(kernels.accumulate(receivers, order, weights))
    for trapezoid in (False, True):
        chi_c, dist_c = kernels.chi(receivers, order, accumulated, weights,
                                    outlets, 1.0, 0.45, trapezoid, 0.0, None)
        chi_py, dist_py = fastops.chi(receivers, py_order, accumulated, weights,
                                      outlets, 1.0, 0.45, trapezoid, 0.0, None)
        assert np.allclose(chi_c, chi_py, equal_nan=True)
        assert np.allclose(dist_c, dist_py, equal_nan=True)

    assert np.allclose(kernels.downstream_distance(receivers, order, weights),
                       fastops.downstream_distance(receivers, py_order, weights))

    carved_c, carved_py = filled.copy(), filled.copy()
    kernels.imposemin(receivers, order, weights, 1e-3, carved_c)
    fastops.imposemin(receivers, py_order, weights, 1e-3, carved_py)
    assert np.allclose(carved_c, carved_py)

    streams = (accumulated > 5).astype(np.uint8)
    for kind in ('strahler', 'shreve'):
        assert np.array_equal(kernels.stream_order(receivers, order, streams, kind),
                              fastops.stream_order(receivers, py_order, streams, kind))


# ---------------------------------------------------------------------------
# Priority-Flood properties (Barnes et al. 2014)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('seed', [10, 11, 12])
@pytest.mark.parametrize('mode,epsilon', [('flat', 0.0), ('epsilon', 1e-6)])
def test_improved_algorithm_matches_the_reference(seed, mode, epsilon):
    """Algorithm 2 must return exactly what Algorithm 1 returns."""
    z = random_surface(seed)
    improved, reference = z.copy(), z.copy()
    kernels.priority_flood(improved, mode=mode, epsilon=epsilon)
    kernels.priority_flood(reference, mode=mode, epsilon=epsilon,
                           reference_algorithm=True)
    assert np.array_equal(improved, reference)


@pytest.mark.parametrize('seed', [20, 21])
def test_filling_never_lowers_a_cell(seed):
    z = random_surface(seed)
    filled = z.copy()
    kernels.priority_flood(filled, mode='flat')
    assert np.all(filled >= z - 1e-12)


def test_flat_fill_produces_the_spill_elevation(single_pit):
    """Every interior cell of the bowl ends up at the 10 m rim."""
    filled = np.ascontiguousarray(single_pit._griddata, dtype=float)
    report = kernels.priority_flood(filled, mode='flat')
    assert np.all(filled == 10.0)
    assert report['cells_filled'] == 9
    assert report['max_fill_depth'] == pytest.approx(9.0)


def test_epsilon_fill_leaves_no_flats(single_pit):
    """Priority-Flood+epsilon must leave a strictly descending path."""
    filled = np.ascontiguousarray(single_pit._griddata, dtype=float)
    kernels.priority_flood(filled, mode='epsilon', epsilon=1e-3, cellsize=1.0)
    interior = filled[1:-1, 1:-1]
    assert np.all(interior > 10.0)
    # The deepest cell must end up highest: it is furthest from the spill.
    assert filled[2, 2] == interior.max()


def test_max_pit_depth_leaves_deep_depressions_alone(single_pit):
    filled = np.ascontiguousarray(single_pit._griddata, dtype=float)
    kernels.priority_flood(filled, mode='flat', max_pit_depth=2.0)
    # The 9 m deep centre stays put ...
    assert filled[2, 2] == pytest.approx(1.0)
    # ... but the flood still reached the far side of it.
    assert np.all(np.isfinite(filled))


def test_nodata_acts_as_an_outlet():
    """A depression touching no-data drains into it and is not filled."""
    z = np.array([
        [10, 10, 10, 10, 10],
        [10, 5, 5, 5, 10],
        [np.nan, 5, 1, 5, 10],
        [10, 5, 5, 5, 10],
        [10, 10, 10, 10, 10]], dtype=float)
    filled = z.copy()
    kernels.priority_flood(filled, mode='flat')
    # (2,1) sits next to the NaN, so it keeps its own elevation ...
    assert filled[2, 1] == pytest.approx(5.0)
    # ... and the pit only rises to the level it must cross to get there.
    assert filled[2, 2] == pytest.approx(5.0)
    assert np.isnan(filled[2, 0])


def test_default_seeds_are_the_perimeter_and_nodata_fringe():
    z = np.arange(25, dtype=float).reshape(5, 5)
    z[2, 2] = np.nan
    seeds = set(kernels.default_seeds(z).tolist())

    perimeter = {i * 5 + j for i in range(5) for j in range(5)
                 if i in (0, 4) or j in (0, 4)}
    nan_fringe = {(2 + di) * 5 + (2 + dj)
                  for di in (-1, 0, 1) for dj in (-1, 0, 1)} - {2 * 5 + 2}

    assert seeds == perimeter | nan_fringe


# ---------------------------------------------------------------------------
# Flow routing
# ---------------------------------------------------------------------------


def test_d8_codes_follow_the_arcgis_convention():
    """Each of the eight codes must point at the neighbour it names."""
    #  32 64 128
    #  16  X   1
    #   8  4   2
    expected = {(0, 1): 1, (1, 1): 2, (1, 0): 4, (1, -1): 8,
                (0, -1): 16, (-1, -1): 32, (-1, 0): 64, (-1, 1): 128}
    for (di, dj), code in expected.items():
        z = np.ones((3, 3))
        z[1 + di, 1 + dj] = 0.0
        codes = kernels.flow_directions(z, cellsize=1.0)
        assert codes[1, 1] == code, 'offset {0} should give code {1}'.format((di, dj), code)


def test_diagonals_compete_on_gradient_not_on_drop():
    """A diagonal neighbour only wins if its *gradient* is steeper."""
    z = np.ones((3, 3)) * 10.0
    z[1, 2] = 9.0          # east: drop 1 over 1  -> gradient 1.00
    z[2, 2] = 8.6          # south-east: drop 1.4 over sqrt(2) -> gradient 0.99
    assert kernels.flow_directions(z, cellsize=1.0)[1, 1] == 1

    z[2, 2] = 8.5          # gradient 1.06, now steeper than east
    assert kernels.flow_directions(z, cellsize=1.0)[1, 1] == 2


def test_accumulation_conserves_weight(rough_dome):
    z = np.ascontiguousarray(rough_dome._griddata, dtype=float)
    _, _, receivers, order = build_network(z, cellsize=30.0)
    weights = np.ones_like(z)
    accumulated = np.asarray(kernels.accumulate(receivers, order, weights))

    outlets = receivers.reshape(-1) < 0
    # Everything that enters the grid must leave through some outlet.
    assert accumulated.reshape(-1)[outlets].sum() == pytest.approx(weights.sum())
    assert accumulated.min() >= 1.0


def test_accumulation_ignores_the_order_it_is_given(rough_dome):
    """Any valid topological order must give the same accumulation."""
    z = np.ascontiguousarray(rough_dome._griddata, dtype=float)
    _, _, receivers, order = build_network(z, cellsize=30.0)
    other_order, _ = fastops.topological_order(receivers)
    weights = np.ones_like(z)
    assert np.allclose(kernels.accumulate(receivers, order, weights),
                       kernels.accumulate(receivers, other_order, weights))


def test_flow_length_tie_break_is_deterministic():
    """Equal-length paths must resolve the same way whatever the order."""
    z = np.array([
        [2.0, 3.0, 2.0],
        [3.0, 1.0, 3.0],
        [2.0, 3.0, 2.0]])
    _, _, receivers, order = build_network(z)
    other_order, _ = fastops.topological_order(receivers)
    step = np.ones_like(z)
    _, codes_a = kernels.flow_length(receivers, order, step)
    _, codes_b = kernels.flow_length(receivers, other_order, step)
    assert np.array_equal(codes_a, codes_b)


def test_upstream_mask_covers_the_whole_basin(rough_dome):
    z = np.ascontiguousarray(rough_dome._griddata, dtype=float)
    _, _, receivers, order = build_network(z, cellsize=30.0)
    outlets = np.flatnonzero(receivers.reshape(-1) < 0).astype(np.int64)
    mask = np.asarray(kernels.upstream_mask(receivers, outlets))
    # Every cell drains to some outlet, so the union covers the grid.
    assert mask.all()


def test_chi_is_zero_at_the_outlet_and_grows_upstream(v_valley):
    z = np.ascontiguousarray(v_valley._griddata, dtype=float)
    filled, _, receivers, order = build_network(z, cellsize=10.0)
    step = np.full(z.shape, 10.0)
    area = np.asarray(kernels.accumulate(receivers, order, np.full(z.shape, 100.0)))

    outlets = np.flatnonzero(receivers.reshape(-1) < 0).astype(np.int64)
    chi, distance = kernels.chi(receivers, order, area, step, outlets,
                                1e4, 0.45, True, 0.0, None)
    chi = np.asarray(chi)

    assert np.nanmin(chi) == 0.0
    # chi must increase monotonically upstream: every cell exceeds its receiver.
    flat_chi = chi.reshape(-1)
    for k, r in enumerate(receivers.reshape(-1)):
        if r >= 0 and np.isfinite(flat_chi[k]) and np.isfinite(flat_chi[r]):
            assert flat_chi[k] >= flat_chi[r] - 1e-12


def test_trapezoid_and_left_endpoint_chi_differ_predictably(v_valley):
    """The trapezoid rule must sit between the two endpoint rules."""
    z = np.ascontiguousarray(v_valley._griddata, dtype=float)
    _, _, receivers, order = build_network(z, cellsize=10.0)
    step = np.full(z.shape, 10.0)
    area = np.asarray(kernels.accumulate(receivers, order, np.full(z.shape, 100.0)))
    outlets = np.flatnonzero(receivers.reshape(-1) < 0).astype(np.int64)

    left, _ = kernels.chi(receivers, order, area, step, outlets, 1e4, 0.45, False, 0.0, None)
    trap, _ = kernels.chi(receivers, order, area, step, outlets, 1e4, 0.45, True, 0.0, None)
    valid = np.isfinite(left) & np.isfinite(trap)
    # Area grows downstream, so the integrand shrinks downstream and the
    # left-endpoint (upstream) value is the larger of the two.
    assert np.all(np.asarray(trap)[valid] <= np.asarray(left)[valid] + 1e-9)


def test_imposemin_only_lowers_and_enforces_the_gradient(rough_dome):
    z = np.ascontiguousarray(rough_dome._griddata, dtype=float)
    filled, _, receivers, order = build_network(z, cellsize=30.0)
    step = np.full(z.shape, 30.0)

    carved = filled.copy()
    kernels.imposemin(receivers, order, step, 1e-3, carved)

    assert np.all(carved <= filled + 1e-12)
    flat = carved.reshape(-1)
    for k, r in enumerate(receivers.reshape(-1)):
        if r >= 0:
            assert flat[k] - flat[r] >= 1e-3 * 30.0 - 1e-9


def test_strahler_order_of_a_simple_confluence():
    """Two order-1 streams meeting make an order-2 stream."""
    z = np.array([[3.0, 9.0, 3.0],
                  [9.0, 2.0, 9.0],
                  [9.0, 1.0, 9.0]])
    _, _, receivers, order = build_network(z)
    streams = np.array([[1, 0, 1], [0, 1, 0], [0, 1, 0]], dtype=np.uint8)

    strahler = np.asarray(kernels.stream_order(receivers, order, streams, 'strahler'))
    assert strahler[0, 0] == 1 and strahler[0, 2] == 1
    assert strahler[1, 1] == 2 and strahler[2, 1] == 2

    shreve = np.asarray(kernels.stream_order(receivers, order, streams, 'shreve'))
    assert shreve[1, 1] == 2 and shreve[2, 1] == 2


def test_drainage_basins_partition_the_grid(rough_dome):
    z = np.ascontiguousarray(rough_dome._griddata, dtype=float)
    _, _, receivers, order = build_network(z, cellsize=30.0)
    labels = np.asarray(kernels.drainage_basins(receivers, order))

    assert labels.min() >= 1, 'every cell drains somewhere on this surface'
    # Each cell must carry the same label as the cell it drains into.
    flat = labels.reshape(-1)
    for k, r in enumerate(receivers.reshape(-1)):
        if r >= 0:
            assert flat[k] == flat[r]


def test_topological_order_reports_cycles():
    """A hand-made loop must be reported, not silently mis-ordered."""
    # Two cells pointing at each other: 1 -> east, 16 -> west.
    codes = np.array([[1, 16]], dtype=np.uint8)
    receivers = kernels.receivers(codes)
    order, cycles = kernels.topological_order(receivers)
    assert cycles == 2
    assert sorted(np.asarray(order).tolist()) == [0, 1]
