"""Pure-NumPy implementations of the TopoAnalysis grid kernels.

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

Conventions
-----------
Grids are 2-D, row-major, with row 0 at the north edge.  ``NaN`` is no-data.
Flow-direction codes follow the ArcGIS convention::

    | 32  64 128 |
    | 16   X   1 |
    |  8   4   2 |
"""

from __future__ import annotations

import heapq
import math

import numpy as np

__all__ = [
    "D8_DI",
    "D8_DJ",
    "D8_DIST",
    "D8_CODES",
    "default_seeds",
    "priority_flood",
    "flow_directions",
    "receivers",
    "topological_order",
    "accumulate",
    "flow_length",
    "propagate_along_main_stem",
    "upstream_mask",
    "chi",
    "drainage_basins",
]

# Neighbour order E, SE, S, SW, W, NW, N, NE -- so that the ArcGIS code for
# direction ``d`` is simply ``1 << d``.
D8_DI = np.array([0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int64)
D8_DJ = np.array([1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)
D8_DIST = np.array([1.0, math.sqrt(2.0)] * 4, dtype=np.float64)
D8_CODES = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.uint8)

_SQRT2 = math.sqrt(2.0)


def _as_2d_float(grid, name):
    arr = np.asarray(grid, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError("{0} must be a 2-D array".format(name))
    return arr


def _shifted(grid, di, dj, fill):
    """``grid`` translated so that element (i, j) holds grid[i+di, j+dj]."""
    out = np.full_like(grid, fill)
    ny, nx = grid.shape
    src_i = slice(max(di, 0), ny + min(di, 0))
    src_j = slice(max(dj, 0), nx + min(dj, 0))
    dst_i = slice(max(-di, 0), ny + min(-di, 0))
    dst_j = slice(max(-dj, 0), nx + min(-dj, 0))
    out[dst_i, dst_j] = grid[src_i, src_j]
    return out


# ---------------------------------------------------------------------------
# Depression filling
# ---------------------------------------------------------------------------


def default_seeds(elevations):
    """Flat indices of the cells a flood should start from.

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
    """
    z = _as_2d_float(elevations, "elevations")
    ny, nx = z.shape
    if ny == 0 or nx == 0:
        return np.empty(0, dtype=np.int64)
    nodata = np.isnan(z)

    seed = np.zeros((ny, nx), dtype=bool)
    seed[0, :] = True
    seed[-1, :] = True
    seed[:, 0] = True
    seed[:, -1] = True
    for di, dj in zip(D8_DI, D8_DJ):
        seed |= _shifted(nodata.astype(np.float64), int(di), int(dj), 0.0) > 0.5
    seed &= ~nodata
    return np.flatnonzero(seed).astype(np.int64)


def priority_flood(
    elevations,
    closed=None,
    seeds=None,
    mode="flat",
    epsilon=0.0,
    cellsize=1.0,
    max_pit_depth=0.0,
    track_visited=False,
    reference_algorithm=False,
):
    """Fill depressions in ``elevations`` in place (Barnes et al. 2014).

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
    """
    if mode not in ("flat", "epsilon"):
        raise ValueError("mode must be 'flat' or 'epsilon'")

    z = elevations
    if not isinstance(z, np.ndarray) or z.dtype != np.float64 or z.ndim != 2:
        raise TypeError("elevations must be a 2-D float64 ndarray (modified in place)")
    if not z.flags['C_CONTIGUOUS']:
        # reshape(-1) would return a *copy* of a non-contiguous array, so the
        # fill would land in the copy and the caller's grid would come back
        # untouched while the report claimed success.  The compiled kernel
        # rejects the same input, so both backends agree.
        raise TypeError(
            "elevations must be C-contiguous to be filled in place; "
            "pass np.ascontiguousarray(elevations)")
    ny, nx = z.shape
    n = ny * nx
    flat = z.reshape(-1)

    if closed is None:
        closed_flat = np.zeros(n, dtype=bool)
    else:
        closed_flat = np.asarray(closed, dtype=bool).reshape(-1).copy()
        if closed_flat.size != n:
            raise ValueError("closed must have the same shape as elevations")
    closed_flat |= np.isnan(flat)

    if seeds is None:
        seed_idx = default_seeds(z)
        keep = seed_idx[~closed_flat[seed_idx]]
        if closed is not None and keep.size == 0:
            keep = _mask_boundary_seeds(closed_flat.reshape(ny, nx))
        seed_idx = keep
    else:
        seed_idx = np.asarray(seeds, dtype=np.int64).reshape(-1)

    use_epsilon = mode == "epsilon"
    limit_depth = max_pit_depth > 0.0
    # The FIFO shortcut is only sound when every cell it holds sits at exactly
    # the current flood level.  Two options break that invariant: an epsilon
    # fill raises each cell by its own increment, and a depth limit re-queues
    # unfilled cells *below* the current level.  In either case everything goes
    # through the priority queue, i.e. the run degrades to Algorithm 1.
    use_pit = not use_epsilon and not limit_depth and not reference_algorithm

    visited = np.zeros(n, dtype=np.uint8) if track_visited else None
    unfilled = np.zeros(n, dtype=np.uint8) if limit_depth else None
    cells_filled = 0
    max_fill_depth = 0.0

    # ``order`` makes the heap ordering total so ties resolve first-in
    # first-out, exactly as the compiled kernel does.
    order = 0
    open_heap = []
    for k in seed_idx:
        k = int(k)
        if k < 0 or k >= n or closed_flat[k]:
            continue
        closed_flat[k] = True
        heapq.heappush(open_heap, (flat[k], order, k))
        order += 1

    from collections import deque

    pit = deque()

    while open_heap or pit:
        if pit and open_heap and open_heap[0][0] <= flat[pit[0]]:
            cur_z, _, cur = heapq.heappop(open_heap)
        elif pit:
            cur = pit.popleft()
            cur_z = flat[cur]
        else:
            cur_z, _, cur = heapq.heappop(open_heap)

        if visited is not None:
            visited[cur] = 1

        i, j = divmod(cur, nx)
        for d in range(8):
            ni = i + int(D8_DI[d])
            nj = j + int(D8_DJ[d])
            if ni < 0 or nj < 0 or ni >= ny or nj >= nx:
                continue
            nk = ni * nx + nj
            if closed_flat[nk]:
                continue
            closed_flat[nk] = True
            if np.isnan(flat[nk]):
                continue

            if use_epsilon:
                if epsilon > 0.0:
                    target = cur_z + epsilon * cellsize * D8_DIST[d]
                    if target <= cur_z:
                        target = np.nextafter(cur_z, np.inf)
                else:
                    target = np.nextafter(cur_z, np.inf)
            else:
                target = cur_z

            if flat[nk] < target:
                depth = float(target) - float(flat[nk])
                if limit_depth and depth > max_pit_depth:
                    unfilled[nk] = 1
                    heapq.heappush(open_heap, (flat[nk], order, nk))
                    order += 1
                    continue
                flat[nk] = target
                cells_filled += 1
                if depth > max_fill_depth:
                    max_fill_depth = depth
                if use_pit:
                    pit.append(nk)
                else:
                    heapq.heappush(open_heap, (flat[nk], order, nk))
                    order += 1
            else:
                heapq.heappush(open_heap, (flat[nk], order, nk))
                order += 1

    return {
        "visited": None if visited is None else visited.reshape(ny, nx),
        "unfilled": None if unfilled is None else unfilled.reshape(ny, nx),
        "cells_filled": cells_filled,
        "max_fill_depth": max_fill_depth,
    }


def _mask_boundary_seeds(closed_2d):
    """Seeds on the edge of the open region when the mask excludes the rim."""
    open_region = ~closed_2d
    ny, nx = closed_2d.shape
    boundary = np.zeros((ny, nx), dtype=bool)
    boundary[0, :] = True
    boundary[-1, :] = True
    boundary[:, 0] = True
    boundary[:, -1] = True
    for di, dj in zip(D8_DI, D8_DJ):
        boundary |= _shifted(closed_2d.astype(np.float64), int(di), int(dj), 0.0) > 0.5
    return np.flatnonzero((boundary & open_region).reshape(-1)).astype(np.int64)


# ---------------------------------------------------------------------------
# Flow directions
# ---------------------------------------------------------------------------


def flow_directions(elevations, cellsize_grid=None, cellsize=1.0):
    """Steepest-descent D8 flow directions as ArcGIS codes.

    The drop to each neighbour is divided by the centre-to-centre distance, so
    a diagonal neighbour only wins if its *gradient* is steeper.  Ties go to
    the first direction in the order E, SE, S, SW, W, NW, N, NE.

    Cells with no lower neighbour -- pits, and no-data cells -- get code 0.
    """
    z = _as_2d_float(elevations, "elevations")
    if cellsize_grid is None:
        de = np.full(z.shape, float(cellsize))
    else:
        de = _as_2d_float(cellsize_grid, "cellsize_grid")
        if de.shape != z.shape:
            raise ValueError("cellsize_grid must match the elevation grid")

    best_slope = np.zeros(z.shape)
    codes = np.zeros(z.shape, dtype=np.uint8)

    for d in range(8):
        neighbour = _shifted(z, int(D8_DI[d]), int(D8_DJ[d]), np.nan)
        slope = (z - neighbour) / (de * D8_DIST[d])
        better = slope > best_slope
        # NaN comparisons are False, so no-data neighbours never win.
        best_slope = np.where(better, slope, best_slope)
        codes = np.where(better, D8_CODES[d], codes).astype(np.uint8)

    codes[np.isnan(z)] = 0
    return codes


def receivers(codes):
    """Flat index of the cell each cell drains to, or -1 if there is none."""
    codes = np.asarray(codes, dtype=np.uint8)
    if codes.ndim != 2:
        raise ValueError("codes must be a 2-D array")
    ny, nx = codes.shape
    ii, jj = np.mgrid[0:ny, 0:nx]

    out = np.full(codes.shape, -1, dtype=np.int64)
    for d in range(8):
        sel = codes == D8_CODES[d]
        if not sel.any():
            continue
        ni = ii + int(D8_DI[d])
        nj = jj + int(D8_DJ[d])
        inside = sel & (ni >= 0) & (nj >= 0) & (ni < ny) & (nj < nx)
        out[inside] = (ni * nx + nj)[inside]
    return out


# ---------------------------------------------------------------------------
# Topological ordering and accumulation
# ---------------------------------------------------------------------------


def _levels(recv_flat, n):
    """Kahn's algorithm, one frontier at a time.

    Returns ``(chunks, n_on_cycles)`` where each chunk is an array of cell
    indices that can be processed simultaneously: no cell in a chunk drains
    into another cell of the same chunk.
    """
    has_recv = (recv_flat >= 0) & (recv_flat < n) & (recv_flat != np.arange(n))
    indegree = np.bincount(recv_flat[has_recv], minlength=n)

    frontier = np.flatnonzero(indegree == 0)
    chunks = []
    seen = 0
    while frontier.size:
        chunks.append(frontier)
        seen += frontier.size
        moving = frontier[has_recv[frontier]]
        if moving.size == 0:
            break
        targets = recv_flat[moving]
        indegree -= np.bincount(targets, minlength=n)
        candidates = np.unique(targets)
        frontier = candidates[indegree[candidates] == 0]

    return chunks, n - seen


def topological_order(recv):
    """Order the cells donors-first.

    Returns ``(order, n_cells_on_cycles)``.  Cells that lie on a cycle -- only
    possible in a hand-edited or externally supplied flow-direction grid --
    cannot be ordered and are appended at the end in index order.
    """
    recv_flat = np.asarray(recv, dtype=np.int64).reshape(-1)
    n = recv_flat.size
    chunks, n_cycles = _levels(recv_flat, n)
    if chunks:
        order = np.concatenate(chunks)
    else:
        order = np.empty(0, dtype=np.int64)
    if n_cycles:
        emitted = np.zeros(n, dtype=bool)
        emitted[order] = True
        order = np.concatenate([order, np.flatnonzero(~emitted)])
    return order.astype(np.int64), int(n_cycles)


def accumulate(recv, order, weights, gate=None):
    """Accumulate ``weights`` downstream along the D8 network.

    ``order`` is accepted for signature compatibility with the compiled
    kernel; the result does not depend on which valid topological order is
    used, so this implementation derives its own level structure -- which
    lets it run vectorised rather than one cell at a time.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    shape = np.asarray(weights).shape
    recv_flat = recv_arr.reshape(-1)
    n = recv_flat.size

    acc = np.asarray(weights, dtype=np.float64).reshape(-1).copy()
    if gate is not None:
        gate_flat = np.asarray(gate).reshape(-1).astype(bool)
        acc[~gate_flat] = 0.0
    else:
        gate_flat = None

    chunks, _ = _levels(recv_flat, n)
    for chunk in chunks:
        if gate_flat is not None:
            chunk = chunk[gate_flat[chunk]]
            if chunk.size == 0:
                continue
        r = recv_flat[chunk]
        valid = (r >= 0) & (r < n) & (r != chunk)
        if not valid.any():
            continue
        np.add.at(acc, r[valid], acc[chunk[valid]])

    return acc.reshape(shape)


def flow_length(recv, order, step_length):
    """Longest upstream flow distance, and the main-stem donor direction.

    Returns ``(length, from_codes)``.  ``from_codes[k]`` is the direction code
    pointing from ``k`` back at the donor that supplies its longest path.
    Ties are settled in favour of the donor with the lower flat index, so the
    main stem does not depend on which topological order was used.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    ny, nx = recv_arr.shape
    n = ny * nx
    recv_flat = recv_arr.reshape(-1)
    step = np.asarray(step_length, dtype=np.float64).reshape(-1)
    order = np.asarray(order, dtype=np.int64).reshape(-1)

    length = np.zeros(n, dtype=np.float64)
    from_index = np.full(n, -1, dtype=np.int64)

    for k in order:
        k = int(k)
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k:
            continue
        candidate = length[k] + step[k]
        incumbent = int(from_index[r])
        if candidate > length[r] or (
            candidate == length[r] and (incumbent < 0 or k < incumbent)
        ):
            length[r] = candidate
            from_index[r] = k

    code_for_offset = {
        (int(D8_DI[d]), int(D8_DJ[d])): int(D8_CODES[d]) for d in range(8)
    }
    from_codes = np.zeros(n, dtype=np.uint8)
    contributing = np.flatnonzero(from_index >= 0)
    for r in contributing:
        r = int(r)
        k = int(from_index[r])
        offset = (k // nx - r // nx, k % nx - r % nx)
        from_codes[r] = code_for_offset.get(offset, 0)

    return length.reshape(ny, nx), from_codes.reshape(ny, nx)


def propagate_along_main_stem(recv, order, main_stem_from_codes, gate, values, mode):
    """Carry or sum ``values`` downstream, but only along longest-flow paths.

    ``mode='carry'`` overwrites the receiver with the donor's value (used by
    :class:`~TopoAnalysis.dem.Relief`); ``mode='sum'`` adds it (used by
    :class:`~TopoAnalysis.dem.Ksi`).
    """
    if mode not in ("carry", "sum"):
        raise ValueError("mode must be 'carry' or 'sum'")

    recv_arr = np.asarray(recv, dtype=np.int64)
    ny, nx = recv_arr.shape
    n = ny * nx
    recv_flat = recv_arr.reshape(-1)
    codes_flat = np.asarray(main_stem_from_codes, dtype=np.uint8).reshape(-1)
    out = np.asarray(values, dtype=np.float64).reshape(-1).copy()
    gate_flat = None if gate is None else np.asarray(gate).reshape(-1).astype(bool)
    order = np.asarray(order, dtype=np.int64).reshape(-1)

    code_for_offset = {
        (int(D8_DI[d]), int(D8_DJ[d])): int(D8_CODES[d]) for d in range(8)
    }

    for k in order:
        k = int(k)
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k:
            continue
        offset = (k // nx - r // nx, k % nx - r % nx)
        code = code_for_offset.get(offset, 0)
        if code == 0 or codes_flat[r] != code:
            continue
        if gate_flat is not None and not gate_flat[k]:
            continue
        if mode == "sum":
            out[r] += out[k]
        else:
            out[r] = out[k]

    return out.reshape(ny, nx)


# ---------------------------------------------------------------------------
# Traversal
# ---------------------------------------------------------------------------


def upstream_mask(recv, outlets):
    """Mask of every cell draining to any of ``outlets`` (flat indices)."""
    recv_arr = np.asarray(recv, dtype=np.int64)
    ny, nx = recv_arr.shape
    n = ny * nx
    recv_flat = recv_arr.reshape(-1)

    valid = (recv_flat >= 0) & (recv_flat < n) & (recv_flat != np.arange(n))
    donor_of = np.flatnonzero(valid)
    target_of = recv_flat[valid]
    sort = np.argsort(target_of, kind="stable")
    donors_sorted = donor_of[sort]
    starts = np.searchsorted(target_of[sort], np.arange(n), side="left")
    ends = np.searchsorted(target_of[sort], np.arange(n), side="right")

    mask = np.zeros(n, dtype=np.uint8)
    stack = []
    for k in np.asarray(outlets, dtype=np.int64).reshape(-1):
        k = int(k)
        if k < 0 or k >= n or mask[k]:
            continue
        mask[k] = 1
        stack.append(k)
    while stack:
        k = stack.pop()
        for donor in donors_sorted[starts[k]:ends[k]]:
            donor = int(donor)
            if mask[donor]:
                continue
            mask[donor] = 1
            stack.append(donor)
    return mask.reshape(ny, nx)


def chi(
    recv,
    order,
    area,
    step_length,
    outlets,
    A0,
    theta,
    trapezoid=False,
    max_length=0.0,
    mask=None,
):
    """Integrate the chi coordinate upstream from ``outlets``.

    ``chi(cell) = chi(receiver) + (A0 / A)^theta * step_length``

    With ``trapezoid=True`` the integrand is averaged between the cell and its
    receiver, which is the rule TopoToolbox's ``chitransform`` uses.  Cells
    outside the outlets' drainage areas stay ``NaN``.

    Returns ``(chi, distance_from_outlet)``.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    ny, nx = recv_arr.shape
    n = ny * nx
    recv_flat = recv_arr.reshape(-1)
    area_flat = np.asarray(area, dtype=np.float64).reshape(-1)
    step = np.asarray(step_length, dtype=np.float64).reshape(-1)
    order = np.asarray(order, dtype=np.int64).reshape(-1)
    mask_flat = None if mask is None else np.asarray(mask).reshape(-1).astype(bool)

    chi_flat = np.full(n, np.nan)
    dist = np.full(n, np.nan)
    for k in np.asarray(outlets, dtype=np.int64).reshape(-1):
        k = int(k)
        if k < 0 or k >= n:
            continue
        if mask_flat is not None and not mask_flat[k]:
            continue
        chi_flat[k] = 0.0
        dist[k] = 0.0

    # Receivers before donors: walk the topological order backwards.
    for k in order[::-1]:
        k = int(k)
        if not np.isnan(chi_flat[k]):
            continue
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k or np.isnan(chi_flat[r]):
            continue
        if mask_flat is not None and not mask_flat[k]:
            continue
        a = area_flat[k]
        if not a > 0.0:
            continue
        dl = step[k]
        if max_length > 0.0 and dist[r] + dl > max_length:
            continue
        integrand = (A0 / a) ** theta
        if trapezoid:
            ar = area_flat[r]
            integrand = 0.5 * (integrand + ((A0 / ar) ** theta if ar > 0.0 else integrand))
        chi_flat[k] = chi_flat[r] + integrand * dl
        dist[k] = dist[r] + dl

    return chi_flat.reshape(ny, nx), dist.reshape(ny, nx)


def downstream_distance(recv, order, step_length):
    """Along-flow distance from each cell to the end of its flow path.

    Zero at outlets and at cells with no receiver, increasing upstream --
    the same quantity as TopoToolbox's ``flowdistance(FD, 'upstream')``.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    shape = recv_arr.shape
    n = recv_arr.size
    recv_flat = recv_arr.reshape(-1)
    step = np.asarray(step_length, dtype=np.float64).reshape(-1)
    order = np.asarray(order, dtype=np.int64).reshape(-1)

    distance = np.zeros(n, dtype=np.float64)
    for k in order[::-1]:
        k = int(k)
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k:
            continue
        distance[k] = distance[r] + step[k]
    return distance.reshape(shape)


def imposemin(recv, order, step_length, sl, elevations):
    """Carve ``elevations`` in place so flow paths descend at gradient >= ``sl``.

    Only lowers cells, never raises them -- the complement of depression
    filling, and what TopoToolbox's ``imposemin`` does.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    n = recv_arr.size
    recv_flat = recv_arr.reshape(-1)
    step = np.asarray(step_length, dtype=np.float64).reshape(-1)
    order = np.asarray(order, dtype=np.int64).reshape(-1)
    if not isinstance(elevations, np.ndarray) or not elevations.flags['C_CONTIGUOUS']:
        # As in priority_flood: reshape(-1) copies a non-contiguous array, so
        # the carving would be thrown away.
        raise TypeError(
            "elevations must be a C-contiguous ndarray to be carved in place; "
            "pass np.ascontiguousarray(elevations)")
    z = elevations.reshape(-1)

    for k in order:
        k = int(k)
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k:
            continue
        candidate = z[k] - sl * step[k]
        if candidate < z[r]:
            z[r] = candidate


def stream_order(recv, order, is_stream, kind="strahler"):
    """Strahler or Shreve stream order over the cells marked in ``is_stream``."""
    recv_arr = np.asarray(recv, dtype=np.int64)
    shape = recv_arr.shape
    n = recv_arr.size
    recv_flat = recv_arr.reshape(-1)
    stream = np.asarray(is_stream).reshape(-1).astype(bool)
    order = np.asarray(order, dtype=np.int64).reshape(-1)

    result = np.zeros(n, dtype=np.int32)

    if kind == "shreve":
        has_donor = np.zeros(n, dtype=bool)
        for k in range(n):
            if not stream[k]:
                continue
            r = int(recv_flat[k])
            if 0 <= r < n and r != k and stream[r]:
                has_donor[r] = True
        result[stream & ~has_donor] = 1
        for k in order:
            k = int(k)
            if not stream[k]:
                continue
            r = int(recv_flat[k])
            if r < 0 or r >= n or r == k or not stream[r]:
                continue
            result[r] += result[k]
        return result.reshape(shape)

    # Each junction is resolved as a whole -- the largest incoming order and
    # how many donors carry it -- rather than donor by donor.  The
    # incremental form depends on the sequence the donors arrive in, which
    # would make the answer depend on which valid topological order was used.
    best = np.zeros(n, dtype=np.int32)
    ties = np.zeros(n, dtype=np.int32)

    for k in order:
        k = int(k)
        if not stream[k]:
            continue
        result[k] = best[k] + 1 if ties[k] >= 2 else max(int(best[k]), 1)

        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k or not stream[r]:
            continue
        if result[k] > best[r]:
            best[r] = result[k]
            ties[r] = 1
        elif result[k] == best[r]:
            ties[r] += 1

    return result.reshape(shape)


def drainage_basins(recv, order, outlets=None, valid=None):
    """Label drainage basins.

    With ``outlets`` given, each cell is labelled with the 1-based position of
    the outlet it drains to (0 if it reaches none).  Without, every terminal
    cell -- pit or edge outlet -- seeds its own basin, matching TopoToolbox's
    ``drainagebasins(FD)``.
    """
    recv_arr = np.asarray(recv, dtype=np.int64)
    ny, nx = recv_arr.shape
    n = ny * nx
    recv_flat = recv_arr.reshape(-1)
    order = np.asarray(order, dtype=np.int64).reshape(-1)
    valid_flat = None if valid is None else np.asarray(valid).reshape(-1).astype(bool)

    labels = np.zeros(n, dtype=np.int32)
    if outlets is None:
        # Only a terminal cell that something actually drains into seeds a
        # basin; an isolated cell stays 0, as in TopoToolbox.
        has_donor = np.zeros(n, dtype=bool)
        for k in range(n):
            if valid_flat is not None and not valid_flat[k]:
                continue
            r = int(recv_flat[k])
            if 0 <= r < n and r != k:
                has_donor[r] = True
        next_label = 0
        for k in range(n):
            if valid_flat is not None and not valid_flat[k]:
                continue
            r = recv_flat[k]
            if (r < 0 or r >= n or r == k) and has_donor[k]:
                next_label += 1
                labels[k] = next_label
    else:
        for position, k in enumerate(np.asarray(outlets, dtype=np.int64).reshape(-1), start=1):
            k = int(k)
            if 0 <= k < n and (valid_flat is None or valid_flat[k]):
                labels[k] = position

    for k in order[::-1]:
        k = int(k)
        if labels[k] != 0:
            continue
        if valid_flat is not None and not valid_flat[k]:
            continue
        r = int(recv_flat[k])
        if r < 0 or r >= n or r == k:
            continue
        labels[k] = labels[r]

    return labels.reshape(ny, nx)
