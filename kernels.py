"""Backend dispatch for the TopoAnalysis grid kernels.

Import the kernels from here rather than from :mod:`TopoAnalysis.fastops` or
the compiled extension directly::

    from TopoAnalysis import kernels
    kernels.priority_flood(z, mode="flat")

If the C++ extension built successfully it is used; otherwise the pure-NumPy
implementations in :mod:`TopoAnalysis.fastops` take over.  Both produce the
same numbers -- the test suite asserts it -- so the choice is purely about
speed.

Set the environment variable ``TOPOANALYSIS_PURE_PYTHON=1`` to force the
NumPy path, which is useful when debugging or comparing the two.
"""

from __future__ import annotations

import importlib.machinery
import importlib.util
import os
import warnings

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import fastops
except ImportError:  # pragma: no cover
    import fastops

__all__ = [
    "backend",
    "have_extension",
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
    "downstream_distance",
    "imposemin",
    "stream_order",
    "default_seeds",
    "D8_DI",
    "D8_DJ",
    "D8_DIST",
    "D8_CODES",
    "PRIORITY_FLOOD_REFERENCE",
]

PRIORITY_FLOOD_REFERENCE = (
    "Barnes, R., Lehman, C., Mulla, D. (2014). Priority-flood: An optimal "
    "depression-filling and watershed-labeling algorithm for digital elevation "
    "models. Computers & Geosciences 62, 117-127. "
    "doi:10.1016/j.cageo.2013.04.024"
)

# Geometry constants are always taken from the NumPy module: both backends use
# the same neighbour ordering, and having one definition avoids them drifting.
D8_DI = fastops.D8_DI
D8_DJ = fastops.D8_DJ
D8_DIST = fastops.D8_DIST
D8_CODES = fastops.D8_CODES

_FORCE_PYTHON = os.environ.get("TOPOANALYSIS_PURE_PYTHON", "").strip() not in ("", "0", "false", "False")

def _load_extension():
    """Import the compiled kernels, or return None if they were not built."""
    try:
        from . import _topoanalysis  # type: ignore[attr-defined]
        return _topoanalysis
    except ImportError:
        pass

    # Flat import (`import kernels` with this directory on sys.path): load the
    # extension from *this file's* directory by path.  A bare
    # `import _topoanalysis` would search all of sys.path and could pick up an
    # unrelated -- or stale, or ABI-incompatible -- module of the same name
    # left lying around somewhere else.
    here = os.path.dirname(os.path.abspath(__file__))
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        candidate = os.path.join(here, "_topoanalysis" + suffix)
        if not os.path.exists(candidate):
            continue
        spec = importlib.util.spec_from_file_location("_topoanalysis", candidate)
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        except ImportError:  # pragma: no cover - a broken build
            return None
        return module
    return None


_ext = None if _FORCE_PYTHON else _load_extension()


def have_extension() -> bool:
    """True when the compiled kernels are in use."""
    return _ext is not None


def backend() -> str:
    """``'c++'`` or ``'python'``."""
    return "c++" if _ext is not None else "python"


def warn_if_slow(n_cells: int, threshold: int = 250_000) -> None:
    """Warn once when a large grid is about to be processed in pure Python."""
    if _ext is None and n_cells >= threshold:
        warnings.warn(
            "TopoAnalysis is using its pure-Python kernels on a grid of "
            "{0:,} cells; this will be slow. Build the C++ extension with "
            "`pip install -e .` from the TopoAnalysis directory for a large "
            "speed-up.".format(n_cells),
            RuntimeWarning,
            stacklevel=3,
        )


default_seeds = fastops.default_seeds


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
    """Fill depressions in place; see :func:`TopoAnalysis.fastops.priority_flood`."""
    if _ext is None:
        return fastops.priority_flood(
            elevations,
            closed=closed,
            seeds=seeds,
            mode=mode,
            epsilon=epsilon,
            cellsize=cellsize,
            max_pit_depth=max_pit_depth,
            track_visited=track_visited,
            reference_algorithm=reference_algorithm,
        )
    return _ext.priority_flood(
        elevations,
        closed,
        seeds,
        mode,
        epsilon,
        cellsize,
        max_pit_depth,
        track_visited,
        reference_algorithm,
    )


def flow_directions(elevations, cellsize_grid=None, cellsize=1.0):
    """Steepest-descent D8 codes; see :func:`TopoAnalysis.fastops.flow_directions`."""
    if _ext is None:
        return fastops.flow_directions(elevations, cellsize_grid, cellsize)
    return _ext.flow_directions(elevations, cellsize_grid, cellsize)


def receivers(codes):
    """Receiver index per cell; see :func:`TopoAnalysis.fastops.receivers`."""
    if _ext is None:
        return fastops.receivers(codes)
    return _ext.receivers(codes)


def topological_order(recv):
    """``(order, n_cells_on_cycles)``; see :func:`TopoAnalysis.fastops.topological_order`."""
    if _ext is None:
        return fastops.topological_order(recv)
    return _ext.topological_order(recv)


def accumulate(recv, order, weights, gate=None):
    """Downstream accumulation; see :func:`TopoAnalysis.fastops.accumulate`."""
    if _ext is None:
        return fastops.accumulate(recv, order, weights, gate)
    return _ext.accumulate(recv, order, weights, gate)


def flow_length(recv, order, step_length):
    """``(length, from_codes)``; see :func:`TopoAnalysis.fastops.flow_length`."""
    if _ext is None:
        return fastops.flow_length(recv, order, step_length)
    return _ext.flow_length(recv, order, step_length)


def propagate_along_main_stem(recv, order, main_stem_from_codes, gate, values, mode):
    """See :func:`TopoAnalysis.fastops.propagate_along_main_stem`."""
    if _ext is None:
        return fastops.propagate_along_main_stem(
            recv, order, main_stem_from_codes, gate, values, mode
        )
    return _ext.propagate_along_main_stem(
        recv, order, main_stem_from_codes, gate, values, mode
    )


def upstream_mask(recv, outlets):
    """See :func:`TopoAnalysis.fastops.upstream_mask`."""
    if _ext is None:
        return fastops.upstream_mask(recv, outlets)
    return _ext.upstream_mask(recv, outlets)


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
    """``(chi, distance)``; see :func:`TopoAnalysis.fastops.chi`."""
    if _ext is None:
        return fastops.chi(
            recv, order, area, step_length, outlets, A0, theta, trapezoid, max_length, mask
        )
    return _ext.chi(
        recv, order, area, step_length, outlets, A0, theta, trapezoid, max_length, mask
    )


def drainage_basins(recv, order, outlets=None, valid=None):
    """See :func:`TopoAnalysis.fastops.drainage_basins`."""
    if _ext is None:
        return fastops.drainage_basins(recv, order, outlets, valid)
    return _ext.drainage_basins(recv, order, outlets, valid)


def downstream_distance(recv, order, step_length):
    """See :func:`TopoAnalysis.fastops.downstream_distance`."""
    if _ext is None:
        return fastops.downstream_distance(recv, order, step_length)
    return _ext.downstream_distance(recv, order, step_length)


def imposemin(recv, order, step_length, sl, elevations):
    """Carve in place; see :func:`TopoAnalysis.fastops.imposemin`."""
    if _ext is None:
        return fastops.imposemin(recv, order, step_length, sl, elevations)
    return _ext.imposemin(recv, order, step_length, sl, elevations)


def stream_order(recv, order, is_stream, kind="strahler"):
    """See :func:`TopoAnalysis.fastops.stream_order`."""
    if _ext is None:
        return fastops.stream_order(recv, order, is_stream, kind)
    return _ext.stream_order(recv, order, is_stream, kind)
