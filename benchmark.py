"""Compare the compiled and pure-NumPy kernels on a synthetic DEM.

Run it with::

    python -m TopoAnalysis.benchmark            # 512 x 512
    python -m TopoAnalysis.benchmark 2048       # a bigger grid
    python -m TopoAnalysis.benchmark 2048 --skip-python
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import fastops, kernels
except ImportError:  # pragma: no cover
    import fastops
    import kernels


def synthetic_dem(n, seed=1):
    """A noisy dome: a radial network with plenty of small depressions."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n] / float(n)
    z = (500.0 * np.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.1)
         + 25.0 * np.sin(12 * np.pi * x) * np.sin(12 * np.pi * y))
    return z + rng.normal(0.0, 2.0, (n, n))


def _time(label, function, *args, **kwargs):
    start = time.perf_counter()
    result = function(*args, **kwargs)
    elapsed = time.perf_counter() - start
    print('    {0:<22s} {1:9.1f} ms'.format(label, elapsed * 1000.0))
    return result, elapsed


def run(module, z, cellsize):
    """Time the full pipeline with one backend; returns the totals."""
    total = 0.0

    filled = np.ascontiguousarray(z, dtype=np.float64)
    _, dt = _time('fill (epsilon)', module.priority_flood, filled,
                  mode='epsilon', epsilon=1e-9, cellsize=cellsize)
    total += dt

    codes, dt = _time('flow directions', module.flow_directions, filled,
                      cellsize=cellsize)
    total += dt

    receivers, dt = _time('receivers', module.receivers, codes)
    total += dt

    (order, _), dt = _time('topological order', module.topological_order, receivers)
    total += dt

    weights = np.full(z.shape, cellsize ** 2)
    area, dt = _time('accumulate', module.accumulate, receivers, order, weights)
    total += dt

    step = np.full(z.shape, cellsize)
    (length, main_stem), dt = _time('flow length', module.flow_length,
                                    receivers, order, step)
    total += dt

    outlets = np.flatnonzero(np.asarray(receivers).reshape(-1) < 0).astype(np.int64)
    _, dt = _time('chi', module.chi, receivers, order, np.asarray(area), step,
                  outlets, 1e6, 0.45, True, 0.0, None)
    total += dt

    print('    {0:<22s} {1:9.1f} ms'.format('TOTAL', total * 1000.0))
    return total


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='python -m TopoAnalysis.benchmark',
        description='Time the TopoAnalysis kernels on a synthetic DEM.')
    parser.add_argument('size', nargs='?', type=int, default=512,
                        help='grid side length in cells (default: 512)')
    parser.add_argument('--cellsize', type=float, default=30.0)
    parser.add_argument('--skip-python', action='store_true',
                        help='do not time the pure-NumPy kernels')
    args = parser.parse_args(argv)

    n = args.size
    z = synthetic_dem(n)
    print('{0} x {1} cells ({2:,} total), dx = {3}'.format(
        n, n, n * n, args.cellsize))

    if not kernels.have_extension():
        print('\nThe compiled extension is not built; only the NumPy kernels '
              'are available.\nBuild it with `pip install -e .` from the '
              'TopoAnalysis directory.')
        print('\nNumPy kernels:')
        run(fastops, z, args.cellsize)
        return 0

    print('\nC++ kernels:')
    cpp_total = run(kernels, z, args.cellsize)

    if args.skip_python:
        return 0
    if n > 1024:
        print('\nSkipping the NumPy kernels: they would take minutes at this '
              'size. Pass a smaller size to compare them.')
        return 0

    print('\nNumPy kernels:')
    python_total = run(fastops, z, args.cellsize)

    print('\nSpeed-up: {0:.0f}x'.format(python_total / cpp_total))
    return 0


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
