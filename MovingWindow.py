"""Moving-window filters over a grid.

A :class:`MovingWindow` applies a reducing function to the cells inside a
window centred on every cell.  Subclasses choose the window shape; a
concrete class supplies the reduction as :attr:`~MovingWindow.function`::

    class WindowMean(CircularMovingWindow):
        function = staticmethod(np.mean)

    grid.apply_moving_window(WindowMean(window_radius=200.0))

Cells outside the grid, and no-data cells, are dropped from the window
rather than filled, so the function always sees real data.

For the common reductions (min, max, mean, median, standard deviation)
:func:`TopoAnalysis.topotoolbox.localtopography` is far faster, because it
uses SciPy's separable filters instead of a Python loop.
"""

from __future__ import annotations

import numpy as np

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import error
except ImportError:  # pragma: no cover
    import error


class MovingWindow(object):
    """Abstract base for moving-window filters.

    Parameters
    ----------
    window_radius : float
        Radius of the window in map units.  ``window_dimension`` is accepted
        as a synonym for backward compatibility.
    """

    #: The reduction applied to the values inside each window.
    function = None

    def __init__(self, *args, **kwargs):
        radius = kwargs.get('window_radius', kwargs.get('window_dimension'))
        if radius is None:
            raise error.InputError('Window radius',
                                   'window_radius is a required parameter')
        # Both names are kept: window_dimension is what the original
        # constructor took, window_radius is what _build_search_kernel reads.
        self.window_radius = radius
        self.window_dimension = radius

        if self.__class__ is MovingWindow:
            raise error.Error('MovingWindow is an abstract base class')

    def _build_search_kernel(self, dx):
        """Relative ``(row, col)`` offsets of the window, as flat arrays."""
        raise NotImplementedError

    def _adjust_kernel(self, row, col, grid, searchKernelRows, searchKernelCols):
        """Clip the window to the grid and drop no-data cells."""
        theseRows = searchKernelRows + row
        theseCols = searchKernelCols + col

        inside = ((theseCols < grid.shape[1]) & (theseRows < grid.shape[0])
                  & (theseRows >= 0) & (theseCols >= 0))

        goodRows = theseRows[inside]
        goodCols = theseCols[inside]

        valid = np.logical_not(np.isnan(grid[goodRows, goodCols]))

        return goodRows[valid], goodCols[valid]

    def apply_moving_window(self, grid, dx, dtype):
        """Run the window over every cell of ``grid``.

        Parameters
        ----------
        grid : ndarray
            2-D data; ``NaN`` is no-data.
        dx : float
            Cell size, used to convert the window radius to cells.
        dtype : numpy dtype
            Type of the output array.

        Returns
        -------
        ndarray
        """
        if self.function is None:
            raise error.Error(
                '{0} has no bound function'.format(self.__class__.__name__))

        outgrid = np.zeros_like(grid, dtype=dtype)

        # The kernel depends on the cell size, not on the data; the original
        # passed the grid itself here, so window_radius/dx divided a float by
        # an array.
        search_kernel_rows, search_kernel_cols = self._build_search_kernel(dx)

        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                # The original called self.__adjustKernel, which does not
                # exist -- the method is __adjust_kernel.
                these_rows, these_cols = self._adjust_kernel(
                    i, j, grid, search_kernel_rows, search_kernel_cols)
                if these_rows.size == 0:
                    outgrid[i, j] = np.nan if np.issubdtype(dtype, np.floating) else 0
                    continue
                outgrid[i, j] = self.function(grid[these_rows, these_cols])

        return outgrid


class RectangularMovingWindow(MovingWindow):
    """A square window of side ``2 * round(window_radius / dx) + 1`` cells."""

    def __init__(self, *args, **kwargs):
        super(RectangularMovingWindow, self).__init__(*args, **kwargs)
        if self.__class__ is RectangularMovingWindow:
            raise error.Error('RectangularMovingWindow has no bound function')

    def _build_search_kernel(self, dx):
        pxlRadius = int(round(self.window_radius / dx))
        relCoords = np.arange(1 + 2 * pxlRadius) - pxlRadius
        searchKernelRow, searchKernelCol = np.meshgrid(relCoords, relCoords)
        return searchKernelRow.flatten(), searchKernelCol.flatten()


class CircularMovingWindow(MovingWindow):
    """A disc-shaped window of radius ``round(window_radius / dx)`` cells."""

    def __init__(self, *args, **kwargs):
        super(CircularMovingWindow, self).__init__(*args, **kwargs)
        if self.__class__ is CircularMovingWindow:
            raise error.Error('CircularMovingWindow has no bound function')

    def _build_search_kernel(self, dx):
        pxlRadius = int(round(self.window_radius / dx))
        relCoords = np.arange(1 + 2 * pxlRadius) - pxlRadius
        searchKernelRow, searchKernelCol = np.meshgrid(relCoords, relCoords)
        # The original squared searchKernelRow twice, which made the window a
        # vertical band rather than a disc.
        dists = np.sqrt(searchKernelRow ** 2 + searchKernelCol ** 2)
        inside = dists <= pxlRadius

        return searchKernelRow[inside].flatten(), searchKernelCol[inside].flatten()
