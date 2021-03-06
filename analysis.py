import numpy as np
import rasterio

import matplotlib.pyplot as plt

from itertools import product


class Quadrats(object):

    def __init__(self, filename=None, data=None, band=1, dx=None, dy=None):
        if filename is not None:
            self.load_data(filename, band=band)
        elif data is not None:
            self.data = data
        else:
            self.data = None

        if dx is not None:
            self.make_quadrats(dx, dy)
        else:
            self.quadrats = []

    def load_data(self, filename, band=1):
        with rasterio.open(filename, 'r') as src:
            self.data = src.read(band)

    def make_quadrats(self, dx, dy=None):
        if dy is None:
            dy = dx
        ny, nx = self.data.shape
        rows = np.arange(0, ny-dy+1, step=dy)
        cols = np.arange(0, nx-dx+1, step=dx)

        # TODO: Implement offset quadrats
        # rows += row_offset
        # cols += col_offset
        # rows = rows[rows <= ny]
        # cols = cols[cols <= nx]

        self.quadrats = [self.data[r:r+dy, c:c+dx] for r in rows for c in cols]

    def map_quadrats(self, func, **kwargs):
        return [func(q, **kwargs) for q in self.quadrats]

    def plot(self, values, ax=None, **kwargs):
        ny, nx = self.data.shape
        dy, dx = self.quadrats[0].shape
        rows = np.arange(0, ny-dy+1, step=dy) + dy / 2
        cols = np.arange(0, nx-dx+1, step=dx) + dx / 2
        y = [p[0] for p in product(rows, cols)]
        x = [p[1] for p in product(rows, cols)]

        if ax is None:
            plt.figure()
            ax = plt.gca()

        ax.scatter(x, y, c=values, **kwargs)
        ax.set_xlim([0, nx])
        ax.set_ylim([0, ny])
        ax.invert_yaxis()

    def quiver(self, u, v, ax=None, **kwargs):
        ny, nx = self.data.shape
        dy, dx = self.quadrats[0].shape
        rows = np.arange(0, ny-dy+1, step=dy) + dy / 2
        cols = np.arange(0, nx-dx+1, step=dx) + dx / 2
        y = [p[0] for p in product(rows, cols)]
        x = [p[1] for p in product(rows, cols)]

        if ax is None:
            plt.figure()
            ax = plt.gca()

        ax.quiver(x, y, u, v, **kwargs)
        ax.set_xlim([0, nx])
        ax.set_ylim([0, ny])
        ax.invert_yaxis()
