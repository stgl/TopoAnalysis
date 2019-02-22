import os

import numpy as np
import rasterio

import matplotlib.pyplot as plt

from dem import BaseSpatialGrid, ScarpWavelet


class Quadrats(object):


    def __init__(self):
        self.data = None
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
        
        rows += row_oofset
        cols += col_offset
        rows = rows[rows <= ny]
        cols = cols[cols <= nx]
        
        self.quadrats = [data[r:r+dy, c:c+dx] for r in rows for c in cols]


    def map_quadrats(self, func, **kwargs):
        return [func(q, **kwargs) for q in self.quadrats]


    def plot_map(self, attr, **kwargs):
        pass
