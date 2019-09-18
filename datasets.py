"""
Utilities for generating synthetic data and loading examples
"""

import numpy as np

from scipy.signal import sawtooth

from dem import Elevation


def triangle_grid(ny, nx, width, amp=1, sig=0, slope_y=None):
    """
    Returns a synthetic landscape with triangular ridges and valleys

    Grid spacing is 1.

    Parameters
    ----------
        ny : int
            Number of rows
        nx : int
            Number of columns
        width : int
            Spacing between valleys in pixels
        amp : float, optional
            Relief between ridge and valleys. Default is 1.
        sig : float, optional
            Amplitude of Gaussian Noise applied to elevations. Default is 0.
        slope_y : float, optional
            Slope of grid in y direction. Default is None (no slope).

    Returns
    -------
        out_obj : Elevation
            Elevation object holding synthetic data
    """
    f = nx / width
    t = np.linspace(0, 1, nx)
    y = amp * sawtooth(2 * np.pi * f * t, width=0.5)
    triangle = np.tile(y, (ny, 1))
    triangle += sig * np.random.randn(ny, nx)
    triangle += np.abs(np.min(triangle))

    if slope_y:
        y = np.linspace(0, ny, num=ny).reshape(ny, 1)
        X = np.tile(y, (1, nx))
        tilt = slope_y * X
        triangle *= tilt

    out_obj = Elevation(nx = nx, ny = ny, dx = 1.0)
    out_obj._griddata = triangle
    return out_obj


def sinusoid_grid(ny, nx, width, amp=1, sig=0, slope_y=None):
    """
    Returns a synthetic landscape with sinusoidal ridges and valleys

    Grid spacing is 1.

    Parameters
    ----------
        ny : int
            Number of rows
        nx : int
            Number of columns
        width : int
            Spacing between valleys in pixels
        amp : float, optional
            Relief between ridge and valleys. Default is 1.
        sig : float, optional
            Amplitude of Gaussian Noise applied to elevations. Default is 0.
        slope_y : float, optional
            Slope of grid in y direction. Default is None (no slope).

    Returns
    -------
        out_obj : Elevation
            Elevation object holding synthetic data
    """
    f = nx / width
    t = np.linspace(0, 1, nx)
    y = amp * np.sin(2 * np.pi * f * t)
    sinusoid = np.tile(y, (ny, 1))
    sinusoid += sig * np.random.randn(ny, nx)
    sinusoid += np.abs(np.min(sinusoid))

    if slope_y:
        y = np.linspace(0, ny, num=ny).reshape(ny, 1)
        X = np.tile(y, (1, nx))
        tilt = slope_y * X
        sinusoid *= tilt

    out_obj = Elevation(nx = nx, ny = ny, dx = 1.0)
    out_obj._griddata = sinusoid
    return out_obj
