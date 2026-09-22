"""Convenience wrappers for fitting channel steepness at a single outlet."""

import numpy as np

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import demRecursionTools
except ImportError:  # pragma: no cover
    import demRecursionTools


def calc_ks_for_outlet(outlet, theta, **kwargs):
    """Fit channel steepness to the basin draining to ``outlet``.

    Parameters
    ----------
    outlet : (float, float)
        Map coordinates of the outlet.
    theta : float
        Reference concavity.
    flow_direction : FlowDirectionD8
        Required, passed by keyword.
    xo : float, optional
        Reference length; the reference area is ``xo**2``.  Default 500.
    plot : float, optional
        If equal to ``theta``, also draw the chi-elevation plot.
    **kwargs
        Any other keyword is a grid to sample along the profiles, e.g.
        ``elevation=dem, area=area``.

    Returns
    -------
    (ndarray, float)
        ``(ks, r_squared)``.
    """
    kwargs = dict(kwargs)
    fd = kwargs.pop('flow_direction')
    # kwargs.pop('xo') without a default raised KeyError whenever xo was
    # left at its default.
    xo = kwargs.pop('xo', None)
    if xo is None:
        xo = 500.0
    plot = kwargs.pop('plot', None)

    de = fd._mean_pixel_dimension()
    ld_list = fd.map_values_to_recursive_list(outlet, **kwargs)

    ret = demRecursionTools.best_ks_with_r2_list(ld_list, de, np.array([theta]), xo=xo)

    if plot is not None and plot == theta:
        import matplotlib.pylab as plt

        e, c = demRecursionTools.chi_elevation(ld_list, de, np.array([theta]), xo=xo)
        plt.figure()
        plt.plot(c, e, 'k.')
        # Matplotlib 3.0 removed the explicit hold call; overlaying is now
        # the default, so no equivalent is needed here.
        (ks, _) = ret
        plt.plot([0, np.max(c)], [0, ks * np.max(c)], 'k-')

    return ret
