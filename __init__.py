"""TopoAnalysis -- digital elevation model analysis for tectonic geomorphology.

A typical analysis chains grid classes together, each one built from the
last::

    from TopoAnalysis import Elevation, FilledElevation, FlowDirectionD8, Area

    dem    = Elevation(gdal_filename='srtm.tif')
    filled = FilledElevation(elevation=dem)
    d8     = FlowDirectionD8(flooded_dem=filled)
    area   = Area(flow_direction=d8)

Grids are constructed by keyword: which keywords you pass decides how the
grid is built.  ``Elevation(gdal_filename=...)`` reads a file,
``Elevation(dx=..., grid=...)`` wraps an array already in memory, and
``Area(flow_direction=...)`` derives one grid from another.  Each class
lists its accepted combinations in ``required_inputs_and_actions``.

Sub-modules
-----------
:mod:`TopoAnalysis.dem`
    The grid classes.  Everything public is re-exported here.
:mod:`TopoAnalysis.topotoolbox`
    Functions that reproduce TopoToolbox 2 (MATLAB) results exactly.
:mod:`TopoAnalysis.kernels`
    The compute kernels, in C++ where available and NumPy otherwise.
:mod:`TopoAnalysis.plotting`, :mod:`TopoAnalysis.demRecursionTools`
    Longitudinal and chi profiles, and steepness fitting along them.
:mod:`TopoAnalysis.datasets`
    Synthetic landscapes for testing.

Depression filling uses the Priority-Flood algorithms of Barnes, Lehman and
Mulla (2014); see :mod:`TopoAnalysis.kernels`.
"""

from __future__ import annotations

__version__ = "1.0.0"

# The package is importable both as ``TopoAnalysis`` and, with this directory
# on sys.path, flat -- the way the original scripts in this repository load
# it.  The same shim appears in every module here.
try:  # pragma: no cover - exercised by whichever import style is in use
    from . import error, fastops, gdal_setup, kernels  # noqa: F401
    from . import dem
except ImportError:  # pragma: no cover
    # The relative form failed, so this file was executed without a parent
    # package.  That happens whenever something imports it by path under a
    # name other than ``TopoAnalysis`` -- notably pytest, which names the
    # module after the checkout *directory*, and so only gets a valid package
    # name when the directory happens to be called TopoAnalysis.
    #
    # Falling back to flat imports only works if this directory is on the
    # path, which used to mean "only when it is also the working directory".
    # Put it there explicitly instead.  Appended rather than prepended: these
    # names only need to be findable, and nothing here should ever take
    # precedence over the standard library or an installed package.
    import os as _os
    import sys as _sys

    _here = _os.path.dirname(_os.path.abspath(__file__))
    if _here not in _sys.path:
        _sys.path.append(_here)
    del _os, _sys, _here

    import error, fastops, gdal_setup, kernels  # noqa: F401
    import dem

#: Everything :mod:`TopoAnalysis.dem` exports through this package.  Listed
#: explicitly rather than filtered out of ``dir(dem)`` so that a name
#: disappearing from ``dem`` breaks the import here, loudly, instead of
#: leaving ``__all__`` promising something that is not there.
_FROM_DEM = (
    "Area",
    "AlongFlowSmoothing",
    "BaseSpatialGrid",
    "BaseSpatialShape",
    "CalculationMixin",
    "ChannelDownSlopeWithSmoothing",
    "ChannelSlope",
    "ChannelSlopeWithSmoothing",
    "ChannelUpSlopeWithSmoothing",
    "Chi",
    "ChiScaledRelief",
    "CrossDivideDChi",
    "Deflection",
    "DiscreteFlowAccumulation",
    "Elevation",
    "FilledElevation",
    "FlowDirection",
    "FlowDirectionD8",
    "FlowLength",
    "GDALMixin",
    "GeographicArea",
    "GeographicChi",
    "GeographicDeflection",
    "GeographicDiscreteFlowAccumulation",
    "GeographicElevation",
    "GeographicFlowLength",
    "GeographicGridMixin",
    "GeographicHillshade",
    "GeographicKsFromChiWithSmoothing",
    "GeographicKsi",
    "GeographicLaplacian",
    "GeographicMainstemValleyArea",
    "GeographicMaxSlope",
    "GeographicRestoredElevation",
    "GeographicThetaFromChiWithSmoothing",
    "GeographicValleyArea",
    "Georef_info",
    "Gradient",
    "Hillshade",
    "Ksi",
    "KsFromChiWithSmoothing",
    "Laplacian",
    "LocalRelief",
    "LogArea",
    "MainstemValleyArea",
    "Mask",
    "MaxFlowLengthTrackingMixin",
    "MaxSlope",
    "MultiscaleCurvatureValleyWidth",
    "NormalizedCrossDivideDChi",
    "PriorityFillGrid",
    "PriorityQueueMixIn",
    "Relief",
    "RestoredElevation",
    "ScaledRelief",
    "ScarpWavelet",
    "ThetaFromChiWithSmoothing",
    "ValleyArea",
    "ValueGrid",
    "gdal_is_available",
    "mosaicFolder",
    "plot",
)

for _name in _FROM_DEM:
    globals()[_name] = getattr(dem, _name)
del _name


def backend() -> str:
    """Which kernel backend is in use: ``'c++'`` or ``'python'``."""
    return kernels.backend()


__all__ = sorted(_FROM_DEM) + [
    "__version__",
    "backend",
    "dem",
    "error",
    "fastops",
    "gdal_setup",
    "kernels",
]
