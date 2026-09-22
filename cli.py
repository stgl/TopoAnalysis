"""Command-line entry point: ``topoanalysis-process``.

Runs the standard pipeline over a DEM and writes the derived grids as
GeoTIFFs beside it::

    topoanalysis-process srtm.tif --outdir results --theta 0.45 --a0 1e6

Use ``topoanalysis-process --info`` to report which compute backend is
active and whether GDAL is available.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

try:  # pragma: no cover - exercised by whichever import style is in use
    from . import dem as _dem
    from . import kernels
    from .__init__ import __version__
except ImportError:  # pragma: no cover
    import dem as _dem
    import kernels
    __version__ = "unknown"


#: Derived grids the pipeline can produce, in dependency order.
PRODUCTS = ('filled', 'flowdir', 'area', 'length', 'slope', 'hillshade',
            'curvature', 'chi', 'ksi', 'relief')


def _build_parser():
    parser = argparse.ArgumentParser(
        prog='topoanalysis-process',
        description='Derive flow routing and channel metrics from a DEM.')
    parser.add_argument('dem', nargs='?', help='input DEM readable by GDAL')
    parser.add_argument('-o', '--outdir', default='.',
                        help='directory for the output grids (default: .)')
    parser.add_argument('-p', '--prefix', default=None,
                        help="output filename prefix (default: the DEM's basename)")
    parser.add_argument('--products', default='filled,flowdir,area,length',
                        help='comma-separated subset of: ' + ','.join(PRODUCTS))
    parser.add_argument('--theta', type=float, default=0.45,
                        help='reference concavity for chi and ksi (default: 0.45)')
    parser.add_argument('--a0', type=float, default=1e6,
                        help='reference drainage area in m^2 (default: 1e6)')
    parser.add_argument('--min-area', type=float, default=1e6,
                        help='drainage area defining the channel network, m^2')
    parser.add_argument('--geographic', action='store_true',
                        help='treat the grid as lat/lon and use spherical cell areas')
    parser.add_argument('--flat-fill', action='store_true',
                        help='fill depressions to a flat surface (as TopoToolbox '
                             'does) instead of adding a small gradient')
    parser.add_argument('--info', action='store_true',
                        help='report the active backend and exit')
    parser.add_argument('--version', action='version',
                        version='TopoAnalysis {0}'.format(__version__))
    return parser


def _report_environment():
    print('TopoAnalysis {0}'.format(__version__))
    print('  compute backend : {0}'.format(kernels.backend()))
    if not kernels.have_extension():
        print('                    (build the C++ extension with `pip install -e .` '
              'for a large speed-up)')
    print('  GDAL available  : {0}'.format(_dem.gdal_is_available()))
    print('  fill algorithm  : {0}'.format(kernels.PRIORITY_FLOOD_REFERENCE))


def main(argv=None):
    """Run the CLI.  Returns a process exit status."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.info:
        _report_environment()
        return 0

    if not args.dem:
        parser.error('a DEM is required unless --info is given')

    requested = [p.strip() for p in args.products.split(',') if p.strip()]
    unknown = [p for p in requested if p not in PRODUCTS]
    if unknown:
        parser.error('unknown product(s): {0}'.format(', '.join(unknown)))

    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or os.path.splitext(os.path.basename(args.dem))[0]

    def out(name):
        return os.path.join(args.outdir, '{0}_{1}.tif'.format(prefix, name))

    elevation_cls = _dem.GeographicElevation if args.geographic else _dem.Elevation
    area_cls = _dem.GeographicArea if args.geographic else _dem.Area
    length_cls = _dem.GeographicFlowLength if args.geographic else _dem.FlowLength
    ksi_cls = _dem.GeographicKsi if args.geographic else _dem.Ksi

    started = time.time()

    def step(label):
        print('[{0:6.1f}s] {1}'.format(time.time() - started, label))
        sys.stdout.flush()

    step('reading {0}'.format(args.dem))
    elevation = elevation_cls(gdal_filename=args.dem)
    print('          {0} x {1} cells, dx = {2}'.format(
        elevation._georef_info.nx, elevation._georef_info.ny, elevation._georef_info.dx))

    step('filling depressions')
    fill_kwargs = {'elevation': elevation}
    if args.flat_fill:
        fill_kwargs['aggradation_slope'] = 0.0
    filled = _dem.FilledElevation(**fill_kwargs)
    if 'filled' in requested:
        filled.save(out('filled'))

    step('routing flow')
    flow_direction = _dem.FlowDirectionD8(flooded_dem=filled)
    if 'flowdir' in requested:
        flow_direction.save(out('flowdir'))

    area = None
    if {'area', 'chi', 'ksi', 'relief'} & set(requested):
        step('accumulating drainage area')
        area = area_cls(flow_direction=flow_direction)
        if 'area' in requested:
            area.save(out('area'))

    flow_length = None
    if {'length', 'ksi', 'relief'} & set(requested):
        step('measuring flow length')
        flow_length = length_cls(flow_direction=flow_direction)
        if 'length' in requested:
            flow_length.save(out('length'))

    if 'slope' in requested:
        step('computing slope')
        slope_cls = _dem.GeographicMaxSlope if args.geographic else _dem.MaxSlope
        slope_cls(elevation=elevation).save(out('slope'))

    if 'hillshade' in requested:
        step('computing hillshade')
        hs_cls = _dem.GeographicHillshade if args.geographic else _dem.Hillshade
        hs_cls(elevation=elevation, azimuth=315, inclination=45).save(out('hillshade'))

    if 'curvature' in requested:
        step('computing curvature')
        lap_cls = _dem.GeographicLaplacian if args.geographic else _dem.Laplacian
        lap_cls(elevation=elevation).save(out('curvature'))

    if 'chi' in requested:
        step('integrating chi')
        outlets = area.areas_greater_than(args.min_area)
        edge_outlets = [xy for xy in outlets
                        if not flow_direction.get_flow_to_cell(
                            *flow_direction._xy_to_rowscols((xy,))[0])[2]]
        _dem.Chi(area=area, flow_direction=flow_direction, theta=args.theta,
                 Ao=args.a0, outlets=edge_outlets or outlets).save(out('chi'))

    if 'ksi' in requested:
        step('integrating ksi')
        ksi_cls(area=area, flow_direction=flow_direction, flow_length=flow_length,
                theta=args.theta, Ao=args.a0).save(out('ksi'))

    if 'relief' in requested:
        step('computing scaled relief')
        _dem.ScaledRelief(flow_direction=flow_direction, elevation=elevation,
                          flow_length=flow_length, Ao=args.a0,
                          theta=args.theta, area=area).save(out('relief'))

    step('done')
    return 0


if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
