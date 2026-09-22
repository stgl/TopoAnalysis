#!/usr/bin/env python3
"""Regenerate ``docs/api.md`` from the live docstrings.

    python docs/generate_api.py

Run it after adding or renaming anything public.  The output is checked in
so the reference is readable on the repository page without building Sphinx.
"""

from __future__ import annotations

import inspect
import os
import pathlib
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

from TopoAnalysis import (MovingWindow, analysis, cli, datasets, dem,  # noqa: E402
                          demRecursionTools, error, fastops, kernels,
                          plotting, topotoolbox, utils)

#: Grid classes, grouped the way the guide introduces them.
GROUPS = [
    ('Core grid types', ['BaseSpatialGrid', 'Georef_info', 'BaseSpatialShape',
                         'ValueGrid', 'Elevation', 'GeographicElevation']),
    ('Depression filling', ['PriorityQueueMixIn', 'FilledElevation', 'PriorityFillGrid']),
    ('Flow routing', ['FlowDirection', 'FlowDirectionD8', 'FlowLength',
                      'GeographicFlowLength', 'Mask', 'DiscreteFlowAccumulation',
                      'GeographicDiscreteFlowAccumulation']),
    ('Drainage area', ['Area', 'GeographicArea', 'LogArea', 'ValleyArea',
                       'GeographicValleyArea', 'MainstemValleyArea',
                       'GeographicMainstemValleyArea']),
    ('Chi and channel steepness', ['Chi', 'GeographicChi', 'Ksi', 'GeographicKsi',
                                   'Relief', 'ScaledRelief', 'ChiScaledRelief',
                                   'ChannelSlope', 'CrossDivideDChi',
                                   'NormalizedCrossDivideDChi', 'RestoredElevation',
                                   'GeographicRestoredElevation']),
    ('Along-flow smoothing', ['AlongFlowSmoothing', 'KsFromChiWithSmoothing',
                              'GeographicKsFromChiWithSmoothing',
                              'ThetaFromChiWithSmoothing',
                              'GeographicThetaFromChiWithSmoothing',
                              'ChannelSlopeWithSmoothing',
                              'ChannelDownSlopeWithSmoothing',
                              'ChannelUpSlopeWithSmoothing']),
    ('Terrain attributes', ['CalculationMixin', 'Hillshade', 'GeographicHillshade',
                            'MaxSlope', 'GeographicMaxSlope', 'Laplacian',
                            'GeographicLaplacian', 'Gradient', 'LocalRelief',
                            'MultiscaleCurvatureValleyWidth', 'ScarpWavelet']),
    ('Isostasy', ['Deflection', 'GeographicDeflection']),
    ('Mixins', ['GDALMixin', 'GeographicGridMixin', 'MaxFlowLengthTrackingMixin']),
]


def first_line(obj):
    doc = inspect.getdoc(obj) or ''
    return doc.strip().split('\n')[0] if doc else ''


def signature(obj):
    """A call signature with the bound first argument removed."""
    try:
        sig = str(inspect.signature(obj))
    except (TypeError, ValueError):
        return '(...)'
    return (sig.replace('self, ', '').replace('(self)', '()')
               .replace('cls, ', '').replace('(cls)', '()'))


def esc(text):
    """Escape the one character that breaks a Markdown table."""
    return text.replace('|', '\\|')


def function_table(module, names=None):
    rows = []
    for name in sorted(names or getattr(module, '__all__', [])):
        obj = getattr(module, name, None)
        if not inspect.isfunction(obj):
            continue
        rows.append('| `{0}{1}` | {2} |'.format(
            name, esc(signature(obj)), esc(first_line(obj))))
    return (['| Function | Summary |', '|---|---|'] + rows + ['']) if rows else []


def class_section(out, cls, name):
    w = out.append
    bases = [b.__name__ for b in cls.__bases__ if b is not object]
    w('#### `{0}`\n'.format(name))
    if bases:
        w('*inherits* `{0}`\n'.format('`, `'.join(bases)))
    doc = inspect.getdoc(cls)
    if doc:
        w(doc + '\n')

    combos = getattr(cls, 'required_inputs_and_actions', None)
    if combos and 'required_inputs_and_actions' in vars(cls):
        w('*Constructor keywords* (any one combination)\n')
        for keys, _action in combos:
            w('- `{0}`'.format('`, `'.join(keys)))
        w('')

    dtype = vars(cls).get('dtype')
    if dtype is not None:
        w('*Stored as* `{0}`\n'.format(np.dtype(dtype).name))

    rows = []
    for mname, member in sorted(vars(cls).items()):
        if mname.startswith('_'):
            continue
        if isinstance(member, property):
            rows.append((mname, '', first_line(member.fget)))
            continue
        target = member.__func__ if isinstance(member, (classmethod, staticmethod)) else member
        if not inspect.isfunction(target):
            continue      # class constants such as dtype, reported above
        kind = ''
        if isinstance(member, classmethod):
            kind = ' *(classmethod)*'
        elif isinstance(member, staticmethod):
            kind = ' *(staticmethod)*'
        rows.append((mname, signature(target), first_line(target) + kind))

    if rows:
        w('| Member | Signature | Summary |')
        w('|---|---|---|')
        for mname, sig, summary in rows:
            w('| `{0}` | `{1}` | {2} |'.format(mname, esc(sig), esc(summary)))
        w('')


def build():
    out = []
    w = out.append

    w('# API reference\n')
    w('Generated from the docstrings by `docs/generate_api.py`, so')
    w('`help(TopoAnalysis.Area)` in a session shows the same text. For how the')
    w('pieces fit together start with [guide.md](guide.md); for what the')
    w('computations actually do see [algorithms.md](algorithms.md).\n')
    w('Attributes prefixed with a single underscore (`_griddata`,')
    w('`_georef_info`) are not private in practice -- they are how you reach')
    w('the data -- but they are not listed here; see')
    w('[guide.md](guide.md#the-grid-model).\n')

    w('## `TopoAnalysis.dem` -- the grid classes\n')
    w('Contents:\n')
    for title, names in GROUPS:
        w('- **{0}**: {1}'.format(title, ', '.join('`{0}`'.format(n) for n in names)))
    w('')

    for title, names in GROUPS:
        w('### {0}\n'.format(title))
        for name in names:
            cls = getattr(dem, name, None)
            if cls is not None:
                class_section(out, cls, name)

    w('### Module-level functions in `TopoAnalysis.dem`\n')
    out.extend(function_table(dem, ['mosaicFolder', 'plot', 'gdal_is_available']))

    w('## `TopoAnalysis.topotoolbox` -- TopoToolbox-equivalent functions\n')
    w(inspect.getdoc(topotoolbox).split('Example')[0].strip() + '\n')
    w('See [topotoolbox_parity.md](topotoolbox_parity.md) for the full')
    w('correspondence table and the places where exact agreement is not')
    w('achievable.\n')
    out.extend(function_table(topotoolbox))
    for name in topotoolbox.__all__:
        obj = getattr(topotoolbox, name)
        if inspect.isfunction(obj):
            w('### `{0}{1}`\n'.format(name, signature(obj)))
            w((inspect.getdoc(obj) or '') + '\n')

    w('## `TopoAnalysis.kernels` -- backend dispatch\n')
    w(inspect.getdoc(kernels) + '\n')
    out.extend(function_table(kernels))

    w('## `TopoAnalysis.fastops` -- the reference NumPy kernels\n')
    w(inspect.getdoc(fastops).split('Conventions')[0].strip() + '\n')
    for name in fastops.__all__:
        obj = getattr(fastops, name)
        if inspect.isfunction(obj):
            w('### `{0}{1}`\n'.format(name, signature(obj)))
            w((inspect.getdoc(obj) or '') + '\n')

    for mod, title in ((demRecursionTools, 'TopoAnalysis.demRecursionTools'),
                       (plotting, 'TopoAnalysis.plotting'),
                       (datasets, 'TopoAnalysis.datasets'),
                       (utils, 'TopoAnalysis.utils'),
                       (analysis, 'TopoAnalysis.analysis'),
                       (cli, 'TopoAnalysis.cli'),
                       (MovingWindow, 'TopoAnalysis.MovingWindow')):
        w('## `{0}`\n'.format(title))
        doc = inspect.getdoc(mod)
        if doc:
            w(doc + '\n')
        names = sorted(n for n, o in vars(mod).items()
                       if not n.startswith('_') and inspect.isfunction(o)
                       and o.__module__ == mod.__name__)
        out.extend(function_table(mod, names))
        for cname in sorted(n for n, o in vars(mod).items()
                            if not n.startswith('_') and inspect.isclass(o)
                            and o.__module__ == mod.__name__):
            cls = getattr(mod, cname)
            w('### `{0}`\n'.format(cname))
            w((inspect.getdoc(cls) or '') + '\n')
            rows = ['| `{0}{1}` | {2} |'.format(mname, esc(signature(member)),
                                                esc(first_line(member)))
                    for mname, member in sorted(vars(cls).items())
                    if not mname.startswith('_') and inspect.isfunction(member)]
            if rows:
                w('| Method | Summary |')
                w('|---|---|')
                out.extend(rows)
                w('')

    w('## `TopoAnalysis.error`\n')
    for cname in ('Error', 'InputError', 'TransitionError'):
        cls = getattr(error, cname)
        w('### `{0}`\n'.format(cname))
        w((inspect.getdoc(cls) or '') + '\n')

    return '\n'.join(out) + '\n'


if __name__ == '__main__':
    target = pathlib.Path(__file__).parent / 'api.md'
    target.write_text(build())
    print('wrote {0} ({1} lines)'.format(target, build().count('\n')))
