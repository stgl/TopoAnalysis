"""Build script for TopoAnalysis.

Almost all configuration lives in ``pyproject.toml``.  This file exists only
to describe the optional C++ extension, which cannot be expressed
declaratively.

The extension is marked ``optional``: if no C++ compiler is available, or the
build fails for any other reason, installation still succeeds and TopoAnalysis
falls back to the pure-NumPy kernels in ``TopoAnalysis.fastops``.  Results are
identical either way; only the speed differs.
"""

from __future__ import annotations

import os
import sys

from setuptools import Extension, setup

HERE = os.path.dirname(os.path.abspath(__file__))


def _pybind11_include_dirs():
    """Locate pybind11's headers, tolerating its absence."""
    try:
        import pybind11
    except ImportError:  # pragma: no cover - only hit on a broken build env
        return None
    return [pybind11.get_include()]


def _extensions():
    include_dirs = _pybind11_include_dirs()
    if include_dirs is None:
        sys.stderr.write(
            "TopoAnalysis: pybind11 not found; building without the compiled "
            "kernels.  Install pybind11 and reinstall for a large speed-up.\n"
        )
        return []

    if sys.platform == "win32":
        extra_compile_args = ["/O2", "/std:c++17", "/EHsc"]
        extra_link_args = []
    else:
        extra_compile_args = ["-O3", "-std=c++17", "-funroll-loops"]
        extra_link_args = []
        if sys.platform == "darwin":
            # Structured bindings and std::optional need a modern libc++.
            extra_compile_args.append("-mmacosx-version-min=10.14")
            extra_link_args.append("-mmacosx-version-min=10.14")

    return [
        Extension(
            "TopoAnalysis._topoanalysis",
            sources=[os.path.join("_cpp", "module.cpp")],
            # Without this a header-only edit would not trigger a rebuild.
            depends=[
                os.path.join("_cpp", "priority_flood.hpp"),
                os.path.join("_cpp", "flow_routing.hpp"),
            ],
            include_dirs=include_dirs + [os.path.join(HERE, "_cpp")],
            language="c++",
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            # Never fail the install because the compiler is unhappy.
            optional=True,
        )
    ]


setup(ext_modules=_extensions())
