# Common tasks.  `make help` lists them.

PYTHON ?= python3
PIP    ?= $(PYTHON) -m pip

.PHONY: help bootstrap gdal gdal-doctor gdal-source install develop build test \
        test-python bench docs clean distclean lint

help:
	@echo "bootstrap    everything from a fresh clone, GDAL included"
	@echo "gdal         install GDAL (C library and Python bindings)"
	@echo "gdal-doctor  report what GDAL is installed and what is wrong"
	@echo "gdal-source  build GDAL from source, no root or package manager"
	@echo "install      install the package and build the C++ kernels"
	@echo "develop      editable install for working on the source"
	@echo "build        compile the C++ kernels in place"
	@echo "test         run the test suite"
	@echo "test-python  run the test suite against the pure-NumPy kernels"
	@echo "bench        compare the two backends on a synthetic DEM"
	@echo "docs         regenerate docs/api.md from the docstrings"
	@echo "clean        remove build artefacts"
	@echo "distclean    also remove the compiled extension"

# The one-command path: creates an environment, installs everything,
# installs GDAL, and checks the result by round-tripping a GeoTIFF.
bootstrap:
	$(PYTHON) install.py

# Run by path rather than as `python -m TopoAnalysis.gdal_setup`: this
# directory *is* the package, so the dotted form only resolves once the
# package is installed, and these targets need to work before that.
gdal:
	$(PYTHON) gdal_setup.py install

gdal-doctor:
	$(PYTHON) gdal_setup.py doctor

gdal-source:
	$(PYTHON) gdal_setup.py install --strategy source

install:
	$(PIP) install .
	@$(PYTHON) -c "import TopoAnalysis; print('backend:', TopoAnalysis.backend())"

develop:
	$(PIP) install -e ".[test]"
	@$(PYTHON) -c "import TopoAnalysis; print('backend:', TopoAnalysis.backend())"

build:
	$(PYTHON) setup.py build_ext --inplace

test:
	$(PYTHON) -m pytest -q

test-python:
	TOPOANALYSIS_PURE_PYTHON=1 $(PYTHON) -m pytest -q

bench:
	$(PYTHON) -m TopoAnalysis.benchmark

docs:
	$(PYTHON) docs/generate_api.py
	@echo "docs/api.md regenerated; the rest of docs/ is hand-written"

lint:
	$(PYTHON) -m compileall -q . >/dev/null && echo "all modules compile"

clean:
	rm -rf build dist *.egg-info .pytest_cache docs/_build
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.py[co]' -delete

distclean: clean
	rm -f _topoanalysis*.so _topoanalysis*.pyd
