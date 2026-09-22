# A reproducible TopoAnalysis image.
#
#     docker build -t topoanalysis .
#     docker run --rm topoanalysis topoanalysis-gdal doctor
#     docker run --rm -v "$PWD:/data" topoanalysis \
#         topoanalysis-process /data/srtm.tif -o /data/out
#
# The base image already carries a matched GDAL C library and Python
# bindings, which is the one thing that is genuinely hard to arrange.
# ``install.py`` therefore finds GDAL already working and does nothing to it;
# everything else it does normally.
FROM ghcr.io/osgeo/gdal:ubuntu-small-3.9.3

# build-essential and python3-dev are for the Priority-Flood kernels; without
# them the package still installs and falls back to its NumPy kernels.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential python3-dev python3-pip python3-venv \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/TopoAnalysis
COPY . .

# --system because the image *is* the environment: there is nothing here to
# protect from a pip install, and a virtualenv would only hide the GDAL that
# the base image went to the trouble of providing.
#
# --break-system-packages because this base image is built on Ubuntu 24.04,
# whose Python carries the PEP 668 "externally managed" marker.  That marker
# exists to stop pip from fighting the host's package manager; inside a
# single-purpose container there is no host to protect, and overriding it is
# the normal and intended thing to do.
RUN python3 install.py --system --break-system-packages --extras test

# Fail the build rather than ship an image whose kernels silently did not
# compile or whose GDAL does not work.
RUN python3 -c "import sys, TopoAnalysis; \
        print('backend:', TopoAnalysis.backend()); \
        sys.exit(0 if TopoAnalysis.backend() == 'c++' else 1)" \
 && topoanalysis-gdal doctor

WORKDIR /data
CMD ["topoanalysis-process", "--info"]
