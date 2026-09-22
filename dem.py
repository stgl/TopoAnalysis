"""Raster grids and DEM analysis for tectonic geomorphology.

Originally developed by Sam Johnstone (January 2015); maintained since by the
Stanford geomorphology group.

Every raster is a :class:`BaseSpatialGrid` or one of its subclasses.  Grids
are constructed by keyword, and the keywords you supply select how the grid is
built -- see :attr:`BaseSpatialGrid.required_inputs_and_actions`::

    elevation = Elevation(gdal_filename='dem.tif')
    filled    = FilledElevation(elevation=elevation)
    d8        = FlowDirectionD8(flooded_dem=filled)
    area      = Area(flow_direction=d8)

Conventions
-----------
* ``_griddata`` is a 2-D NumPy array with row 0 at the **north** edge.
* ``NaN`` marks no-data in floating-point grids.
* D8 flow-direction codes follow the ArcGIS convention::

      | 32  64 128 |
      | 16   X   1 |
      |  8   4   2 |

* Depression filling uses the Priority-Flood algorithms of Barnes, Lehman and
  Mulla (2014); see :mod:`TopoAnalysis.kernels`.
"""

import glob  # Used for finding files that I want to mosaic
import os
import subprocess  # Used to run gdal_merge.py from the command line
import sys
import warnings

import numpy as np
from numpy import float64, uint8
from matplotlib import pyplot as plt

# These modules are importable both as a package (``from TopoAnalysis import
# dem``) and flat (``import dem`` with this directory on sys.path).  The flat
# form is how the original scripts in this repository load them, so it has to
# keep working.
try:  # pragma: no cover - exercised by whichever import style is in use
    from . import error as Error
    from . import kernels
except ImportError:  # pragma: no cover
    import error as Error
    import kernels


class _LazyGDALModule(object):
    """Defer ``from osgeo import ...`` until something actually touches GDAL.

    GDAL is only needed to read and write files.  Importing it lazily keeps
    the rest of the library usable -- and testable -- where GDAL is absent,
    which is common because it cannot be installed reliably from PyPI.  When
    it really is needed and really is missing, the error says what to do
    about it instead of surfacing as a bare ``ImportError`` at import time.
    """

    __slots__ = ("_submodule", "_module")

    def __init__(self, submodule):
        self._submodule = submodule
        self._module = None

    def _load(self):
        if self._module is None:
            # Point GDAL and PROJ at the data files of a prefix TopoAnalysis
            # built for itself, if there is one.  A no-op otherwise, and it
            # never overrides a setting the user made.
            try:
                from . import gdal_setup
            except ImportError:  # pragma: no cover - flat import style
                try:
                    import gdal_setup
                except ImportError:
                    gdal_setup = None
            if gdal_setup is not None:
                gdal_setup.activate_environment()

            try:
                from osgeo import gdal, ogr, osr
            except ImportError as exc:  # pragma: no cover - environment dependent
                raise ImportError(
                    "TopoAnalysis needs the GDAL Python bindings to read and "
                    "write raster files, and they are not installed.\n"
                    "\n"
                    "    topoanalysis-gdal install\n"
                    "\n"
                    "installs both the GDAL C library and the bindings, "
                    "choosing whichever route suits this machine; run "
                    "`topoanalysis-gdal doctor` first to see what it would "
                    "do. Grids built in memory (dx=..., grid=...) do not "
                    "require GDAL at all."
                ) from exc
            # Without this GDAL signals failure by returning None, which turns
            # a bad file path into an AttributeError three frames later.
            gdal.UseExceptions()
            ogr.UseExceptions()
            osr.UseExceptions()
            self._module = {"gdal": gdal, "ogr": ogr, "osr": osr}[self._submodule]
        return self._module

    def __getattr__(self, name):
        return getattr(self._load(), name)


gdal = _LazyGDALModule("gdal")
ogr = _LazyGDALModule("ogr")
osr = _LazyGDALModule("osr")


def gdal_is_available():
    """True when the GDAL Python bindings can be imported."""
    try:
        gdal._load()
    except ImportError:
        return False
    return True


def _require_statsmodels():
    """Import statsmodels on demand; only the regression grids need it."""
    try:
        import statsmodels.api as sm
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Fitting channel steepness needs statsmodels. Install it with "
            "`pip install statsmodels`."
        ) from exc
    return sm


SQRT2 = np.sqrt(2.0)


class _Progress(object):
    """Percent-complete ticker for the long per-cell regression loops.

    Prints ``Percent completion... 10...20...`` on one line and finishes with
    ``100``.  Pass ``verbose=False`` to silence it.
    """

    def __init__(self, total, verbose=True, stream=None):
        self.total = int(total)
        self.verbose = verbose and self.total > 0
        self.stream = stream if stream is not None else sys.stdout
        self.count = 0
        self.next_readout = 10
        if self.verbose:
            self.stream.write('Percent completion...')
            self.stream.flush()

    def tick(self):
        if not self.verbose:
            return
        self.count += 1
        # Integer arithmetic: accumulating 1.0/total in a float drifted and
        # could skip or repeat a readout.
        percent = (100 * self.count) // self.total
        while percent >= self.next_readout and self.next_readout < 100:
            self.stream.write(str(self.next_readout) + '...')
            self.stream.flush()
            self.next_readout += 10

    def done(self):
        if not self.verbose:
            return
        self.stream.write('100\n')
        self.stream.flush()


class GDALMixin(object):
    """Reading and writing rasters through GDAL.

    Mixed into :class:`BaseSpatialGrid`; nothing here is meant to be called
    directly.  Every method reaches GDAL through the module-level lazy
    proxies, so the mixin can be present without GDAL being installed.
    """

    
    def _get_projection_from_EPSG_projection_code(self, EPSGprojectionCode):
        # Get raster projection
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(EPSGprojectionCode)
        return srs.ExportToWkt()

    def _get_gdal_type_for_numpy_type(self, numpy_type):
    
        from numpy import float64, uint8, uint16, int16, uint32, int32, float32, complex64
        from osgeo.gdal import GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32, GDT_Float32, GDT_Float64, GDT_CFloat64, GDT_Unknown
        
        type_map = { uint8: GDT_Byte,
                uint16: GDT_UInt16,
                int16: GDT_Int16,
                uint32: GDT_UInt32,
                int32: GDT_Int32,
                float32: GDT_Float32,
                float64: GDT_Float64,
                complex64: GDT_CFloat64 }
        
        gdal_type = type_map.get(numpy_type)
        
        if gdal_type is None:
            return GDT_Unknown
        else:
            return gdal_type
    
    def _get_numpy_type_for_gdal_type(self, gdal_type):
        
        from numpy import float64, uint8, uint16, int16, uint32, int32, float32, complex64
        from osgeo.gdal import GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32, GDT_Float32, GDT_Float64, GDT_CFloat64
                    
        type_map = { GDT_Byte: uint8,
                GDT_UInt16: uint16,
                GDT_Int16: int16,
                GDT_UInt32: uint32,
                GDT_Int32: int32,
                GDT_Float32: float32,
                GDT_Float64: float64,
                GDT_CFloat64: complex64}
        
        numpy_type = type_map.get(gdal_type)
        
        if numpy_type is None:
            return float64
        else:
            return numpy_type
    
    def _readGDALFile(self, filename, dtype):
        gdal_file = gdal.Open(filename)
        try:
            return self._read_GDAL_dataset(gdal_file, dtype)
        finally:
            # GDAL only flushes and releases the file handle when the last
            # reference goes away, so drop it explicitly rather than waiting
            # for the garbage collector.
            gdal_file = None

    @staticmethod
    def _read_band(band, dtype):
        """Read one band, converting its no-data value to ``NaN``.

        Integer grids have no NaN to convert to, so their no-data value is
        left as-is.
        """
        data = band.ReadAsArray().astype(dtype)
        nodata = band.GetNoDataValue()
        if nodata is not None and np.issubdtype(np.dtype(dtype), np.floating):
            if np.isnan(nodata):
                pass  # already NaN in the array
            else:
                data[data == nodata] = np.nan
        return data

    def _read_GDAL_dataset(self, gdal_dataset, dtype):
        band = gdal_dataset.GetRasterBand(1)
        data = self._read_band(band, dtype)
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        return geoTransform, nx, ny, data
    
    def _getGeoRefInfo(self, gdalDataset):
        #Get info needed to initialize new dataset
        nx = gdalDataset.RasterXSize
        ny = gdalDataset.RasterYSize
    
        #Write geographic information
        geoTransform = gdalDataset.GetGeoTransform()  # Steal the coordinate system from the old dataset
        projection = gdalDataset.GetProjection()  # Steal the Projections from the old dataset
    
        return nx, ny, projection, geoTransform
 
    def getDEMcoords(self, GdalData, dx):
        """Cell-centre x and y coordinate lists for an open GDAL dataset."""
    
        #Get grid size
        nx, ny = GdalData.RasterXSize, GdalData.RasterYSize
    
        #Get information about the spatial reference
        (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = GdalData.GetGeoTransform()
        xllcenter = upper_left_x + dx/2.0  # x coordinate center of lower left pxl
        yllcenter = upper_left_y - (ny-0.5)*dx # y coordinate center of lower left pxl
    
        #Create arrays of the x and y coordinates of each pixel (the axes)
        xcoordinates = [x*dx + xllcenter for x in range(nx)]
        ycoordinates = [y*dx + yllcenter for y in range(ny)][::-1] #Flip the ys so that the first row corresponds to the first entry of this array
    
        return xcoordinates, ycoordinates

    def _create_gdal_representation_from_array(self, georef_info, GDALDRIVERNAME, array_data, dtype,
                                               outfile_path='name', dst_options=None,
                                               multiple_bands=False, verbose=False):
        """Write a NumPy array out as a georeferenced GDAL dataset.

        Driver names are listed at https://gdal.org/drivers/raster/.  The
        returned dataset is still open; let it go out of scope (or set it to
        ``None``) to flush it to disk.
        """
        if dst_options is None:
            dst_options = []

        drvr = gdal.GetDriverByName(GDALDRIVERNAME)
        if drvr is None:
            raise Error.InputError(
                'GDAL driver', "no GDAL driver named '{0}'".format(GDALDRIVERNAME))

        bands = len(array_data) if multiple_bands else 1
        outRaster = drvr.Create(outfile_path, georef_info.nx, georef_info.ny, bands,
                                self._get_gdal_type_for_numpy_type(dtype), dst_options)

        outRaster.SetGeoTransform(georef_info.geoTransform)
        # An unset projection is stored as the integer 0 by Georef_info.
        if georef_info.projection != 0 and georef_info.projection is not None:
            outRaster.SetProjection(georef_info.projection)

        is_float = np.issubdtype(np.dtype(dtype), np.floating)
        band_arrays = array_data if multiple_bands else [array_data]
        for index, band_array in enumerate(band_arrays):
            if verbose:
                print('writing band {0}/{1}'.format(index + 1, bands))
            band = outRaster.GetRasterBand(index + 1)
            if is_float:
                # Recording NaN as the no-data value means other GIS software
                # renders the holes as holes instead of as real elevations.
                band.SetNoDataValue(float('nan'))
            band.WriteArray(np.asarray(band_array))

        return outRaster
   
    def _clipRasterToRaster(self, input_gdal_dataset, clipping_gdal_dataset, dtype):

        # Source
        src_proj = input_gdal_dataset.GetProjection()
        src_geotrans = input_gdal_dataset.GetGeoTransform()
    
        # We want a section of source that matches this:
        match_proj = clipping_gdal_dataset.GetProjection()
        match_geotrans = clipping_gdal_dataset.GetGeoTransform()
        wide = clipping_gdal_dataset.RasterXSize
        high = clipping_gdal_dataset.RasterYSize
    
        # Output / destination
        dst = gdal.GetDriverByName('MEM').Create('name', wide, high, 1, dtype)
        dst.SetGeoTransform( match_geotrans )
        dst.SetProjection( match_proj)
    
        # Do the work
        gdal.ReprojectImage(input_gdal_dataset, dst, src_proj, match_proj, gdal.GRA_Bilinear)
        # gdal.ReprojectImage(src, dst, None, None, GRA_Bilinear)
        # gdal.ReprojectImage(dst, src, None, None, GRA_Bilinear)
    
        #
        return dst
 
    def _clipRasterToShape(self, raster, shape):
        
        # TODO: This needs implementation to write out the raster and shape files, execute the warp, read in the resulting file, and delete the filenames.  What a hack.
        
        pass
        # drvr = gdal.GetDriverByName(GdalDriver)
        # drvr.Create(outputFilename,1,1,1)
    #    warp= 'gdalwarp -cutline \'%s\' -crop_to_cutline -dstalpha \'%s\' \'%s\'' % (shpFilename, srcFilename, outputFilename)
        
        #warp= 'gdalwarp -cutline \'%s\' -crop_to_cutline \'%s\' \'%s\'' % (shpFilename, srcFilename, outputFilename)
    
        #os.system(warp)
    
    def _convertToUTM(self, dataset, dx, utmZone):

        #Get Spatial reference info
        oldRef = osr.SpatialReference()  # Initiate a spatial reference
    
        oldRef.ImportFromWkt(dataset.GetProjectionRef())  # Clone the spatial reference from the dataset
    
        newRef = osr.SpatialReference()
        newRef.SetUTM(abs(utmZone), utmZone > 0)
    
        #Set up the transform
        transform = osr.CoordinateTransformation(oldRef, newRef) # Create the coordinate transform object
        tVect = dataset.GetGeoTransform()  # Get the coordinate transform vector
        nx, ny = dataset.RasterXSize, dataset.RasterYSize  # Size of the original raster
        (ulx, uly, ulz ) = transform.TransformPoint(tVect[0], tVect[3])
        (lrx, lry, lrz ) = transform.TransformPoint(tVect[0] + tVect[1]*nx, tVect[3] + tVect[5]*ny)
        memDrv = gdal.GetDriverByName('MEM')  # Create a gdal driver in memory
        dataOut = memDrv.Create('name', int((lrx - ulx)/dx), int((uly - lry)/dx), 1, gdal.GDT_Float32)
        newtVect = (ulx, dx, tVect[2], uly, tVect[4], -dx)
    
    
        dataOut.SetGeoTransform(newtVect)  # Set the new geotransform
        dataOut.SetProjection(newRef.ExportToWkt())
        # Perform the projection/resampling
        res = gdal.ReprojectImage(dataset, dataOut, oldRef.ExportToWkt(), newRef.ExportToWkt(), gdal.GRA_Cubic)
    
        return dataOut
    
    def _getRasterGeoTransformFromAsciiRaster(self, fileName):
        #Read in the components of the geotransform from the raster, BEWARE! This
        #is specific to how some matlab/ C scripts I have write these rasters. I believe
        #this is that standard arc ascii raster export format, but could be wrong
        georef_data = dict()
        
        with open(fileName, "r") as ascii_file:
            for _ in range(6):
                line = ascii_file.readline()
                (key, value) = (line.split()[0], float(line.split()[-1]))
                georef_data[key.lower()] = value

        required_values = ('ncols', 'nrows', 'cellsize')
        
        if len(set(required_values).difference(set(georef_data.keys()))) != 0:
            raise Error.InputError('A/I ASCII grid error','The following properties are missing: ' + str(set(required_values).difference(set(georef_data.keys()))))
        
        if georef_data.get('xllcorner') is None and georef_data.get('xllcenter') is None:
            raise Error.InputError('A/I ASCII grid error','Neither XLLCorner nor XLLCenter is present.')
        
        if georef_data.get('yllcorner') is None and georef_data.get('yllcenter') is None:
            raise Error.InputError('A/I ASCII grid error','Neither YLLCorner nor YLLCenter is present.')

        dx = georef_data.get('cellsize')
        nx = int(georef_data.get('ncols'))
        ny = int(georef_data.get('nrows'))
        
        if georef_data.get('xllcenter') is not None:
            xUL = georef_data.get('xllcenter') - (dx/2.0)
        else:
            xUL = georef_data.get('xllcorner');
        
        if georef_data.get('yllcenter') is not None:
            yUL = (georef_data.get('yllcenter') - (dx/2.0)) + dx*ny
        else:
            yUL = georef_data.get('yllcorner') + dx*ny
    
        return (xUL, dx, 0, yUL, 0, -dx), nx, ny
    
    def _writeArcAsciiRaster(self, georef_info, outfile_path, np_array_data, nodata_value,
                             format_string):
        """Write an ArcInfo ASCII grid.

        GDAL's AAIGrid driver cannot ``Create()`` -- only ``CreateCopy()`` --
        so the file is written directly.  ``NaN`` cells are replaced by
        ``nodata_value`` on the way out, because most GIS software does not
        recognise the literal text ``nan`` in an ASCII grid.
        """
        data = np.asarray(np_array_data)
        if np.issubdtype(data.dtype, np.floating):
            if nodata_value is None or (isinstance(nodata_value, float) and np.isnan(nodata_value)):
                nodata_value = -9999.0
            data = np.where(np.isnan(data), nodata_value, data)

        header = "ncols     %s\n" % georef_info.nx
        header += "nrows    %s\n" % georef_info.ny
        header += "xllcenter %s\n" % georef_info.xllcenter
        header += "yllcenter %s\n" % georef_info.yllcenter
        header += "cellsize %s\n" % georef_info.dx
        header += "NODATA_value %s" % nodata_value

        np.savetxt(outfile_path, data, header=header, fmt=format_string, comments='')

    def _asciiRasterToMemory(self, fileName):
        """Read an ArcInfo ASCII grid, honouring its ``NODATA_value``."""
        # The geotransform is (xUL, dx, skewX, yUL, skewY, -dy).  It is parsed
        # from the header rather than taken from GDAL so that grids written
        # with xllcenter (rather than xllcorner) land in the right place.
        gt, nx, ny = self._getRasterGeoTransformFromAsciiRaster(fileName)

        ds = gdal.Open(fileName)
        try:
            # Ignoring NODATA_value here is what used to turn -9999 flags into
            # real elevations several thousand metres below sea level.
            data = self._read_band(ds.GetRasterBand(1), self.dtype)
        finally:
            ds = None

        return gt, nx, ny, data

class GeographicGridMixin(object):
    """Spherical geometry for grids in degrees of latitude and longitude.

    Mix it in *before* the grid class to override the cell area and cell
    dimension with their true, latitude-dependent values::

        class GeographicArea(GeographicGridMixin, Area):
            pass

    Every ``Geographic*`` class in this module is exactly that.
    """

    
    def _getUTMZone(self, dataset):
        #Function to get the approximate UTM zone %NOTE: I need to check how east and west are handled...
    
        #Utm zone boundary (zones are numbered in order, 1:60) #NEED TO DOUBLE CHECK THIS
        westBound = np.array([-180 + x*6 for x in range(60)]) #west boundary of 6 degree UTM zone bounds
        eastBound = np.array([-174 + x*6 for x in range(60)]) #east boundary of 6 degree UTM zone bounds
    
        #Midpoint of dataset
        tVect = dataset.GetGeoTransform()  # Get the coordinate transform vector, (ulx, dx, xRot, uly, yRot, -dx)
        nx, ny = dataset.RasterXSize, dataset.RasterYSize #Get the number of colums and rows
        midLat = tVect[3]-tVect[1]*ny/2.0 #half way down the dataset
        midLong = tVect[0]+tVect[1]*nx/2.0 #half way across the dataset
    
        #Convert UTM zone to negative to distinguish it as south (if appropriate)
        southMultiplier = 1
        if midLat < 0:
            southMultiplier = -1
    
        #The utm zone, the index of the boundaries that surround the point incremented to account for pythons 0 indexing
        zone = np.nonzero(np.logical_and(midLong > westBound, midLong < eastBound))[0] + 1
    
        return zone*southMultiplier

    def _getLatsLongsFromGeoTransform(self, geoTransform, nx, ny):
        dLong = geoTransform[1]
        dLat = geoTransform[5]
    
        #Determine the location of the center of the first pixel of the data
        xllcenter = geoTransform[0]+dLong/2.0
        yllcenter = geoTransform[3]+dLat/2.0
    
        #Assign the latitudinal (i.e. y direction) coordinates
        lats = np.zeros(ny)
        for i in range(len(lats)):
            lats[i] = yllcenter + i*dLat
    
        #Assign the longitudinal (i.e. x direction) coordinates
        longs = np.zeros(nx)
        for i in range(len(longs)):
            longs[i] = xllcenter + i*dLong
    
        return lats,longs
    
    def _approximateDxFromGeographicData(self, geoTransform):
        #Function to return the approximate grid spacing in Meters. Will return the closest integer value
        dTheta = geoTransform[1] #Angular grid spacing
        metersPerDegree = 110000 #110 km per degree (approximation)
        return int(dTheta*metersPerDegree) #convert degrees to meters
    
    
    #: Earth radius in metres used for spherical cell areas.  TopoToolbox's
    #: ``cellarea`` instead derives its radius from a total surface area of
    #: 510 072 000 km^2, i.e. 6 371 047 m, which makes its areas larger by
    #: about 1.5e-5 relative.  Override this attribute to match exactly.
    earth_radius = 6371.0 * 1000.0

    def _area_per_pixel(self, *args, **kwargs):
        """Area of each cell in square metres, on a sphere.

        Cells shrink towards the poles, so this varies with row but not with
        column.  The result is cached because the trigonometry is not free
        and several derived grids ask for it.
        """
        if getattr(self, '_GeographicGridMixin__app', None) is not None:
            return self._GeographicGridMixin__app

        re = self.earth_radius

        dLong = np.abs(self._georef_info.geoTransform[1])
        dLat = np.abs(self._georef_info.geoTransform[5])

        # np.float64(range(n)) is a TypeError on modern NumPy -- np.float64 is
        # a scalar type, not an array constructor.
        lats = self._georef_info.yllcenter + np.arange(self._georef_info.ny, dtype=float64) * dLat
        longs = self._georef_info.xllcenter + np.arange(self._georef_info.nx, dtype=float64) * dLong

        [LONG, LAT] = np.meshgrid(longs, lats)
        LAT1 = np.radians(LAT - dLat / 2.0)
        LAT2 = np.radians(LAT + dLat / 2.0)
        LONG1 = np.radians(LONG - dLong / 2.0)
        LONG2 = np.radians(LONG + dLong / 2.0)

        # lats ascends northwards; flip so row 0 is the north edge again.
        initialAreas = np.flipud(np.abs((re**2)*(LONG1 - LONG2)*(np.sin(LAT2) - np.sin(LAT1))))

        self.__app = initialAreas
        return self.__app

    def _mean_pixel_dimension(self, *args, **kwargs):
        """Square root of the cell area: an equivalent isotropic cell size."""
        return np.sqrt(self._area_per_pixel())


class CalculationMixin(object):
    """Finite-difference terrain derivatives.

    Supplies the slope, curvature and boundary-padding helpers that
    :class:`Elevation`, :class:`Hillshade` and :class:`MaxSlope` share.
    """

    
    def _calcFiniteSlopes(self, grid, dx, nx, ny):
        # sx,sy = calcFiniteDiffs(elevGrid,dx)
        # calculates finite differences in X and Y direction using the
        # 2nd order/centered difference method.
        # Applies a boundary condition such that the size and location
        # of the grids in is the same as that out.
    
        # Assign boundary conditions
        
        Zbc = self.assignBCs(grid, nx, ny)
    
        #Compute finite differences
        Sx = (Zbc[1:-1, 2:] - Zbc[1:-1, :-2])/(2*dx)
        Sy = (Zbc[2:,1:-1] - Zbc[:-2, 1:-1])/(2*dx)
    
        return Sx, Sy
    
    def assignBCs(self, grid, nx, ny):
        """Pad a grid by one cell, replicating the edge values.

        Returns an array two rows and two columns larger, so that a centred
        difference can be taken everywhere including the perimeter.
        """
        # Pads the boundaries of a grid
        # Boundary condition pads the boundaries with equivalent values
        # to the data margins, e.g. x[-1,1] = x[1,1]
        # This creates a grid 2 rows and 2 columns larger than the input
    
        Zbc = np.zeros((ny + 2, nx + 2))  # Create boundary condition array
        Zbc[1:-1,1:-1] = grid  # Insert old grid in center
    
        #Assign boundary conditions - sides
        Zbc[0, 1:-1] = grid[0, :]
        Zbc[-1, 1:-1] = grid[-1, :]
        Zbc[1:-1, 0] = grid[:, 0]
        Zbc[1:-1, -1] = grid[:,-1]
    
        #Assign boundary conditions - corners
        Zbc[0, 0] = grid[0, 0]
        Zbc[0, -1] = grid[0, -1]
        Zbc[-1, 0] = grid[-1, 0]
        Zbc[-1, -1] = grid[-1, -1]

        return Zbc

    def calcFiniteCurv(self, grid, dx):
        """Laplacian from a 5-point centred stencil, same shape as the input."""
        #C = calcFiniteCurv(elevGrid, dx)
        #calculates finite differnces in X and Y direction using the centered difference method.
        #Applies a boundary condition such that the size and location of the grids in is the same as that out.
    
        #Assign boundary conditions
        Zbc = self.assignBCs(grid,grid.shape[1],grid.shape[0])
    
        #Compute finite differences
        Cx = (Zbc[1:-1, 2:] - 2*Zbc[1:-1, 1:-1] + Zbc[1:-1, :-2])/dx**2
        Cy = (Zbc[2:, 1:-1] - 2*Zbc[1:-1, 1:-1] + Zbc[:-2, 1:-1])/dx**2
    
        return Cx+Cy
    
    def calcContourCurvature(self, grid,dx):
        """Contour (plan) curvature: the curvature of the contour lines.

        Positive where contours are convex downslope, i.e. on noses; negative
        in hollows.  The one-cell border comes back as ``NaN``.
        """
        # kt = (fxx*fy^2 - 2*fxyfxfy + fyy*fx^2)/((fx^2 + fy^2)*sqrt((fx^2 + fy^2)+1)
    
        #Preallocate
        Kt = np.zeros_like(grid)*np.nan
    
        #First derivatives, 2nd order centered difference
        fx = (grid[1:-1,2:] - grid[1:-1,:-2])/(dx*2)
        fy = (grid[2:,1:-1] - grid[:-2,1:-1])/(dx*2)
    
        #Second derivatives, 2nd order centered differece
        fxx = (grid[1:-1,2:] - 2*grid[1:-1,1:-1] + grid[1:-1,:-2])/(dx**2)
        fyy = (grid[2:,1:-1] - 2*grid[1:-1,1:-1] + grid[:-2,1:-1])/(dx**2);
    
        #Partial derivative
        fxy = (grid[2:,2:] - grid[2:,1:-1] - grid[1:-1,2:] + 2*grid[1:-1,1:-1] - grid[:-2,1:-1] - grid[1:-1,:-2] + grid[:-2,:-2])
        fxy = fxy/(4*dx**2)
    
        #Contour curvature
        Kt[1:-1, 1:-1] = (fxx*fy**2 - 2*fxy*fx*fy + fyy*fx**2)/((fx**2 + fy**2)*np.sqrt((fx**2 + fy**2)+1))
    
        return Kt

    def calcAverageSlopeOfGridSubset(self, gridSubset,dx):
        """Least-squares plane through a block of cells; returns its slopes.

        Returns ``(dz/dx, dz/dy)`` of the best-fit plane ``z = ax + by + c``.
        """
        ## Sx,Sy = calcAverageSlopeOfGridSubset(numpy matrix, dx)
        #Compute the average slope over a subset of a grid (or a whole grid if you're into that),
        #by fitting a plane to the elevation data stored in grid subset
        ny, nx = gridSubset.shape
        xs = (0.5+np.arange(nx))*dx - (nx*dx)/2.0
        ys = (0.5+np.arange(ny))*dx - (ny*dx)/2.0
        X,Y = np.meshgrid(xs,ys)
        #Fit a plane of the form z = ax + by + c, where beta = [a b c]
        M=np.vstack((X.flatten(),Y.flatten(),np.ones(nx*ny))).T
        beta = np.linalg.lstsq(M, gridSubset.flatten(), rcond=None)[0]

        return beta[0], beta[1] #Return the slope in the x and y directions respectively

class Georef_info(object):
    """Where a grid sits on the ground.

    Attributes
    ----------
    geoTransform : tuple
        The GDAL 6-tuple ``(x_origin, dx, 0, y_origin, 0, -dy)``, whose
        origin is the *outer corner* of the north-west cell.
    projection : str or int
        WKT, or the integer 0 when no projection is set.
    dx : float
        Cell size in map units.  Cells are assumed square.
    xllcenter, yllcenter : float
        Centre of the south-west cell.
    nx, ny : int
        Columns and rows.
    """

    def __init__(self):
        self.geoTransform = 0
        self.projection = 0
        self.dx = 0
        self.xllcenter = 0
        self.yllcenter = 0
        self.nx = 0
        self.ny = 0
        
    def __deepcopy__(self, memo):
        import copy
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))
        return result
            
    
class BaseSpatialShape(object):
    """A vector layer, used to cut rasters to a boundary.

    ::

        shape = BaseSpatialShape(shapefile_name='basin.shp')
        mask = shape.createMaskFromShape(dem._georef_info,
                                         dem._georef_info.projection,
                                         gdal.GDT_Byte)
    """

    def __init__(self, *args, **kwargs):
        if kwargs.get('shapefile_name') is None:
            raise Error.InputError('Input Error', 'Inputs not satisfied')
        self.shapedata = ogr.Open(kwargs.get('shapefile_name'))

    def createMaskFromShape(self, geoRefInfo, projection, dtype, noDataValue = 0):
        """Rasterise the shape onto ``geoRefInfo``'s frame, 1 inside and 0 out."""
    
        #Open Shapefile
        source_ds = self.shapedata
        source_layer = source_ds.GetLayer()
    
        maskSrc = gdal.GetDriverByName('MEM').Create('name',geoRefInfo.nx,geoRefInfo.ny, 1, dtype)
        maskSrc.SetGeoTransform(geoRefInfo.geoTransform)
        maskBand = maskSrc.GetRasterBand(1)
        maskBand.SetNoDataValue(noDataValue)
    
        # 5. Rasterize why is the burn value 0... isn't that the same as the background?
        gdal.RasterizeLayer(maskSrc, [1], source_layer, burn_values=[1])
        grid = maskSrc.ReadAsArray().astype(dtype)
        maskSrc = None
        return BaseSpatialGrid(nx = geoRefInfo.nx, ny = geoRefInfo.ny, projection = projection, geo_transform = geoRefInfo.geoTransform, grid= grid)

    
class BaseSpatialGrid(GDALMixin):
    """A georeferenced raster, and the base class of every grid here.

    Construct one by keyword; the combination you pass selects how it is
    built (see :attr:`required_inputs_and_actions`)::

        BaseSpatialGrid(gdal_filename='grid.tif')
        BaseSpatialGrid(dx=30.0, grid=array)
        BaseSpatialGrid(nx=100, ny=80, dx=30.0)          # random values
        BaseSpatialGrid(nx=100, ny=80, projection=wkt, geo_transform=gt)

    The data live in ``_griddata``, a 2-D NumPy array with row 0 at the
    **north** edge and ``NaN`` for no-data; the georeferencing lives in
    ``_georef_info`` (a :class:`Georef_info`).

    Indexing is by ``(row, col)`` and is forgiving: reading outside the grid
    returns ``None`` and writing outside it is ignored, which is what the
    neighbourhood traversals in this module rely on.
    """

    #: Keyword combinations this class accepts, each paired with the method
    #: that builds the grid from them.  Subclasses extend the list; the
    #: entries here are available to every grid class.
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('nx', 'ny', 'dx'), '_create_random_grid'),
                                   (('dx', 'grid'), '_create_from_grid'))
    dtype = float64

    # BaseSpatialGrid is a class for all rasters
    
    def __deepcopy__(self, memo):
        import copy
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))
        return result
        
    def __init__(self, *args, **kwargs):

        super(BaseSpatialGrid,self).__init__()
        self._georef_info = Georef_info()
        self._sorted = False

        if len(kwargs.keys()) == 0:
            return

        evaluative_action = self.__get_evaluative_action(*args, **kwargs)

        if evaluative_action is None:
            raise Error.InputError(
                'Input Error',
                '{0} could not be built from keywords {1}. Supported '
                'combinations are: {2}'.format(
                    self.__class__.__name__,
                    sorted(kwargs),
                    '; '.join(', '.join(keys)
                              for keys, _ in self.required_inputs_and_actions)))

        # Resolve the method by name rather than by building and executing a
        # source string: the action comes from class data, but compiling text
        # to call a method is both slower and a needless way to run code.
        getattr(self, evaluative_action)(*args, **kwargs)

    def __setitem__(self, key, value):
        """Set a cell, ignoring writes that fall outside the grid.

        Out-of-range writes are silently dropped rather than raising: several
        of the neighbourhood traversals in this module rely on being able to
        write past an edge without checking first.
        """
        i, j = key
        if i < 0 or j < 0 or i >= self._georef_info.ny or j >= self._georef_info.nx:
            return
        self._griddata[i, j] = value

    def __getitem__(self, key):
        """Return a cell value, or ``None`` if the index is outside the grid.

        The ``None`` return -- rather than an ``IndexError`` -- is part of the
        interface: neighbourhood code throughout this module tests the result
        against ``None`` instead of bounds-checking every access.
        """
        i, j = key

        i = int(i)
        j = int(j)

        if 0 <= i < self._georef_info.ny and 0 <= j < self._georef_info.nx:
            return self._griddata[i, j]

        return None

    def __get_evaluative_action(self, *args, **kwargs):
        """Pick the constructor that the supplied keywords satisfy."""
        these_kw = set(kwargs.keys())

        # Subclass-specific combinations win, then the generic ones every
        # grid supports.  Falling back to the base list is what lets
        # `Mask(dx=..., grid=...)` and friends be built in memory without
        # every subclass having to restate the common entries.
        for source in (self.required_inputs_and_actions,
                       BaseSpatialGrid.required_inputs_and_actions):
            for required_input_set, evaluative_action in source:
                if set(required_input_set).issubset(these_kw):
                    return evaluative_action

        return None
    
    def __populate_georef_info_using_geoTransform(self, nx, ny):
        """Derive dx and the south-west cell centre from the geotransform."""
        # nx and ny are assigned FIRST: yllcenter depends on ny, and reading
        # self._georef_info.ny before setting it picked up the default of 0,
        # which put the grid ny*dx too far north for anyone building one with
        # nx/ny/projection/geo_transform.
        self._georef_info.nx = nx
        self._georef_info.ny = ny
        self._georef_info.dx = self._georef_info.geoTransform[1]
        self._georef_info.xllcenter = self._georef_info.geoTransform[0] + self._georef_info.dx / 2.0
        self._georef_info.yllcenter = (self._georef_info.geoTransform[3]
                                       - self._georef_info.dx * (ny - 0.5))
        
    def _create(self, *args, **kwargs):
            
        self._georef_info.geoTransform = kwargs.get('geo_transform')
        self._georef_info.projection = kwargs.get('projection')
        self.__populate_georef_info_using_geoTransform(kwargs['nx'], kwargs['ny'])
        if kwargs.get('grid') is None:
            self._griddata = np.zeros(shape = (self._georef_info.ny,self._georef_info.nx), dtype = self.dtype)
        else:
            self._griddata = kwargs.get('grid')

    def _geoTransform_from_corner(self):
        """Derive the GDAL geotransform from dx and the lower-left centre.

        Grids built in memory used to leave ``geoTransform`` at its default
        of 0, which made ``save()`` fail with an unhelpful GDAL error.
        """
        dx = self._georef_info.dx
        return (self._georef_info.xllcenter - 0.5 * dx,
                dx,
                0.0,
                self._georef_info.yllcenter + (self._georef_info.ny - 0.5) * dx,
                0.0,
                -dx)

    def _create_random_grid(self, *args, **kwargs):
        self._georef_info.dx = kwargs['dx']
        self._georef_info.nx = kwargs['nx']
        self._georef_info.ny = kwargs['ny']
        self._georef_info.xllcenter = kwargs.get('xllcenter', 0.0)
        self._georef_info.yllcenter = kwargs.get('yllcenter', 0.0)
        self._griddata = np.zeros((self._georef_info.ny, self._georef_info.nx),
                                  dtype=self.dtype)
        self._randomize_grid_values(*args, **kwargs)
        self._georef_info.geoTransform = self._geoTransform_from_corner()

    def _create_from_grid(self, *args, **kwargs):
        self._georef_info.dx = kwargs['dx']
        grid = np.asarray(kwargs['grid'])
        (ny, nx) = grid.shape
        self._georef_info.nx = nx
        self._georef_info.ny = ny
        self._georef_info.xllcenter = kwargs.get('xllcenter', 0.0)
        self._georef_info.yllcenter = kwargs.get('yllcenter', 0.0)
        self._griddata = grid
        self._georef_info.geoTransform = self._geoTransform_from_corner()
        if kwargs.get('projection') is not None:
            self._georef_info.projection = kwargs['projection']

    def _randomize_grid_values(self, *args, **kwargs):
        if kwargs.get('mask') is not None:
            i = np.where(kwargs.get('mask')._griddata == 1)
            self._griddata[i] = np.random.rand(len(i[0]))
        else:
            self._griddata = np.random.rand(self._georef_info.ny, self._georef_info.nx)
            
    def _read_ai(self, *args, **kwargs):
        
        self._georef_info.projection = self._get_projection_from_EPSG_projection_code(kwargs['EPSGprojectionCode'])
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._asciiRasterToMemory(kwargs['ai_ascii_filename'])
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)
    
    def _read_gdal(self, *args, **kwargs):
        
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._readGDALFile(kwargs['gdal_filename'], self.dtype)
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)
    
    def _copy_info_from_grid(self, grid, set_zeros = False):
        import copy
        self._georef_info = copy.deepcopy(grid._georef_info)
        if not set_zeros:
            self._griddata = grid._griddata.copy()
        else:
            self._griddata = np.zeros_like(grid._griddata, self.dtype)
                
    def _getNeighborIndices(self, row, col):
        #Search kernel for D8 flow routing, the relative indices of each of the 8 points surrounding a pixel
        # |i-1,j-1  i-1,j  i-1,j+1|
        # |i,j-1     i,j     i,j+1|
        # |i+1,j-1  i+1,j  i+1,j+1|
        rowKernel = np.array([1, 1, 1, 0, 0, -1, -1, -1])
        colKernel = np.array([-1, 0, 1, -1, 1, -1, 0, 1])
    
        rt2 = np.sqrt(2)
        dxMults = np.array([rt2, 1.0, rt2, 1.0, 1.0, rt2, 1.0, rt2])  # Unit Distance from pixel to surrounding coordinates
    
        #Find all the surrounding indices
        outRows = (rowKernel + row).astype(int)
        outCols = (colKernel + col).astype(int)
    
        #Determine which indices are out of bounds
        inBounds = (outRows >= 0)*(outRows < self._georef_info.ny)*(outCols >= 0)*(outCols < self._georef_info.nx)
        return (outRows[inBounds], outCols[inBounds], dxMults[inBounds])
    
    def _xy_to_rowscols(self, v):
        l = list()
        for (x,y) in v:

            col = int(round((x-self._georef_info.xllcenter)/self._georef_info.dx))
            row = int((self._georef_info.ny - 1) - round((y-self._georef_info.yllcenter)/self._georef_info.dx))

            # ``>=``, not ``>``: a point one cell past the right or bottom
            # edge used to come back as a valid subscript and then raise an
            # IndexError somewhere else entirely.
            if col >= self._georef_info.nx or row >= self._georef_info.ny or col < 0 or row < 0:
                l.append((None, None))
            else:
                l.append((row,col))
        return tuple(l)
    
    def _rowscols_to_xy(self, l):
        v = list()
        for(row,col) in l:
            x = float64(col)*self._georef_info.dx + self._georef_info.xllcenter
            y = (float64(self._georef_info.ny - 1.0) - float64(row))*self._georef_info.dx + self._georef_info.yllcenter
            v.append((x,y))
        return tuple(v)
    
    def _area_per_pixel(self, *args, **kwargs):
        return self._georef_info.dx**2 * np.ones((self._georef_info.ny, self._georef_info.nx))

    def _mean_pixel_dimension(self, *args, **kwargs):
        return self._georef_info.dx * np.ones_like(self._griddata, dtype = float64)

    def coordinate_vectors(self):
        """``(x, y)`` coordinates of the cell centres, west-to-east and
        south-to-north.

        Returns
        -------
        (ndarray, ndarray)
            ``x`` has ``nx`` entries and ``y`` has ``ny`` entries.  Note that
            ``y`` ascends northwards, i.e. ``y[0]`` is the *last* row of
            ``_griddata``.
        """
        info = self._georef_info
        x = info.xllcenter + np.arange(info.nx) * info.dx
        y = info.yllcenter + np.arange(info.ny) * info.dx
        return x, y

    def get_XY_matricies(self):
        """Meshgrid of cell-centre coordinates.

        The name's spelling is preserved for backward compatibility;
        :meth:`get_XY_matrices` is the same method.
        """
        # np.arange(start, start + (n-1)*dx, dx) yields n-1 points, not n, so
        # the returned grids used to be one column and one row short.
        x, y = self.coordinate_vectors()
        return np.meshgrid(x, y)

    #: Correctly spelled alias for :meth:`get_XY_matricies`.
    get_XY_matrices = get_XY_matricies

    def resample(self, de, interpolation='quintic'):
        """Resample onto a grid of spacing ``de``.

        Parameters
        ----------
        de : float
            Target cell size, in the grid's own units.
        interpolation : str
            ``'nearest'``, ``'linear'``, ``'cubic'`` or ``'quintic'``.  The
            last two use spline interpolation of order 3 and 5.

        Returns
        -------
        BaseSpatialGrid
            A new grid of the same class.
        """
        import copy
        from scipy.interpolate import RegularGridInterpolator

        order = {'nearest': 0, 'linear': 1, 'cubic': 3, 'quintic': 5}.get(interpolation)
        if order is None:
            raise Error.InputError(
                'resample', "interpolation must be one of 'nearest', 'linear', "
                            "'cubic' or 'quintic'; got {0!r}".format(interpolation))

        info = self._georef_info
        georef = copy.deepcopy(info)
        georef.dx = de
        georef.nx = max(int(np.round(info.nx * info.dx / de)), 1)
        georef.ny = max(int(np.round(info.ny * info.dx / de)), 1)
        # Keep the south-west cell centre fixed, then rebuild the transform
        # from it -- the old code kept the source transform's origin and
        # merely swapped in the new spacing, which shifted the whole grid.
        georef.xllcenter = info.xllcenter
        georef.yllcenter = info.yllcenter
        georef.geoTransform = (georef.xllcenter - 0.5 * de,
                               de,
                               info.geoTransform[2] if info.geoTransform != 0 else 0.0,
                               georef.yllcenter + (georef.ny - 0.5) * de,
                               info.geoTransform[4] if info.geoTransform != 0 else 0.0,
                               -de)

        # Interpolate in south-to-north order so the y axis is increasing,
        # which is what RegularGridInterpolator requires.
        source = np.flipud(np.asarray(self._griddata, dtype=float64))
        xo = info.xllcenter + np.arange(info.nx) * info.dx
        yo = info.yllcenter + np.arange(info.ny) * info.dx

        method = {0: 'nearest', 1: 'linear', 3: 'cubic', 5: 'quintic'}[order]
        try:
            interpolator = RegularGridInterpolator(
                (yo, xo), source, method=method, bounds_error=False, fill_value=None)
        except ValueError:
            # 'cubic' and 'quintic' need SciPy >= 1.10 and enough points per
            # axis.  Fall back rather than fail, but say so: silently
            # changing the interpolation order changes the answer.
            warnings.warn(
                "SciPy cannot do {0!r} interpolation on this grid (it needs "
                "SciPy >= 1.10 and at least {1} points per axis); falling "
                "back to 'linear'.".format(interpolation, order + 1),
                RuntimeWarning, stacklevel=2)
            interpolator = RegularGridInterpolator(
                (yo, xo), source, method='linear', bounds_error=False, fill_value=None)

        xf = georef.xllcenter + np.arange(georef.nx) * de
        yf = georef.yllcenter + np.arange(georef.ny) * de
        YF, XF = np.meshgrid(yf, xf, indexing='ij')
        resampled = interpolator(np.column_stack((YF.ravel(), XF.ravel())))
        resampled = resampled.reshape(georef.ny, georef.nx)

        return_grid = self.__class__()
        return_grid._georef_info = georef
        return_grid._griddata = np.flipud(resampled).astype(self.dtype)
        return return_grid

    def location_in_grid(self, xo):
        """True when ``xo = (x, y)`` lands on a cell holding real data.

        Real means: inside the grid, not ``NaN``, and not zero.
        """
        index = self._xy_to_rowscols((xo,))[0]
        return index[0] is not None and index[1] is not None and self[index[0], index[1]] is not None and self[index[0], index[1]] != 0 and not np.isnan(self[index[0], index[1]])
    
    def tile(self, tile_xdim = 400, tile_ydim=400, tile_xpadding = 10, tile_ypadding = 10):
        """Split the grid into overlapping tiles.

        Interior tiles overlap their neighbours by the padding width; tiles on
        the outside of the grid are mirrored outwards instead, so every tile
        has the same amount of context around its core.  Analyses that need a
        neighbourhood (filling, curvature at a scale) can then be run
        tile-by-tile and reassembled with :meth:`remove_padding` followed by
        :meth:`mosaic`.

        Returns
        -------
        list
            Tiles, each an independent grid of this class.
        """
        dx = self._georef_info.dx
        grid_nx, grid_ny = self._georef_info.nx, self._georef_info.ny

        num_xtiles = int(np.ceil(float(grid_nx) / float(tile_xdim)))
        num_ytiles = int(np.ceil(float(grid_ny) / float(tile_ydim)))

        tiles = []

        for x_left in range(num_xtiles):
            for y_top in range(num_ytiles):

                left = max(x_left * tile_xdim - tile_xpadding, 0)
                right = min((x_left + 1) * tile_xdim + tile_xpadding, grid_nx)
                top = max(y_top * tile_ydim - tile_ypadding, 0)
                bottom = min((y_top + 1) * tile_ydim + tile_ypadding, grid_ny)

                needs_left_padding = (x_left == 0)
                needs_right_padding = (x_left == (num_xtiles - 1))
                needs_top_padding = (y_top == 0)
                needs_bottom_padding = (y_top == (num_ytiles - 1))

                griddata = self._griddata[top:bottom, left:right].copy()

                xllcenter = self._georef_info.xllcenter + left * dx
                yllcenter = self._georef_info.yllcenter + (grid_ny - bottom) * dx
                nx = right - left
                ny = bottom - top

                # Each `and tile_?padding` guard matters: with zero padding
                # `griddata[:, -0:]` is the whole array, so the edge tiles
                # came back with their data duplicated.
                if needs_left_padding and tile_xpadding:
                    griddata = np.concatenate((np.fliplr(griddata[:, 0:tile_xpadding]), griddata), axis=1)
                    xllcenter -= tile_xpadding * dx
                    nx += tile_xpadding

                if needs_right_padding and tile_xpadding:
                    griddata = np.concatenate((griddata, np.fliplr(griddata[:, -tile_xpadding:])), axis=1)
                    nx += tile_xpadding

                if needs_top_padding and tile_ypadding:
                    griddata = np.concatenate((np.flipud(griddata[0:tile_ypadding, :]), griddata), axis=0)
                    ny += tile_ypadding

                if needs_bottom_padding and tile_ypadding:
                    griddata = np.concatenate((griddata, np.flipud(griddata[-tile_ypadding:, :])), axis=0)
                    ny += tile_ypadding
                    yllcenter -= tile_ypadding * dx

                new_tile = self.__class__()
                # (ny - 0.5), not (ny + 0.5): the north edge of the transform
                # has to invert yllcenter = gt[3] - dx*(ny - 0.5), otherwise
                # every tile is georeferenced one cell too far north.
                new_tile._georef_info.geoTransform = (
                    xllcenter - 0.5 * dx, dx, 0, yllcenter + (ny - 0.5) * dx, 0, -dx)
                new_tile._georef_info.projection = self._georef_info.projection
                new_tile._georef_info.xllcenter = xllcenter
                new_tile._georef_info.yllcenter = yllcenter
                new_tile._georef_info.dx = dx
                new_tile._georef_info.nx = nx
                new_tile._georef_info.ny = ny
                new_tile._griddata = griddata

                tiles.append(new_tile)

        return tiles

    def remove_padding(self, xpadding = 10, ypadding = 10):
        """Strip the overlap that :meth:`tile` added, ready for :meth:`mosaic`."""
        import copy
        unpadded = copy.deepcopy(self)
        # `[-0:]`-style slices select everything, so guard the zero case.
        rows = slice(ypadding, -ypadding if ypadding else None)
        cols = slice(xpadding, -xpadding if xpadding else None)
        unpadded._griddata = unpadded._griddata[rows, cols]
        unpadded._georef_info.xllcenter += xpadding*unpadded._georef_info.dx
        unpadded._georef_info.yllcenter += ypadding*unpadded._georef_info.dx
        unpadded._georef_info.nx -= 2*xpadding
        unpadded._georef_info.ny -= 2*ypadding
        unpadded._georef_info.geoTransform = (unpadded._georef_info.xllcenter - 0.5*unpadded._georef_info.dx, unpadded._georef_info.dx, 0, unpadded._georef_info.yllcenter + (float(unpadded._georef_info.ny-0.5))*self._georef_info.dx, 0, -self._georef_info.dx)
        
        return unpadded
        
    def sort(self, reverse=True, force = False, mask = None):
        """Flat indices of the cells, ordered by value.

        Parameters
        ----------
        reverse : bool
            ``True`` (the default) returns highest-value first.
        force : bool
            Recompute even if a sort is already cached.
        mask : BaseSpatialGrid, optional
            Keep only cells where the mask equals 1.

        Returns
        -------
        ndarray
            Flat (C-order) indices into ``_griddata``.

        Notes
        -----
        The result is cached, so passing a different ``mask`` on a later call
        has no effect unless ``force=True`` is also given.
        """
        if force:
            self._sorted = False
        if not self._sorted:
            sort_indexes = self._griddata.argsort(axis=None)
            if mask is not None:
                # Select the sorted positions whose cell is inside the mask.
                # The previous version fed positions-in-the-sorted-array to
                # ravel_multi_index as if they were (row, col) subscripts,
                # which produced indices unrelated to the mask.
                mask_flat = np.asarray(mask._griddata).reshape(-1)
                sort_indexes = sort_indexes[mask_flat[sort_indexes] == 1]
            self._sort_indexes = sort_indexes
            self._sorted = True
        if reverse:
            return self._sort_indexes[::-1]
        else:
            return self._sort_indexes

    def apply_moving_window(self, moving_window):
        """Apply a :class:`~TopoAnalysis.MovingWindow.MovingWindow` to the grid."""
        out_grid = self.__class__()
        out_grid._copy_info_from_grid(self, set_zeros=True)
        out_grid._griddata = moving_window.apply_moving_window(
            self._griddata, self._georef_info.dx, self.dtype)
        return out_grid

    def average_over_distance(self, distance, grid=None):
        """Mean of the grid within ``distance`` of each cell.

        The averaging kernel is a disc of radius ``distance`` in map units,
        normalised to unit sum, applied by FFT convolution.  Edges wrap, as
        they always do for an FFT convolution.
        """
        info = self._georef_info
        dx = info.dx
        if grid is None:
            grid = self._griddata
        grid = np.asarray(grid, dtype=float64)
        ny, nx = grid.shape

        # Build the kernel about index (0, 0) with wrap-around distances.
        # That is exactly the origin a circular convolution assumes, so no
        # fftshift is needed afterwards -- and unlike a centre-of-array
        # kernel it stays correct for even grid dimensions.  (Building the
        # distances with np.arange(xmin, xmax, dx) could also produce
        # nx +/- 1 samples and make the FFT multiply below fail outright.)
        di = np.minimum(np.arange(ny), ny - np.arange(ny)) * dx
        dj = np.minimum(np.arange(nx), nx - np.arange(nx)) * dx
        D = np.hypot(di[:, None], dj[None, :])
        template = (D <= distance).astype(float64)
        total = template.sum()
        if total == 0:
            raise Error.InputError(
                'average_over_distance',
                'distance {0} is smaller than one cell ({1})'.format(distance, dx))
        template /= total

        from numpy.fft import fft2, ifft2

        return np.real(ifft2(fft2(grid) * fft2(template)))

    def clip_to_bounds(self, bounds):
        """Clip to ``((xmin, xmax), (ymin, ymax))``."""
        extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
        return self.clip_to_extent(extent)

    def clip_to_extent(self, extent):
        """Clip to ``(xmin, xmax, ymin, ymax)``, keeping every cell inside it.

        Returns
        -------
        BaseSpatialGrid
            A new grid.  The original is untouched.
        """
        import copy

        xmin, xmax, ymin, ymax = extent
        info = self._georef_info
        if xmax < xmin or ymax < ymin:
            raise Error.InputError('clip_to_extent',
                                   'extent must be (xmin, xmax, ymin, ymax)')

        # Work in index space and clamp, rather than round-tripping through
        # _xy_to_rowscols (which returns None outside the grid) and then
        # nudging the result by a cell in each direction.
        col_min = int(np.ceil((xmin - info.xllcenter) / info.dx - 1e-9))
        col_max = int(np.floor((xmax - info.xllcenter) / info.dx + 1e-9))
        row_from_y_min = int(np.ceil((ymin - info.yllcenter) / info.dx - 1e-9))
        row_from_y_max = int(np.floor((ymax - info.yllcenter) / info.dx + 1e-9))

        col_min = max(col_min, 0)
        col_max = min(col_max, info.nx - 1)
        row_from_y_min = max(row_from_y_min, 0)
        row_from_y_max = min(row_from_y_max, info.ny - 1)

        if col_max < col_min or row_from_y_max < row_from_y_min:
            raise Error.InputError('clip_to_extent',
                                   'extent does not overlap the grid')

        # Row 0 is the north edge, so the largest y is the smallest row.
        top = (info.ny - 1) - row_from_y_max
        bottom = (info.ny - 1) - row_from_y_min

        return_grid = self.__class__()
        return_grid._georef_info = copy.deepcopy(info)
        return_grid._griddata = self._griddata[top:bottom + 1, col_min:col_max + 1].copy()
        return_grid._georef_info.nx = return_grid._griddata.shape[1]
        return_grid._georef_info.ny = return_grid._griddata.shape[0]
        return_grid._georef_info.xllcenter = info.xllcenter + col_min * info.dx
        return_grid._georef_info.yllcenter = info.yllcenter + row_from_y_min * info.dx
        return_grid._georef_info.geoTransform = return_grid._geoTransform_from_corner()

        return return_grid

    def extent_of_data(self):
        """``(xmin, xmax, ymin, ymax)`` of the smallest box holding all valid data."""
        valid = ~np.isnan(np.asarray(self._griddata, dtype=float64))
        if not valid.any():
            raise Error.InputError('extent_of_data', 'the grid is entirely no-data')

        rows = np.flatnonzero(valid.any(axis=1))
        # Columns need any() along axis 0; the original indexed rows for both,
        # so a grid with empty left columns reported the wrong extent (and
        # ran off the end for non-square grids).
        cols = np.flatnonzero(valid.any(axis=0))
        top, bottom = int(rows[0]), int(rows[-1])
        left, right = int(cols[0]), int(cols[-1])

        (ll, ur) = self._rowscols_to_xy(((bottom, left), (top, right)))
        return (ll[0], ur[0], ll[1], ur[1])


    def clip_to_mask_grid(self, mask_grid):
        """Reproject and resample this grid onto ``mask_grid``'s frame, in place."""
        gdal_source = self._create_gdal_representation_from_array(self._georef_info, 'MEM', self._griddata, self.dtype)
        gdal_mask = mask_grid._create_gdal_representation_from_array(mask_grid._georef_info, 'MEM', mask_grid._griddata, mask_grid.dtype)
        gdal_clip = self._clipRasterToRaster(gdal_source, gdal_mask, self.dtype)
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._read_GDAL_dataset(gdal_clip, self.dtype)
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)

    def clip_to_shapefile(self, shapefile):
        """Cut the grid to a vector boundary, in place.

        .. warning::
           ``_clipRasterToShape`` is a stub, so this currently does nothing.
        """
        gdal_source = self._create_gdal_representation_from_array(self._georef_info, 'MEM', self._griddata, self.dtype)
        gdal_clip = shapefile.shapedata
        gdal_result = self._clipRasterToShape(gdal_source, gdal_clip)
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._read_GDAL_dataset(gdal_result, self.dtype)
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)
    
    def calculate_gradient_over_length_scale(self, length_scale):
        """Centred finite-difference gradient measured over ``length_scale``.

        The stencil reaches ``N = ceil(length_scale / dx)`` cells either side
        of the centre, so the two sampled cells are ``2*N*dx`` apart.  A band
        ``N`` cells wide around the grid edge comes back as ``NaN``.

        Returns
        -------
        (ndarray, ndarray)
            ``(dz/dx, dz/dy)`` with y increasing northwards.
        """
        elevGrid = np.asarray(self._griddata, dtype=float64)
        dx = self._georef_info.dx
        # int(), not a bare np.ceil(): NumPy returns a float and Python 3
        # refuses to slice with one.
        N = int(np.ceil(length_scale / dx))
        if N < 1:
            raise Error.InputError('calculate_gradient_over_length_scale',
                                   'length_scale must be at least one cell')
        if 2 * N >= min(elevGrid.shape):
            raise Error.InputError('calculate_gradient_over_length_scale',
                                   'length_scale is larger than the grid')

        # The separation between the sampled cells is 2*N*dx, not (2*N+1)*dx.
        span = 2.0 * N * dx
        Sx = (elevGrid[N:-N, (2 * N):] - elevGrid[N:-N, :-(2 * N)]) / span
        # Row index increases southwards, so negate to get dz/dy in map space.
        Sy = -(elevGrid[(2 * N):, N:-N] - elevGrid[:-(2 * N), N:-N]) / span

        SxPadded = np.full(elevGrid.shape, np.nan)
        SyPadded = np.full(elevGrid.shape, np.nan)
        SxPadded[N:-N, N:-N] = Sx
        SyPadded[N:-N, N:-N] = Sy

        return SxPadded, SyPadded

    def calculate_laplacian_over_length_scale(self, length_scale):
        """Laplacian measured over ``length_scale`` with a centred stencil.

        A band ``N = ceil(length_scale / dx)`` cells wide around the grid edge
        comes back as ``NaN``.
        """
        dx = self._georef_info.dx
        grid = np.asarray(self._griddata, dtype=float64)
        N = int(np.ceil(length_scale / dx))
        if N < 1:
            raise Error.InputError('calculate_laplacian_over_length_scale',
                                   'length_scale must be at least one cell')
        if 2 * N >= min(grid.shape):
            raise Error.InputError('calculate_laplacian_over_length_scale',
                                   'length_scale is larger than the grid')

        Curv = np.full(grid.shape, np.nan)
        h2 = (N * dx) ** 2
        centre = grid[N:-N, N:-N]
        Cx = (grid[N:-N, (2 * N):] - 2 * centre + grid[N:-N, :-(2 * N)]) / h2
        Cy = (grid[(2 * N):, N:-N] - 2 * centre + grid[:-(2 * N), N:-N]) / h2

        Curv[N:-N, N:-N] = Cx + Cy
        return Curv

    def principal_curvatures(self):
        """Maximum and minimum principal curvature of the surface.

        Returns
        -------
        (BaseSpatialGrid, BaseSpatialGrid)
            ``(k1, k2)`` with ``k1 >= k2``.  Positive curvature is convex-up.
        """
        dx = self._georef_info.dx
        grid = np.asarray(self._griddata, dtype=float64)
        # np.gradient returns the derivative along axis 0 (rows) first, and
        # rows increase southwards, so flip its sign to get dz/dy.
        dz_drow, Zx = np.gradient(grid, dx)
        Zy = -dz_drow
        _, Zxx = np.gradient(Zx, dx)
        dZy_drow, Zxy = np.gradient(Zy, dx)
        Zyy = -dZy_drow

        p = Zx ** 2 + Zy ** 2
        H = (Zx ** 2 + 1) * Zyy - 2 * Zx * Zy * Zxy + (Zy ** 2 + 1) * Zxx
        H = -H / (2 * (p + 1) ** 1.5)

        K = (Zxx * Zyy - Zxy ** 2) / (1 + p) ** 2

        # H**2 - K is non-negative in exact arithmetic; clip so that rounding
        # near an umbilic point does not turn both curvatures into NaN.
        discriminant = np.sqrt(np.maximum(H ** 2 - K, 0.0))
        k1 = H + discriminant
        k2 = H - discriminant

        k1re = BaseSpatialGrid()
        k1re._copy_info_from_grid(self, True)
        k2re = BaseSpatialGrid()
        k2re._copy_info_from_grid(self, True)
        k1re._griddata = k1
        k2re._griddata = k2

        return (k1re, k2re)
        
    def extent(self):
        """``(xmin, xmax, ymin, ymax)`` of the grid's outer pixel edges.

        These are the bounds matplotlib's ``imshow(extent=...)`` expects.
        """
        info = self._georef_info
        half = 0.5 * info.dx
        return (info.xllcenter - half,
                info.xllcenter + (info.nx - 0.5) * info.dx,
                info.yllcenter - half,
                info.yllcenter + (info.ny - 0.5) * info.dx)

    def plot(self, **kwargs):
        """Display the grid with ``imshow``, georeferenced in map coordinates.

        Parameters
        ----------
        interactive : bool
            Leave the figure open without blocking (default ``True``).
        colorbar : bool
            Draw a colour bar (default ``True``).
        ax : matplotlib Axes, optional
            Draw into an existing axes instead of a new figure.
        **kwargs
            Passed through to ``imshow`` (``cmap``, ``vmin``, ``alpha``, ...).

        Returns
        -------
        matplotlib.figure.Figure
        """
        interactive = kwargs.pop('interactive', True)
        colorbar = kwargs.pop('colorbar', True)
        ax = kwargs.pop('ax', None)
        kwargs.setdefault('extent', self.extent())

        if ax is None:
            fig = plt.figure()
            ax = fig.gca()
        else:
            fig = ax.figure

        im = ax.imshow(self._griddata, **kwargs)
        if colorbar:
            fig.colorbar(im, ax=ax)
        if interactive:
            plt.ion()
            plt.show(block=False)
            plt.pause(0.001)
        return fig

    def _neighbourhood_indices(self, index, pixel_radius):
        """Row/column indices of the in-bounds cells within ``pixel_radius``.

        Returns ``(rows, cols)`` as flat integer arrays.  The disc is
        inclusive of the radius and of the far edge -- the previous
        ``range(i - r, i + r)`` stopped one cell short and let negative
        indices wrap around to the opposite side of the grid.
        """
        i, j = int(index[0]), int(index[1])
        r = int(np.ceil(pixel_radius))
        ny, nx = self._griddata.shape

        rows = np.arange(max(i - r, 0), min(i + r, ny - 1) + 1)
        cols = np.arange(max(j - r, 0), min(j + r, nx - 1) + 1)
        if rows.size == 0 or cols.size == 0:
            return np.empty(0, dtype=int), np.empty(0, dtype=int)

        RR, CC = np.meshgrid(rows, cols, indexing='ij')
        inside = np.hypot(RR - i, CC - j) <= pixel_radius
        return RR[inside], CC[inside]

    def find_nearest_cell_with_value(self, index, value, pixel_radius):
        """Index of the cell within ``pixel_radius`` whose value is closest to ``value``.

        Returns ``(None, None)`` when the neighbourhood holds no valid data.
        """
        rows, cols = self._neighbourhood_indices(index, pixel_radius)
        if rows.size == 0:
            return None, None
        difference = np.abs(value - self._griddata[rows, cols])
        if np.all(np.isnan(difference)):
            return None, None
        best = int(np.nanargmin(difference))
        return int(rows[best]), int(cols[best])

    def find_nearest_cell_with_value_greater_than(self, index, value, pixel_radius):
        """Index of the cell within ``pixel_radius`` whose value is the smallest
        one at or above ``value``.

        If no cell reaches ``value``, falls back to the cell closest to it
        from below.  Returns ``(None, None)`` when there is no valid data.
        """
        rows, cols = self._neighbourhood_indices(index, pixel_radius)
        if rows.size == 0:
            return None, None
        difference = self._griddata[rows, cols] - value

        above = difference >= 0
        if np.any(above & ~np.isnan(difference)):
            candidates = np.where(above, difference, np.nan)
        else:
            candidates = np.where(~np.isnan(difference), -difference, np.nan)
        if np.all(np.isnan(candidates)):
            return None, None
        best = int(np.nanargmin(candidates))
        return int(rows[best]), int(cols[best])

    def find_nearest_cell_with_greatest_value(self, index, pixel_radius = 5):
        """Index of the highest-valued cell within ``pixel_radius``.

        Returns ``(None, None)`` when the neighbourhood holds no valid data.
        """
        rows, cols = self._neighbourhood_indices(index, pixel_radius)
        if rows.size == 0:
            return None, None
        values = self._griddata[rows, cols]
        if np.all(np.isnan(values)):
            return None, None
        best = int(np.nanargmax(values))
        return int(rows[best]), int(cols[best])


    def snap_locations_to_greatest_value(self, v, pixel_radius = 5):
        """Move each ``(x, y)`` to the highest-valued cell within ``pixel_radius``.

        Locations for which no valid cell is found are dropped.
        """
        snap_idxs = list()
        for idx in self._xy_to_rowscols(v):
            if idx[0] is None:
                continue
            i, j = self.find_nearest_cell_with_greatest_value(idx, pixel_radius)
            if i is not None:
                snap_idxs.append((i, j))
        return self._rowscols_to_xy(snap_idxs)

    def snap_locations_to_closest_value(self, v, value, pixel_radius = 5):
        """Move each ``(x, y)`` to the nearby cell whose value best matches.

        ``value`` is a sequence with one target value per location.
        Locations for which no valid cell is found are dropped.
        """
        snap_idxs = list()
        for (idx, target) in zip(self._xy_to_rowscols(v), value):
            if idx[0] is None:
                continue
            i, j = self.find_nearest_cell_with_value(idx, target, pixel_radius)
            if i is not None:
                snap_idxs.append((i, j))
        return self._rowscols_to_xy(snap_idxs)

    def export_arc_e00_grid(self, filename):
        """Write an ArcInfo E00 interchange grid."""
        self._create_gdal_representation_from_array(self._georef_info, 'E00GRID', self._griddata, self.dtype, filename)

    def set_value_at_rowscols(self, value, rowscols):
        """Set ``value`` at each ``(row, col)`` in ``rowscols``."""
        rowscols = list(rowscols)
        if not rowscols:
            return
        rows, cols = zip(*rowscols)
        self._griddata[np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)] = value

    def save(self, filename):
        """Write the grid to an LZW-compressed GeoTIFF."""
        dataset = self._create_gdal_representation_from_array(
            self._georef_info, 'GTiff', self._griddata, self.dtype, filename,
            ['COMPRESS=LZW'])
        # Dropping the reference is what flushes GDAL's write cache; without
        # it the file can be truncated if the process exits soon after.
        dataset.FlushCache()
        del dataset

    def write_to_ai(self, filename, nodata_value=-9999.0):
        """Write the grid as an ArcInfo ASCII grid."""
        self._writeArcAsciiRaster(self._georef_info, filename, self._griddata,
                                  nodata_value, '%10.2f')

    def vectorize(self, filename):
        """Polygonize the grid into an ESRI shapefile called ``filename.shp``.

        Each polygon carries the raster value in a ``RASTER_VAL`` field.
        """
        import tempfile

        handle, tmpfilename = tempfile.mkstemp(suffix='.tif')
        os.close(handle)
        try:
            self.save(tmpfilename)
            gdal_dataset = gdal.Open(tmpfilename)
            try:
                srcband = gdal_dataset.GetRasterBand(1)
                drv = ogr.GetDriverByName("ESRI Shapefile")
                dst_ds = drv.CreateDataSource(filename + ".shp")
                try:
                    dst_layer = dst_ds.CreateLayer(os.path.basename(filename), srs=None)
                    dst_layer.CreateField(ogr.FieldDefn('RASTER_VAL', ogr.OFTReal))
                    gdal.Polygonize(srcband, None, dst_layer, 0, [], callback=None)
                finally:
                    # The shapefile is only written when the data source is
                    # released; the original code left it dangling.
                    dst_ds = None
            finally:
                gdal_dataset = None
        finally:
            if os.path.exists(tmpfilename):
                os.remove(tmpfilename)

    @classmethod
    def _georef_from_dataset(cls, return_object, gdal_dataset):
        """Populate a grid's georeferencing from an open GDAL dataset."""
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        info = return_object._georef_info
        info.geoTransform = gdal_dataset.GetGeoTransform()
        info.projection = gdal_dataset.GetProjection()
        info.dx = info.geoTransform[1]
        info.xllcenter = info.geoTransform[0] + info.dx / 2.0
        info.yllcenter = info.geoTransform[3] - (info.dx * (ny - 0.5))
        info.nx = nx
        info.ny = ny
        return return_object

    @classmethod
    def load(cls, filename):
        """Read a grid previously written by :meth:`save`."""
        return_object = cls()
        gdal_dataset = gdal.Open(filename)
        try:
            return_object._griddata = cls._read_band(gdal_dataset.GetRasterBand(1), cls.dtype)
            cls._georef_from_dataset(return_object, gdal_dataset)
        finally:
            gdal_dataset = None
        return return_object

    @classmethod
    def _mosaic_frame(cls, tiles):
        """Georeferencing for the smallest grid covering every tile."""
        if not tiles:
            raise Error.InputError('mosaic', 'no tiles to mosaic')

        dx = tiles[0]._georef_info.dx
        xmin = min(t._georef_info.xllcenter for t in tiles)
        ymin = min(t._georef_info.yllcenter for t in tiles)
        xmax = max(t._georef_info.xllcenter + (t._georef_info.nx - 1) * t._georef_info.dx
                   for t in tiles)
        ymax = max(t._georef_info.yllcenter + (t._georef_info.ny - 1) * t._georef_info.dx
                   for t in tiles)

        ny = int(round((ymax - ymin) / dx)) + 1
        nx = int(round((xmax - xmin) / dx)) + 1

        info = Georef_info()
        info.nx = nx
        info.ny = ny
        info.dx = dx
        info.xllcenter = xmin
        info.yllcenter = ymin
        info.projection = tiles[0]._georef_info.projection
        info.geoTransform = (xmin - 0.5 * dx, dx, 0,
                             ymin + (ny - 0.5) * dx, 0, -dx)
        return info

    @staticmethod
    def _tile_slice(info, tile):
        """Where ``tile`` belongs inside a mosaic described by ``info``."""
        j_min = int(round((tile._georef_info.xllcenter - info.xllcenter) / info.dx))
        i_min = int(round((info.yllcenter + (info.ny - 1) * info.dx
                           - tile._georef_info.yllcenter
                           - (tile._georef_info.ny - 1) * tile._georef_info.dx) / info.dx))
        return (slice(i_min, i_min + tile._georef_info.ny),
                slice(j_min, j_min + tile._georef_info.nx))

    @classmethod
    def mosaic(cls, tiles):
        """Reassemble tiles produced by :meth:`tile` into a single grid.

        Tiles must share a cell size and be aligned on a common grid.  Where
        tiles overlap, the last one in ``tiles`` wins.
        """
        tiles = list(tiles)
        info = cls._mosaic_frame(tiles)

        return_object = cls()
        return_object._georef_info = info
        return_object._griddata = np.full((info.ny, info.nx), np.nan)

        for tile in tiles:
            rows, cols = cls._tile_slice(info, tile)
            return_object._griddata[rows, cols] = tile._griddata

        return return_object

class ValueGrid(BaseSpatialGrid):
    """A plain raster of values, with a helper for scattered writes."""


    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('dx', 'grid'), '_create_from_grid'))

    def set_value_at_indexes(self, indexes, value):
        """Set ``value`` at each ``(row, col)`` pair in ``indexes``."""
        indexes = list(indexes)
        if not indexes:
            return
        # zip() is a lazy iterator in Python 3; using one directly as a NumPy
        # subscript raises rather than doing what it did in Python 2.
        rows, cols = zip(*indexes)
        self._griddata[np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)] = value


class FlowDirection(BaseSpatialGrid):
    """Base class for flow-direction grids."""
    pass

class FlowDirectionD8(FlowDirection):
    """Single-direction (D8) flow directions.

    Values are ArcGIS direction codes; 0 means the cell has no downstream
    neighbour (a pit, an outlet at the grid edge, or no-data)::

        | 32  64 128 |
        | 16   X   1 |
        |  8   4   2 |

    Build one from a filled DEM::

        d8 = FlowDirectionD8(flooded_dem=FilledElevation(elevation=dem))

    Direction is chosen by steepest *gradient*, so the drop to a diagonal
    neighbour is divided by ``sqrt(2) * dx`` before being compared with the
    cardinal ones.  Ties go to the first direction in the order E, SE, S, SW,
    W, NW, N, NE.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('flooded_dem',), '_create_from_flooded_dem'),
                                   (('elevation',), '_create_from_elevation'))

    dtype = uint8

    def _create_from_elevation(self, *args, **kwargs):
        """Route an already-hydrologically-correct DEM.

        Alias for :meth:`_create_from_flooded_dem`; use it when the DEM has
        been conditioned elsewhere.
        """
        kwargs = dict(kwargs)
        kwargs['flooded_dem'] = kwargs.pop('elevation')
        self._create_from_flooded_dem(*args, **kwargs)

    def _create_from_flooded_dem(self, *args, **kwargs):
        flooded_dem = kwargs['flooded_dem']
        self._copy_info_from_grid(flooded_dem)

        elevations = np.ascontiguousarray(flooded_dem._griddata, dtype=float64)
        kernels.warn_if_slow(elevations.size)

        # A geographic grid's cell size varies with latitude, so ask the grid
        # itself rather than assuming a single dx.
        cellsize_grid = np.ascontiguousarray(
            flooded_dem._mean_pixel_dimension(*args, **kwargs), dtype=float64)

        # The old hand-rolled version compared each cell only with a window
        # that stopped two columns/rows short, so the last row and column
        # never received a direction; it also scaled diagonals by 1.41 rather
        # than sqrt(2), and resolved ties towards NE rather than
        # deterministically.
        self._griddata = np.asarray(
            kernels.flow_directions(elevations, cellsize_grid), dtype=self.dtype)
        self._griddata[np.isnan(elevations)] = 0

        mask = kwargs.get('mask')
        self._sort_indexes = flooded_dem.sort(reverse=False, force=True, mask=mask)
        self._sorted = True
        self._invalidate_network()

    # -- flow network -----------------------------------------------------
    #
    # ``receivers`` and ``topological_order`` are derived from the direction
    # codes and cached.  Everything that walks the network downstream uses
    # them, which is what lets accumulation, flow length and chi run in a
    # single linear pass instead of a per-cell Python loop.

    def _invalidate_network(self):
        self._receivers = None
        self._topological_order = None
        self._cycle_count = 0

    @property
    def receivers(self):
        """Flat index of the cell each cell drains to, or -1 for none."""
        if getattr(self, '_receivers', None) is None:
            self._receivers = np.asarray(kernels.receivers(self._griddata), dtype=np.int64)
        return self._receivers

    @property
    def topological_order(self):
        """Flat cell indices ordered so every donor precedes its receiver."""
        if getattr(self, '_topological_order', None) is None:
            order, cycles = kernels.topological_order(self.receivers)
            self._topological_order = np.asarray(order, dtype=np.int64)
            self._cycle_count = int(cycles)
            if cycles:
                warnings.warn(
                    "{0} cells of this flow-direction grid lie on circular "
                    "flow paths and cannot be ordered; results that depend on "
                    "routing through them will be wrong. This usually means "
                    "the DEM was not filled before routing.".format(cycles),
                    RuntimeWarning, stacklevel=2)
        return self._topological_order

    @property
    def cycle_count(self):
        """How many cells lie on circular flow paths (0 for a valid network)."""
        self.topological_order  # populate the cache
        return self._cycle_count

    def step_lengths(self, pixel_dimension=None):
        """Distance from each cell to its receiver, in map units."""
        if pixel_dimension is None:
            pixel_dimension = self._mean_pixel_dimension()
        return np.asarray(pixel_dimension, dtype=float64) * self.pixel_scale(float64)

    # -- keeping the cache honest -----------------------------------------
    #
    # Everything downstream now reads `receivers`/`topological_order` rather
    # than re-deriving routing from `_griddata` on every call, so any write
    # to a flow code has to drop the cache.  Editing a code and having the
    # edit silently ignored would be far worse than the cost of rebuilding.

    def __setitem__(self, key, value):
        super(FlowDirectionD8, self).__setitem__(key, value)
        self._invalidate_network()

    def set_value_at_rowscols(self, value, rowscols):
        super(FlowDirectionD8, self).set_value_at_rowscols(value, rowscols)
        self._invalidate_network()

    def invalidate_network(self):
        """Drop the cached receivers and ordering.

        Call this after modifying ``_griddata`` directly; the indexed write
        paths do it for you.
        """
        self._invalidate_network()

    def __flow_code_for_position(self, flooded_dem, i, j):
        """Steepest-descent code for a single cell.

        Kept for :meth:`update_flow_codes_in_mask`, which re-routes a handful
        of cells after their elevations change.
        """
        z = flooded_dem[i, j]
        if z is None or np.isnan(z):
            return 0

        best_code = 0
        best_slope = 0.0
        de = self._georef_info.dx
        for d, (di, dj) in enumerate(zip(kernels.D8_DI, kernels.D8_DJ)):
            neighbour = flooded_dem[i + int(di), j + int(dj)]
            if neighbour is None or np.isnan(neighbour):
                continue
            # The scaling used to be inverted -- cardinal drops were divided
            # by 1.41 and diagonal ones were not -- which biased routing
            # towards the diagonals.
            slope = (z - neighbour) / (de * kernels.D8_DIST[d])
            if slope > best_slope:
                best_slope = slope
                best_code = int(kernels.D8_CODES[d])
        return best_code

    def get_flow_to_cell(self,i,j):
        """Index of the cell that ``(i, j)`` drains into.

        Returns
        -------
        (int or None, int or None, bool)
            Row, column and a flag that is ``False`` when the cell has no
            downstream neighbour inside the grid.
        """
        code = self._griddata[i, j]
        if code == 0:
            return None, None, False

        for d in range(8):
            if code != kernels.D8_CODES[d]:
                continue
            i_out = i + int(kernels.D8_DI[d])
            j_out = j + int(kernels.D8_DJ[d])
            if 0 <= i_out < self._georef_info.ny and 0 <= j_out < self._georef_info.nx:
                return i_out, j_out, True
            break

        return None, None, False

    def get_upstream_cell_indexes(self, i, j):
        """Indices of the immediate neighbours that drain into ``(i, j)``."""
        options = list()
        for d in range(8):
            i_up = i + int(kernels.D8_DI[d])
            j_up = j + int(kernels.D8_DJ[d])
            # A neighbour in direction d drains here if it points back the
            # opposite way, i.e. four bit positions round the compass.
            opposite = int(kernels.D8_CODES[(d + 4) % 8])
            if self[i_up, j_up] == opposite:
                options.append((i_up, j_up))
        return options

    def __get_flow_from_cell(self, i, j, max_recursion_depth=None, depth=0):
        """Every cell upstream of ``(i, j)``, including ``(i, j)`` itself.

        Returns ``(rows, cols)`` as plain lists.  ``max_recursion_depth``
        limits how many steps upstream the search goes; ``depth`` is accepted
        for backward compatibility and gives the starting depth.

        The traversal is iterative.  The original recursive version needed
        ``sys.setrecursionlimit(1000000)`` and still crashed the interpreter
        on basins more than a few tens of thousands of cells long.
        """
        if max_recursion_depth is None and depth == 0:
            # No depth limit: let the compiled kernel do the whole basin in
            # one linear pass over the flow network.
            flat = np.int64(i) * self._georef_info.nx + np.int64(j)
            mask = kernels.upstream_mask(self.receivers, np.array([flat], dtype=np.int64))
            rows, cols = np.nonzero(mask)
            return rows.tolist(), cols.tolist()

        ny, nx = self._georef_info.ny, self._georef_info.nx
        limit = None if max_recursion_depth is None else max_recursion_depth - depth

        i_source = []
        j_source = []
        seen = set()
        stack = [((int(i), int(j)), 0)]
        while stack:
            (ci, cj), d = stack.pop()
            if (ci, cj) in seen:
                continue
            seen.add((ci, cj))
            i_source.append(ci)
            j_source.append(cj)
            if limit is not None and d >= limit:
                continue
            for k in range(8):
                ui = ci + int(kernels.D8_DI[k])
                uj = cj + int(kernels.D8_DJ[k])
                if not (0 <= ui < ny and 0 <= uj < nx):
                    continue
                if self._griddata[ui, uj] == kernels.D8_CODES[(k + 4) % 8]:
                    stack.append(((ui, uj), d + 1))

        return i_source, j_source

    def __map_flow_from_cell(self, index, **kwargs):
        """Build the nested upstream dictionary rooted at ``index``.

        Each node is a dict with:

        ``index``
            the ``(row, col)`` of the cell;
        ``distance_scale``
            1 or sqrt(2), the step length from this node's *parent* in cell
            widths -- 1 at the root;
        ``next``
            the list of upstream child nodes, absent at a headwater;
        plus one entry per keyword grid, holding that grid's value at the cell.

        ``assume_valid_routing=False`` guards against cycles in a
        flow-direction grid that was not derived from a filled DEM.

        Notes
        -----
        ``distance_scale`` used to be written onto the *parent* node and
        overwritten by each child in turn, so it ended up describing whichever
        child happened to be visited last. Chi profiles built from these
        dictionaries were wrong wherever a confluence had a diagonal and a
        cardinal tributary.
        """
        if len(index) == 1:
            root = tuple(index[0])
        else:
            root = tuple(index)

        assume_valid_routing = kwargs.get('assume_valid_routing', True)
        value_grids = {name: grid for name, grid in kwargs.items()
                       if name != 'assume_valid_routing'}
        ny, nx = self._georef_info.ny, self._georef_info.nx

        def make_node(cell, scale):
            i, j = cell
            node = {'index': cell, 'distance_scale': scale}
            for name, grid in value_grids.items():
                node[name] = grid[i, j]
            return node

        # Every node carries `index` as a plain (row, col) pair, the root
        # included. It used to keep the nested form the caller passed, so
        # consumers had to special-case the root -- and those that did not
        # failed to unpack it.
        root_node = make_node(root, 1.0)
        # Explicit stack of (node, cell) pairs instead of recursion.
        stack = [(root_node, root)]
        while stack:
            node, (i, j) = stack.pop()
            children = []
            for k in range(8):
                ui = i + int(kernels.D8_DI[k])
                uj = j + int(kernels.D8_DJ[k])
                if not (0 <= ui < ny and 0 <= uj < nx):
                    continue
                if self._griddata[ui, uj] != kernels.D8_CODES[(k + 4) % 8]:
                    continue
                if not assume_valid_routing:
                    if self.visited[ui, uj]:
                        continue
                    self.visited[ui, uj] = True
                child = make_node((ui, uj), float(kernels.D8_DIST[k]))
                children.append(child)
                stack.append((child, (ui, uj)))
            if children:
                node['next'] = children

        return root_node
    
    def update_flow_codes_in_mask(self, *args, **kwargs):
        """Re-derive flow codes inside a mask after the DEM has changed."""
        flooded_dem = args[0]
        mask = args[1]

        indexes = np.where(mask._griddata == 1)
        for (i,j) in zip(indexes[0],indexes[1]):
            flow_code = self.__flow_code_for_position(flooded_dem, i, j)
            self._griddata[i,j] = flow_code
        self._invalidate_network()

    def basin_mask(self, outlets):
        """Boolean array marking every cell draining to any of ``outlets``.

        ``outlets`` is a sequence of ``(x, y)`` map coordinates.
        """
        flat = []
        for (row, col) in self._xy_to_rowscols(outlets):
            if row is None:
                continue
            flat.append(row * self._georef_info.nx + col)
        mask = kernels.upstream_mask(self.receivers,
                                     np.array(flat, dtype=np.int64))
        return np.asarray(mask).astype(bool)

    def divides_for_outlets(self, outlet1, outlet2):
        """Cells on each side of the divide shared by two basins.

        Returns
        -------
        (tuple, tuple)
            ``(cells_in_basin1, cells_in_basin2)``, each a tuple of
            ``(row, col)`` pairs lying against the shared divide.
        """
        from scipy.ndimage import binary_dilation, generate_binary_structure

        in_basin1 = self.basin_mask((outlet1,))
        in_basin2 = self.basin_mask((outlet2,))

        # An 8-connected structure, so cells meeting only at a corner still
        # count as facing each other across the divide.
        structure = generate_binary_structure(2, 2)
        ind1 = np.argwhere(in_basin1 & binary_dilation(in_basin2, structure=structure))
        ind2 = np.argwhere(in_basin2 & binary_dilation(in_basin1, structure=structure))

        # Materialise the pairs: the old version returned zip objects, which
        # are single-use and silently came back empty on a second read.
        return (tuple(map(tuple, ind1)), tuple(map(tuple, ind2)))

    def locations_of_paired_hollows(self, outlet1, outlet2, area, Ao=1E5):
        """Channel heads reached by flowing down from a shared divide.

        Walks downstream from each divide cell until drainage area first
        exceeds ``Ao``, and returns the distinct cells where that happens on
        each side.
        """
        def map_down_to_hollow(ind):
            """Follow flow down until area first reaches Ao; None if never."""
            (i, j) = ind
            seen = set()
            while True:
                if (i, j) in seen:
                    return None  # circular routing
                seen.add((i, j))
                (next_i, next_j, good) = self.get_flow_to_cell(i, j)
                if not good:
                    return None
                here = area[i, j]
                there = area[next_i, next_j]
                if here is None or there is None:
                    return None
                if here < Ao <= there:
                    return next_i, next_j
                (i, j) = (next_i, next_j)

        ind1, ind2 = self.divides_for_outlets(outlet1, outlet2)
        h1 = {h for h in (map_down_to_hollow(ind) for ind in ind1) if h is not None}
        h2 = {h for h in (map_down_to_hollow(ind) for ind in ind2) if h is not None}

        return tuple(sorted(h1)), tuple(sorted(h2))

    def get_indexes_of_upstream_cells(self, i, j):
        """Tuple of ``(row, col)`` for every cell upstream of ``(i, j)``."""
        rows, cols = self.__get_flow_from_cell(i, j)
        # A tuple rather than a zip object: callers iterate the result more
        # than once, and a zip is exhausted after the first pass.
        return tuple(zip(rows, cols))

    def get_indexes_of_upstream_cells_for_location(self, x, y):
        """As :meth:`get_indexes_of_upstream_cells`, from a map coordinate."""
        v = ((x, y), )
        ((i, j),) = self._xy_to_rowscols(v)
        if i is None:
            return tuple()
        return self.get_indexes_of_upstream_cells(i, j)

    def search_down_flow_direction_with_length(self, start, search_length = np.inf):
        """Follow flow downstream from ``start``, recording distance travelled.

        Stops when flow leaves the grid, at a pit, once ``search_length`` is
        covered, or if the path revisits a cell.  The terminal cell is
        included in the result.

        Returns
        -------
        (tuple, tuple)
            ``(cells, distances)``, where ``cells`` are ``(row, col)`` pairs
            and ``distances`` are the along-path distances at each cell.
        """
        cells = list()
        length_list = list()
        ((row, col),) = self._xy_to_rowscols((start,))
        if row is None:
            return tuple(), tuple()

        # A set for the revisit test: the original scanned the whole path
        # list at every step, making long flow paths quadratic.  The edge
        # test is now `no in-grid receiver` rather than `on the perimeter`,
        # so a path that merely starts on the edge still gets followed.
        seen = set()
        length = 0.0
        while True:
            if (row, col) in seen or length >= search_length:
                break
            seen.add((row, col))
            cells.append((row, col))
            length_list.append(length)
            row_n, col_n, inBounds = self.get_flow_to_cell(row, col)
            if not inBounds:
                break
            diagonal = (row_n != row) and (col_n != col)
            length += self._georef_info.dx * (SQRT2 if diagonal else 1.0)
            (row, col) = (row_n, col_n)

        return tuple(cells), tuple(length_list)

    def search_down_flow_direction(self, start, search_length = np.inf):
        """The cells along the flow path from ``start``, without distances."""
        return self.search_down_flow_direction_with_length(start, search_length=search_length)[0]
    
    def search_down_flow_direction_from_rowscols_location(self, start, return_rowscols = False, search_length = np.inf):
        """Follow flow downstream from a ``(row, col)`` start.

        Returns map coordinates unless ``return_rowscols`` is true.
        """
        rcs = self.search_down_flow_direction(self._rowscols_to_xy(((start),))[0], search_length = search_length)
        return rcs if return_rowscols else self._rowscols_to_xy(rcs)

    def search_down_flow_direction_from_xy_location(self, start, return_rowscols=False, search_length=np.inf):
        """Follow flow downstream from an ``(x, y)`` start.

        Returns map coordinates unless ``return_rowscols`` is true.
        """
        rcs = self.search_down_flow_direction(start, search_length=search_length)
        return rcs if return_rowscols else self._rowscols_to_xy(rcs)

    def convert_rivertools_directions_to_arc(self, no_data=0):
        """Convert RiverTools flow codes in place to the ArcGIS convention.

        RiverTools numbers the same eight directions but starts from the top
        right, so every code is rotated by one position::

            ArcGIS      | 32 64 128 |
                        | 16  X   1 |
                        |  8  4   2 |
        """
        codes = np.asarray(self._griddata)
        valid = codes > 0
        # np.log2 of an array cannot be wrapped in int(); use integer
        # arithmetic on the exponent instead, and take the no-data value as an
        # argument rather than reading an attribute that never existed.
        exponent = np.zeros(codes.shape, dtype=np.int16)
        exponent[valid] = np.round(np.log2(codes[valid].astype(float64))).astype(np.int16)
        exponent -= 1
        exponent[exponent == -1] = 7
        converted = np.zeros(codes.shape, dtype=self.dtype)
        converted[valid] = (2 ** exponent[valid]).astype(self.dtype)
        converted[codes == no_data] = no_data
        self._griddata = converted
        self._invalidate_network()

    def pixel_scale(self, dtype = np.float32):
        """Step-length multiplier for each cell: 1 cardinal, sqrt(2) diagonal."""
        dim = np.ones_like(self._griddata, dtype = dtype)
        diagonal = np.isin(self._griddata, [2, 8, 32, 128])
        dim[diagonal] = SQRT2
        return dim

    def map_values_to_recursive_list(self, outlet, **kwargs):
        """Nested dictionary of the basin upstream of ``outlet``.

        Pass grids as keywords to sample them at every cell, e.g.
        ``map_values_to_recursive_list(outlet, elevation=dem, area=area)``.
        See :meth:`__map_flow_from_cell` for the node layout.
        """
        v = (outlet, )
        (ij_outlet, ) = self._xy_to_rowscols(v)
        if ij_outlet[0] is None:
            raise Error.InputError('map_values_to_recursive_list',
                                   'outlet {0} lies outside the grid'.format(outlet))
        if kwargs.get('assume_valid_routing', True) is False:
            self.visited = np.zeros_like(self._griddata).astype(bool)
        return self.__map_flow_from_cell((ij_outlet,), **kwargs)

    def bounds_of_basin_for_outlet(self, outlet):
        """``((xmin, xmax), (ymin, ymax))`` of the basin draining to ``outlet``.

        ``outlet`` is given as ``(y, x)`` -- latitude then longitude -- which
        is the order the rest of this method's callers use.
        """
        (lat, lon) = outlet
        indexes = self.get_indexes_of_upstream_cells_for_location(lon, lat)
        xy = self._rowscols_to_xy(indexes)
        if not xy:
            raise Error.InputError('bounds_of_basin_for_outlet',
                                   'no cells drain to {0}'.format(outlet))
        # Comparing a float with None raises in Python 3, and `x > None or
        # None is None` evaluates the comparison first, so the original could
        # not survive its own first iteration.
        lons = [p[0] for p in xy]
        lats = [p[1] for p in xy]
        return ((min(lons), max(lons)), (min(lats), max(lats)))

    def divides(self):
        """Grid marking drainage divides: cells with no upstream neighbour."""
        divides = BaseSpatialGrid()
        divides._copy_info_from_grid(self, True)

        codes = self._griddata
        # A cell is a divide when none of its eight neighbours drains into it.
        is_divide = np.ones(codes.shape, dtype=bool)
        for d in range(8):
            di, dj = int(kernels.D8_DI[d]), int(kernels.D8_DJ[d])
            opposite = kernels.D8_CODES[(d + 4) % 8]
            neighbour = np.zeros(codes.shape, dtype=codes.dtype)
            src_i = slice(max(di, 0), codes.shape[0] + min(di, 0))
            src_j = slice(max(dj, 0), codes.shape[1] + min(dj, 0))
            dst_i = slice(max(-di, 0), codes.shape[0] + min(-di, 0))
            dst_j = slice(max(-dj, 0), codes.shape[1] + min(-dj, 0))
            neighbour[dst_i, dst_j] = codes[src_i, src_j]
            is_divide &= (neighbour != opposite)

        divides._griddata = is_divide.astype(divides.dtype)
        return divides

    def paired_divides(self, mask = None):
        """Pair each divide cell with the cell directly across the divide.

        The partner is the reflection of the divide cell's downstream
        neighbour through the cell itself, i.e. the cell one step *up* the
        opposite hillslope.

        Returns
        -------
        tuple
            ``((xy_from, xy_to), ...)`` in map coordinates.
        """
        divides = self.divides()
        if mask is not None:
            selected = np.where((divides._griddata == 1) & (mask._griddata > 0))
        else:
            selected = np.where(divides._griddata == 1)

        pairs = list()
        for (i, j) in zip(selected[0].tolist(), selected[1].tolist()):

            (i_next, j_next, good) = self.get_flow_to_cell(i, j)

            if good:
                i_opposite = 2*i - i_next
                j_opposite = 2*j - j_next

                if (0 <= i_opposite < self._georef_info.ny
                        and 0 <= j_opposite < self._georef_info.nx):
                    (xy, ) = self._rowscols_to_xy(((i, j),))
                    (xy_next, ) = self._rowscols_to_xy(((i_opposite, j_opposite),))
                    pairs.append((xy, xy_next))

        return tuple(pairs)

class Elevation(CalculationMixin, BaseSpatialGrid):
    """A DEM.

    Adds the terrain-derivative helpers of :class:`CalculationMixin` to the
    base raster, plus the edge detection that depression filling starts from.
    """

    def findDEMedge(self):
        """Cells where water can leave the DEM.

        Those are the grid perimeter plus every valid cell touching a no-data
        cell -- the same seed set TopoToolbox's ``fillsinks`` uses.  Priority
        flooding starts here.

        Returns
        -------
        (ndarray, ndarray)
            ``(rows, cols)``, as :func:`numpy.where` returns.
        """
        # Delegates to the kernel so that the cells the flood is seeded with
        # are, by construction, exactly the cells reported here.
        flat = kernels.default_seeds(np.asarray(self._griddata, dtype=float64))
        return np.unravel_index(flat, (self._georef_info.ny, self._georef_info.nx))

    def outlets_at_coastlines(self, iterations=3):
        """Map coordinates of the cells forming the land/sea boundary.

        Land is taken to be every cell with a positive elevation.  The
        coastline is the band left over when the land mask is dilated and
        eroded by ``iterations`` cells and the two are compared.
        """
        from scipy.ndimage import binary_dilation, binary_erosion

        land = np.asarray(self._griddata) > 0
        dilated = binary_dilation(land, iterations=iterations)
        eroded = binary_erosion(land, iterations=iterations)
        # Exactly one of the two is set only in the fringe between them.
        coastline = dilated.astype(np.int16) + eroded.astype(np.int16) == 1

        rows, cols = np.nonzero(coastline)
        return self._rowscols_to_xy(list(zip(rows.tolist(), cols.tolist())))

    def track_flow_downhill(self, starting_point, maximum_pit_depth = 20):
        """Walk downhill from ``starting_point``, hopping over shallow pits.

        At each step the lowest unvisited neighbour is taken, provided it is
        no more than ``maximum_pit_depth`` above the current cell.  That lets
        the path climb out of small closed depressions instead of stopping in
        them, without filling the DEM.

        Returns
        -------
        (tuple, tuple, tuple)
            ``(xy, distance, elevation)`` along the path.
        """
        adjust = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]

        (ij, ) = self._xy_to_rowscols((starting_point,))
        if ij[0] is None:
            raise Error.InputError('track_flow_downhill',
                                   'starting point lies outside the grid')

        def neighbour_elevations(cell):
            """Elevations of the eight neighbours; NaN where off-grid."""
            # Building this with `self[...] or None` produced an object-dtype
            # array, and `None not in array` then silently never matched.
            values = np.empty(8, dtype=float64)
            for k, (di, dj) in enumerate(adjust):
                value = self[cell[0] + di, cell[1] + dj]
                values[k] = np.nan if value is None else value
            return values

        visited = [tuple(ij)]
        seen = {tuple(ij)}
        e_a = neighbour_elevations(ij)

        while not np.any(np.isnan(e_a)):
            here = self[ij[0], ij[1]]
            moved = False
            for sl in np.argsort(e_a):
                this_ij = (ij[0] + adjust[sl][0], ij[1] + adjust[sl][1])
                if this_ij in seen:
                    continue
                if (e_a[sl] - here) > maximum_pit_depth:
                    continue
                visited.append(this_ij)
                seen.add(this_ij)
                ij = this_ij
                moved = True
                break
            if not moved:
                break
            e_a = neighbour_elevations(ij)

        xy = self._rowscols_to_xy(visited)
        elevations = tuple(self[i, j] for (i, j) in visited)
        distances = [0.0]
        for previous, current in zip(xy[:-1], xy[1:]):
            step = np.hypot(current[0] - previous[0], current[1] - previous[1])
            distances.append(distances[-1] + step)
        return xy, tuple(distances), elevations

class GeographicElevation(GeographicGridMixin, Elevation):
    """A DEM on a latitude/longitude grid."""

    pass

class Gradient(BaseSpatialGrid):
    """The two components of the topographic gradient.

    ::

        gradient = Gradient(elevation=dem)
        gradient._gx, gradient._gy

    Stored as two bands rather than one, so direction survives ``save()``.
    :meth:`average_gradient` smooths both components over a length scale,
    which is how a regional aspect is obtained.
    """

    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('elevation',), '_create_from_elevation'))
     
    def _create_from_elevation(self, *args, **kwargs):
        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation,True)
        dx = elevation._georef_info.dx
        [self._gy, self._gx] = np.gradient(elevation._griddata,dx)
    
    def average_gradient(self, distance):
        """Both gradient components averaged over a disc of radius ``distance``."""
        outgrid = Gradient()
        outgrid._copy_info_from_grid(self, True)
        outgrid._gy = self.average_over_distance(distance, grid = self._gy)
        outgrid._gx = self.average_over_distance(distance, grid = self._gx)
        return outgrid
    
    def save(self, filename):    
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', [self._gx, self._gy], self.dtype, filename, ['COMPRESS=LZW'], multiple_bands=True)
    
    @classmethod
    def load(cls, filename):
        """Read back the multi-band grid written by :meth:`save`."""
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.nan
            return grid
        
        return_object = cls()
        gdal_dataset = gdal.Open(filename)
        
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        
        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[0]+return_object._georef_info.dx/2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
        return_object._gx = get_band(gdal_dataset, 1)
        return_object._gy = get_band(gdal_dataset, 2)        
            
        gdal_file = None
        return return_object  
    
    def plot(self, **kwargs):

        interactive = kwargs.pop('interactive', True)
        extent = [self._georef_info.xllcenter, self._georef_info.xllcenter+(self._georef_info.nx-0.5)*self._georef_info.dx, self._georef_info.yllcenter, self._georef_info.yllcenter+(self._georef_info.ny-0.5)*self._georef_info.dx]
        mag = np.sqrt(np.power(self._gx, 2) + np.power(self._gy, 2))
        direction = np.arctan2(self._gy,self._gx)*180 / np.pi
        if kwargs.get('azimuth'):
            direction = 90.0 - direction
            direction = (direction >= 180)*(direction - 360) + (direction < 180) * (direction)
            kwargs.pop('azimuth')
        if kwargs.get('reflect'):
            direction = (direction > 90)*(direction - 180) + (direction < -90)*(direction + 180) + ((direction <= 90) & (direction >= -90))*direction
            kwargs.pop('reflect')
        plt.figure()
        plt.imshow(mag, extent = extent, **kwargs)
        plt.figure()
        plt.imshow(direction, extent = extent, **kwargs)
        if interactive:
            plt.ion()
            plt.show(block=False)
        else:
            plt.ioff()
            plt.show(block=True)

class ScarpWavelet(BaseSpatialGrid):
    """Template-matched scarp amplitude, age, orientation and signal-to-noise.

    Holds the four-band output of a scarp template search.  Load one with
    :meth:`load`; the bands become ``_A``, ``_kt``, ``_orientation`` and
    ``_SNR``.  :meth:`template_window` needs the separate ``scarplet``
    package.
    """

    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'))
        
    @classmethod
    def load(cls, filename):
        """Read back the multi-band grid written by :meth:`save`."""
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.nan
            return grid
        
        return_object = cls()
        gdal_dataset = gdal.Open(filename)
        
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        
        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[0]+return_object._georef_info.dx/2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
        return_object._griddata = np.zeros((ny,nx))
        return_object._A = get_band(gdal_dataset, 1)
        return_object._kt = get_band(gdal_dataset, 2)
        return_object._orientation = get_band(gdal_dataset, 3)
        return_object._SNR = get_band(gdal_dataset, 4)    
            
        gdal_file = None
        return return_object

    def load_elevation(self, filename):
        """Attach the DEM the templates were matched against."""
        self.elevation = Elevation.load(filename)

    def plot_orientations(self, *args, **kwargs):
        """Draw scarp orientations over a hillshade, faded by signal-to-noise.

        Pass ``elevation`` and ``distance`` to have the regional slope
        direction removed first, so what is drawn is the scarp orientation
        *relative to* the hillslope it sits on.
        """
        # pop inputs, create those that are needed:
        
        elevation = kwargs.pop('elevation', None)
        distance = kwargs.pop('distance', None)
        orientation = kwargs.pop('orientation', None)
        hillshade = kwargs.pop('hillshade', None)
        interactive = kwargs.pop('interactive', True)
        adjust_orientations = False
        
        if elevation is not None and distance is not None:
            gradient = Gradient(elevation = elevation)
            average_gradient = gradient.average_gradient(distance)
            orientation = BaseSpatialGrid()
            orientation._copy_info_from_grid(self, True)
            orientation._griddata = np.arctan2(average_gradient._gy,average_gradient._gx)*180 / np.pi
            if hillshade is None:
                hillshade = Hillshade(elevation = elevation, azimuth = 320, inclination = 20)
            adjust_orientations = True
        elif orientation is not None:
            adjust_orientations = True
        
        adjusted_orientations = BaseSpatialGrid()
        adjusted_orientations._copy_info_from_grid(self, True)
        adjusted_orientations._griddata = self._orientation

        if adjust_orientations:
            # Align orientations:
            orientation._griddata = -orientation._griddata
            orientation._griddata = (orientation._griddata < -90.0).astype(float)*(orientation._griddata + 180.0) + (orientation._griddata >= 90.0).astype(float)*(orientation._griddata - 180.0) + ((orientation._griddata > -90.0) & (orientation._griddata < 90.0))*orientation._griddata
            adjusted_orientations._griddata = adjusted_orientations._griddata - orientation._griddata

        adjusted_orientations._griddata = (adjusted_orientations._griddata < -90.0).astype(float)*(adjusted_orientations._griddata + 180.0) + (adjusted_orientations._griddata >= 90.0).astype(float)*(adjusted_orientations._griddata - 180.0) + ((adjusted_orientations._griddata > -90.0) & (adjusted_orientations._griddata < 90.0)).astype(float)*adjusted_orientations._griddata

        extent = kwargs.pop('extent', None)
        if extent is None:
            extent = [self._georef_info.xllcenter, self._georef_info.xllcenter+(self._georef_info.nx-0.5)*self._georef_info.dx, self._georef_info.yllcenter, self._georef_info.yllcenter+(self._georef_info.ny-0.5)*self._georef_info.dx]

        plt.figure()
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', -90)
        vmax = kwargs.pop('vmax', 90)
        title = kwargs.pop('title', '')

        if hillshade is not None:
            from matplotlib import cm
            plt.imshow(hillshade._griddata, extent = extent, cmap = cm.gray, **kwargs)
        self._adjusted_orientations = adjusted_orientations._griddata
        
        from matplotlib import colors
        norm = colors.Normalize(vmin = vmin, vmax = vmax, clip = True)
        normalized_data = norm(adjusted_orientations._griddata)
        adjusted_orientations_rgba = cm.ScalarMappable(cmap = cmap, norm = norm).to_rgba(adjusted_orientations._griddata)
        norm = colors.LogNorm()
        adjusted_orientations_rgba[:,:,3] = (~normalized_data.mask).astype(float)*norm(self._SNR)
        plt.imshow(adjusted_orientations_rgba, extent = extent, **kwargs)
        plt.title(title)
        if interactive:
            plt.ion()
            plt.show(block=False)
        else:
            plt.ioff()
            plt.show(block=True)


    def calculate_valid_orientations(self,
                                     window_size,
                                     age,
                                     num=46,
                                     discard_max_min=True):
        """Range of template orientations that stay inside the valid data.

        Returns ``(max_orientation, min_orientation)`` in radians, ``NaN``
        where no orientation fits.  Cells within ``window_size`` of an edge
        are always ``NaN``.
        """
        data = self.valid_data()
        max = -np.inf * np.ones_like(data)
        min = np.inf * np.ones_like(data)

        orientations = np.linspace(-np.pi / 2, np.pi / 2, num=num)
        for orientation in orientations:
            # XXX: This uses orientation to conform to N-E coordinates
            # i.e. N = 0, E = -90, W = +90
            this = self._convolve_mask(data, window_size, age, orientation)
            mask = (this != this.max()) * (max < orientation)
            max[mask] = orientation
            mask = (this != this.max()) * (min > orientation)
            min[mask] = orientation

        max[np.isinf(max)] = np.nan
        min[np.isinf(min)] = np.nan
        max[np.isnan(data)] = np.nan
        min[np.isnan(data)] = np.nan

        max[0:window_size, :] = np.nan
        max[-window_size:, :] = np.nan
        max[:, 0:window_size] = np.nan
        max[:, -window_size:] = np.nan
        min[0:window_size, :] = np.nan
        min[-window_size:, :] = np.nan
        min[:, 0:window_size] = np.nan
        min[:, -window_size:] = np.nan

        if discard_max_min:
            max[max == orientations.max()] = np.nan
            min[min == orientations.min()] = np.nan

        return max, min

    def _convolve_mask(self, data, window_size, age, orientation):
        """
        Returns the number of valid pixels in specified window at each pixel
        """
        from numpy.fft import fft2, ifft2, fftshift
        mask = ~np.isnan(data)
        window = self.template_window(window_size,
                                      age,
                                      orientation,
                                      use_pixels=True)
        count = np.round(np.real(fftshift(ifft2(fft2(window) * fft2(mask)))))
        return count

    def template_window(self, window_size, age, orientation, use_pixels=False):
        """The scarp template mask at one size, age and orientation.

        Needs the separate ``scarplet`` package.
        """
        nx = self._georef_info.nx
        ny = self._georef_info.ny
        if use_pixels:
            dx = 1
        else:
            dx = self._georef_info.dx

        from scarplet.WindowedTemplate import Scarp
        template = Scarp(window_size, age, orientation, nx, ny, dx)
        window = template.get_mask()

        return window

    def valid_data(self):
        """1 where the attached DEM has usable, above-sea-level data; NaN elsewhere."""
        if self.elevation is None:
            raise AttributeError('No elevation data! Use load_elevation first')
        valid = np.ones_like(self.elevation._griddata)
        valid[np.isnan(self.elevation._griddata)] = np.nan
        valid[self.elevation._griddata <= 0] = np.nan
        return valid

class LocalRelief(BaseSpatialGrid):
    """Elevation range within a circular window.

    ::

        relief = LocalRelief(elevation=dem, pixel_radius=25)
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('elevation','pixel_radius'), '_create_from_elevation_and_radius'))

    def _create_from_elevation_and_radius(self, *args, **kwargs):
        from scipy.ndimage import maximum_filter, minimum_filter

        elevation = kwargs['elevation']
        pixel_radius = int(kwargs['pixel_radius'])
        self._copy_info_from_grid(elevation, True)

        # np.ogrid[-r:r] yields 2r samples, so the window used to be one cell
        # short on the far side and therefore off-centre.
        y, x = np.ogrid[-pixel_radius:pixel_radius + 1, -pixel_radius:pixel_radius + 1]
        footprint = x * x + y * y <= pixel_radius * pixel_radius

        data = np.asarray(elevation._griddata, dtype=float64)
        maximum = maximum_filter(data, footprint=footprint)
        minimum = minimum_filter(data, footprint=footprint)
        self._griddata = maximum - minimum

class Mask(BaseSpatialGrid):
    """0/1 grid marking the cells that drain to a set of outlets.

    ::

        mask = Mask(flow_direction=d8, outlets=[(x, y), ...])
    """

    dtype = uint8

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('flow_direction','outlets'), '_create_from_flow_direction_and_outlets'))

    def _create_from_flow_direction_and_outlets(self, *args, **kwargs):
        flow_direction = kwargs['flow_direction']
        outlets = kwargs['outlets']
        self._copy_info_from_grid(flow_direction, True)
        # One linear pass over the flow network for all the outlets at once,
        # rather than a separate upstream walk per outlet.
        self._griddata = flow_direction.basin_mask(outlets).astype(self.dtype)

    def perform_opening(self, structure = None, iterations = 1):
        """Morphological opening: erode then dilate, removing specks."""
        from scipy.ndimage import binary_opening
        self._griddata = binary_opening(
            self._griddata, iterations=iterations, structure=structure).astype(self.dtype)

    def perform_erosion(self, structure = None, iterations = 1):
        """Morphological erosion: shrink the mask by ``iterations`` cells."""
        from scipy.ndimage import binary_erosion
        self._griddata = binary_erosion(
            self._griddata, iterations=iterations, structure=structure).astype(self.dtype)

class LogArea(BaseSpatialGrid):
    """Base-10 logarithm of drainage area."""

    # float64, not uint8: log10(area) is a real number, and declaring it as a
    # byte meant `save()` wrote a grid of 6s and 7s.
    dtype = float64

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('area',), '_create_from_area'))

    def _create_from_area(self, *args, **kwargs):
        area = kwargs['area']
        self._copy_info_from_grid(area, True)
        data = np.asarray(area._griddata, dtype=float64)
        with np.errstate(divide='ignore', invalid='ignore'):
            self._griddata = np.where(data > 0, np.log10(np.where(data > 0, data, 1.0)), np.nan)

class Hillshade(CalculationMixin, BaseSpatialGrid):
    """Shaded relief, following the ESRI hillshade definition.

    ::

        hs = Hillshade(elevation=dem, azimuth=315, inclination=45)

    ``azimuth`` is the compass bearing of the light source in degrees
    (0 = north, 90 = east) and ``inclination`` its height above the horizon.
    Values run 0-255.
    """

    dtype = uint8

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('elevation','azimuth', 'inclination'), '_create_from_elevation'))

    def _create_from_elevation(self, *args, **kwargs):

        elevation = kwargs['elevation']
        az = kwargs['azimuth']
        elev = kwargs['inclination']

        self._copy_info_from_grid(elevation)
        self.calcHillshade(az, elev)

    def calcHillshade(self, az, elev, z_factor=1.0):
        """Fill the grid with an ESRI-style hillshade.

        Parameters
        ----------
        az : float
            Illumination azimuth, degrees clockwise from north.
        elev : float
            Illumination altitude above the horizon, in degrees.
        z_factor : float
            Vertical exaggeration.
        """
        azRad = (360.0 - az + 90.0) * np.pi / 180.0
        zenithRad = (90.0 - elev) * np.pi / 180.0

        # Sx is dz/dx (positive east); Sy is dz/drow (positive south), which
        # is the same sign convention ESRI's [dz/dy] uses.
        Sx, Sy = self._calcFiniteSlopes(self._griddata, self._mean_pixel_dimension(),
                                        self._georef_info.nx, self._georef_info.ny)
        Sx = Sx * z_factor
        Sy = Sy * z_factor

        # ESRI's aspect is atan2(dz/dy, -dz/dx). Without the negation the
        # illumination came out mirrored east-west.
        AspectRad = np.arctan2(Sy, -Sx)
        SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))

        shade = 255.0 * ((np.cos(zenithRad) * np.cos(SmagRad))
                         + (np.sin(zenithRad) * np.sin(SmagRad)
                            * np.cos(azRad - AspectRad)))
        # Slopes facing away from the light give a negative value; ESRI
        # clamps those to 0 rather than letting them wrap round the byte.
        self._griddata = np.clip(np.nan_to_num(shade), 0.0, 255.0).astype(self.dtype)

class GeographicHillshade(GeographicGridMixin, Hillshade):
    """Hillshade on a latitude/longitude grid."""

    pass

class MaxSlope(CalculationMixin, BaseSpatialGrid):
    """Magnitude of the topographic gradient, from centred differences.

    ::

        slope = MaxSlope(elevation=dem)

    Values are dimensionless (rise over run); take ``arctan`` for degrees.
    """

    dtype = float64

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('elevation',), '_create_from_elevation'))

    def _create_from_elevation(self, *args, **kwargs):

        elevation = kwargs['elevation']

        self._copy_info_from_grid(elevation)
        self.calcSlope()

    def calcSlope(self):
        """Replace the grid with the gradient magnitude of its own values."""
        Sx, Sy = self._calcFiniteSlopes(self._griddata, self._mean_pixel_dimension(),
                                        self._georef_info.nx, self._georef_info.ny)
        self._griddata = np.sqrt(Sx**2 + Sy**2)

class GeographicMaxSlope(GeographicGridMixin, MaxSlope):
    """Gradient magnitude on a latitude/longitude grid."""

    pass

class Laplacian(CalculationMixin, BaseSpatialGrid):
    """Laplacian of elevation -- positive in valleys, negative on ridges.

    ::

        curvature = Laplacian(elevation=dem)

    Uses SciPy's 5-point stencil, scaled by ``1 / dx**2``, so the units are
    inverse length.
    """

    dtype = float64

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('elevation',), '_create_from_elevation'))

    def _create_from_elevation(self, *args, **kwargs):

        from scipy.ndimage import laplace
        elevation = kwargs['elevation']

        self._copy_info_from_grid(elevation)
        # _mean_pixel_dimension(), not _georef_info.dx: on a geographic grid
        # dx is in DEGREES, and dividing by it gave a curvature ten orders of
        # magnitude too large.  The mixin overrides this with the true metric
        # cell size; for a projected grid the two are the same number.
        cell = self._mean_pixel_dimension(*args, **kwargs)
        self._griddata = (laplace(np.asarray(elevation._griddata, dtype=float64))
                          / np.power(cell, 2))


class GeographicLaplacian(GeographicGridMixin, Laplacian):
    """Laplacian of elevation on a latitude/longitude grid."""

    pass

class PriorityQueueMixIn(object):
    """Depression filling by Priority-Flood.

    Mix into a grid class to give it :meth:`_flood`, which removes closed
    depressions so that every cell has a downhill path to an outlet.

    The algorithms are those of

        Barnes, R., Lehman, C., Mulla, D. (2014). "Priority-flood: An optimal
        depression-filling and watershed-labeling algorithm for digital
        elevation models." *Computers & Geosciences* 62, 117-127.

    and run in :mod:`TopoAnalysis.kernels`, in C++ where the extension is
    built and in NumPy otherwise.

    Attributes
    ----------
    aggradation_slope : float
        Gradient of the increment added while filling.  The default of 1e-12
        makes filled surfaces very slightly convergent, which gives D8 routing
        something to follow across what would otherwise be a dead flat
        plateau.  Set it to 0 for the flat fill that TopoToolbox's
        ``fillsinks`` produces.
    """

    aggradation_slope = 1E-12

    class priorityQueue:
        """Stably-ordered priority queue.

        Retained because external code constructs it directly.  The filling
        routines no longer use it -- they call the compiled kernel instead --
        so it is kept purely for compatibility.
        """

        def __init__(self):
            import heapq
            self._heapq = heapq
            self.__pq = []
            self.__counter = 0
            self.__nItems = 0

        def get(self):
            """Pop the lowest-priority item; returns ``(priority, item)``."""
            priority, count, item = self._heapq.heappop(self.__pq)
            self.__nItems -= 1
            return priority, item

        def put(self, priority, item):
            """Insert ``item``.  Equal priorities come back out in FIFO order."""
            self.__counter += 1
            self.__nItems += 1
            self._heapq.heappush(self.__pq, [priority, self.__counter, item])

        def isEmpty(self):
            return self.__nItems == 0

    def randomize_subbasins_with_mask(self, *args, **kwargs):
        """Replace elevations inside a mask with random values and re-flood.

        Used to generate synthetic drainage networks with the same outlets as
        a real one.
        """
        mask = args[0]
        outlets = args[1]

        # `self.__flood` inside this class mangles to _PriorityQueueMixIn__flood,
        # which does not exist -- the method is _flood.
        self._flood(mask=mask, outlets=outlets, randomize=True)

    def _flood(self, *args, **kwargs):
        """Fill depressions in ``self._griddata``, in place.

        Keyword arguments
        -----------------
        mask : BaseSpatialGrid, optional
            Only cells where the mask equals 1 take part.
        outlets : sequence of (x, y), optional
            Flood inwards from these points instead of from the DEM edge.
        randomize : bool, optional
            Randomise elevations inside the mask before flooding.
        maximum_pit_depth : float, optional
            Leave depressions deeper than this unfilled.
        binary_result : bool, optional
            Replace the grid with a 0/1 map of the cells the flood reached.
        clip_to_fill : bool, optional
            Set cells the flood never reached to ``NaN``.

        Returns
        -------
        dict
            ``{'cells_filled': int, 'max_fill_depth': float,
            'visited': ndarray|None}``.
        """
        elevations = np.ascontiguousarray(self._griddata, dtype=float64)
        # Kept for the clip_to_fill/maximum_pit_depth comparison below.
        original = elevations.copy() if kwargs.get('clip_to_fill') else None

        mask = kwargs.get('mask')
        closed = None
        if mask is not None:
            closed = (np.asarray(mask._griddata) != 1).astype(np.uint8)

        if kwargs.get('randomize') is True:
            self._griddata = elevations
            self._randomize_grid_values(mask=mask)
            elevations = np.ascontiguousarray(self._griddata, dtype=float64)
            if original is not None:
                original = elevations.copy()

        seeds = None
        outlets = kwargs.get('outlets')
        if outlets is not None:
            flat = []
            for (row, col) in self._xy_to_rowscols(outlets):
                if row is None:
                    continue
                flat.append(row * self._georef_info.nx + col)
            if not flat:
                raise Error.InputError('_flood', 'none of the outlets lie inside the grid')
            seeds = np.array(flat, dtype=np.int64)

        want_visited = bool(kwargs.get('binary_result') or kwargs.get('clip_to_fill'))
        max_pit_depth = kwargs.get('maximum_pit_depth') or 0.0

        kernels.warn_if_slow(elevations.size)
        result = kernels.priority_flood(
            elevations,
            closed=closed,
            seeds=seeds,
            # A non-zero aggradation slope is Barnes et al.'s Algorithm 4
            # (Priority-Flood+epsilon); zero gives the flat fill of
            # Algorithm 2, which is what TopoToolbox's fillsinks returns.
            mode='epsilon' if self.aggradation_slope > 0 else 'flat',
            epsilon=self.aggradation_slope,
            cellsize=self._georef_info.dx,
            max_pit_depth=max_pit_depth,
            track_visited=want_visited,
        )

        self._griddata = elevations.astype(self.dtype, copy=False)

        visited = result['visited']
        if kwargs.get('binary_result'):
            self._griddata = np.asarray(visited, dtype=self.dtype)
        elif kwargs.get('clip_to_fill') is True:
            self._griddata = self._griddata.astype(float64, copy=True)
            self._griddata[np.asarray(visited) == 0] = np.nan
            if max_pit_depth:
                # A depth limit leaves whole depressions unfilled that the
                # flood nonetheless traverses, so `visited` no longer marks
                # them.  Compare against an unlimited fill: the cells that
                # would have risen but did not are exactly the depressions
                # that were declined, which is what the old implementation
                # excised here.
                unlimited = np.ascontiguousarray(
                    np.array(original, dtype=float64, copy=True))
                kernels.priority_flood(
                    unlimited,
                    closed=None if closed is None else closed.copy(),
                    seeds=seeds,
                    mode='epsilon' if self.aggradation_slope > 0 else 'flat',
                    epsilon=self.aggradation_slope,
                    cellsize=self._georef_info.dx)
                declined = unlimited > self._griddata + 1e-12
                self._griddata[declined] = np.nan

        return result

class PriorityFillGrid(PriorityQueueMixIn, BaseSpatialGrid):
    """Cells of a mask that are connected to a set of outlets.

    Floods outwards from ``outlets`` through the cells where ``mask`` is 1
    and returns a 0/1 grid of everything reached.  Used to grow a valley
    network from a curvature mask, keeping only the parts that connect to a
    real channel.
    """

    dtype = uint8

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('outlets', 'mask'), '_create_from_outlets_and_mask'))

    def _create_from_outlets_and_mask(self, *args, **kwargs):

        mask = kwargs['mask']
        outlets = kwargs['outlets']
        kwargs = dict(kwargs)
        kwargs['randomize'] = False
        kwargs['binary_result'] = True

        self._copy_info_from_grid(mask, True)
        # Flooding needs a surface to descend; give the outlets a high value
        # so the flood spreads from them through the rest of the mask.
        self._griddata = np.zeros((self._georef_info.ny, self._georef_info.nx),
                                  dtype=float64)
        for (row, col) in self._xy_to_rowscols(outlets):
            if row is not None:
                self._griddata[row, col] = 100.0
        self._griddata[np.asarray(mask._griddata) != 1] = 0.0

        self._flood(*args, **kwargs)
        self._griddata = np.asarray(self._griddata, dtype=self.dtype)

class FilledElevation(PriorityQueueMixIn, Elevation):
    """A DEM with its closed depressions removed.

    ::

        filled = FilledElevation(elevation=dem)

    Filling uses Priority-Flood (Barnes et al. 2014).  By default the filled
    surface carries a gradient of
    :attr:`PriorityQueueMixIn.aggradation_slope` across former depressions so
    that D8 routing has a direction to follow; pass
    ``aggradation_slope=0`` -- or set the class attribute -- for the flat fill
    that TopoToolbox's ``fillsinks`` produces.

    Other keywords (``mask``, ``outlets``, ``maximum_pit_depth``,
    ``clip_to_fill``) are passed to :meth:`PriorityQueueMixIn._flood`.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('elevation',), '_create_from_elevation'))

    def _create_from_elevation(self, *args, **kwargs):

        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation)
        if 'aggradation_slope' in kwargs:
            # Per-instance override, so one call can ask for a flat fill
            # without changing the default for every other grid.
            self.aggradation_slope = kwargs['aggradation_slope']
        self.fill_report = self._flood(*args, **kwargs)


class Area(BaseSpatialGrid):
    """D8 contributing (drainage) area.

    ::

        area = Area(flow_direction=d8)

    Values are in squared map units: each cell contributes its own area
    (``dx**2``, or the true spherical cell area for a
    :class:`GeographicArea`) plus everything upstream of it.

    To count contributing *cells* instead -- which is what TopoToolbox's
    ``flowacc`` returns -- divide by the per-cell area, or pass
    ``weights=np.ones(...)``.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('flow_direction',), '_create_from_flow_direction'),
                                   (('dx', 'grid'), '_create_from_grid'))

    def _create_from_flow_direction(self, *args, **kwargs):

        flow_dir = kwargs['flow_direction']
        self._copy_info_from_grid(flow_dir, True)
        self._calc_d8_area(*args, **kwargs)

    def _calc_d8_area(self, *args, **kwargs):
        """Accumulate per-cell weights down the D8 network.

        Keyword arguments
        -----------------
        flow_direction : FlowDirectionD8
        weights : ndarray, optional
            Per-cell contribution. Defaults to the cell area.
        mask, evaluate_at : BaseSpatialGrid, optional
            Cells where these are 0 contribute nothing and pass nothing on.
        """
        flow_dir = kwargs['flow_direction']
        kernels.warn_if_slow(flow_dir._griddata.size)

        weights = kwargs.get('weights')
        if weights is None:
            weights = self._area_per_pixel(*args, **kwargs)
        weights = np.ascontiguousarray(weights, dtype=float64)

        gate = None
        for gate_name in ('evaluate_at', 'mask'):
            gate_grid = kwargs.get(gate_name)
            if gate_grid is None:
                continue
            # `mask[i, j] is not None` was always true -- __getitem__ returns
            # the value, not None, for a masked-out cell -- so a supplied
            # mask never actually excluded anything.
            this_gate = (np.asarray(gate_grid._griddata) != 0).astype(np.uint8)
            gate = this_gate if gate is None else (gate & this_gate)

        # Ordering comes from the flow network itself rather than from an
        # elevation sort, so it stays correct for a flow-direction grid
        # loaded from disk, where `sort()` would have ordered cells by their
        # direction *code*.
        receivers = flow_dir.receivers
        order = flow_dir.topological_order
        self._griddata = np.asarray(
            kernels.accumulate(receivers, order, weights, gate))

    # Preserved under its original name-mangled spelling for any external
    # caller that reached in for it.
    _Area__calcD8Area = _calc_d8_area

    def areas_greater_than(self, min_area):
        """Map coordinates of every cell whose area is at least ``min_area``."""
        rows, cols = np.where(self._griddata >= min_area)
        return self._rowscols_to_xy(list(zip(rows.tolist(), cols.tolist())))

    def areas_between(self, fd, min_area, max_area):
        """Channel heads: cells in an area band whose donors are all below it.

        These are the points where drainage area first crosses ``min_area``,
        so they mark the upstream ends of the channel network.
        """
        ij_out = []
        rows, cols = np.where((self._griddata >= min_area) & (self._griddata <= max_area))
        for (i, j) in zip(rows.tolist(), cols.tolist()):
            options = fd.get_upstream_cell_indexes(i, j)
            # No upstream neighbour at all means this cell *is* a source, so
            # it qualifies. The original compared None < min_area and raised.
            if not options:
                ij_out.append((i, j))
                continue
            max_a = max(self._griddata[i_p, j_p] for (i_p, j_p) in options)
            if max_a < min_area:
                ij_out.append((i, j))
        return self._rowscols_to_xy(ij_out)


class GeographicArea(GeographicGridMixin, Area):
    """Drainage area on a latitude/longitude grid, using true cell areas."""

    pass

class _ValleyMaskMixin(object):
    """Shared construction of a connected valley-floor mask.

    Cells whose Laplacian exceeds ``valley_laplace_value`` are convergent
    (valley-like).  Those that are connected to a cell with more than
    ``min_area_value`` of drainage area form the valley network; isolated
    convergent patches on hillslopes are discarded.
    """

    def _valley_mask(self, *args, **kwargs):
        from scipy.ndimage import binary_dilation, binary_erosion

        mask = Mask()
        mask._copy_info_from_grid(kwargs['laplace'], True)
        mask._griddata[kwargs['laplace']._griddata >= kwargs['valley_laplace_value']] = 1

        outlets = kwargs['area'].areas_greater_than(kwargs['min_area_value'])
        pfg = PriorityFillGrid(mask=mask, outlets=outlets)

        iterations = kwargs.get('iterations') or 0
        if iterations:
            # Closing: dilate then erode, to bridge one-cell gaps in the
            # valley network without growing it overall.
            pfg._griddata = binary_dilation(pfg._griddata, iterations=iterations)
            pfg._griddata = binary_erosion(
                pfg._griddata, iterations=iterations).astype(pfg.dtype)
        return pfg

    def _create_from_flow_direction_and_valley_slope_and_area(self, *args, **kwargs):
        kwargs = dict(kwargs)
        kwargs['evaluate_at'] = self._valley_mask(*args, **kwargs)
        self._create_from_flow_direction(*args, **kwargs)


class ValleyArea(_ValleyMaskMixin, Area):
    """Drainage area accumulated only over valley-floor cells.

    ::

        va = ValleyArea(flow_direction=d8, area=area, laplace=laplacian,
                        valley_laplace_value=0.01, min_area_value=1e6)
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('flow_direction', 'area', 'laplace', 'valley_laplace_value', 'min_area_value'), '_create_from_flow_direction_and_valley_slope_and_area'))


class GeographicValleyArea(GeographicGridMixin, ValleyArea):
    """Valley-floor area on a latitude/longitude grid."""

    pass


class MainstemValleyArea(_ValleyMaskMixin, Area):
    """Valley-floor area accumulated along the main stem only.

    Where :class:`ValleyArea` sums every tributary, this follows only the
    longest flow path into each cell, so the result is the valley area of the
    trunk stream rather than of the whole basin.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('area', 'laplace', 'valley_laplace_value', 'min_area_value', 'flow_direction'), '_create_from_flow_direction_and_valley_slope_and_area'))

    def _calc_d8_area(self, *args, **kwargs):
        """Accumulate down the longest-flow-length path only.

        Notes
        -----
        The original version of this method was named ``__calcD8Area``, so
        Python mangled it to ``_MainstemValleyArea__calcD8Area`` while the
        inherited caller looked for ``_Area__calcD8Area``.  It was therefore
        never called, and ``MainstemValleyArea`` silently behaved exactly like
        :class:`ValleyArea`.  It also compared a whole grid object against a
        scalar, which would have raised as soon as it did run.
        """
        flow_dir = kwargs['flow_direction']
        kernels.warn_if_slow(flow_dir._griddata.size)

        dA = np.ascontiguousarray(self._area_per_pixel(*args, **kwargs), dtype=float64)
        step = flow_dir.step_lengths(self._mean_pixel_dimension(*args, **kwargs))

        gate = None
        for gate_name in ('evaluate_at', 'mask'):
            gate_grid = kwargs.get(gate_name)
            if gate_grid is not None:
                this_gate = (np.asarray(gate_grid._griddata) != 0).astype(np.uint8)
                gate = this_gate if gate is None else (gate & this_gate)
        if gate is not None:
            dA = np.where(gate != 0, dA, 0.0)

        receivers = flow_dir.receivers
        order = flow_dir.topological_order
        _, main_stem = kernels.flow_length(receivers, order,
                                           np.ascontiguousarray(step, dtype=float64))

        self._griddata = np.asarray(
            kernels.propagate_along_main_stem(receivers, order, main_stem, gate, dA, 'sum'))

    _Area__calcD8Area = _calc_d8_area


class GeographicMainstemValleyArea(GeographicGridMixin, MainstemValleyArea):
    """Main-stem valley area on a latitude/longitude grid."""

    pass

def _chi_profile_along_path(points, area_grid, elevation_grid, de, theta):
    """Chi and elevation sampled along an ordered list of cells.

    ``points`` is a 2 x n integer array of rows and columns running
    downstream to upstream.  The integral uses the trapezoidal rule over the
    along-path distance, with diagonal steps counted as sqrt(2) cell widths.

    Returns
    -------
    (ndarray, ndarray)
        ``(chi, elevation)``, both length n, with ``chi[0] == 0``.
    """
    n = points.shape[1]
    adjustment = np.ones(n)
    if n > 2:
        diagonal = ((points[0, 1:-1] != points[0, 2:])
                    & (points[1, 1:-1] != points[1, 2:]))
        adjustment[np.flatnonzero(diagonal) + 1] += (SQRT2 - 1.0)

    area_profile = area_grid[points[0], points[1]]
    elevation_profile = elevation_grid[points[0], points[1]]
    de_profile = de[points[0], points[1]]

    chi_profile = np.zeros_like(elevation_profile, dtype=float64)
    if n > 1:
        integrand = 0.5 * (np.power(area_profile[1:], -theta)
                           + np.power(area_profile[:-1], -theta))
        step = 0.5 * (de_profile[1:] + de_profile[:-1]) * adjustment[1:]
        chi_profile[1:] = np.cumsum(integrand * step)
    return chi_profile, elevation_profile


class AlongFlowSmoothing(object):
    """Collect the cells within a window along each flow path.

    Mixed into the ``*WithSmoothing`` classes, which fit a model over that
    window at every cell.  The window is set either by drop
    (``vertical_interval``) or by distance (``horizontal_interval``), and
    extends both upstream and downstream of the cell being evaluated.
    """


    def _find_points_along_path(self, de, **kwargs):
        area = kwargs['area']
        flow_direction = kwargs['flow_direction']
        elevation = kwargs['elevation']
        area_threshold = kwargs.get('area_threshold', 0)

        import time
        t1 = time.time()

        upstream_i = (np.ones_like(area._griddata) * -1.0).astype(int)
        upstream_j = (np.ones_like(area._griddata) * -1.0).astype(int)
        downstream_i = (np.ones_like(area._griddata) * -1.0).astype(int)
        downstream_j = (np.ones_like(area._griddata) * -1.0).astype(int)
        shape = upstream_i.shape
        indexes = area.sort(reverse=False)
        (i_s, j_s) = np.unravel_index(indexes, shape)

        def set_usdsindexes(ij):
            i, j = ij
            ds_i, ds_j, good = flow_direction.get_flow_to_cell(i, j)
            if good:
                (downstream_i[i, j], downstream_j[i, j]) = (ds_i, ds_j)
                (upstream_i[ds_i, ds_j], upstream_j[ds_i, ds_j]) = (i, j)

        list(map(set_usdsindexes, zip(i_s, j_s)))

        t2 = time.time()
        print('completed flow graph in: ' + str(t2 - t1) + " s")

        vertical_interval = kwargs.get('vertical_interval', None)
        if vertical_interval is not None:
            def find_points_along_path(this_i, this_j):
                ret = list()
                ret += [(this_i, this_j)]
                (ups_i, ups_j) = (upstream_i[this_i, this_j], upstream_j[this_i, this_j])
                (ds_i, ds_j) = (downstream_i[this_i, this_j], downstream_j[this_i, this_j])
                if (ups_i < 0) or (ds_i < 0):
                    return None
                ret = [(ds_i, ds_j)] + ret + [(ups_i, ups_j)]
                delta_e = elevation._griddata[ups_i, ups_j] - elevation._griddata[ds_i, ds_j]

                while (delta_e < vertical_interval) & (ups_i >= 0) & (area._griddata[ups_i, ups_j] > area_threshold):
                    (ups_i, ups_j) = (upstream_i[ups_i, ups_j], upstream_j[ups_i, ups_j])
                    (ds_i, ds_j) = (downstream_i[ds_i, ds_j], downstream_j[ds_i, ds_j])
                    if (ups_i < 0) or (ds_i < 0) or ((ups_i == ds_i) and (ups_j == ds_j)):
                        return None
                    ret = [(ds_i, ds_j)] + ret + [(ups_i, ups_j)]
                    delta_e = elevation._griddata[ups_i, ups_j] - elevation._griddata[ds_i, ds_j]
                return ret
            return find_points_along_path
        else:
            horizontal_interval = kwargs['horizontal_interval']

            def find_points_along_path(this_i, this_j):
                horizontal_distance = 0
                ret = list()
                ret += [(this_i, this_j)]
                (ups_i, ups_j) = (upstream_i[this_i, this_j], upstream_j[this_i, this_j])
                (ds_i, ds_j) = (downstream_i[this_i, this_j], downstream_j[this_i, this_j])
                if (ups_i < 0) or (ds_i < 0):
                    return None
                ret = [(ds_i, ds_j)] + ret + [(ups_i, ups_j)]
                horizontal_distance += ((1.0 if ((ds_i == this_i) or (ds_j == this_j)) else 1.414) + (
                    1.0 if ((ups_i == this_i) or (ups_j == this_j)) else 1.414)) * de[this_i, this_j]
                while (horizontal_distance < horizontal_interval) & (ups_i >= 0) & (area._griddata[ups_i, ups_j] > area_threshold):
                    (ups_i, ups_j) = (upstream_i[ups_i, ups_j], upstream_j[ups_i, ups_j])
                    (ds_i, ds_j) = (downstream_i[ds_i, ds_j], downstream_j[ds_i, ds_j])
                    if (ups_i < 0) or (ds_i < 0) or ((ups_i == ds_i) and (ups_j == ds_j)):
                        return None
                    horizontal_distance += (1.0 if ((ds_i == ret[0][0]) or (ds_j == ret[0][1])) else 1.414) * de[
                        ds_i, ds_j] + (1.0 if ((ups_i == ret[-1][0]) or (ups_j == ret[-1][1])) else 1.414) * de[
                                               ups_i, ups_j]
                    ret = [(ds_i, ds_j)] + ret + [(ups_i, ups_j)]
                return ret
            return find_points_along_path
    
class KsFromChiWithSmoothing(BaseSpatialGrid, AlongFlowSmoothing):
    """Channel steepness fitted to chi--elevation in a moving window.

    ::

        ks = KsFromChiWithSmoothing(elevation=dem, area=area,
                                    flow_direction=d8, theta=0.45,
                                    vertical_interval=50.0)

    At every channel cell, chi and elevation are extracted over a window
    along the flow path and a line through the origin is fitted; its slope
    is k_s.  ``save()`` writes seven bands: ks, the number of profiles
    crossing each cell, mean squared error, sum of squares, r-squared,
    p-value, and the number of points in the regression.

    Needs ``statsmodels``.  Pass ``verbose=False`` to silence the progress
    ticker.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('elevation', 'area', 'flow_direction', 'theta', 'vertical_interval'), '_create_from_elevation_area_flow_direction'),
                                   (('elevation', 'area', 'flow_direction', 'theta', 'horizontal_interval'), '_create_from_elevation_area_flow_direction'),
                                   )
            
    def _create_from_elevation_area_flow_direction(self, *args, **kwargs):

        sm = _require_statsmodels()

        elevation = kwargs['elevation']
        area = kwargs['area']
        theta = kwargs['theta']
        area_threshold = kwargs.get('area_threshold', 0)
        de = area._mean_pixel_dimension()

        self._copy_info_from_grid(elevation)
        self._griddata = np.zeros_like(elevation._griddata)
        self._n = np.zeros_like(self._griddata).astype(int)
        self._n_regression = np.zeros_like(self._griddata).astype(int)
        self._mse = np.zeros_like(self._griddata)
        self._ss = np.zeros_like(self._griddata)
        self._r2 = np.zeros_like(self._griddata)
        self._pval = np.zeros_like(self._griddata)
        self._griddata[:] = np.nan

        find_points_along_path = self._find_points_along_path(de, **kwargs)

        def calc_ks(i, j):
            points = find_points_along_path(i, j)
            if points is None:
                return np.nan, np.nan, np.nan, np.nan, [[], []], np.nan, 0

            points = np.array(list(zip(*points))).astype(int)
            chi_profile, elevation_profile = _chi_profile_along_path(
                points, area._griddata, elevation._griddata, de, theta)

            X = np.array(chi_profile)
            y = np.array(elevation_profile) - elevation_profile[0]
            res = sm.OLS(y, X).fit()
            SS = res.ssr
            return (res.params[0], SS / float(len(chi_profile)), SS, res.rsquared,
                    points, res.pvalues[0], len(chi_profile))

        rows, cols = np.where(~np.isnan(area._griddata) & ~np.isnan(elevation._griddata)
                              & (area._griddata > area_threshold))
        ij = list(zip(rows, cols))
        progress = _Progress(len(ij), kwargs.get('verbose', True))
        for (i, j) in ij:
            (self._griddata[i, j], self._mse[i, j], self._ss[i, j], self._r2[i, j],
             pts, self._pval[i, j], self._n_regression[i, j]) = calc_ks(i, j)
            # np.add.at, not `+=`: with repeated subscripts a fancy-indexed
            # `+=` applies each increment only once, so cells crossed by
            # several profiles were undercounted.
            if len(pts[0]):
                np.add.at(self._n, (pts[0], pts[1]), 1)
            progress.tick()
        progress.done()


    def save(self, filename):
        
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', [self._griddata, self._n, self._mse, self._ss, self._r2, self._pval, self._n_regression], self.dtype, filename, ['COMPRESS=LZW', 'BIGTIFF=YES'], multiple_bands=True)
            
    @classmethod
    def load(cls, filename):
        """Read back the multi-band grid written by :meth:`save`."""
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.nan
            return grid
        
        return_object = cls()
        gdal_dataset = gdal.Open(filename)
        
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        
        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[0]+return_object._georef_info.dx/2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
        return_object._griddata = get_band(gdal_dataset, 1)
        return_object._n = get_band(gdal_dataset, 2)
        return_object._mse = get_band(gdal_dataset, 3)
        return_object._ss = get_band(gdal_dataset, 4)
        return_object._r2 = get_band(gdal_dataset, 5)
        return_object._pval = get_band(gdal_dataset, 6)
        return_object._n_regression = get_band(gdal_dataset, 7)
        
            
        gdal_file = None
        return return_object


class ThetaFromChiWithSmoothing(BaseSpatialGrid, AlongFlowSmoothing):
    """Concavity fitted to chi--elevation in a moving window.

    As :class:`KsFromChiWithSmoothing`, but ``theta`` is optimised at each
    cell rather than held fixed: the value returned is the one that makes
    the chi--elevation relation most nearly linear.  The search is bounded
    at +/-10.

    Needs ``statsmodels``.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction', 'vertical_interval', 'min_area'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation', 'area', 'flow_direction', 'horizontal_interval', 'min_area'),
                                    '_create_from_elevation_area_flow_direction'),
                                   )

    def _create_from_elevation_area_flow_direction(self, *args, **kwargs):

        from scipy.optimize import fmin

        sm = _require_statsmodels()

        elevation = kwargs['elevation']
        area = kwargs['area']
        min_area = kwargs['min_area']
        de = area._mean_pixel_dimension()

        self._copy_info_from_grid(elevation)
        self._griddata = np.zeros_like(elevation._griddata)
        self._n = np.zeros_like(self._griddata).astype(int)
        self._n_regression = np.zeros_like(self._griddata).astype(int)
        self._mse = np.zeros_like(self._griddata)
        self._ss = np.zeros_like(self._griddata)
        self._r2 = np.zeros_like(self._griddata)
        self._pval = np.zeros_like(self._griddata)
        self._griddata[:] = np.nan

        find_points_along_path = self._find_points_along_path(de, **kwargs)
        empty = (np.nan, np.nan, np.nan, np.nan, [[], []], np.nan, 0)

        def calc_theta(i, j):

            points = find_points_along_path(i, j)
            if points is None:
                return empty

            points = np.array(list(zip(*points))).astype(int)

            def fit_for_theta(theta):
                chi_profile, elevation_profile = _chi_profile_along_path(
                    points, area._griddata, elevation._griddata, de, theta)
                X = np.array(chi_profile)
                y = np.array(elevation_profile) - elevation_profile[0]
                return sm.OLS(y, X).fit(), chi_profile

            def ssr_for_theta(theta):
                theta = float(np.atleast_1d(theta)[0])
                if theta > 10 or theta < -10:
                    return np.inf
                return fit_for_theta(theta)[0].ssr

            try:
                theta_bf, funval, iterations, funcalls, warnflag = fmin(
                    ssr_for_theta, np.array([0.5]), xtol=1E-5, ftol=1E-5,
                    maxiter=100, maxfun=200, full_output=True, disp=False)
                res, chi_profile = fit_for_theta(float(theta_bf[0]))
                SS = res.ssr
                return (float(theta_bf[0]), SS / float(len(chi_profile)), SS,
                        res.rsquared, points, res.pvalues[0], len(chi_profile))
            except (ValueError, np.linalg.LinAlgError, ZeroDivisionError, FloatingPointError):
                # A bare `except:` here also swallowed KeyboardInterrupt and
                # genuine programming errors.
                return empty

        rows, cols = np.where((area._griddata != 0) & ~np.isnan(area._griddata)
                              & ~np.isnan(elevation._griddata) & (area._griddata > min_area))
        ij = list(zip(rows, cols))
        progress = _Progress(len(ij), kwargs.get('verbose', True))
        for (i, j) in ij:
            (self._griddata[i, j], self._mse[i, j], self._ss[i, j], self._r2[i, j],
             pts, self._pval[i, j], self._n_regression[i, j]) = calc_theta(i, j)
            if len(pts[0]):
                np.add.at(self._n, (pts[0], pts[1]), 1)
            progress.tick()
        progress.done()

    def save(self, filename):

        self._create_gdal_representation_from_array(self._georef_info, 'GTiff',
                                                    [self._griddata, self._n, self._mse, self._ss, self._r2, self._pval,
                                                     self._n_regression], self.dtype, filename,
                                                    ['COMPRESS=LZW', 'BIGTIFF=YES'], multiple_bands=True)

    @classmethod
    def load(cls, filename):
        """Read back the multi-band grid written by :meth:`save`."""
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.nan
            return grid

        return_object = cls()
        gdal_dataset = gdal.Open(filename)

        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize

        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[
                                                   0] + return_object._georef_info.dx / 2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3] - (
                    return_object._georef_info.dx * (ny - 0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny

        return_object._griddata = get_band(gdal_dataset, 1)
        return_object._n = get_band(gdal_dataset, 2)
        return_object._mse = get_band(gdal_dataset, 3)
        return_object._ss = get_band(gdal_dataset, 4)
        return_object._r2 = get_band(gdal_dataset, 5)
        return_object._pval = get_band(gdal_dataset, 6)
        return_object._n_regression = get_band(gdal_dataset, 7)

        gdal_file = None
        return return_object





class ChannelSlopeWithSmoothing(BaseSpatialGrid, AlongFlowSmoothing):
    """Channel gradient measured over a window along the flow path.

    Far less noisy than a cell-to-cell slope, which on a real DEM is mostly
    vertical quantisation.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   )

    scale_factor = 1.0

    def points_along_path(self, points, i, j):
        """Which part of the window to measure over; the whole of it here."""
        return points

    def calc_channel_slope(self, i, j, elevation, de, find_points_along_path):
        """Gradient across the window centred on ``(i, j)``.

        Returns ``NaN`` where the window runs off the end of the network.
        """
        points = find_points_along_path(i, j)

        if points is None or len(points) == 0:
            return np.nan

        points = self.points_along_path(points, i, j)

        if points is not None:

            pts = list(zip(*(points)))
            points = np.array(pts).astype(int)
            adjustment = np.ones((len(points[0])))
            ind = np.where((points[0, 1:-1] != points[0, 2:]) & (points[1, 1:-1] != points[1, 2:]))
            adjustment[ind[0] + 1] += 0.414
            elevation_profile = elevation._griddata[points[0], points[1]]
            de_profile = de[points[0], points[1]]
            dx = np.sum(de_profile * adjustment)
            dy = elevation_profile[-1] - elevation_profile[0]
            return dy / dx

        else:

            return np.nan

    def _create_from_elevation_area_flow_direction(self, *args, **kwargs):

        elevation = kwargs['elevation']
        area = kwargs['area']
        if kwargs.get('horizontal_interval') is not None:
            kwargs['horizontal_interval'] = kwargs['horizontal_interval']*self.scale_factor
        if kwargs.get('vertical_interval') is not None:
            kwargs['vertical_interval'] = kwargs['vertical_interval']*self.scale_factor
        min_area = kwargs.get('min_area', 1E6)

        de = elevation._mean_pixel_dimension()

        self._copy_info_from_grid(elevation)
        self._griddata = np.zeros_like(elevation._griddata)
        self._griddata[:] = np.nan

        find_points_along_path = self._find_points_along_path(de, **kwargs)

        rows, cols = np.where((area._griddata != 0) & ~np.isnan(area._griddata)
                              & ~np.isnan(elevation._griddata) & (area._griddata > min_area))
        ij = list(zip(rows, cols))
        progress = _Progress(len(ij), kwargs.get('verbose', True))
        for (i, j) in ij:
            self._griddata[i, j] = self.calc_channel_slope(i, j, elevation, de, find_points_along_path)
            progress.tick()
        progress.done()

class ChannelDownSlopeWithSmoothing(ChannelSlopeWithSmoothing):
    """Channel gradient over the downstream half of the window."""

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   )
    scale_factor = 2.0

    def points_along_path(self, points, i, j):
        """Only the downstream half of the window."""
        position_of_center = list([ind[0] for ind in zip(range(len(points)),points) if (ind[1][0] == i and ind[1][1] == j)])[0]
        return points[0:position_of_center] if len(points[0:position_of_center]) > 0 else None

class ChannelUpSlopeWithSmoothing(ChannelSlopeWithSmoothing):
    """Channel gradient over the upstream half of the window."""

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   )
    scale_factor = 2.0

    def points_along_path(self, points, i, j):
        """Only the upstream half of the window."""
        position_of_center = list([ind[0] for ind in zip(range(len(points)),points) if (ind[1][0] == i and ind[1][1] == j)])[0]
        return points[position_of_center:] if len(points[position_of_center:]) > 0 else None

class GeographicKsFromChiWithSmoothing(GeographicGridMixin, KsFromChiWithSmoothing):
    """Windowed steepness fit on a latitude/longitude grid."""

    pass

class GeographicThetaFromChiWithSmoothing(GeographicGridMixin, ThetaFromChiWithSmoothing):
    """Windowed concavity fit on a latitude/longitude grid."""

    pass

class MultiscaleCurvatureValleyWidth(BaseSpatialGrid):
    """Valley width from the scale at which curvature is most negative.

    A local quadratic is fitted at a range of window sizes; the width
    recorded at each cell is the window that minimises the smaller principal
    curvature, i.e. the scale at which the valley is best resolved.  The
    curvature itself is kept in ``_minC``.

    ::

        width = MultiscaleCurvatureValleyWidth(
            elevation=dem, area=area, area_cutoff=1e6,
            min_width=30.0, max_width=600.0)

    Pass ``use_dask=True`` to evaluate the scales in parallel.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('elevation', 'area', 'area_cutoff', 'max_width', 'min_width'), '_create_from_inputs'),)

    class Utilities(object):
        """The quadratic-fit machinery, as classmethods so dask can pickle it."""

        @classmethod
        def _calc_inv_G_for_kernel(cls, X, Y, N, fix_center=False):
            x4 = np.sum(np.power(X, 4)[:])
            x2y2 = np.sum((np.power(X, 2)[:] * np.power(Y, 2))[:])
            x2 = np.sum(np.power(X, 2)[:])

            if fix_center:
                G = np.asarray([[x4, x2y2, 0, 0, 0],
                                [x2y2, x4, 0, 0, 0],
                                [0, 0, x2y2, 0, 0],
                                [0, 0, 0, x2, 0],
                                [0, 0, 0, 0, x2]])
            else:
                G = np.asarray([[x4, x2y2, 0, 0, 0, x2],
                                [x2y2, x4, 0, 0, 0, x2],
                                [0, 0, x2y2, 0, 0, 0],
                                [0, 0, 0, x2, 0, 0],
                                [0, 0, 0, 0, x2, 0],
                                [x2, x2, 0, 0, 0, N]])

            from numpy.linalg import inv

            return inv(G)

        @classmethod
        def _convolve(cls, X, Y, Z, K, fix_center=False):
            """Cross-correlate the elevation grid with the five moment kernels.

            Each returned array is the windowed sum of ``z * basis`` with the
            centre cell's own contribution removed, which is what the normal
            equations for the local quadratic fit need.
            """
            from numpy.fft import fft2, ifft2, fftshift

            Xt = X
            Yt = Y
            data = np.asarray(Z._griddata, dtype=float64)

            # The kernels are symmetric about the grid centre, so correlation
            # and convolution coincide and the elevation grid can be used as
            # is.  The original wrote fliplr(fliplr(...)), which is the
            # identity -- presumably a typo for a 180-degree rotation, and
            # harmless only because of that symmetry.
            FZ = fft2(data)
            FX1 = fft2(np.power(Xt, 2))
            FX2 = fft2(np.power(Yt, 2))
            FX3 = fft2(Xt * Yt)
            FX4 = fft2(Xt)
            FX5 = fft2(Yt)
            FX6 = fft2(K)

            # fftshift, not ifftshift: the kernels are centred in the array,
            # so the convolution output must be rolled forward by half the
            # grid.  For the odd-sized grids this class constructs the two
            # differ by one cell, which displaced every fitted coefficient.
            def correlate(kernel_hat):
                return np.real(fftshift(ifft2(FZ * kernel_hat)))

            g = correlate(FX1) - np.sum(np.power(Xt, 2)) * data
            h = correlate(FX2) - np.sum(np.power(Yt, 2)) * data
            i = correlate(FX3) - np.sum(Xt * Yt) * data
            j = correlate(FX4) - np.sum(Xt) * data
            k = correlate(FX5) - np.sum(Yt) * data
            if not fix_center:
                l = correlate(FX6) - np.sum(K) * data
            else:
                l = None

            return g, h, i, j, k, l

        @classmethod
        def _Cmin(cls, H, g, h, i, j, k, l, fix_center=False):

            if fix_center:
                a = H[0, 0] * g + H[0, 1] * h + H[0, 2] * i + H[0, 3] * j + H[0, 4] * k
                b = H[1, 0] * g + H[1, 1] * h + H[1, 2] * i + H[1, 3] * j + H[1, 4] * k
                c = H[2, 0] * g + H[2, 1] * h + H[2, 2] * i + H[2, 3] * j + H[2, 4] * k
            else:
                a = H[0, 0] * g + H[0, 1] * h + H[0, 2] * i + H[0, 3] * j + H[0, 4] * k + H[0, 5] * l
                b = H[1, 0] * g + H[1, 1] * h + H[1, 2] * i + H[1, 3] * j + H[1, 4] * k + H[1, 5] * l
                c = H[2, 0] * g + H[2, 1] * h + H[2, 2] * i + H[2, 3] * j + H[2, 4] * k + H[2, 5] * l

            return -a - b - np.sqrt(np.power((a - b), 2) + np.power(c, 2))

        @classmethod
        def _calc_coefficients(cls, H, g, h, i, j, k, l, fix_center=False):
            if fix_center:
                a = H[0, 0] * g + H[0, 1] * h + H[0, 2] * i + H[0, 3] * j + H[0, 4] * k
                b = H[1, 0] * g + H[1, 1] * h + H[1, 2] * i + H[1, 3] * j + H[1, 4] * k
                c = H[2, 0] * g + H[2, 1] * h + H[2, 2] * i + H[2, 3] * j + H[2, 4] * k
                d = H[3, 0] * g + H[3, 1] * h + H[3, 2] * i + H[3, 3] * j + H[3, 4] * k
                e = H[4, 0] * g + H[4, 1] * h + H[4, 2] * i + H[4, 3] * j + H[4, 4] * k
                f = 0
            else:
                a = H[0, 0] * g + H[0, 1] * h + H[0, 2] * i + H[0, 3] * j + H[0, 4] * k + H[0, 5] * l
                b = H[1, 0] * g + H[1, 1] * h + H[1, 2] * i + H[1, 3] * j + H[1, 4] * k + H[1, 5] * l
                c = H[2, 0] * g + H[2, 1] * h + H[2, 2] * i + H[2, 3] * j + H[2, 4] * k + H[2, 5] * l
                d = H[3, 0] * g + H[3, 1] * h + H[3, 2] * i + H[3, 3] * j + H[3, 4] * k + H[3, 5] * l
                e = H[4, 0] * g + H[4, 1] * h + H[4, 2] * i + H[4, 3] * j + H[4, 4] * k + H[4, 5] * l
                f = H[5, 0] * g + H[5, 1] * h + H[5, 2] * i + H[5, 3] * j + H[5, 4] * k + H[5, 5] * l
            return a, b, c, d, e, f

        @classmethod
        def _calc_coefficients_for_scale(cls, Z, de, fix_center=False):
            X, Y, N, K = cls._create_kernel(Z, de)
            H = cls._calc_inv_G_for_kernel(X, Y, N, fix_center)
            g, h, i, j, k, l = cls._convolve(X, Y, Z, K, fix_center)
            return cls._calc_coefficients(H, g, h, i, j, k, l, fix_center)

        @classmethod
        def _create_kernel(cls, Z, de):

            center = (Z._georef_info.xllcenter + (Z._georef_info.dx / 2) * (Z._georef_info.nx - 1),
                      Z._georef_info.yllcenter + (Z._georef_info.dx / 2) * (Z._georef_info.ny - 1))
            x = np.arange(Z._georef_info.nx) * Z._georef_info.dx + Z._georef_info.xllcenter - center[0]
            y = np.arange(Z._georef_info.ny) * Z._georef_info.dx + Z._georef_info.yllcenter - center[1]

            X, Y = np.meshgrid(x, y)
            K = ((np.abs(X) <= (de)) & (np.abs(Y) <= (de))).astype(float)

            N = np.sum(K[:])
            return X * K, Y * K, N, K

        @classmethod
        def _Cmin_for_scale(cls, Z, de, fix_center=False, verbose=True):
            """Smaller principal curvature of the quadratic fitted at one scale."""
            X, Y, N, K = cls._create_kernel(Z, de)
            H = cls._calc_inv_G_for_kernel(X, Y, N, fix_center)
            g, h, i, j, k, l = cls._convolve(X, Y, Z, K, fix_center)
            Cmin = cls._Cmin(H, g, h, i, j, k, l, fix_center)
            if verbose:
                sys.stdout.write('scale ' + str(de) + '\n')
                sys.stdout.flush()
            return Cmin, de

    @classmethod
    def _elevation_fit_for_location(cls, x, y, elevation, de, fix_center = False):

        Z = Elevation()
        Z._copy_info_from_grid(elevation, set_zeros=False)

        needs_reshaping_x = ((Z._georef_info.nx % 2) != 1)
        needs_reshaping_y = ((Z._georef_info.ny % 2) != 1)

        if needs_reshaping_x:
            nx = Z._georef_info.nx + 1
            xllcenter = Z._georef_info.xllcenter - Z._georef_info.dx
            Z._georef_info.nx = nx
            Z._georef_info.xllcenter = xllcenter
            Z._griddata = np.concatenate(((np.reshape(Z._griddata[:, 0], (Z._georef_info.ny, 1)), Z._griddata)), axis=1)
        if needs_reshaping_y:
            ny = Z._georef_info.ny + 1
            yllcenter = Z._georef_info.yllcenter - Z._georef_info.dx
            Z._georef_info.ny = ny
            Z._georef_info.yllcenter = yllcenter
            Z._griddata = np.concatenate(((np.reshape(Z._griddata[0, :], (1, Z._georef_info.nx)), Z._griddata)), axis=0)

        # Determine location in index space:
        ((i,j),) = elevation._xy_to_rowscols(((x,y),))
        ((xa, ya),) = elevation._rowscols_to_xy(((i,j),))
        a,b,c,d,e,f = cls.Utilities._calc_coefficients_for_scale(Z, de)
        elevation_center = Z[i,j]
        a = a[i,j]
        b = b[i,j]
        c = c[i,j]
        d = d[i,j]
        e = e[i,j]
        f = f[i,j] + elevation_center
        print('Window size = ' + str(de) + '\n' + 'a = ' + str(a) + ', b = ' + str(b) + ', c = ' + str(c)  + ', d = ' + str(d) + ', e = ' + str(e) + ', f = ' + str(f))

        (nx, ny, dx) = (Z._georef_info.nx, Z._georef_info.ny, Z._georef_info.dx)
        (xllcenter, yllcenter) = (Z._georef_info.xllcenter, Z._georef_info.yllcenter)

        # arange over a half-open interval of exactly n*dx yields n-1 samples
        # as often as n; build the axes from the cell count instead.
        [X, Y] = np.meshgrid(xllcenter + dx / 2 + np.arange(nx) * dx,
                             yllcenter + dx / 2 + np.arange(ny) * dx)
        X -= xa
        Y -= ya

        Z = Elevation()
        Z._copy_info_from_grid(elevation, set_zeros=True)
        Z._griddata = a*np.power(X,2) + b*np.power(Y,2) + c*X*Y + d*X + e*Y + f

        start_x_index = 1 if needs_reshaping_x else 0
        start_y_index = 1 if needs_reshaping_y else 0
        Z._griddata = Z._griddata[start_y_index:, start_x_index:]

        return Z


    def _create_from_inputs(self, *args, **kwargs):
        
        # Condition inputs to ensure that grids produce square convolution matrices:
        
        (area_cutoff, max_width, min_width, normalize, fix_center, use_dask) = (
            kwargs['area_cutoff'], kwargs['max_width'], kwargs['min_width'],
            kwargs.get('normalize', False), kwargs.get('fix_center', False),
            kwargs.get('use_dask', False))
        verbose = kwargs.get('verbose', True)
        
        if use_dask:
            from dask import compute, delayed
        
        Z = Elevation()
        A = Area()
        Z._copy_info_from_grid(kwargs['elevation'], set_zeros = False)
        A._copy_info_from_grid(kwargs['area'], set_zeros = False)
        
        needs_reshaping_x = ((Z._georef_info.nx % 2) != 1) 
        needs_reshaping_y = ((Z._georef_info.ny % 2) != 1)
            
        if needs_reshaping_x:
            nx = Z._georef_info.nx + 1
            xllcenter = Z._georef_info.xllcenter - Z._georef_info.dx
            Z._georef_info.nx = nx
            Z._georef_info.xllcenter = xllcenter
            Z._griddata = np.concatenate(((np.reshape(Z._griddata[:,0],(Z._georef_info.ny, 1)), Z._griddata)), axis=1)
            A._georef_info.nx = nx
            A._georef_info.xllcenter = xllcenter
            A._griddata = np.concatenate((np.zeros((Z._georef_info.ny, 1)), A._griddata), axis=1)
        if needs_reshaping_y:
            ny = Z._georef_info.ny + 1
            yllcenter = Z._georef_info.yllcenter - Z._georef_info.dx
            Z._georef_info.ny = ny
            Z._georef_info.yllcenter = yllcenter
            Z._griddata = np.concatenate(((np.reshape(Z._griddata[0,:], (1, Z._georef_info.nx)), Z._griddata)), axis=0)
            A._georef_info.ny = ny
            A._georef_info.yllcenter = yllcenter
            A._griddata = np.concatenate((np.zeros((1, Z._georef_info.nx)), A._griddata), axis=0)
        self._copy_info_from_grid(kwargs['elevation'], set_zeros = True)
        de = Z._georef_info.dx
        scales = np.arange(min_width, max_width, de)
        g_minC = np.zeros_like(Z._griddata)
        g_w = np.zeros_like(Z._griddata)
        ind = 1
        
        if use_dask:
            from functools import partial
            wrapper = partial(self.Utilities._Cmin_for_scale, Z,
                              fix_center=fix_center, verbose=verbose)
            tasks = [delayed(wrapper)(s) for s in scales]
            results = compute(*tasks)
            for result in results:
                minC, scale = result
                if normalize:
                    minC *= scale
                i = np.where(minC < g_minC)
                g_w[i] = np.ones(i[0].shape)*scale
                g_minC[i] = minC[i]
                ind += 1
        else:
            for scale in scales:   
                minC, _ = self.Utilities._Cmin_for_scale(Z, scale, fix_center,
                                                         verbose=verbose)
                if normalize:
                    minC *= scale
                i = np.where(minC < g_minC)
                g_w[i] = np.ones(i[0].shape)*scale
                g_minC[i] = minC[i]
                ind += 1
                
        i = np.where(A._griddata < area_cutoff)
        g_minC[i] = np.nan
        g_w[i] = np.nan
        start_x_index = 1 if needs_reshaping_x else 0
        start_y_index = 1 if needs_reshaping_y else 0
        self._griddata = g_w[start_y_index:,start_x_index:]
        self._minC = g_minC[start_y_index:, start_x_index:]
    
    def remove_padding(self, xpadding = 10, ypadding = 10):
        """Strip tile overlap from both the widths and the curvatures."""
        import copy
        unpadded = copy.deepcopy(self)
        rows = slice(ypadding, -ypadding if ypadding else None)
        cols = slice(xpadding, -xpadding if xpadding else None)
        unpadded._griddata = unpadded._griddata[rows, cols]
        unpadded._minC = unpadded._minC[rows, cols]
        unpadded._georef_info.xllcenter += xpadding*unpadded._georef_info.dx
        unpadded._georef_info.yllcenter += ypadding*unpadded._georef_info.dx
        unpadded._georef_info.nx -= 2*xpadding
        unpadded._georef_info.ny -= 2*ypadding
        unpadded._georef_info.geoTransform = (unpadded._georef_info.xllcenter - 0.5*unpadded._georef_info.dx, unpadded._georef_info.dx, 0, unpadded._georef_info.yllcenter + (float(unpadded._georef_info.ny-0.5))*self._georef_info.dx, 0, -self._georef_info.dx)
        
        return unpadded
    
    def save(self, filename):
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', [self._griddata, self._minC], self.dtype, filename, ['COMPRESS=LZW'], multiple_bands=True)
    
    @classmethod
    def load(cls, filename):
        """Read back the multi-band grid written by :meth:`save`."""
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.nan
            return grid
        
        return_object = cls()
        gdal_dataset = gdal.Open(filename)
        
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        
        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[0]+return_object._georef_info.dx/2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
        return_object._griddata = get_band(gdal_dataset, 1)
        return_object._minC = get_band(gdal_dataset, 2)        
            
        gdal_file = None
        return return_object
    
    @classmethod
    def mosaic(cls, tiles):
        """Reassemble tiles, carrying the curvature band along with the widths."""
        xmin = np.nan
        xmax = np.nan
        ymin = np.nan
        ymax = np.nan
        
        for tile in tiles:
            
            xmin = tile._georef_info.xllcenter if (np.isnan(xmin) or (xmin > tile._georef_info.xllcenter)) else xmin
            ymin = tile._georef_info.yllcenter if (np.isnan(ymin) or (ymin > tile._georef_info.yllcenter)) else ymin
            xmax = tile._georef_info.xllcenter+(tile._georef_info.nx-1)*tile._georef_info.dx if (np.isnan(xmax) or (xmax < tile._georef_info.xllcenter+(tile._georef_info.nx-1)*tile._georef_info.dx)) else xmax
            ymax = tile._georef_info.yllcenter+(tile._georef_info.ny-1)*tile._georef_info.dx if (np.isnan(ymax) or (ymax < tile._georef_info.yllcenter+(tile._georef_info.ny-1)*tile._georef_info.dx)) else ymax
        
        ny = int(round((ymax - ymin) / tiles[0]._georef_info.dx)) + 1
        nx = int(round((xmax - xmin) / tiles[0]._georef_info.dx)) + 1
        
        return_object = cls()
        return_object._griddata = np.zeros((ny, nx))
        return_object._griddata[:] = np.nan
        return_object._minC = np.zeros((ny, nx))
        return_object._minC[:] = np.nan
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        return_object._georef_info.dx = tiles[0]._georef_info.dx
        return_object._georef_info.xllcenter = xmin
        return_object._georef_info.yllcenter = ymin
        return_object._georef_info.geoTransform = (return_object._georef_info.xllcenter - 0.5*return_object._georef_info.dx, return_object._georef_info.dx, 0, return_object._georef_info.yllcenter + (float(return_object._georef_info.ny-0.5))*return_object._georef_info.dx, 0, -return_object._georef_info.dx)
        
        for tile in tiles:
            j_min = int(round((tile._georef_info.xllcenter - return_object._georef_info.xllcenter) / return_object._georef_info.dx))
            j_max = j_min + tile._georef_info.nx
            i_min = int(round((return_object._georef_info.yllcenter + (return_object._georef_info.ny-1)*return_object._georef_info.dx - tile._georef_info.yllcenter - (tile._georef_info.ny-1)*tile._georef_info.dx) / return_object._georef_info.dx))
            i_max = i_min + tile._georef_info.ny
                        
            return_object._griddata[i_min:i_max, j_min:j_max] = tile._griddata
            return_object._minC[i_min:i_max, j_min:j_max] = tile._minC
            
        return return_object                   
     
class DiscreteFlowAccumulation(BaseSpatialGrid):
    """Area drained by a single steepest-descent walk from each outlet.

    Unlike :class:`Area`, which accumulates over the whole D8 network, this
    follows one path per outlet, stepping to the lowest unvisited neighbour
    each time, and records the running area along it.  It is used to trace
    individual flow paths across surfaces (debris flows, lava) where the
    network abstraction does not apply.

    ::

        dfa = DiscreteFlowAccumulation(elevation=dem, outlets=[(x, y)])

    Keyword arguments
    -----------------
    mask : BaseSpatialGrid
        Confine the walk to cells where the mask is 1.
    terminations_only : bool
        Record the total area only at the end of each path.
    display_output : bool
        Report progress outlet by outlet.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'outlets'), '_create_from_elevation_outlets'),
                                   (('elevation', ), '_create_from_elevation'),)

    def _create_from_elevation(self, *args, **kwargs):
        """Start a walk from every valid cell."""
        kwargs = dict(kwargs)
        elevation = Elevation()
        elevation._copy_info_from_grid(kwargs['elevation'])
        if kwargs.get('mask') is not None:
            outside = np.where(kwargs['mask']._griddata != 1)
            elevation._griddata[outside] = np.nan
        kwargs['elevation'] = elevation
        rows, cols = np.where(~np.isnan(elevation._griddata))
        kwargs['outlets'] = elevation._rowscols_to_xy(
            list(zip(rows.tolist(), cols.tolist())))
        self._create_from_elevation_outlets(*args, **kwargs)

    def _create_from_elevation_outlets(self, *args, **kwargs):

        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation, True)
        adjust = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]

        per_pixel_area = self._area_per_pixel()
        terminations_only = kwargs.get('terminations_only') is True
        outlets = list(kwargs['outlets'])

        def neighbour_elevations(cell):
            """Neighbour elevations, NaN off-grid, as a plain float array."""
            values = np.empty(8, dtype=float64)
            for k, (di, dj) in enumerate(adjust):
                value = elevation[cell[0] + di, cell[1] + dj]
                values[k] = np.nan if value is None else value
            return values

        for counter, outlet in enumerate(outlets, start=1):
            if kwargs.get('display_output') is True:
                print('Evaluating outlet {0} / {1}'.format(counter, len(outlets)))

            (ij, ) = elevation._xy_to_rowscols((outlet,))
            if ij[0] is None:
                continue

            seen = {tuple(ij)}
            e_a = neighbour_elevations(ij)
            area = 0.0

            # `None not in numpy_array` never matched, so the original loop
            # only ever stopped on NaN.
            while not np.any(np.isnan(e_a)):
                area += per_pixel_area[ij[0], ij[1]]
                if not terminations_only:
                    self._griddata[ij[0], ij[1]] = area
                moved = False
                for sl in np.argsort(e_a):
                    this_ij = (ij[0] + adjust[sl][0], ij[1] + adjust[sl][1])
                    if this_ij in seen:
                        continue
                    seen.add(this_ij)
                    ij = this_ij
                    moved = True
                    break
                if not moved:
                    break
                e_a = neighbour_elevations(ij)

            if terminations_only:
                self._griddata[ij[0], ij[1]] = area

    def _area_per_pixel(self, *args, **kwargs):
        return self._georef_info.dx**2 * np.ones((self._georef_info.ny, self._georef_info.nx))

    def _mean_pixel_dimension(self, *args, **kwargs):
        # Reading the shape from kwargs['elevation'] raised KeyError whenever
        # this was called without one; the grid knows its own shape.
        return self._georef_info.dx * np.ones(
            (self._georef_info.ny, self._georef_info.nx), dtype=float64)

class GeographicDiscreteFlowAccumulation(GeographicGridMixin, DiscreteFlowAccumulation):
    """Discrete flow accumulation on a latitude/longitude grid."""

    pass
    
class FlowLength(BaseSpatialGrid):
    """Longest upstream flow distance reaching each cell.

    ::

        length = FlowLength(flow_direction=d8)

    Alongside the distances, the grid records which upstream neighbour
    supplies each cell's longest path.  That "main stem" network is what
    :class:`Relief`, :class:`ScaledRelief` and :class:`Ksi` follow, and what
    :meth:`locations_along_flow_path_from_outlet` walks.

    .. note::
       The main-stem codes use the same ArcGIS convention as
       :class:`FlowDirectionD8`, pointing *from* a cell *at* its upstream
       donor.  Before version 1.0 they used a row-flipped convention that
       was internally consistent but disagreed with every other
       flow-direction grid in the library, so ``*_directions`` side-files
       written by older versions must be regenerated.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('flow_direction',), '_create_from_flow_direction_and_sorted_indexes'))

    def _create_from_flow_direction_and_sorted_indexes(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        self.__calculate_flow_length(*args, **kwargs)

    def __calculate_flow_length(self, *args, **kwargs):

        flow_dir = kwargs['flow_direction']
        kernels.warn_if_slow(flow_dir._griddata.size)

        step = flow_dir.step_lengths(self._mean_pixel_dimension(*args, **kwargs))
        length, main_stem = kernels.flow_length(
            flow_dir.receivers,
            flow_dir.topological_order,
            np.ascontiguousarray(step, dtype=float64))

        self._griddata = np.asarray(length, dtype=self.dtype)
        self.__flow_directions = np.asarray(main_stem, dtype=np.uint8)

    @property
    def main_stem_directions(self):
        """Codes pointing from each cell at its longest-path donor."""
        return self.__flow_directions

    @staticmethod
    def __flow_direction_for_length(from_index, to_index):
        """Code stored at ``to_index`` that points back at ``from_index``."""
        d_i = from_index[0] - to_index[0]
        d_j = from_index[1] - to_index[1]
        for d in range(8):
            if int(kernels.D8_DI[d]) == d_i and int(kernels.D8_DJ[d]) == d_j:
                return int(kernels.D8_CODES[d])
        return 0

    def points_with_length(self, length, fd):
        """Map coordinates where the flow path first exceeds ``length``.

        Used to pick a consistent set of basin outlets: every returned point
        has ``length`` of channel upstream of it.
        """
        tolerance = np.nanmin(self._mean_pixel_dimension())
        min_length = length - tolerance * 2.0

        indexes_of_locations = list()
        ind = np.where((self._griddata >= min_length) & (self._griddata <= length))
        for (this_row, this_col) in zip(ind[0], ind[1]):
            (next_row, next_col, is_good) = fd.get_flow_to_cell(this_row, this_col)
            if not is_good:
                continue
            if self[this_row, this_col] <= length <= self[next_row, next_col]:
                indexes_of_locations.append((this_row, this_col))

        return self._rowscols_to_xy(indexes_of_locations)

    def indexes_along_flow_path_from_outlet(self, outlet):
        """``(row, col)`` along the longest flow path upstream of ``outlet``."""
        ij = self._xy_to_rowscols((outlet,))[0]
        if ij[0] is None:
            return []
        return self.__get_upstream_indexes((ij,))

    def locations_along_flow_path_from_outlet(self, outlet):
        """As above, in map coordinates."""
        return self._rowscols_to_xy(self.indexes_along_flow_path_from_outlet(outlet))

    def locations_along_flow_path_from_outlets(self, outlets):
        """Concatenated main-stem paths for several outlets."""
        points = tuple()
        for outlet in outlets:
            points += tuple(self.locations_along_flow_path_from_outlet(outlet))
        return points

    def __step_upstream(self, i, j):
        """The single main-stem donor of ``(i, j)``, or ``None``."""
        code = self.__flow_directions[i, j]
        if code == 0:
            return None
        for d in range(8):
            if code == kernels.D8_CODES[d]:
                up_i = i + int(kernels.D8_DI[d])
                up_j = j + int(kernels.D8_DJ[d])
                if 0 <= up_i < self._georef_info.ny and 0 <= up_j < self._georef_info.nx:
                    return (up_i, up_j)
                return None
        return None

    def __get_upstream_indexes(self, index):
        """Walk the main stem upstream.  Iterative, so long rivers are fine."""
        (i, j) = index[0]
        indexes = [(i, j)]
        seen = {(i, j)}
        while True:
            nxt = self.__step_upstream(i, j)
            if nxt is None or nxt in seen:
                break
            indexes.append(nxt)
            seen.add(nxt)
            (i, j) = nxt
        return indexes

    def __map_flow_from_cell(self, index, **kwargs):
        """Nested dictionary along the main stem upstream of ``index``.

        Each node carries ``index``, ``distance`` (the flow length at that
        cell), ``distance_scale`` (the step length from its parent, in cell
        widths) and one entry per keyword grid.
        """
        (i, j) = index[0]

        def make_node(cell, scale):
            node = {'index': cell,
                    'distance': self._griddata[cell[0], cell[1]],
                    'distance_scale': scale}
            for name, grid in kwargs.items():
                node[name] = grid[cell[0], cell[1]]
            return node

        root = make_node((i, j), 1.0)
        node = root
        seen = {(i, j)}
        while True:
            nxt = self.__step_upstream(i, j)
            if nxt is None or nxt in seen:
                break
            diagonal = (nxt[0] != i) and (nxt[1] != j)
            child = make_node(nxt, SQRT2 if diagonal else 1.0)
            node['next'] = [child]
            node = child
            seen.add(nxt)
            (i, j) = nxt

        return root

    def is_along_flow_length(self, from_index, to_index):
        """True when ``from_index`` supplies ``to_index``'s longest path."""
        (i_to, j_to) = to_index
        flow_code = self.__flow_directions[i_to, j_to]

        return flow_code == self.__flow_direction_for_length(from_index, to_index)

    def clip_to_extent(self, extent):
        """Clip both the lengths and the main-stem codes to ``extent``."""
        return_grid = super(FlowLength, self).clip_to_extent(extent)

        # Clip the direction grid the same way as the data, through a
        # temporary grid -- the old version assigned the clipped *lengths*
        # into the direction array.
        helper = BaseSpatialGrid()
        helper._georef_info = self._georef_info
        helper._griddata = self.__flow_directions
        return_grid.__flow_directions = helper.clip_to_extent(extent)._griddata.astype(np.uint8)

        return return_grid

    def map_values_to_recursive_list(self, outlet, **kwargs):
        """Nested main-stem dictionary upstream of ``outlet``."""
        v = (outlet, )
        (ij_outlet, ) = self._xy_to_rowscols(v)
        if ij_outlet[0] is None:
            raise Error.InputError('map_values_to_recursive_list',
                                   'outlet {0} lies outside the grid'.format(outlet))
        return self.__map_flow_from_cell((ij_outlet,), **kwargs)

    def save(self, filename):
        """Write the lengths, plus the main-stem codes to ``filename_directions``."""
        super(FlowLength, self).save(filename)
        flow_dir_name = filename + "_directions"
        dataset = self._create_gdal_representation_from_array(
            self._georef_info, 'GTiff', self.__flow_directions, np.uint8,
            flow_dir_name, ['COMPRESS=LZW'])
        dataset.FlushCache()
        del dataset

    @classmethod
    def load(cls, filename):
        """Read back a grid written by :meth:`save`."""
        return_object_bsp = BaseSpatialGrid.load(filename)
        return_object = cls()
        return_object._georef_info = return_object_bsp._georef_info
        return_object._griddata = return_object_bsp._griddata
        flow_dir_filename = filename + "_directions"
        gdal_dataset = gdal.Open(flow_dir_filename)
        try:
            band = gdal_dataset.GetRasterBand(1)
            return_object.__flow_directions = band.ReadAsArray().astype(np.uint8)
        finally:
            gdal_dataset = None
        return return_object

class GeographicFlowLength(GeographicGridMixin, FlowLength):
    """Flow length on a latitude/longitude grid."""

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('flow_direction',), '_create_from_flow_direction_and_sorted_indexes'))
    
class MaxFlowLengthTrackingMixin(object):
    """Propagate a value downstream along longest-flow-length paths only.

    Subclasses set :attr:`_propagate_mode` to ``'carry'`` (the receiver takes
    the donor's value) or ``'sum'`` (the receiver adds it) and may define
    :meth:`_propagation_gate` to suppress the transfer for some cells.
    """

    #: ``'carry'`` or ``'sum'``.
    _propagate_mode = 'carry'

    def _propagation_gate(self, *args, **kwargs):
        """0/1 grid of cells allowed to pass their value on; ``None`` = all."""
        return None

    def _calculate_by_tracking_down_max_flow_length(self, *args, **kwargs):
        flow_dir = kwargs['flow_direction']
        flow_length = kwargs['flow_length']
        kernels.warn_if_slow(flow_dir._griddata.size)

        self._griddata = np.asarray(kernels.propagate_along_main_stem(
            flow_dir.receivers,
            flow_dir.topological_order,
            np.ascontiguousarray(flow_length.main_stem_directions, dtype=np.uint8),
            self._propagation_gate(*args, **kwargs),
            np.ascontiguousarray(self._griddata, dtype=float64),
            self._propagate_mode))


class Relief(BaseSpatialGrid, MaxFlowLengthTrackingMixin):
    """Height of the head of the longest flow path above each cell.

    ::

        relief = Relief(flow_direction=d8, elevation=dem, flow_length=length)

    Pass ``area`` and ``Ao`` together to stop the relief signal being carried
    down from cells with less than ``Ao`` of drainage area, which keeps
    hillslope noise out of the channel network.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('flow_direction', 'elevation', 'flow_length'), '_create_from_flow_direction_sorted_indexes_and_elevation'))

    _propagate_mode = 'carry'

    def _propagation_gate(self, *args, **kwargs):
        area = kwargs.get('area')
        Ao = kwargs.get('Ao')
        if area is None or Ao is None:
            return None
        return (np.asarray(area._griddata) > Ao).astype(np.uint8)

    def _create_from_flow_direction_sorted_indexes_and_elevation(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['elevation'])
        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)
        self._griddata = self._griddata - kwargs['elevation']._griddata

    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        """Per-cell transfer rule.

        Retained for subclasses that override it; the bulk path uses the
        compiled kernel instead.
        """
        (i, j) = pos
        (i_next, j_next) = next_pos

        if kwargs.get('area') is not None and kwargs.get('Ao') is not None:
            if kwargs['area'][i, j] <= kwargs['Ao']:
                return

        self._griddata[i_next, j_next] = self._griddata[i, j]


class ScaledRelief(Relief):
    """:class:`Relief` scaled by ``Ao**theta``.

    That scaling puts relief in the same units as :class:`Ksi`, so the ratio
    of the two is a channel steepness.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('flow_direction', 'flooded_dem', 'elevation', 'flow_length', 'Ao', 'theta'), '_create_scaled_from_flow_direction_flooded_dem_and_elevation'),
                           (('flow_direction', 'elevation', 'flow_length', 'Ao', 'theta'), '_create_scaled_from_flow_direction_sorted_indexes_and_elevation'))

    def _create_scaled_from_flow_direction_flooded_dem_and_elevation(self, *args, **kwargs):
        # The original called a `_create_from_flow_direction_flooded_dem_and_elevation`
        # that was never defined, so supplying `flooded_dem` raised
        # AttributeError. The flooded DEM plays no part in the calculation
        # beyond being an elevation surface, so route it to the same code.
        self._create_from_flow_direction_sorted_indexes_and_elevation(*args, **kwargs)
        self._griddata = self._griddata * (np.power(kwargs['Ao'], kwargs['theta']))

    def _create_scaled_from_flow_direction_sorted_indexes_and_elevation(self, *args, **kwargs):
        self._create_from_flow_direction_sorted_indexes_and_elevation(*args, **kwargs)
        self._griddata = self._griddata * (np.power(kwargs['Ao'], kwargs['theta']))


class Ksi(BaseSpatialGrid, MaxFlowLengthTrackingMixin):
    r"""The chi coordinate integrated along longest-flow-length paths.

    .. math::
        \chi = \int \left(\frac{A_0}{A}\right)^{\theta} \mathrm{d}x

    ::

        ksi = Ksi(area=area, flow_direction=d8, flow_length=length,
                  theta=0.45, Ao=1e6)

    Cells with less than ``Ao`` of drainage area are left out of the
    integral, so chi is defined only on the channel network.

    .. note::
       Before version 1.0 the integrand was ``(Ao / (A - Ao))**theta``, which
       diverges as ``A`` approaches ``Ao`` and does not match the definition
       of chi in the literature. Values from older runs are not comparable
       with these.

    See also :class:`Chi`, which integrates upstream from chosen outlets over
    the whole network rather than along main stems.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('area','flow_direction','theta', 'Ao', 'flow_length'), '_create_from_inputs'))

    _propagate_mode = 'sum'

    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)

        area = np.asarray(kwargs['area']._griddata, dtype=float64)
        Ao = kwargs['Ao']
        theta = kwargs['theta']

        step = (self._mean_pixel_dimension(*args, **kwargs)
                * kwargs['flow_direction'].pixel_scale(float64))

        self._griddata = np.zeros(area.shape, dtype=float64)
        on_network = area > Ao
        self._griddata[on_network] = (
            (Ao / area[on_network]) ** theta) * step[on_network]

        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)

    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        """Per-cell transfer rule; the bulk path uses the compiled kernel."""
        (i, j) = pos
        (i_next, j_next) = next_pos
        self._griddata[i_next, j_next] += self._griddata[i, j]

class GeographicKsi(GeographicGridMixin, Ksi):
    """Chi along longest-flow-length paths on a latitude/longitude grid."""

    pass

class RestoredElevation(BaseSpatialGrid):
    r"""Topography implied by a uniform channel steepness.

    Starting from the given outlets, elevations are rebuilt upstream from

    .. math::  \frac{\mathrm{d}z}{\mathrm{d}x} = k_s A^{-\theta}

    and the divides are then allowed to migrate towards whichever side ends
    up lower.  Repeating that (``iterations`` times) converges on the
    landscape a steady, spatially uniform ``ks`` would produce, which can be
    compared against the real DEM to find where it is out of steady state.

    ::

        restored = RestoredElevation(flow_direction=d8, elevation=dem,
                                     area=area, theta=0.45, ks=50,
                                     outlets=[(x, y)], iterations=5)

    Keyword arguments
    -----------------
    randomize : bool
        Randomise elevations inside the basin first, then re-route, so the
        result does not inherit the real network's planform.
    fix_external_outlets : bool
        Hold the outer boundary of the basin fixed while divides migrate.
    verbose : bool
        Report convergence each iteration (default ``True``).
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('flow_direction','elevation','area','theta','ks','outlets', 'iterations'), '_fill_dem'))

    def _fill_dem(self, *args, **kwargs):
        import copy
        self._copy_info_from_grid(kwargs['elevation'], False)
        outlets = kwargs['outlets']
        area = copy.deepcopy(kwargs['area'])
        flow_direction = copy.deepcopy(kwargs['flow_direction'])
        randomize = kwargs.get('randomize')
        iterations = kwargs['iterations']
        ks = kwargs['ks']
        theta = kwargs['theta']
        verbose = kwargs.get('verbose', True)
        v = [rc for rc in self._xy_to_rowscols(outlets) if rc[0] is not None]
        pixel_dimension = self._mean_pixel_dimension(*args, **kwargs)

        mask = Mask(flow_direction=flow_direction, outlets=outlets)
        external_divides = copy.deepcopy(mask)

        if kwargs.get('fix_external_outlets') is True:
            from scipy.ndimage import binary_erosion as erosion
            external_divides._griddata = (
                external_divides._griddata
                - erosion(external_divides._griddata, np.ones((3, 3)),
                          border_value=0).astype(int))
        else:
            external_divides._griddata = np.zeros_like(external_divides._griddata, int)

        if randomize:
            filled = FilledElevation(elevation=self, mask=mask, randomize=randomize,
                                     outlets=outlets)
            flow_direction.update_flow_codes_in_mask(filled, mask)
            inside = np.where(mask._griddata == 1)
            self._griddata[inside] = filled._griddata[inside]
            area = self.__recalculate_area(area, flow_direction, pixel_dimension, outlets)

        self.convergence = []
        for iteration in range(iterations):
            last_grid = self._griddata.copy()
            divides = self.__fill_outlets(area, flow_direction, pixel_dimension, v, ks, theta)
            flow_direction = self.__migrate_divides(flow_direction, divides, external_divides,
                                                    verbose=verbose)
            area = self.__recalculate_area(area, flow_direction, pixel_dimension, outlets)
            change_in_elevation = float(np.nanmean((self._griddata - last_grid) ** 2))
            self.convergence.append(change_in_elevation)
            if verbose:
                print('Iteration {0}: mean squared change {1:g}'.format(
                    iteration, change_in_elevation))

        self.flow_direction = flow_direction
        self.area = area

    def __is_unit_area(self, idx, flow_direction):
        """True when no neighbour drains into ``idx`` (it is a divide cell)."""
        return not flow_direction.get_upstream_cell_indexes(idx[0], idx[1])

    def __migrate_divides(self, *args, **kwargs):
        """Redirect divide-adjacent cells towards whichever side is lower."""
        flow_direction = args[0]
        divides = args[1]
        external_divides = args[2]
        verbose = kwargs.get('verbose', True)
        migrated = 0

        ny, nx = self._georef_info.ny, self._georef_info.nx

        for (i, j) in divides:
            if external_divides[i, j] == 1:
                continue
            here = self[i, j]
            if here is None:
                continue
            for d in range(8):
                ni = i + int(kernels.D8_DI[d])
                nj = j + int(kernels.D8_DJ[d])
                # Explicit bounds checks instead of eight try/except blocks:
                # the exception handlers also swallowed real errors, and
                # negative indices wrapped instead of being rejected.
                if not (0 <= ni < ny and 0 <= nj < nx):
                    continue
                if external_divides[ni, nj] == 1:
                    continue
                neighbour = self[ni, nj]
                if neighbour is None or not neighbour > here:
                    continue
                flow_direction[ni, nj] = int(kernels.D8_CODES[(d + 4) % 8])
                migrated += 1

        flow_direction._invalidate_network()
        if verbose:
            print("Migrated " + str(migrated) + " divides.")
        return flow_direction

    def __fill_outlets(self, area, flow_direction, pixel_dimension, outlet_indexes, ks, theta):
        divides = list()
        for ind in outlet_indexes:
            divides += self.__fill_upstream_points(ind, ks, theta, area, flow_direction,
                                                   pixel_dimension)
        return divides

    def __fill_upstream_points(self, ind, ks, theta, area, flow_direction, pixel_dimension):
        """Integrate dz/dx = ks * A**-theta upstream from ``ind``.

        Returns the divide cells reached.  Iterative: the recursive original
        needed a million-frame recursion limit and still overflowed on real
        basins.
        """
        ny, nx = self._georef_info.ny, self._georef_info.nx
        divides = list()
        stack = [tuple(ind)]
        seen = {tuple(ind)}

        while stack:
            (i, j) = stack.pop()
            elevation_for_cell = self._griddata[i, j]
            a = area[i, j]
            is_divide = True

            for d in range(8):
                ui = i + int(kernels.D8_DI[d])
                uj = j + int(kernels.D8_DJ[d])
                if not (0 <= ui < ny and 0 <= uj < nx):
                    continue
                if flow_direction[ui, uj] != kernels.D8_CODES[(d + 4) % 8]:
                    continue
                is_divide = False
                if (ui, uj) in seen:
                    continue
                seen.add((ui, uj))
                step = pixel_dimension[i, j] * kernels.D8_DIST[d]
                self._griddata[ui, uj] = elevation_for_cell + step * ks * a ** (-theta)
                stack.append((ui, uj))

            if is_divide:
                divides.append((i, j))

        return divides

    def __recalculate_area(self, area, flow_direction, pixel_dimension, outlets):
        """Re-accumulate drainage area after the network has been rewired."""
        weights = np.asarray(pixel_dimension, dtype=float64) ** 2
        accumulated = np.asarray(kernels.accumulate(
            flow_direction.receivers,
            flow_direction.topological_order,
            np.ascontiguousarray(weights)))

        inside = self.__basin_mask(flow_direction, outlets)
        area._griddata = np.where(inside, accumulated, area._griddata)
        return area

    @staticmethod
    def __basin_mask(flow_direction, outlets):
        return flow_direction.basin_mask(outlets)

class GeographicRestoredElevation(GeographicGridMixin, RestoredElevation):
    """Restored topography on a latitude/longitude grid."""

    pass

class ChiScaledRelief(BaseSpatialGrid):
    """Height above the basin outlet, scaled by ``Ao**theta``.

    Plotted against :class:`Chi`, the slope of this quantity is the channel
    steepness index.

    ::

        relief = ChiScaledRelief(elevation=dem, flow_direction=d8,
                                 theta=0.45, Ao=1e6, outlets=[(x, y)])
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('elevation', 'flow_direction', 'theta', 'Ao', 'outlets'), '_create_from_inputs'),
                           (('elevation', 'flow_direction', 'flow_length', 'theta', 'Ao', 'basin_length'), '_create_from_basin_length'))

    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['elevation'], True)
        elevation = kwargs['elevation']
        flow_direction = kwargs['flow_direction']
        scale = np.power(kwargs['Ao'], kwargs['theta'])

        outlet_indexes = [rc for rc in self._xy_to_rowscols(kwargs['outlets'])
                          if rc[0] is not None]
        for outlet_number, outlet_index in enumerate(outlet_indexes, start=1):
            # One vectorised basin extraction per outlet, rather than an
            # element-by-element write through __setitem__.
            basin = flow_direction.basin_mask(
                self._rowscols_to_xy((outlet_index,)))
            elevation_of_outlet = elevation[outlet_index[0], outlet_index[1]]
            self._griddata[basin] = (
                np.asarray(elevation._griddata)[basin] - elevation_of_outlet) * scale
            if kwargs.get('output_flag', False):
                print('Outlet {0}/{1} completed.'.format(outlet_number, len(outlet_indexes)))

    def _create_from_basin_length(self, *args, **kwargs):
        """Use every point with ``basin_length`` of channel upstream as an outlet."""
        kwargs = dict(kwargs)
        kwargs['outlets'] = kwargs['flow_length'].points_with_length(
            kwargs['basin_length'], kwargs['flow_direction'])
        kwargs['output_flag'] = True
        return self._create_from_inputs(*args, **kwargs)

class Chi(BaseSpatialGrid):
    r"""The chi coordinate, integrated upstream from chosen outlets.

    .. math::
        \chi(x) = \int_{x_{\mathrm{outlet}}}^{x}
                  \left(\frac{A_0}{A(x')}\right)^{\theta}\,\mathrm{d}x'

    ::

        chi = Chi(area=area, flow_direction=d8, theta=0.45, Ao=1e6,
                  outlets=[(x, y)])

    Chi is 0 at each outlet and increases upstream.  Cells outside the
    outlets' drainage areas are left as 0.

    Keyword arguments
    -----------------
    trapezoid : bool
        Average the integrand across each step instead of taking it at the
        upstream cell.  This is the rule TopoToolbox's ``chitransform`` uses;
        the default (``False``) reproduces the original TopoAnalysis result.
    maximum_length : float
        Stop integrating beyond this flow distance from the outlet.
    mask : BaseSpatialGrid
        Restrict the integration to cells where the mask is non-zero.

    .. note::
       Before version 1.0, chi at the outlet itself was one cell-step of the
       integrand rather than 0.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('area','flow_direction', 'theta', 'Ao', 'outlets'), '_create_from_inputs'),
                           (('area', 'flow_direction', 'flow_length', 'theta', 'Ao', 'basin_length'), '_create_from_basin_length'))

    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        self.__calculate_chi(*args, **kwargs)

    def _create_from_basin_length(self, *args, **kwargs):
        """Use every point with ``basin_length`` of channel upstream as an outlet."""
        kwargs = dict(kwargs)
        kwargs['outlets'] = kwargs['flow_length'].points_with_length(
            kwargs['basin_length'], kwargs['flow_direction'])
        kwargs['output_flag'] = True
        return self._create_from_inputs(*args, **kwargs)

    def __calculate_chi(self, *args, **kwargs):
        area = kwargs['area']
        flow_direction = kwargs['flow_direction']
        kernels.warn_if_slow(flow_direction._griddata.size)

        step = flow_direction.step_lengths(self._mean_pixel_dimension(*args, **kwargs))

        outlets = []
        for (row, col) in self._xy_to_rowscols(kwargs['outlets']):
            if row is not None:
                outlets.append(row * self._georef_info.nx + col)

        mask = kwargs.get('mask')
        mask_array = None if mask is None else (np.asarray(mask._griddata) != 0).astype(np.uint8)

        chi, distance = kernels.chi(
            flow_direction.receivers,
            flow_direction.topological_order,
            np.ascontiguousarray(area._griddata, dtype=float64),
            np.ascontiguousarray(step, dtype=float64),
            np.array(outlets, dtype=np.int64),
            kwargs['Ao'],
            kwargs['theta'],
            bool(kwargs.get('trapezoid', False)),
            float(kwargs.get('maximum_length') or 0.0),
            mask_array)

        chi = np.asarray(chi)
        self.distance_from_outlet = np.asarray(distance)
        # Cells the integration never reached come back as NaN; the grid's
        # documented empty value is 0.
        self._griddata = np.where(np.isnan(chi), 0.0, chi)

        if kwargs.get('output_flag', False):
            print('{0} outlets completed.'.format(len(outlets)))


class GeographicChi(GeographicGridMixin, Chi):
    """Chi on a latitude/longitude grid."""

    pass


class _CrossDivideMixin(object):
    """Shared machinery for the across-divide chi contrast grids."""

    def _divide_pairs(self, *args, **kwargs):
        """``(from_rows, from_cols, to_rows, to_cols)`` for each divide pair."""
        paired_divides = self._paired_divides(*args, **kwargs)
        if not paired_divides:
            return None

        # zip() is a one-shot iterator in Python 3 and is not subscriptable,
        # so the original `zip(*pairs)[0]` raised immediately.
        from_point, to_point = zip(*paired_divides)

        from_idx = self._xy_to_rowscols(from_point)
        to_idx = self._xy_to_rowscols(to_point)

        keep = [k for k in range(len(from_idx))
                if from_idx[k][0] is not None and to_idx[k][0] is not None]
        if not keep:
            return None

        from_rows = np.array([from_idx[k][0] for k in keep], dtype=int)
        from_cols = np.array([from_idx[k][1] for k in keep], dtype=int)
        to_rows = np.array([to_idx[k][0] for k in keep], dtype=int)
        to_cols = np.array([to_idx[k][1] for k in keep], dtype=int)
        return from_rows, from_cols, to_rows, to_cols


class CrossDivideDChi(_CrossDivideMixin, BaseSpatialGrid):
    """Difference in chi across each drainage divide.

    A large value means the two sides of the divide are far from
    steady state, so the divide is expected to migrate towards the
    higher-chi side (Willett et al. 2014).
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('chi','flow_direction'), '_create_from_inputs'))

    def _paired_divides(self, *args, **kwargs):
        return kwargs['flow_direction'].paired_divides()

    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        pairs = self._divide_pairs(*args, **kwargs)
        if pairs is None:
            return
        from_rows, from_cols, to_rows, to_cols = pairs

        chi = kwargs['chi']._griddata
        a = chi[to_rows, to_cols]
        b = chi[from_rows, from_cols]
        minchi = np.minimum(a, b)
        maxchi = np.maximum(a, b)

        self._griddata[from_rows, from_cols] = (maxchi - minchi) * (minchi != 0).astype(float)


class NormalizedCrossDivideDChi(_CrossDivideMixin, BaseSpatialGrid):
    """Chi contrast across each divide, divided by the mean chi there.

    Normalising makes divides in high- and low-chi parts of a landscape
    comparable.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'),
                           (('chi','flow_direction'), '_create_from_inputs'))

    def _paired_divides(self, *args, **kwargs):
        return kwargs['flow_direction'].paired_divides(mask=kwargs['chi'])

    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        pairs = self._divide_pairs(*args, **kwargs)
        if pairs is None:
            return
        from_rows, from_cols, to_rows, to_cols = pairs

        chi = kwargs['chi']._griddata
        a = chi[to_rows, to_cols]
        b = chi[from_rows, from_cols]
        minchi = np.minimum(a, b)
        maxchi = np.maximum(a, b)

        rangechi = (maxchi - minchi) * (minchi != 0).astype(float)
        meanchi = (maxchi + minchi) / 2.0
        with np.errstate(divide='ignore', invalid='ignore'):
            self._griddata[from_rows, from_cols] = np.where(
                meanchi != 0, rangechi / meanchi, 0.0)

class Deflection(BaseSpatialGrid):
    r"""Flexural deflection of an elastic plate under a topographic load.

    Solves the thin-plate equation in the Fourier domain:

    .. math::
        w(k) = \frac{-\rho_c\,g\,h(k)}
                    {(\rho_m - \rho_c)\,g + D\,(2\pi k)^4}

    ::

        w = Deflection(elevation=dem, D=1e23, rho_m=3300, rho_c=2700, g=9.81)

    Pass ``restored_elevation`` as well to get the *change* in deflection
    between two surfaces, which is the isostatic response to the erosion
    between them.

    Because the solution is spectral, the load is implicitly periodic; pad
    the DEM if edge effects matter.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('elevation', 'D', 'rho_m', 'rho_c', 'g'), '_deflect'))

    def _deflect(self, *args, **kwargs):
        if kwargs.get('restored_elevation'):
            self._griddata = (self.__deflect_elevation(kwargs['elevation'], **kwargs)
                              - self.__deflect_elevation(kwargs['restored_elevation'], **kwargs))
        else:
            self._griddata = self.__deflect_elevation(kwargs['elevation'], **kwargs)

    def __deflect_elevation(self, *args, **kwargs):
        elevation = args[0]
        self._copy_info_from_grid(elevation, False)
        kwargs = dict(kwargs)
        kwargs['flow_direction'] = elevation
        dx = float(np.mean(self._mean_pixel_dimension(*args, **kwargs)))

        load = -np.nan_to_num(np.asarray(self._griddata, dtype=float64)) \
            * kwargs['rho_c'] * kwargs['g']
        load_hat = np.fft.fft2(load)

        # np.fft.fftfreq gives exactly one wavenumber per sample, in the
        # order fft2 expects. Building the axes by hand from half-width
        # ranges produced arrays whose length did not have to match the grid,
        # so the multiply below could raise a broadcast error.
        wn_x = np.fft.fftfreq(self._georef_info.nx, d=dx)
        wn_y = np.fft.fftfreq(self._georef_info.ny, d=dx)
        WN_x, WN_y = np.meshgrid(wn_x, wn_y)

        # (2*pi)**4, not (2*3.141)**4 -- the truncated constant was 0.2% low
        # in the flexural term.
        response = 1.0 / ((kwargs['rho_m'] - kwargs['rho_c']) * kwargs['g']
                          + (2 * np.pi) ** 4 * kwargs['D']
                          * (WN_x ** 2 + WN_y ** 2) ** 2)

        return np.real(np.fft.ifft2(load_hat * response))

class GeographicDeflection(GeographicGridMixin, Deflection):
    """Flexural deflection on a latitude/longitude grid."""

    pass


class ChannelSlope(BaseSpatialGrid):
    """Downstream gradient of each cell along its D8 flow direction.

    ::

        slope = ChannelSlope(flow_direction=d8, elevation=dem)

    Positive values are downhill.  Cells with no downstream neighbour stay 0.
    """

    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'),
                               (('flow_direction','elevation'), '_create_from_flow_direction_and_elevation'))

    def _create_from_flow_direction_and_elevation(self, *args, **kwargs):

        flow_dir = kwargs['flow_direction']
        elevation = kwargs['elevation']
        self._copy_info_from_grid(flow_dir, True)
        # isinstance, not an exact class comparison, so GeographicFlowDirection
        # and any other subclass is handled too.
        if isinstance(flow_dir, FlowDirectionD8):
            self.__calc_D8_slope(flow_dir, elevation, *args, **kwargs)

    def __calc_D8_slope(self, fd, dem, *args, **kwargs):

        receivers = fd.receivers.reshape(-1)
        has_receiver = receivers >= 0

        z = np.asarray(dem._griddata, dtype=float64).reshape(-1)
        step = fd.step_lengths(self._mean_pixel_dimension(*args, **kwargs)).reshape(-1)

        slope = np.zeros(z.shape, dtype=float64)
        donors = np.flatnonzero(has_receiver)
        targets = receivers[donors]
        # Vectorised: the original walked every cell in Python, which cost
        # minutes on a large DEM for what is a single subtraction per cell.
        slope[donors] = (z[donors] - z[targets]) / step[donors]

        self._griddata = slope.reshape(self._georef_info.ny, self._georef_info.nx)


def mosaicFolder(folderPath, fileSuffix, outfile):
    """Merge every ``*fileSuffix`` raster in a folder into ``outfile``.

    Uses GDAL's Python API directly rather than shelling out to
    ``gdal_merge.py``, which is not on the path of most installations.
    """
    files = sorted(glob.glob1(folderPath, '*' + fileSuffix))
    if not files:
        raise Error.InputError(
            'mosaicFolder',
            'no files matching *{0} in {1}'.format(fileSuffix, folderPath))

    paths = [os.path.join(folderPath, name) for name in files]
    vrt = gdal.BuildVRT('', paths)
    try:
        gdal.Translate(outfile, vrt, creationOptions=['COMPRESS=LZW'])
    finally:
        vrt = None


def plot(*args, **kwargs):
    """Scatter one grid's values against another's, cell by cell.

    ::

        plot(area, slope, xlabel='Drainage area', ylabel='Slope',
             decimation_factor=50)

    Parameters
    ----------
    xlabel, ylabel : str
        Axis labels.
    symbol : str
        Matplotlib format string; ``'k.'`` by default.
    indexes : tuple of ndarray
        Restrict to these cells instead of every cell valid in both grids.
    decimation_factor : int
        Plot only every n-th point.
    interactive : bool
        Leave the figure open without blocking.

    Returns
    -------
    matplotlib.axes.Axes
    """
    grid1 = args[0]._griddata
    grid2 = args[1]._griddata

    interactive = kwargs.pop('interactive', True)
    # Reading these with .get() into a local only inside the `is None`
    # branch meant that actually *passing* xlabel left the name unbound and
    # raised NameError two lines later.
    xlabel = kwargs.pop('xlabel', None) or 'Grid 1'
    ylabel = kwargs.pop('ylabel', None) or 'Grid 2'
    symbol = kwargs.pop('symbol', None) or 'k.'
    valid_indexes = kwargs.pop('indexes', None)
    decimation_factor = kwargs.pop('decimation_factor', None)

    if valid_indexes is None:
        valid_indexes = np.where(~np.isnan(grid1 + grid2))

    if decimation_factor is not None:
        valid_indexes = tuple(axis[::decimation_factor] for axis in valid_indexes)

    plt.plot(grid1[valid_indexes], grid2[valid_indexes], symbol, **kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if interactive:
        plt.ion()
        plt.show(block=False)
    else:
        plt.ioff()
        plt.show(block=True)

    return plt.gca()
