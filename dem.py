#Functions for demMethods calculations typical of geomorphology
#Developed by Sam Johnstone, January 2015, samuelj@stanford.edu , johnsamstone@gmail.com
#Not all functions are benchmarked
#If you somehow have this, and have never spoken to me, please reach out
#- I'd be curious to hear what you are doing (maybe I can even help with something!)

from osgeo import gdal, ogr, osr # Used to load gis in files
import os # Used to join file paths, iteract with os in other ways..
import glob  # Used for finding files that I want to mosaic (by allowing wildcard searches of the filesystem)
import heapq # Used for constructing priority queue, which is used for filling dems
import numpy as np # Used for tons o stuff, keeping most data stored as numpy arrays
import subprocess # Used to run gdal_merge.py from the command line
from numpy import uint8, int8, float64
from matplotlib import pyplot as plt
import sys
from scipy import stats
from . import error as Error

try:
    import statsmodels.api as sm
except:
    print('Warning: no statsmodels present.  If you want to compute steepness, you will need to install this.')

sys.setrecursionlimit(1000000)

        
class GDALMixin(object):
    
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
        geoTransform, nx, ny, data = self._read_GDAL_dataset(gdal_file, dtype)
        gdal_file = None
        return geoTransform, nx, ny, data
        
    def _read_GDAL_dataset(self, gdal_dataset, dtype):
        band = gdal_dataset.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        data = band.ReadAsArray().astype(dtype)
        nodata_elements = np.where(data == nodata)
        from numpy import uint8
        if dtype is not uint8:
            data[nodata_elements] = np.NAN
        
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

    def _create_gdal_representation_from_array(self, georef_info, GDALDRIVERNAME, array_data, dtype, outfile_path='name', dst_options = [], multiple_bands = False):
        #A function to write the data in the numpy array arrayData into a georeferenced dataset of type
        #  specified by GDALDRIVERNAME, a string, options here: http://www.gdal.org/formats_list.html
        #  This is accomplished by copying the georeferencing information from an existing GDAL dataset,
        #  provided by createDataSetFromArray
    
        #Initialize new data
        drvr = gdal.GetDriverByName(GDALDRIVERNAME)  #  Get the desired driver
        bands = 1 if multiple_bands is False else len(array_data)
        outRaster = drvr.Create(outfile_path, georef_info.nx, georef_info.ny, bands , self._get_gdal_type_for_numpy_type(dtype), dst_options)  # Open the file
    
        #Write geographic information
        outRaster.SetGeoTransform(georef_info.geoTransform)  # Steal the coordinate system from the old dataset
        if georef_info.projection != 0:
            outRaster.SetProjection(georef_info.projection)   # Steal the Projections from the old dataset
    
        #Write the array
        if not multiple_bands:
            outRaster.GetRasterBand(1).WriteArray(array_data)   # Writes my array to the raster
        else:
            for i in range(len(array_data)):
                print('writing band ' + str(i) + "/" + str(bands))
                outRaster.GetRasterBand(i+1).WriteArray(array_data[i])
                
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
    
    def _writeArcAsciiRaster(self, georef_info, outfile_path, np_array_data, nodata_value, format_string):
        #A function to write the data stored in the numpy array npArrayData to a ArcInfo Text grid. Gdal doesn't
        #allow creation of these types of data for whatever reason
    
        header = "ncols     %s\n" % georef_info.nx
        header += "nrows    %s\n" % georef_info.ny
        header += "xllcenter %s\n" % georef_info.xllcenter
        header += "yllcenter %s\n" % georef_info.yllcenter
        header += "cellsize %s\n" % georef_info.dx
        header += "NODATA_value %s" % nodata_value
    
        np.savetxt(outfile_path, np_array_data, header=header, fmt=format_string, comments='')
    
    def _asciiRasterToMemory(self, fileName):

        # the geotransfrom structured as (xUL, dx, skewX, yUL, scewY, -dy)
        gt, nx, ny = self._getRasterGeoTransformFromAsciiRaster(fileName)
    
        #Open gdal dataset
        ds = gdal.Open(fileName)
    
        #Prep output
        memDrv = gdal.GetDriverByName('MEM')  # Create a gdal driver in memory
        dataOut = memDrv.CreateCopy('name', ds, 0) #Copy data to gdal driver
    
        data = dataOut.ReadAsArray()
        dataOut = None
        
        data = self.dtype(data)
        return gt, nx, ny, data

class GeographicGridMixin(object):
    
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
    
    
    def _area_per_pixel(self, *args, **kwargs):

        if hasattr(self, '_GeographicGridMixin__app'):
            return self._GeographicGridMixin__app
        
        #Returns a grid of the area of each pixel within the domain specified by the gdalDataset
    
        re = 6371.0 * 1000.0 #radius of the earth in meters
    
        dLong = np.abs(self._georef_info.geoTransform[1])
        dLat = np.abs(self._georef_info.geoTransform[5])
        
        lats = self._georef_info.yllcenter + np.float64(range(self._georef_info.ny))*dLat
        longs = self._georef_info.xllcenter + np.float64(range(self._georef_info.nx))*dLong
            
        #Get the size of pixels in the lat and long direction

        [LONG, LAT] = np.meshgrid(longs, lats)
        LAT1 = np.radians(LAT - dLat / 2.0)
        LAT2 = np.radians(LAT + dLat / 2.0)
        LONG1 = np.radians(LONG - dLong / 2.0)
        LONG2 = np.radians(LONG + dLong / 2.0)
        
        initialAreas = np.flipud(np.abs((re**2)*(LONG1 - LONG2)*(np.sin(LAT2) - np.sin(LAT1))))
        
        self.__app = initialAreas
        return self.__app
    
    def _mean_pixel_dimension(self, *args, **kwargs):
        
        return np.sqrt(self._area_per_pixel())


class CalculationMixin(object):
    
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
        Zbc[-1, -1] = grid[-1, 0]
    
        return Zbc

    def calcFiniteCurv(self, grid, dx):
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
        ## Sx,Sy = calcAverageSlopeOfGridSubset(numpy matrix, dx)
        #Compute the average slope over a subset of a grid (or a whole grid if you're into that),
        #by fitting a plane to the elevation data stored in grid subset
        nx,ny = len(gridSubset[1,:]), len(gridSubset[:,1])
        xs = (0.5+np.arange(nx))*dx - (nx*dx)/2.0
        ys = (0.5+np.arange(ny))*dx - (ny*dx)/2.0
        X,Y = np.meshgrid(xs,ys)
        #Fit a plane of the form z = ax + by + c, where beta = [a b c]
        M=np.vstack((X.flatten(),Y.flatten(),np.ones((1,nx*ny)))).T
        beta = np.linalg.lstsq(M,gridSubset.flatten())[0]
    
        return beta[0], beta[1] #Return the slope in the x and y directions respecticely

class Georef_info(object):
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
    # Wrapper for GDAL shapes.
    def __init__(self, *args, **kwargs):
        if kwargs.get('shapefile_name') is None:
            raise Error.InputError('Input Error', 'Inputs not satisfied')
        self.shapedata = ogr.Open(kwargs.get('shapefile_name'))

    def createMaskFromShape(self, geoRefInfo, projection, dtype, noDataValue = 0):
    
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
        
        if evaluative_action == None:
            raise Error.InputError('Input Error', 'Inputs not satisfied')
        
        eval('self.' + evaluative_action + '(*args, **kwargs)')
    
    def __setitem__(self, key, value):
        i, j = key
        if i < 0 or j < 0 or i > self._georef_info.ny or j > self._georef_info.nx:
            return
        self._griddata[i, j] = value
    
    def __getitem__(self, key):
        i, j = key

        i = int(i)
        j = int(j)
        
        if i >= 0 and j >= 0 and i < self._georef_info.ny and j < self._georef_info.nx:
            try:
                return self._griddata[i, j]
            except:
                print('error: '+ str(i) + str(j))
        
        return None
        
    def __get_evaluative_action(self, *args, **kwargs):
                
        for required_input_set, evaluative_action in self.required_inputs_and_actions:
            these_kw = set(kwargs.keys())
            required_kw = set(required_input_set)
            if required_kw.issubset(these_kw):
                return evaluative_action 
        
        return None
    
    def __populate_georef_info_using_geoTransform(self, nx, ny):
        self._georef_info.dx = self._georef_info.geoTransform[1]
        self._georef_info.xllcenter = self._georef_info.geoTransform[0]+self._georef_info.dx/2.0
        self._georef_info.yllcenter = self._georef_info.geoTransform[3]-(self._georef_info.dx*(self._georef_info.ny-0.5))
        self._georef_info.nx = nx
        self._georef_info.ny = ny
        
    def _create(self, *args, **kwargs):
            
        self._georef_info.geoTransform = kwargs.get('geo_transform')
        self._georef_info.projection = kwargs.get('projection')
        self.__populate_georef_info_using_geoTransform(kwargs['nx'], kwargs['ny'])
        if kwargs.get('grid') is None:
            self._griddata = np.zeros(shape = (self._georef_info.ny,self._georef_info.nx), dtype = self.dtype)
        else:
            self._griddata = kwargs.get('grid')

    def _create_random_grid(self, *args, **kwargs):
        self._georef_info.dx = kwargs['dx']
        self._georef_info.nx = kwargs['nx']
        self._georef_info.ny = kwargs['ny']
        self._georef_info.xllcenter = 0
        self._georef_info.yllcenter = 0
        self._randomize_grid_values(*args, **kwargs)

    def _create_from_grid(self, *args, **kwargs):
        self._georef_info.dx = kwargs['dx']
        (ny, nx) = kwargs['grid'].shape
        self._georef_info.nx = nx
        self._georef_info.ny = ny
        self._georef_info.xllcenter = 0
        self._georef_info.yllcenter = 0
        self._griddata = kwargs['grid']
    
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

            if col > self._georef_info.nx or row > self._georef_info.ny or col < 0 or row < 0:
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

    def get_XY_matricies(self):
        xllc, yllc, nx, ny, dx = (self._georef_info.xllcenter, self._georef_info.yllcenter, self._georef_info.nx, self._georef_info.ny, self._georef_info.dx)
        x = np.arange(xllc, xllc+(nx-1)*dx, dx)
        y = np.arange(yllc, yllc+(ny-1)*dx, dx)
        return np.meshgrid(x,y)

    def resample(self, de):
        return_grid = self.__class__()
        return_grid._copy_info_from_grid(self, set_zeros=True)
        import copy
        geoTransform = (self._georef_info.geoTransform[0], 
                        np.sign(self._georef_info.geoTransform[1])*de, 
                        self._georef_info.geoTransform[2], 
                        self._georef_info.geoTransform[3], 
                        self._georef_info.geoTransform[4], 
                        np.sign(self._georef_info.geoTransform[5])*de)
        
        georef = copy.deepcopy(self._georef_info)
        georef.geoTransform = geoTransform
        georef.dx = de
        georef.nx = int(np.round(self._georef_info.nx * self._georef_info.dx / de))
        georef.ny = int(np.round(self._georef_info.ny * self._georef_info.dx / de))
        return_grid._georef_info = georef
        xo = np.linspace(self._georef_info.xllcenter,(self._georef_info.nx-1)*(self._georef_info.dx), num = self._georef_info.nx)
        yo = np.linspace(self._georef_info.yllcenter,(self._georef_info.ny-1)*(self._georef_info.dx), num = self._georef_info.ny)
        xf = np.linspace(georef.xllcenter,(georef.nx-1)*georef.dx, num = georef.nx)
        yf = np.linspace(georef.yllcenter,(georef.ny-1)*georef.dx, num = georef.ny)
        from scipy.interpolate import interp2d
        interp_grid = copy.deepcopy(self._griddata)
        if geoTransform[1] > 0:
            interp_grid = np.fliplr(interp_grid)
        if geoTransform[5] < 0:
            interp_grid = np.flipud(interp_grid)
        f = interp2d(xo, yo, interp_grid, kind='quintic')
        return_grid._griddata = f(xf, yf)
        return return_grid
    
    def location_in_grid(self, xo):
        index = self._xy_to_rowscols((xo,))[0]
        return index[0] is not None and index[1] is not None and self[index[0], index[1]] is not None and self[index[0], index[1]] != 0 and not np.isnan(self[index[0], index[1]])
    
    def tile(self, tile_xdim = 400, tile_ydim=400, tile_xpadding = 10, tile_ypadding = 10):
        
        num_xtiles = int(np.ceil(float(self._georef_info.nx) / float(tile_xdim)))
        num_ytiles = int(np.ceil(float(self._georef_info.ny) / float(tile_ydim)))
        
        tiles = []
        import copy
        
        for x_left in range(num_xtiles):
            for y_top in range(num_ytiles):
                
                left = x_left*tile_xdim-tile_xpadding if (x_left*tile_xdim-tile_xpadding) >= 0 else 0
                right = (x_left+1)*tile_xdim+tile_xpadding if ((x_left+1)*tile_xdim+tile_xpadding) <= (self._georef_info.nx-1) else self._georef_info.nx
                top = y_top*tile_ydim-tile_ypadding if (y_top*tile_ydim-tile_ypadding) >= 0 else 0
                bottom = (y_top+1)*tile_ydim+tile_ypadding if ((y_top+1)*tile_ydim+tile_ypadding) <= (self._georef_info.ny-1) else self._georef_info.ny
                
                new_tile = self.__class__()
                
                needs_left_padding = (x_left == 0)
                needs_right_padding = (x_left == (num_xtiles-1))
                
                needs_top_padding = (y_top == 0)
                needs_bottom_padding = (y_top == (num_ytiles-1))
                
                griddata = self._griddata[top:bottom,left:right]
                
                xllcenter = self._georef_info.xllcenter + left*self._georef_info.dx
                yllcenter = self._georef_info.yllcenter + (self._georef_info.ny - bottom)*self._georef_info.dx
                nx = right-left
                ny = bottom - top
                
                if needs_left_padding:
                    griddata = np.concatenate((np.fliplr(griddata[:,0:tile_xpadding]), griddata), axis = 1)
                    xllcenter -= tile_xpadding*self._georef_info.dx
                    nx += tile_xpadding
                    
                if needs_right_padding:
                    griddata = np.concatenate((griddata, np.fliplr(griddata[:,-tile_xpadding:])), axis = 1)
                    nx += tile_xpadding
                
                if needs_top_padding:
                    griddata = np.concatenate((np.flipud(griddata[0:tile_ypadding,:]), griddata), axis = 0)
                    ny += tile_ypadding
                
                if needs_bottom_padding:
                    griddata = np.concatenate((griddata, np.flipud(griddata[-tile_ypadding:,:])), axis = 0)
                    ny += tile_ypadding
                    yllcenter -= tile_ypadding*self._georef_info.dx
    
                new_tile._georef_info.geoTransform = (xllcenter - 0.5*self._georef_info.dx, self._georef_info.dx, 0, yllcenter + (ny + 0.5)*self._georef_info.dx, 0, -self._georef_info.dx)
                new_tile._georef_info.xllcenter = xllcenter
                new_tile._georef_info.yllcenter = yllcenter
                new_tile._georef_info.dx = self._georef_info.dx
                new_tile._georef_info.nx = nx
                new_tile._georef_info.ny = ny
                new_tile._griddata = griddata
                
                tiles += [copy.deepcopy(new_tile)]
        
        return tiles
    
    def remove_padding(self, xpadding = 10, ypadding = 10):
        
        import copy
        unpadded = copy.deepcopy(self)
        unpadded._griddata =  unpadded._griddata[ypadding:-ypadding,xpadding:-xpadding]
        unpadded._georef_info.xllcenter += xpadding*unpadded._georef_info.dx
        unpadded._georef_info.yllcenter += ypadding*unpadded._georef_info.dx
        unpadded._georef_info.nx -= 2*xpadding
        unpadded._georef_info.ny -= 2*ypadding
        unpadded._georef_info.geoTransform = (unpadded._georef_info.xllcenter - 0.5*unpadded._georef_info.dx, unpadded._georef_info.dx, 0, unpadded._georef_info.yllcenter + (float(unpadded._georef_info.ny-0.5))*self._georef_info.dx, 0, -self._georef_info.dx)
        
        return unpadded
        
    def sort(self, reverse=True, force = False, mask = None):
        if force:
            self._sorted = False
        if not self._sorted:
            if mask is not None:
                sort_indexes = self._griddata.argsort(axis = None)
                unraveled_sort_indexes = np.unravel_index(sort_indexes, mask.shape)
                query_indexes = np.where(mask._griddata[unraveled_sort_indexes] == 1)
                query_indexes_raveled = np.ravel_multi_index(query_indexes,mask.shape)
                self._sort_indexes = sort_indexes[query_indexes_raveled]
            else:
                self._sort_indexes = self._griddata.argsort(axis = None)
            self._sorted = True
        if reverse:
            return self._sort_indexes[::-1]
        else:
            return self._sort_indexes
    
    def apply_moving_window(self, moving_window):
        out_grid = None
        out_grid._georef_info = self._georef_info
        out_grid._griddata = moving_window.apply_moving_window(self._griddata, self._georef_info.dx, self.dtype)
        return out_grid
    
    def average_over_distance(self, distance, grid = None):
        (nx, ny) = (self._georef_info.nx, self._georef_info.ny)
        (dx, dy) = (self._georef_info.dx, self._georef_info.dx)
        (xmin, xmax) = (self._georef_info.xllcenter, nx*dx+self._georef_info.xllcenter)
        (ymin, ymax) = (self._georef_info.yllcenter, ny*dy+self._georef_info.yllcenter)
        (xcenter, ycenter) = ((xmin+xmax) / 2.0, (ymin+ymax) / 2.0)
        [x,y] = np.meshgrid(np.arange(xmin, xmax, dx),np.arange(ymin, ymax, dy))
        D = np.sqrt(np.power((x - xcenter),2) + np.power((y - ycenter), 2))
        template = (D <= distance).astype(float)
        template /= np.sum(template[:])
        if grid is None:
            grid = self._griddata
        
        from numpy.fft import fft2 as fft2
        from numpy.fft import ifft2 as ifft2
        from numpy.fft import ifftshift
        
        return np.real(ifftshift(ifft2(fft2(grid)*fft2(template))))
        
    def clip_to_bounds(self, bounds):
        extent = (bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
        return self.clip_to_extent(extent)
    
    def clip_to_extent(self, extent):
        import copy
        return_grid = copy.deepcopy(self)
        return_grid._georef_info = copy.deepcopy(self._georef_info)
        lower_left = (extent[0]-self._georef_info.dx, extent[2]-self._georef_info.dx)
        upper_right = (extent[1]+self._georef_info.dx, extent[3]+self._georef_info.dx)
        idx = self._xy_to_rowscols((lower_left, upper_right))
        return_grid._griddata = return_grid._griddata[idx[1][0]+1:idx[0][0]+1,idx[0][1]:idx[1][1]]
        return_grid._georef_info.nx = return_grid._griddata.shape[1]
        return_grid._georef_info.ny = return_grid._griddata.shape[0]
        return_grid._georef_info.xllcenter = self._rowscols_to_xy((idx[0],))[0][0]
        return_grid._georef_info.yllcenter = self._rowscols_to_xy((idx[0],))[0][1]
        return_grid._georef_info.geoTransform = (return_grid._georef_info.xllcenter, return_grid._georef_info.geoTransform[1], return_grid._georef_info.geoTransform[2], return_grid._georef_info.yllcenter + (return_grid._georef_info.dx*(return_grid._georef_info.ny-1)), return_grid._georef_info.geoTransform[4], return_grid._georef_info.geoTransform[5],)
        
        return return_grid
    
    def extent_of_data(self):
        
        top = 0
        bottom = self._georef_info.ny - 1
        left = 0
        right = self._georef_info.nx - 1
        
        while np.sum(np.isnan(self._griddata[top,:])) == self._georef_info.nx:
            top = top + 1
        
        while np.sum(np.isnan(self._griddata[bottom,:])) == self._georef_info.nx:
            bottom = bottom - 1
        
        while np.sum(np.isnan(self._griddata[left,:])) == self._georef_info.ny:
            left = left + 1
            
        while np.sum(np.isnan(self._griddata[right,:])) == self._georef_info.ny:
            right = right - 1
        
        (ll, ur) = self._rowscols_to_xy(((bottom, left), (top, right)))   
        return (ll[0], ur[0], ll[1], ur[1])
    
    
    def clip_to_mask_grid(self, mask_grid):
        gdal_source = self._create_gdal_representation_from_array(self._georef_info, 'MEM', self._griddata, self.dtype)
        gdal_mask = mask_grid._create_gdal_representation_from_array(mask_grid._georef_info, 'MEM', mask_grid._griddata, mask_grid.dtype)
        gdal_clip = self._clipRasterToRaster(gdal_source, gdal_mask, self.dtype)
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._read_GDAL_dataset(gdal_clip, self.dtype)
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)

    def clip_to_shapefile(self, shapefile):
        gdal_source = self._create_gdal_representation_from_array(self._georef_info, 'MEM', self._griddata, self.dtype)
        gdal_clip = shapefile.shapedata
        gdal_result = self._clipRasterToShape(gdal_source, gdal_clip)
        self._georef_info.geoTransform, self._georef_info.nx, self._georef_info.ny, self._griddata = self._read_GDAL_dataset(gdal_result, self.dtype)
        self.__populate_georef_info_using_geoTransform(self._georef_info.nx, self._georef_info.ny)
    
    def calculate_gradient_over_length_scale(self,length_scale):
        # sx,sy = calcFiniteDiffs(elevGrid,dx)
        # calculates finite differences in X and Y direction using a finite difference
        # kernel that extends N cells from the center. The width of the kernel is then
        # 2N + 1 by 2N + 1. Applies a boundary condition such that the size and location
        # of the grids in is the same as that out. However, the larger N is, the more NoData
        #Will be around the edges .
    
        #Compute finite differences
        elevGrid = self._griddata
        dx = self._georef_info.dx
        N = np.ceil(length_scale / dx)
        
        Sx = (elevGrid[N:-N, (2*N):] - elevGrid[N:-N, :-(2*N)])/(((2*N)+1)*dx)
        Sy = (elevGrid[(2*N):,N:-N] - elevGrid[:-(2*N), N:-N])/(((2*N)+1)*dx)
    
        #Create two new arrays of the original DEMs size
        SxPadded = np.empty(elevGrid.shape)
        SxPadded[:] = np.NAN
        SyPadded = np.empty(elevGrid.shape)
        SyPadded[:] = np.NAN
    
        SyPadded[N:-N, N:-N] = Sy
        SxPadded[N:-N, N:-N] = Sx
    
        return SxPadded, SyPadded
    
    def calculate_laplacian_over_length_scale(self, length_scale):
        #C = calcFiniteCurv(elevGrid, dx)
        #calculates finite differnces in X and Y direction using the centered difference method.
        #Applies a boundary condition such that the size and location of the grids in is the same as that out.
    
        dx = self._georef_info.dx
        grid = self._griddata
        winRad = np.ceil(length_scale / dx)
        #Assign boundary conditions
        Curv = np.zeros_like(grid)*np.nan
    
        #Compute finite differences
        Cx = (grid[winRad:-winRad, (2*winRad):] - 2*grid[winRad:-winRad, winRad:-winRad] + grid[winRad:-winRad, :-(2*winRad)])/(2*dx*winRad)**2
        Cy = (grid[(2*winRad):, winRad:-winRad] - 2*grid[winRad:-winRad, winRad:-winRad] + grid[:-(2*winRad), winRad:-winRad])/(2*dx*winRad)**2
    
        Curv[winRad:-winRad,winRad:-winRad] = Cx+Cy
        return Curv

    def principal_curvatures(self):
        Zy, Zx  = np.gradient(self._griddata, self._georef_info.dx)
        Zxy, Zxx = np.gradient(Zx, self._georef_info.dx)
        Zyy, _ = np.gradient(Zy, self._georef_info.dx)
        
        H = (Zx**2 + 1)*Zyy - 2*Zx*Zy*Zxy + (Zy**2 + 1)*Zxx
        H = -H/(2*(Zx**2 + Zy**2 + 1)**(1.5))
        
        K = (Zxx * Zyy - (Zxy ** 2)) /  (1 + (Zx ** 2) + (Zy **2)) ** 2
        
        k1 = H + np.sqrt(np.power(H,2) - K)
        k2 = H - np.sqrt(np.power(H,2) - K)
        
        k1re = BaseSpatialGrid()
        k1re._copy_info_from_grid(self, True)
        k2re = BaseSpatialGrid()
        k2re._copy_info_from_grid(self, True)
        k1re._griddata = k1
        k2re._griddata = k2
        
        return (k1re, k2re)
        
    def plot(self, **kwargs):

        interactive = kwargs.pop('interactive', True)
        colorbar = kwargs.pop('colorbar', True)
        extent = [self._georef_info.xllcenter, self._georef_info.xllcenter+(self._georef_info.nx-0.5)*self._georef_info.dx, self._georef_info.yllcenter, self._georef_info.yllcenter+(self._georef_info.ny-0.5)*self._georef_info.dx]
        fig = plt.figure()
        plt.axis(extent)
        if interactive:
            plt.ion()
            plt.show(block = False)
        else:
            plt.ioff()
        im = plt.imshow(self._griddata, extent = extent, **kwargs)
        if colorbar:
            plt.colorbar(im)
        plt.draw()
        plt.pause(0.001)
        if interactive:
            plt.show()
        return fig

    
    def find_nearest_cell_with_value(self, index, value, pixel_radius):

        (i, j) = index
        nearest = tuple()
        minimum_difference = np.nan
        for y in range(int(i-pixel_radius),int(i+pixel_radius)):
            for x in range(int(j-pixel_radius),int(j+pixel_radius)):
                if x is not None and y is not None and np.sqrt( (y-i)**2 + (x-j)**2) <= pixel_radius and x < self._griddata.shape[1] and y < self._griddata.shape[0]:
                    value_difference = np.abs(value - self._griddata[y,x])
                    if np.isnan(minimum_difference) or minimum_difference > value_difference:
                        minimum_difference = value_difference
                        nearest = (y,x)
        (y,x) = nearest
        return y, x
    
    def find_nearest_cell_with_value_greater_than(self, index, value, pixel_radius):

        (i, j) = index
        nearest = tuple()
        minimum_difference = np.nan
        for y in range(int(i-pixel_radius),int(i+pixel_radius)):
            for x in range(int(j-pixel_radius),int(j+pixel_radius)):
                if x is not None and y is not None and np.sqrt( (y-i)**2 + (x-j)**2) <= pixel_radius and x < self._griddata.shape[1] and y < self._griddata.shape[0]:
                    value_difference = self._griddata[y,x] - value
                    if (np.isnan(minimum_difference) or (minimum_difference > value_difference)) and (value_difference >= 0):
                        minimum_difference = value_difference
                        nearest = (y,x)
        if np.isnan(minimum_difference):
            for y in range(int(i-pixel_radius),int(i+pixel_radius)):
                for x in range(int(j-pixel_radius),int(j+pixel_radius)):
                    if x is not None and y is not None and np.sqrt( (y-i)**2 + (x-j)**2) <= pixel_radius and x < self._griddata.shape[1] and y < self._griddata.shape[0]:
                        value_difference = value - self._griddata[y,x]
                        if (np.isnan(minimum_difference) or (minimum_difference > value_difference)) :
                            minimum_difference = value_difference
                            nearest = (y,x)
        (y,x) = nearest
        return y, x
    
    def find_nearest_cell_with_greatest_value(self, index, pixel_radius = 5):
        
        (i, j) = index
        nearest = tuple()
        maximum_value = np.nan
        for y in range(int(i-pixel_radius),int(i+pixel_radius)):
            for x in range(int(j-pixel_radius),int(j+pixel_radius)):
                if np.sqrt( (y-i)**2 + (x-j)**2) <= pixel_radius:
                    if np.isnan(maximum_value) or maximum_value < self._griddata[y,x]:
                        maximum_value = self._griddata[y,x]
                        nearest = (y,x)
        
        (y,x) = nearest
        return y, x
    

    def snap_locations_to_greatest_value(self, v, pixel_radius = 5):
        
        idxs = self._xy_to_rowscols(v)
        snap_idxs = list()
        for idx in idxs:
            i, j = self.find_nearest_cell_with_greatest_value(idx, pixel_radius)
            snap_idxs.append((i,j))
        return self._rowscols_to_xy(snap_idxs)
    
    def snap_locations_to_closest_value(self, v, value, pixel_radius = 5):
        idxs = self._xy_to_rowscols(v)
        snap_idxs = list()
        values_with_idx = zip(idxs, value)
        for (idx, value) in values_with_idx:
            i, j = self.find_nearest_cell_with_value(idx, value, pixel_radius)
            snap_idxs.append((i,j))
        return self._rowscols_to_xy(snap_idxs)
    
    def export_arc_e00_grid(self, filename):

        self._create_gdal_representation_from_array(self._georef_info, 'E00GRID', self._griddata, self.dtype, filename)

    def set_value_at_rowscols(self, value, rowscols):
        rowscols_array = (np.array(list(zip(*rowscols))).T[:,0],np.array(list(zip(*rowscols))).T[:,1])
        self._griddata[rowscols_array] = value
        
    def save(self, filename):
        
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', self._griddata, self.dtype, filename, ['COMPRESS=LZW'])
    
    def write_to_ai(self, filename):
        self._writeArcAsciiRaster(self._georef_info, filename, self._griddata, np.NAN, '%10.2f')
    
    def vectorize(self, filename):
        
        import uuid
        tmpfilename = str(uuid.uuid4())
        self.save(tmpfilename)
        gdal_dataset = gdal.Open(tmpfilename)
        srcband = gdal_dataset.GetRasterBand(1)
        
        dst_layername = filename
        drv = ogr.GetDriverByName("ESRI Shapefile")
        dst_ds = drv.CreateDataSource( dst_layername + ".shp" )
        dst_layer = dst_ds.CreateLayer(dst_layername, srs = None )
        newField = ogr.FieldDefn('RASTER_VAL', ogr.OFTReal)
        dst_layer.CreateField(newField)
        gdal.Polygonize( srcband, None, dst_layer, 0, [], callback=None )
        os.remove(tmpfilename)
        
    @classmethod
    def load(cls, filename):
        
        return_object = cls()
        
        gdal_dataset = gdal.Open(filename)
        band = gdal_dataset.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        return_object._griddata = band.ReadAsArray().astype(cls.dtype)
        if nodata is not None:
            nodata_elements = np.where(return_object._griddata == nodata)
            from numpy import uint8
            if cls.dtype is not uint8:
                return_object._griddata[nodata_elements] = np.NAN
        
        
        geoTransform = gdal_dataset.GetGeoTransform()
        nx = gdal_dataset.RasterXSize
        ny = gdal_dataset.RasterYSize
        
        return_object._georef_info.geoTransform = geoTransform
        return_object._georef_info.dx = return_object._georef_info.geoTransform[1]
        return_object._georef_info.xllcenter = return_object._georef_info.geoTransform[0]+return_object._georef_info.dx/2.0
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
            
        gdal_file = None
        return return_object

    @classmethod
    def mosaic(cls, tiles):
        
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
            
        return return_object
    
class ValueGrid(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('dx', 'grid'), '_create_from_grid'))
    
    def set_value_at_indexes(self, indexes, value):
        ij = zip(*indexes)
        self._griddata[ij] = value
        
    
class FlowDirection(BaseSpatialGrid):
    pass

class FlowDirectionD8(FlowDirection):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('flooded_dem',), '_create_from_flooded_dem'))
    
    from numpy import uint8
    dtype = uint8

    def _create_from_flooded_dem(self, *args, **kwargs):
        flooded_dem = kwargs['flooded_dem']
        self._copy_info_from_grid(flooded_dem)
        self._griddata = np.zeros_like(flooded_dem._griddata, dtype = self.dtype)
        
        flooded_dem_ijp1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_ijp1[:,0:-2] = (flooded_dem._griddata[:,0:-2] - flooded_dem._griddata[:, 1:-1])
        
        flooded_dem_ip1jp1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_ip1jp1[0:-2,0:-2] = (flooded_dem._griddata[0:-2,0:-2] - flooded_dem._griddata[1:-1, 1:-1]) / 1.41
        
        flooded_dem_ip1j = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_ip1j[0:-2,:] = (flooded_dem._griddata[0:-2,:] - flooded_dem._griddata[1:-1, :])
        
        flooded_dem_ip1jm1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_ip1jm1[0:-2,1:-1] = (flooded_dem._griddata[0:-2,1:-1] - flooded_dem._griddata[1:-1, 0:-2]) / 1.41
        
        flooded_dem_ijm1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_ijm1[:,1:-1] = (flooded_dem._griddata[:,1:-1] - flooded_dem._griddata[:, 0:-2])
        
        flooded_dem_im1jm1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_im1jm1[1:-1,1:-1] = (flooded_dem._griddata[1:-1,1:-1] - flooded_dem._griddata[0:-2, 0:-2]) / 1.41
        
        flooded_dem_im1j = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_im1j[1:-1,:] = (flooded_dem._griddata[1:-1,:] - flooded_dem._griddata[0:-2, :])
        
        flooded_dem_im1jp1 = -np.ones_like(flooded_dem._griddata) * 1E22
        flooded_dem_im1jp1[1:-1,0:-2] = (flooded_dem._griddata[1:-1,0:-2] - flooded_dem._griddata[0:-2, 1:-1])  / 1.41
        
        min_sl = np.fmax(np.fmax(np.fmax(np.fmax(np.fmax(np.fmax(np.fmax(flooded_dem_ijp1, flooded_dem_ip1jp1), flooded_dem_ip1j), flooded_dem_ip1jm1), flooded_dem_ijm1), flooded_dem_im1jm1), flooded_dem_im1j), flooded_dem_im1jp1)
        
        self._griddata[flooded_dem_ijp1 == min_sl] = 1
        self._griddata[flooded_dem_ip1jp1 == min_sl] = 2
        self._griddata[flooded_dem_ip1j == min_sl] = 4
        self._griddata[flooded_dem_ip1jm1 == min_sl] = 8
        self._griddata[flooded_dem_ijm1 == min_sl] = 16
        self._griddata[flooded_dem_im1jm1 == min_sl] = 32
        self._griddata[flooded_dem_im1j == min_sl] = 64
        self._griddata[flooded_dem_im1jp1 == min_sl] = 128
        
        self._griddata[np.isnan(flooded_dem._griddata)] = 0
        
        mask = kwargs.get('mask')
        self._sort_indexes = flooded_dem.sort(reverse = False, force = True, mask = mask)
        self._sorted = True
        
    def __flow_code_for_position(self, flooded_dem, i, j):
        
        flow_code = 0
        min_diff = 0
                
        if flooded_dem[i,j+1] is not None and ((flooded_dem[i, j] - flooded_dem[i,j+1])) / 1.41 > min_diff:
            flow_code = 1
            min_diff = (flooded_dem[i, j] - flooded_dem[i,j+1]) / 1.41
        
        if flooded_dem[i+1,j+1] is not None and ((flooded_dem[i, j] - flooded_dem[i+1,j+1])) > min_diff:
            flow_code = 2
            min_diff = flooded_dem[i, j] - flooded_dem[i+1,j+1]
        
        if flooded_dem[i+1,j] is not None and ((flooded_dem[i, j] - flooded_dem[i+1,j])) / 1.41 > min_diff:
            flow_code = 4
            min_diff = (flooded_dem[i, j] - flooded_dem[i+1,j]) / 1.41
        
        if flooded_dem[i+1,j-1] is not None and ((flooded_dem[i, j] - flooded_dem[i+1,j-1])) > min_diff:
            flow_code = 8
            min_diff = flooded_dem[i, j] - flooded_dem[i+1,j-1]
        
        if flooded_dem[i,j-1] is not None and ((flooded_dem[i, j] - flooded_dem[i,j-1])) / 1.41 > min_diff:
            flow_code = 16
            min_diff = (flooded_dem[i, j] - flooded_dem[i,j-1]) / 1.41
            
        if flooded_dem[i-1,j-1] is not None and ((flooded_dem[i, j] - flooded_dem[i-1,j-1])) > min_diff:
            flow_code = 32
            min_diff = flooded_dem[i, j] - flooded_dem[i-1,j-1]
            
        if flooded_dem[i-1,j] is not None and ((flooded_dem[i, j] - flooded_dem[i-1,j]) / 1.41) > min_diff:
            flow_code = 64
            min_diff = (flooded_dem[i, j] - flooded_dem[i-1,j]) / 1.41
            
        if flooded_dem[i-1,j+1] is not None and ((flooded_dem[i, j] - flooded_dem[i-1,j+1])) > min_diff:
            flow_code = 128
            min_diff = flooded_dem[i, j] - flooded_dem[i-1,j+1]
            
        return flow_code
    
    def get_flow_to_cell(self,i,j):
        #Function to get the indices of the cell that is drained to based on the flow direction specified in fd
            
        iOut = None
        jOut = None
        isGood = False

        if self._griddata[i,j] == 1 and j+1 < self._georef_info.nx:
            iOut = i
            jOut = j+1
        elif self._griddata[i,j] == 2 and i+1 < self._georef_info.ny and j+1 < self._georef_info.nx:
            iOut = i+1
            jOut = j+1
        elif self._griddata[i,j] == 4 and i+1 < self._georef_info.ny:
            iOut = i+1
            jOut = j
        elif self._griddata[i,j] == 8 and i+1 < self._georef_info.ny and j-1 >= 0:
            iOut = i+1
            jOut = j-1
        elif self._griddata[i,j] == 16 and j-1 >= 0:
            iOut = i
            jOut = j-1
        elif self._griddata[i,j] == 32 and i-1 >= 0 and j-1 >= 0:
            iOut = i-1
            jOut = j-1
        elif self._griddata[i,j] == 64 and i-1 >= 0:
            iOut = i-1
            jOut = j
        elif self._griddata[i,j] == 128 and i-1 >= 0 and j+1 < self._georef_info.nx:
            iOut = i-1
            jOut = j+1
    
        if not(iOut is None):
            isGood = True
    
        return iOut, jOut, isGood

    def get_upstream_cell_indexes(self, i, j):
        
        options = list()
        
        if self[i, j+1] == 16:
            options += [(i, j+1)]
        
        if self[i+1, j+1] == 32:
            options += [(i+1, j+1)]
            
        if self[i+1, j] == 64:
            options += [(i+1, j)]
                    
        if self[i+1, j-1] == 128:
            options += [(i+1, j-1)]
                
        if self[i, j-1] == 1:
            options += [(i, j-1)]
         
        if self[i-1, j-1] == 2:
            options += [(i-1, j-1)]
            
        if self[i-1, j] == 4:
            options += [(i-1, j)]
            
        if self[i-1, j+1] == 8:
            options += [(i-1, j+1)]

        return options
    
    def __get_flow_from_cell(self, i, j, max_recursion_depth=None, depth=0):
        
        i_source = [i]
        j_source = [j]

        if max_recursion_depth is not None and depth > max_recursion_depth:
            return i_source, j_source

        recursion_depth = depth + 1


        if self[i, j+1] == 16:
            [i_append, j_append] = self.__get_flow_from_cell(i,j+1,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
        
        if self[i+1, j+1] == 32:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j+1,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i+1, j] == 64:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
        
        if self[i+1, j-1] == 128:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j-1,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
                
        if self[i, j-1] == 1:
            [i_append, j_append] = self.__get_flow_from_cell(i,j-1,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
         
        if self[i-1, j-1] == 2:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j-1,max_recursion_depth=max_recursion_depth,depth=recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i-1, j] == 4:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j,max_recursion_depth = max_recursion_depth,depth = recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i-1, j+1] == 8:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j+1,max_recursion_depth = max_recursion_depth,depth = recursion_depth)
            i_source = i_source + i_append
            j_source = j_source + j_append
               
        return i_source, j_source

    def __map_flow_from_cell(self, index, **kwargs):

        if len(index) == 1:
            (i, j) = index[0]
        else:
            (i, j) = index

        return_dict = dict()
        
        return_dict['index'] = index
        for arg in kwargs:
            if arg != 'assume_valid_routing':
                return_dict[arg] = kwargs[arg][i,j]
        
        return_dict['next'] = []
        return_dict['distance_scale'] = 1.0
        
        if self[i, j+1] == 16:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i, j+1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i, j+1] = True
                return_dict['distance_scale'] = 1.0
                child_dict = self.__map_flow_from_cell((i,j+1), **kwargs)
                return_dict['next'].append(child_dict)
        
        if self[i+1, j+1] == 32:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i+1, j+1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i+1, j+1] = True
                return_dict['distance_scale'] = 1.4142135623730951
                child_dict = self.__map_flow_from_cell((i+1,j+1), **kwargs)
                return_dict['next'].append(child_dict)

        if self[i+1, j] == 64:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i+1, j]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i+1, j] = True
                return_dict['distance_scale'] = 1.0
                child_dict = self.__map_flow_from_cell((i+1,j), **kwargs)
                return_dict['next'].append(child_dict)
        
        if self[i+1, j-1] == 128:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i+1, j-1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i+1, j-1] = True
                return_dict['distance_scale'] = 1.4142135623730951
                child_dict = self.__map_flow_from_cell((i+1,j-1), **kwargs)
                return_dict['next'].append(child_dict)
                
        if self[i, j-1] == 1:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i, j-1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i, j-1] = True
                return_dict['distance_scale'] = 1.0
                child_dict = self.__map_flow_from_cell((i,j-1), **kwargs)
                return_dict['next'].append(child_dict)

        if self[i-1, j-1] == 2:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i-1, j-1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i-1, j-1] = True
                return_dict['distance_scale'] = 1.4142135623730951
                child_dict = self.__map_flow_from_cell((i-1,j-1), **kwargs)
                return_dict['next'].append(child_dict)
            
        if self[i-1, j] == 4:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i-1, j]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i-1, j] = True
                return_dict['distance_scale'] = 1.0
                child_dict = self.__map_flow_from_cell((i-1,j), **kwargs)
                return_dict['next'].append(child_dict)
            
        if self[i-1, j+1] == 8:
            if kwargs.get('assume_valid_routing', True) or not self.visited[i-1, j+1]:
                if not kwargs.get('assume_valid_routing', True):
                    self.visited[i-1, j+1] = True
                return_dict['distance_scale'] = 1.4142135623730951
                child_dict = self.__map_flow_from_cell((i-1,j+1), **kwargs)
                return_dict['next'].append(child_dict)
        
        if len(return_dict.get('next', 0)) == 0:
            return_dict.pop('next')
                   
        return return_dict
    
    def update_flow_codes_in_mask(self, *args, **kwargs):
        
        flooded_dem = args[0]
        mask = args[1]
        
        indexes = np.where(mask._griddata == 1)
        for (i,j) in zip(indexes[0],indexes[1]):
            flow_code = self.__flow_code_for_position(flooded_dem, i, j)
            self._griddata[i,j] = flow_code
    
    def divides_for_outlets(self, outlet1, outlet2):
        basin1 = BaseSpatialGrid()
        basin1._copy_info_from_grid(self, True)
        
        ind = self.get_indexes_of_upstream_cells_for_location(outlet1[0], outlet1[1])
        i,j = zip(*ind)
        basin1._griddata[i,j] = 1
        
        basin2 = BaseSpatialGrid()
        basin2._copy_info_from_grid(self, True)
        
        ind = self.get_indexes_of_upstream_cells_for_location(outlet2[0], outlet2[1])
        i,j = zip(*ind)
        basin2._griddata[i,j] = 1
        
        from scipy.ndimage.morphology import binary_dilation as dilate
        
        basin3 = BaseSpatialGrid()
        basin3._copy_info_from_grid(self, True)
        basin3._griddata = basin1._griddata + dilate(basin2._griddata)
        
        basin4 = BaseSpatialGrid()
        basin4._copy_info_from_grid(self, True)
        basin4._griddata = basin2._griddata + dilate(basin1._griddata)
        
        (i, j) = np.where(basin3._griddata == 2)
        ind1 = zip(i,j)
        
        (i,j) = np.where(basin4._griddata == 2)
        ind2 = zip(i,j)
        
        return ind1, ind2
    
    def locations_of_paired_hollows(self, outlet1, outlet2, area, Ao=1E5):
        
        def map_down_to_hollow(ind, fd):
            (i,j) = ind
            (next_i, next_j, _) = fd.get_flow_to_cell(i,j)
            if((area[i,j] < Ao) and (area[next_i, next_j] >= Ao)):
                return next_i, next_j
            else:
                return map_down_to_hollow((next_i, next_j), fd)
            
        ind1, ind2 = self.divides_for_outlets(outlet1, outlet2)
        h1 = set()
        for ind in ind1:
            h1.add(map_down_to_hollow(ind, self))
        h2 = set()
        for ind in ind2:
            h2.add(map_down_to_hollow(ind, self))
        
        return tuple(h1), tuple(h2)
        
        
    def get_indexes_of_upstream_cells(self, i, j):
        
        i,j = self.__get_flow_from_cell(i,j)
        return zip(i,j)
    
    def get_indexes_of_upstream_cells_for_location(self, x, y):
        
        v = ((x, y), )
        ((i, j),) = self._xy_to_rowscols(v)
        return self.get_indexes_of_upstream_cells(i, j)

    def search_down_flow_direction_with_length(self, start, search_length = np.inf):
        l = list()
        length_list = list()
        ((row,col),) = self._xy_to_rowscols((start,))
        length = 0.0
        while not (row == 0 or row == self._georef_info.ny-1 or col == 0 or col == self._georef_info.nx - 1 or \
                   (row,col) in l) and (length < search_length):
            l.append((row, col))
            length_list.append(length)
            row,col,inBounds = self.get_flow_to_cell(row, col)
            length += self._georef_info.dx * (1.0 if (row == col) else 1.414)
            if not inBounds:
                break
        return tuple(l), tuple(length_list)

    def search_down_flow_direction(self, start, search_length = np.inf):
        return self.search_down_flow_direction_with_length(start, search_length=search_length)[0]
    
    def search_down_flow_direction_from_rowscols_location(self, start, return_rowscols = False, search_length = np.inf):
        
        rcs = self.search_down_flow_direction(self._rowscols_to_xy(((start),))[0], search_length = search_length)
        return rcs if return_rowscols else self._rowscols_to_xy(rcs)

    def search_down_flow_direction_from_xy_location(self, start, return_rowscols=False, search_length=np.inf):

        rcs = self.search_down_flow_direction(start, search_length=search_length)
        return rcs if return_rowscols else self._rowscols_to_xy(rcs)

    def convert_rivertools_directions_to_arc(self):
        # Function to convert river tools flow directions to arcGisFlowDirections
        # ArcGIS convention
        # |i-1,j-1  i-1,j  i-1,j+1|  |32 64 128|
        # |i,j-1     i,j     i,j+1|  |16  X  1 |
        # |i+1,j-1  i+1,j  i+1,j+1|  |8   4  2 |
        # In river tools convention the 1 is in the top right
                
        convertedFlowDir = int(np.log2(self._griddata))
        convertedFlowDir -= 1
        convertedFlowDir[convertedFlowDir == -1] = 7
        convertedFlowDir = 2**convertedFlowDir
        convertedFlowDir[self._griddata == self.noData] = self.noData
        self._griddata = convertedFlowDir
    
    def pixel_scale(self, dtype = np.float32):
        
        dim = np.ones_like(self._griddata, dtype = dtype)
        dim[np.where((self._griddata == 32) | (self._griddata == 128) | (self._griddata == 2) | (self._griddata == 8))] = 1.41421356
        return dim
    
    def map_values_to_recursive_list(self, outlet, **kwargs):
        
        v = (outlet, )
        (ij_outlet, ) = self._xy_to_rowscols(v)
        if kwargs.get('assume_valid_routing', True) is False:
            self.visited = np.zeros_like(self._griddata).astype(bool)
        return self.__map_flow_from_cell((ij_outlet,), **kwargs)


    def bounds_of_basin_for_outlet(self, outlet):
        
        (lat, lon) = outlet
        indexes = self.get_indexes_of_upstream_cells_for_location(lon, lat)
        xy = self._rowscols_to_xy(indexes) 
        bounds = [[None, None], [None, None]]
        for (lon, lat) in xy:
            if lon > bounds[0][1] or bounds[0][1] is None:
                bounds[0][1] = lon
            if lon < bounds[0][0] or bounds[0][0] is None:
                bounds[0][0] = lon
            if lat > bounds[1][1] or bounds[1][1] is None:
                bounds[1][1] = lat
            if lat < bounds[1][0] or bounds[1][0] is None:
                bounds[1][0] = lat
        return ((bounds[0][0], bounds[0][1]), (bounds[1][0], bounds[1][1]))
    
    def divides(self):
        
        divides = BaseSpatialGrid()
        divides._copy_info_from_grid(self, True)
        
        ip1j = self._griddata[2:-1,1:-2]
        im1j = self._griddata[0:-3,1:-2]
        ijp1 = self._griddata[1:-2,2:-1]
        ijm1 = self._griddata[1:-2,0:-3]
        
        ip1jp1 = self._griddata[2:-1, 2:-1]
        ip1jm1 = self._griddata[2:-1, 0:-3]
        im1jp1 = self._griddata[0:-3, 2:-1]
        im1jm1 = self._griddata[0:-3, 0:-3]
        
        divides._griddata[1:-2,1:-2] = ((ijm1 != 1) & (im1jm1 != 2) & (im1j != 4) & (im1jp1 != 8) & (ijp1 != 16) & (ip1jp1 != 32) & (ip1j != 64) & (ip1jm1 != 128))
        return divides
    
    def paired_divides(self, mask = None):
        
        divides = self.divides()
        if mask is not None:
            i = np.where((divides._griddata == 1) & (mask._griddata > 0))
        else:
            i = np.where(divides._griddata == 1)

        rcs = zip(i[0].tolist(),i[1].tolist())
        
        pairs = list()
        for (i,j) in rcs:
            
            (i_next, j_next, good) = self.get_flow_to_cell(i,j)
            
            if good and i_next is not None and j_next is not None:
                i_next = 2*i - i_next
                j_next = 2*j - j_next
                
                if i_next >= 0 and i_next < self._georef_info.ny and j_next >= 0 and j_next < self._georef_info.nx:
                    (xy, ) = self._rowscols_to_xy(((i,j),))
                    (xy_next, ) = self._rowscols_to_xy(((i_next, j_next),))
                    pairs.append( (xy, xy_next))
                
        return tuple(pairs)
                
class Elevation(CalculationMixin, BaseSpatialGrid):

    def findDEMedge(self):
        #Pad the data so that we can take a windowed max
        original = np.zeros((self._georef_info.ny, self._georef_info.nx))
        i = np.where(np.isnan(self._griddata) | (self._griddata == 0))
        original[i] = 1
        
        from scipy.ndimage.morphology import binary_dilation as dilation
        from scipy.ndimage import generate_binary_structure as generate_binary_structure
        structure = generate_binary_structure(2, 2)
        dilated = dilation(original, structure = structure, border_value = 1).astype(int)
        edges = dilated - original
        
        return np.where(edges)
    
    def outlets_at_coastlines(self):
        
        import scipy.ndimage.morphology as morph
        
        mask1 = BaseSpatialGrid()
        mask2 = BaseSpatialGrid()
        mask3 = BaseSpatialGrid()
        mask1._copy_info_from_grid(self, True)
        mask1.dtype = np.uint8
        mask1._griddata = mask1._griddata.astype(np.uint8)
        mask1._griddata[(self._griddata > 0)] = 1
        mask2._copy_info_from_grid(mask1, False)
        mask3._copy_info_from_grid(mask1, False)
        mask2._griddata = morph.binary_dilation(mask2._griddata, iterations = 3).astype(uint8)
        mask3._griddata = morph.binary_erosion(mask3._griddata, iterations = 3).astype(uint8)
        mask3._griddata = (mask3._griddata.astype(float) + mask2._griddata.astype(float)).astype(uint8)
        i = np.where(mask3._griddata == 1)
        rc = zip(i[0].tolist(),i[1].tolist())
        return self._rowscols_to_xy(rc)
    
    def track_flow_downhill(self, starting_point, maximum_pit_depth = 20):
        
        adjust = [(-1, -1), (0, -1), (1,-1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
        
        (ij, ) = self._xy_to_rowscols((starting_point,))
        visited = (ij, )
        ij_a = [(ij[0] + x[0], ij[1] + x[1]) for x in adjust]
        e_a = np.array([self[i[0], i[1]] for i in ij_a])
        new = True
        while None not in e_a and np.sum(np.isnan(e_a)) == 0 and new:
            sls = np.argsort(e_a)
            new = False
            for sl in sls:
                this_ij = (ij[0] + adjust[sl][0], ij[1] + adjust[sl][1])
                if (this_ij not in visited) and ( (e_a[sl] - self[ij[0], ij[1]]) <= maximum_pit_depth):
                    visited += (this_ij, )
                    ij = this_ij
                    new = True
                    break
            if new:
                ij_a = [(ij[0] + x[0], ij[1] + x[1]) for x in adjust]
                e_a = np.array([self[i[0], i[1]] for i in ij_a]) 
        xy = self._rowscols_to_xy(visited)
        this_xy = xy[0]
        l = tuple()
        e = tuple()
        last_l = 0
        for txy in xy:
            (ij, ) = self._xy_to_rowscols((txy,))
            e += (self[ij[0], ij[1]], )
            this_dl = np.sqrt(np.power(txy[0] - this_xy[0], 2) + np.power(txy[1] - this_xy[1], 2))
            last_l = last_l + this_dl
            l += (last_l, )
            this_xy = txy
        return xy, l, e

class Gradient(BaseSpatialGrid):
    
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
        outgrid = Gradient()
        outgrid._copy_info_from_grid(self, True)
        outgrid._gy = self.average_over_distance(distance, grid = self._gy)
        outgrid._gx = self.average_over_distance(distance, grid = self._gx)
        return outgrid
    
    def save(self, filename):    
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', [self._gx, self._gy], self.dtype, filename, ['COMPRESS=LZW'], multiple_bands=True)
    
    @classmethod
    def load(cls, filename):
        
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.NAN
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
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'))
        
    @classmethod
    def load(cls, filename):
        
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.NAN
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
        self.elevation = Elevation.load(filename)

    def plot_orientations(self, *args, **kwargs):
        
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
        if self.elevation is None:
            raise AttributeError('No elevation data! Use load_elevation first')
        valid = np.ones_like(self.elevation._griddata)
        valid[np.isnan(self.elevation._griddata)] = np.nan
        valid[self.elevation._griddata <= 0] = np.nan
        return valid

class LocalRelief(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('elevation','pixel_radius'), '_create_from_elevation_and_radius'))
    
    def _create_from_elevation_and_radius(self, *args, **kwargs):
        elevation = kwargs['elevation']
        pixel_radius = kwargs['pixel_radius']
        self._copy_info_from_grid(elevation,True)
        
        import scipy.ndimage.filters as filters
        y,x = np.ogrid[-pixel_radius:pixel_radius, -pixel_radius:pixel_radius]
        footprint = x*x + y*y <= pixel_radius*pixel_radius
        maximum = np.zeros_like(elevation._griddata)
        minimum = np.zeros_like(elevation._griddata)
        filters.maximum_filter(elevation._griddata, footprint = footprint, output = maximum)
        filters.minimum_filter(elevation._griddata, footprint = footprint, output = minimum)
        self._griddata = maximum - minimum
        
class Mask(BaseSpatialGrid):
    
    dtype = uint8
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('flow_direction','outlets'), '_create_from_flow_direction_and_outlets'))
    
    def _create_from_flow_direction_and_outlets(self, *args, **kwargs):
        flow_direction = kwargs['flow_direction']
        outlets = kwargs['outlets']
        self._copy_info_from_grid(flow_direction, True)
        for (outlet_x, outlet_y) in outlets:
            indexes = flow_direction.get_indexes_of_upstream_cells_for_location(outlet_x, outlet_y)
            for index in indexes:
                self[index[0],index[1]] = 1

    def perform_opening(self, structure = None, iterations = 1):
        from scipy.ndimage.morphology import binary_opening
        self._griddata = binary_opening(self._griddata, iterations = iterations, structure = structure)

    def perform_erosion(self, structure = None, iterations = 1):
        from scipy.ndimage.morphology import binary_erosion
        self._griddata = binary_erosion(self._griddata, iterations = iterations, structure = structure)

class LogArea(BaseSpatialGrid):
    
    dtype = uint8
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('area',), '_create_from_area'))
    
    def _create_from_area(self, *args, **kwargs):
        area = kwargs['area']
        self._copy_info_from_grid(area, True)
        self._griddata = np.log10(area._griddata)
                        
class Hillshade(CalculationMixin, BaseSpatialGrid):
    
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
    
    def calcHillshade(self,az,elev):
        #Hillshade = calcHillshade(elevGrid,az,elev)
        #Esri calculation for generating a hillshade, elevGrid is expected to be a numpy array
    
        # Convert angular measurements to radians
                
        azRad, elevRad = (360 - az + 90)*np.pi/180, (90-elev)*np.pi/180
        Sx, Sy = self._calcFiniteSlopes(self._griddata, self._mean_pixel_dimension(), self._georef_info.nx, self._georef_info.ny)  # Calculate slope in X and Y directions
    
        AspectRad = np.arctan2(Sy, Sx) # Angle of aspect
        SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))  # magnitude of slope in radians
    
        self._griddata = 255.0 * ((np.cos(elevRad) * np.cos(SmagRad)) + (np.sin(elevRad)* np.sin(SmagRad) * np.cos(azRad - AspectRad)))

class GeographicHillshade(GeographicGridMixin, Hillshade):
    pass

class MaxSlope(CalculationMixin, BaseSpatialGrid):
    
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
 
        Sx, Sy = self._calcFiniteSlopes(self._griddata, self._mean_pixel_dimension(), self._georef_info.nx, self._georef_info.ny)  # Calculate slope in X and Y directions    
        self._griddata = np.sqrt(Sx**2 + Sy**2)
        
class GeographicMaxSlope(GeographicGridMixin, MaxSlope):
    
    pass

class Laplacian(CalculationMixin, BaseSpatialGrid):
    
    dtype = float64
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('elevation',), '_create_from_elevation'))
    
    def _create_from_elevation(self, *args, **kwargs):

        from scipy.ndimage import laplace
        elevation = kwargs['elevation']
        
        self._copy_info_from_grid(elevation)
        self._griddata = laplace(elevation._griddata) / np.power(elevation._georef_info.dx,2)
    

class GeographicLaplacian(GeographicGridMixin, Laplacian):
    
    pass

class PriorityQueueMixIn(object):
    
    aggradation_slope = 1E-12

    class priorityQueue:
        #Implements a priority queue using heapq. Python has a priority queue module built in, but it
        # is not stably sorted (meaning that two items who are tied in priority are treated arbitrarily, as opposed to being
        # returned on a first in first out basis). This circumvents that by keeping a count on the items inserted and using that
        # count as a secondary priority
    
        def __init__(self):
            # A counter and the number of items are stored separately to ensure that items remain stably sorted and to
            # keep track of the size of the queue (so that we can check if its empty, which will be useful will iterating
            # through the queue)
            self.__pq = []
            self.__counter = 0
            self.__nItems = 0
    
        def get(self):
            #Remove an item and its priority from the queue
            priority, count, item = heapq.heappop(self.__pq)
            self.__nItems -= 1
            return priority, item
    
        def put(self, priority, item):
            #Add an item to the priority queue
            self.__counter += 1
            self.__nItems += 1
            entry = [priority, self.__counter, item]
            heapq.heappush(self.__pq, entry)
    
        def isEmpty(self):
            return self.__nItems == 0

    def randomize_subbasins_with_mask(self, *args, **kwargs):
        mask = args[0]
        outlets = args[1]
        
        self.__flood(mask = mask, outlets = outlets, randomize = True)
        
    def _flood(self, *args, **kwargs):
        # dem is a numpy array of elevations to be flooded, aggInc is the minimum amount to increment elevations by moving upstream
        # use priority flood algorithm described in  Barnes et al., 2013
        # Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models
        # NOTE: They have another algorithm to make this more efficient, but to use that and a slope makes things more
        # complicated
        
        #Create a grid to keep track of which cells have been filled        
        closed = np.zeros_like(self._griddata)
        i = np.where(np.isnan(self._griddata))
        closed[i] = 1
        #using_mask = False
                
        if kwargs.get('mask') is not None:
            closed = (kwargs.get('mask')._griddata != 1).astype(int)
        if kwargs.get('outlets') is not None:
            v = self._xy_to_rowscols(kwargs.get('outlets'))
            #using_mask = True
            edgeRows = [a[0] for a in v]
            edgeCols = [a[1] for a in v]
        if kwargs.get('randomize') is True:
            self._randomize_grid_values(mask = kwargs['mask'])
        if kwargs.get('outlets') is None:
            #Add all the edge cells to the priority queue, mark those cells as draining (not closed)
            edgeRows, edgeCols = self.findDEMedge()
        
        should_randomize_priority_queue = False
        
        if kwargs.get('randomize') is not None:
            if kwargs.get('randomize') == True:
                should_randomize_priority_queue = True
                
        priority_queue = FilledElevation.priorityQueue() # priority queue to sort filling operation
        
        for i in range(len(edgeCols)):
            row, col = edgeRows[i], edgeCols[i]
    
            closed[row, col] = True
            priority_queue.put(self._griddata[row, col], (row, col)) # store the indices as a vector of row column, in the priority queue prioritized by the dem value
    
        if kwargs.get('binary_result') is True or kwargs.get('clip_to_fill') is True:
            visited = np.zeros_like(self._griddata, dtype = np.float)
                
        #While there is anything left in the priority queue, continue to fill holes
        while not priority_queue.isEmpty():
                        
            priority, (row, col) = priority_queue.get()
            
            if kwargs.get('binary_result') is True or kwargs.get('clip_to_fill') is True:
                visited[row, col] = 1.0
            
            elevation = self[row,col]
            
            neighborRows, neighborCols, dxMults = self._getNeighborIndices(row, col)
            
            dxs = self._georef_info.dx * dxMults
                
            #Look through the upstream neighbors
            for i in range(len(neighborCols)):
                should_fill = True
                
                if not closed[neighborRows[i], neighborCols[i]]:    
                    #If this was a hole (lower than the cell downstream), fill it
                    if self._griddata[neighborRows[i], neighborCols[i]] <= elevation:

                        if kwargs.get('maximum_pit_depth'):
                            fill_depth = elevation - self._griddata[neighborRows[i], neighborCols[i]]
                            if fill_depth > kwargs.get('maximum_pit_depth'):
                                should_fill = False
                        
                        if should_fill:
                            self._griddata[neighborRows[i], neighborCols[i]] = elevation + self.aggradation_slope*dxs[i]
    
                    closed[neighborRows[i], neighborCols[i]] = True
                    if should_randomize_priority_queue:
                        #if not using_mask or (using_mask and ~edge_of_mask_is_in_kernel):
                        
                        if should_fill:
                            priority_queue.put(np.random.rand(1)[0], [neighborRows[i], neighborCols[i]])
                    else:
                        if should_fill:
                            priority_queue.put(self._griddata[neighborRows[i], neighborCols[i]], [neighborRows[i], neighborCols[i]])
            
        if kwargs.get('binary_result'):
            self._griddata = visited
        if kwargs.get('clip_to_fill') is True:
            self._griddata[visited == 0] = np.NAN
            
        
            
class PriorityFillGrid(PriorityQueueMixIn, BaseSpatialGrid):
    
    dtype = uint8
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('outlets', 'mask'), '_create_from_outlets_and_mask'))
    
    def _create_from_outlets_and_mask(self, *args, **kwargs):
        
        mask = kwargs['mask']
        outlets = kwargs['outlets']
        kwargs['randomize'] = False
        kwargs['binary_result'] = True
        self._copy_info_from_grid(mask,True)
        rcs = self._xy_to_rowscols(outlets)
        for rc in rcs:
            self._griddata[rc[0],rc[1]] = 100
        i = np.where(np.isnan(mask._griddata))
        self._griddata[i] = 0
        i = np.where(mask._griddata != 1)
        self._griddata[i] = 0
        self._flood(*args, **kwargs)
                
class FilledElevation(PriorityQueueMixIn, Elevation):
    
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('elevation',), '_create_from_elevation'))
    
    def _create_from_elevation(self, *args, **kwargs):

        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation)
        self._flood(*args, **kwargs) 
                                
class Area(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('flow_direction',), '_create_from_flow_direction'),
                                   (('dx', 'grid'), '_create_from_grid'))

    def _create_from_flow_direction(self, *args, **kwargs):
        
        flow_dir = kwargs['flow_direction']
        self._copy_info_from_grid(flow_dir, True)
        #if flow_dir.__class__ == FlowDirectionD8:
        self.__calcD8Area(*args, **kwargs)
    
    def __calcD8Area(self, *args, **kwargs):

        # I am returning area and flowDirections but NOTE!I might be in danger the way I handle this...
        # Flow directions are indices of arrays which may be different sizes, this is not truly giving me flow directions
        # Because in the getNeighbor function I don't return values at edges.... might want to think about improving this,
        # although it should work in this application
        # Calculate the D8 drainage area for the numpy array representing a demMethods in elevGrid, assumes that elevGrid has already been filled
        # Assigns BCs to deal with edges - these should be changed to handle flow differently, currently it will not force flow through the edges
    
        flow_dir = kwargs['flow_direction']

        if kwargs.get('sorted_indexes') is not None:
            idcs = kwargs.get('sorted_indexes')
        else:
            idcs = flow_dir.sort()
            
        area = self._area_per_pixel(*args, **kwargs)  # area of a pixel
        
        [ind_i, ind_j] = np.unravel_index(idcs, flow_dir._griddata.shape)
        
        has_mask = kwargs.get('mask') is not None
        has_evaluate_at = kwargs.get('evaluate_at') is not None
        
        if has_evaluate_at:
            area[kwargs['evaluate_at']._griddata == 0] = 0
            
        
        for i, j in zip(ind_i, ind_j):  # Loop through all the data in sorted order    
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)
            if is_good and not has_mask:
                if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                    area[i_next, j_next] += area[i,j]
            elif is_good and kwargs['mask'][i, j] is not None:
                if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                    area[i_next, j_next] += area[i,j]
    
        self._griddata = area # Return non bc version of area
    
    def areas_greater_than(self, min_area):
        ij_cols = np.where(self._griddata >= min_area)
        ij = zip(ij_cols[0].tolist(), ij_cols[1].tolist())
        return self._rowscols_to_xy(ij)
    
    def areas_between(self, fd, min_area, max_area):
        ij_out = []
        ij_cols = np.where((self._griddata >= min_area) & (self._griddata <= max_area))
        for (i,j) in zip(ij_cols[0].tolist(), ij_cols[1].tolist()):
            options = fd.get_upstream_cell_indexes(i,j)
            max_a = None
            for (i_p, j_p) in options:
                if max_a is None or (self._griddata[i_p, j_p] > max_a):
                    max_a = self._griddata[i_p, j_p]
            if max_a < min_area:
                ij_out += [(i,j)]
        return self._rowscols_to_xy(ij_out)
        
    
class GeographicArea(GeographicGridMixin, Area):
    pass

class ValleyArea(Area):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('flow_direction', 'area', 'laplace', 'valley_laplace_value', 'min_area_value'), '_create_from_flow_direction_and_valley_slope_and_area'))

    def _create_from_flow_direction_and_valley_slope_and_area(self, *args, **kwargs):
        
        mask = Mask()
        mask._copy_info_from_grid(kwargs['laplace'], True)
        mask._griddata[kwargs['laplace']._griddata >= kwargs['valley_laplace_value'] ] = 1
        outlets = kwargs['area'].areas_greater_than(kwargs['min_area_value'])
        pfg = PriorityFillGrid(mask = mask, outlets = outlets)
        import scipy.ndimage.morphology as morph
        if kwargs.get('iterations') is None or kwargs.get('iterations') == 0:
            pass 
        else:
            iterations = kwargs['iterations']
            pfg._griddata = morph.binary_dilation(pfg._griddata, iterations = iterations)
            pfg._griddata = morph.binary_erosion(pfg._griddata, iterations = iterations)
        
        kwargs['evaluate_at'] = pfg    
        self._create_from_flow_direction(*args, **kwargs)
    
            
class GeographicValleyArea(GeographicGridMixin, ValleyArea):
    pass    
        
class MainstemValleyArea(Area):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('area', 'laplace', 'valley_laplace_value', 'min_area_value', 'flow_direction'), '_create_from_flow_direction_and_valley_slope_and_area'))

    def _create_from_flow_direction_and_valley_slope_and_area(self, *args, **kwargs):
        
        mask = Mask()
        mask._copy_info_from_grid(kwargs['laplace'], True)
        mask._griddata[kwargs['laplace']._griddata >= kwargs['valley_laplace_value'] ] = 1
        outlets = kwargs['area'].areas_greater_than(kwargs['min_area_value'])
        pfg = PriorityFillGrid(mask = mask, outlets = outlets)
        import scipy.ndimage.morphology as morph
        if kwargs.get('iterations') is None or kwargs.get('iterations') == 0:
            pass 
        else:
            iterations = kwargs['iterations']
            pfg._griddata = morph.binary_dilation(pfg._griddata, iterations = iterations)
            pfg._griddata = morph.binary_erosion(pfg._griddata, iterations = iterations)
        
        kwargs['evaluate_at'] = pfg    
        self._create_from_flow_direction(*args, **kwargs)
    
    def __calcD8Area(self, *args, **kwargs):
    
        flow_dir = kwargs['flow_direction']

        if kwargs.get('sorted_indexes') is not None:
            idcs = kwargs.get('sorted_indexes')
        else:
            idcs = flow_dir.sort()
            
        dA = self._area_per_pixel(*args, **kwargs)
        dx = self._mean_pixel_dimension(*args, **kwargs) * flow_dir.pixel_scale()
        
        length = np.zeros_like(dA)
        
        [ind_i, ind_j] = np.unravel_index(idcs, flow_dir._griddata.shape)
        
        has_mask = kwargs.get('mask') is not None
        has_evaluate_at = kwargs.get('evaluate_at') is not None
        
        if has_evaluate_at:
            dA[kwargs['evaluate_at']._griddata == 0] = 0
                    
        for i, j in zip(ind_i, ind_j):  # Loop through all the data in sorted order    
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)
            
            if is_good:
                next_l = length[i,j] + dx[i,j]
                if next_l > length[i_next, j_next]:
                    length[i_next, j_next] = next_l
                    if not has_mask:
                        if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                            self._griddata[i_next, j_next] = dA[i_next,j_next] + self._griddata[i,j]
                    elif kwargs['mask'][i, j] is not None:
                        if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                            self._griddata[i_next, j_next] = dA[i_next,j_next] + self._griddata[i,j]
                elif kwargs['area'] <= kwargs['min_area_value']:
                    if not has_mask:
                        if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                            self._griddata[i_next, j_next] += self._griddata[i,j]
                    elif kwargs['mask'][i, j] is not None:
                        if not has_evaluate_at or (has_evaluate_at and (kwargs['evaluate_at'][i,j] == 1)):
                            self._griddata[i_next, j_next] += self._griddata[i,j]
                            
class GeographicMainstemValleyArea(GeographicGridMixin, MainstemValleyArea):
    pass

class AlongFlowSmoothing(object):

    def _find_points_along_path(self, de, **kwargs):
        area = kwargs['area']
        flow_direction = kwargs['flow_direction']

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
                while (delta_e < vertical_interval) & (ups_i >= 0):
                    (ups_i, ups_j) = (upstream_i[ups_i, ups_j], upstream_j[ups_i, ups_j])
                    (ds_i, ds_j) = (downstream_i[ds_i, ds_j], downstream_j[ds_i, ds_j])
                    if (ups_i < 0) or (ds_i < 0):
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
                while (horizontal_distance < horizontal_interval) & (ups_i >= 0):
                    (ups_i, ups_j) = (upstream_i[ups_i, ups_j], upstream_j[ups_i, ups_j])
                    (ds_i, ds_j) = (downstream_i[ds_i, ds_j], downstream_j[ds_i, ds_j])
                    if (ups_i < 0) or (ds_i < 0):
                        return None
                    horizontal_distance += (1.0 if ((ds_i == ret[0][0]) or (ds_j == ret[0][1])) else 1.414) * de[
                        ds_i, ds_j] + (1.0 if ((ups_i == ret[-1][0]) or (ups_j == ret[-1][1])) else 1.414) * de[
                                               ups_i, ups_j]
                    ret = [(ds_i, ds_j)] + ret + [(ups_i, ups_j)]
                return ret
            return find_points_along_path
    
class KsFromChiWithSmoothing(BaseSpatialGrid, AlongFlowSmoothing):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('elevation', 'area', 'flow_direction', 'theta', 'vertical_interval'), '_create_from_elevation_area_flow_direction'),
                                   (('elevation', 'area', 'flow_direction', 'theta', 'horizontal_interval'), '_create_from_elevation_area_flow_direction'),
                                   )
            
    def _create_from_elevation_area_flow_direction(self, *args, **kwargs):
        
        elevation = kwargs['elevation']
        area = kwargs['area']
        theta = kwargs['theta']
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

        def calc_ks(i,j):
            points = find_points_along_path(i, j)     
            if points is not None:
                pts = list(zip(*(points)))
                points = np.array(pts).astype(int)
                adjustment = np.ones((len(points[0])))
                i = np.where((points[0,1:-1] != points[0,2:]) & (points[1,1:-1] != points[1,2:]))
                adjustment[i[0]+1] += 0.414
                area_profile = area._griddata[points[0],points[1]]
                elevation_profile = elevation._griddata[points[0], points[1]]
                de_profile = de[points[0], points[1]]
                chi_profile = np.zeros_like(elevation_profile)
                chi_profile[1:] = np.cumsum(0.25*(np.power(area_profile[1:], -theta) + np.power(area_profile[0:-1], -theta)) * (de_profile[1:] + de_profile[0:-1])*adjustment[1:]) 
                X = np.array(chi_profile)
                y = np.array(elevation_profile) - elevation_profile[0]
                model = sm.OLS(y, X)
                res = model.fit()
                SS = res.ssr
                return res.params[0], SS / float(len(chi_profile)), SS, res.rsquared, points, res.pvalues[0], len(chi_profile)
            
            else:
                
                return np.nan, np.nan, np.nan, np.nan, [[],[]], np.nan, 0 
        
        i = np.where((area._griddata != 0) & ~np.isnan(area._griddata) & ~np.isnan(elevation._griddata))
        ij = list(zip(i[0],i[1]))
        totalnumber = len(ij)
        counter = 0.0    
        next_readout = 0.1    
        sys.stdout.write('Percent completion...')    
        sys.stdout.flush()    
        for (i,j) in ij:    
            self._griddata[i,j], self._mse[i,j], self._ss[i,j], self._r2[i,j], pts, self._pval[i,j], self._n_regression[i,j] = calc_ks(i,j)    
            self._n[pts[0], pts[1]] += 1    
            counter += 1.0 / totalnumber    
            if counter > next_readout:
                sys.stdout.write(str(int(next_readout*100)) + "...")
                sys.stdout.flush()
                next_readout += 0.1
                
        sys.stdout.write('Percent completion...')
        sys.stdout.flush()            
        sys.stdout.write('100')
        sys.stdout.flush()
                                            
    def save(self, filename):
        
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', [self._griddata, self._n, self._mse, self._ss, self._r2, self._pval, self._n_regression], self.dtype, filename, ['COMPRESS=LZW', 'BIGTIFF=YES'], multiple_bands=True)
            
    @classmethod
    def load(cls, filename):
        
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.NAN
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
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction', 'vertical_interval', 'min_area'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation', 'area', 'flow_direction', 'horizontal_interval', 'min_area'),
                                    '_create_from_elevation_area_flow_direction'),
                                   )

    def _create_from_elevation_area_flow_direction(self, *args, **kwargs):

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

        def calc_theta(i, j):

            points = find_points_along_path(i, j)
            if points is not None:
                pts = list(zip(*(points)))
                points = np.array(pts).astype(int)
                adjustment = np.ones((len(points[0])))
                i = np.where((points[0, 1:-1] != points[0, 2:]) & (points[1, 1:-1] != points[1, 2:]))
                adjustment[i[0] + 1] += 0.414
                area_profile = area._griddata[points[0], points[1]]
                elevation_profile = elevation._griddata[points[0], points[1]]
                de_profile = de[points[0], points[1]]

                def r2_for_theta(theta):
                    if theta > 10 or theta < -10:
                        return np.inf
                    chi_profile = np.zeros_like(elevation_profile)
                    chi_profile[1:] = np.cumsum(
                        0.25 * (np.power(area_profile[1:], -theta) + np.power(area_profile[0:-1], -theta)) * (
                                de_profile[1:] + de_profile[0:-1]) * adjustment[1:])
                    X = np.array(chi_profile)
                    y = np.array(elevation_profile) - elevation_profile[0]
                    model = sm.OLS(y, X)
                    res = model.fit()
                    return res.ssr
                from scipy.optimize import fmin
                try:
                    (theta_bf, funval, iter, funcalls, warnflag) = fmin(r2_for_theta, np.array([0.5]), (), 1E-5, 1E-5, 100, 200, True, False, 0, None)
                    chi_profile = np.zeros_like(elevation_profile)
                    chi_profile[1:] = np.cumsum(
                        0.25 * (np.power(area_profile[1:], -theta_bf) + np.power(area_profile[0:-1], -theta_bf)) * (
                                de_profile[1:] + de_profile[0:-1]) * adjustment[1:])
                    X = np.array(chi_profile)
                    y = np.array(elevation_profile) - elevation_profile[0]
                    model = sm.OLS(y, X)
                    res = model.fit()
                    SS = res.ssr
                    return theta_bf, SS / float(len(chi_profile)), SS, res.rsquared, points, res.pvalues[0], len(chi_profile)
                except:
                    return np.nan, np.nan, np.nan, np.nan, [[], []], np.nan, 0
            else:

                return np.nan, np.nan, np.nan, np.nan, [[],[]], np.nan, 0

        i = np.where((area._griddata != 0) & ~np.isnan(area._griddata) & ~np.isnan(elevation._griddata) & (area._griddata > min_area))
        ij = list(zip(i[0], i[1]))
        totalnumber = len(ij)
        counter = 0.0
        next_readout = 0.1
        sys.stdout.write('Percent completion...')
        sys.stdout.flush()
        for (i, j) in ij:
            self._griddata[i, j], self._mse[i, j], self._ss[i, j], self._r2[i, j], pts, self._pval[i, j], \
            self._n_regression[i, j] = calc_theta(i, j)
            self._n[pts[0], pts[1]] += 1
            counter += 1.0 / totalnumber
            if counter > next_readout:
                sys.stdout.write(str(int(next_readout * 100)) + "...")
                sys.stdout.flush()
                next_readout += 0.1

        sys.stdout.write('Percent completion...')
        sys.stdout.flush()
        sys.stdout.write('100')
        sys.stdout.flush()

    def save(self, filename):

        self._create_gdal_representation_from_array(self._georef_info, 'GTiff',
                                                    [self._griddata, self._n, self._mse, self._ss, self._r2, self._pval,
                                                     self._n_regression], self.dtype, filename,
                                                    ['COMPRESS=LZW', 'BIGTIFF=YES'], multiple_bands=True)

    @classmethod
    def load(cls, filename):

        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.NAN
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
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_flow_direction'),
                                   )

    scale_factor = 1.0

    def points_along_path(self, points, i, j):
        return points

    def calc_channel_slope(self, i, j, elevation, de, find_points_along_path):

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

    def _create_from_elevation_flow_direction(self, *args, **kwargs):

        elevation = kwargs['elevation']
        if kwargs.get('horizontal_interval') is not None:
            kwargs['horizontal_interval'] = kwargs['horizontal_interval']*self.scale_factor
        if kwargs.get('vertical_interval') is not None:
            kwargs['vertical_interval'] = kwargs['vertical_interval']*self.scale_factor

        de = elevation._mean_pixel_dimension()

        self._copy_info_from_grid(elevation)
        self._griddata = np.zeros_like(elevation._griddata)
        self._griddata[:] = np.nan

        find_points_along_path = self._find_points_along_path(de, **kwargs)

        i = np.where(~np.isnan(elevation._griddata))
        ij = list(zip(i[0], i[1]))
        totalnumber = len(ij)
        counter = 0.0
        next_readout = 0.1
        sys.stdout.write('Percent completion...')
        sys.stdout.flush()
        for (i, j) in ij:
            self._griddata[i, j] = self.calc_channel_slope(i, j, elevation, de, find_points_along_path)
            counter += 1.0 / totalnumber
            if counter > next_readout:
                sys.stdout.write(str(int(next_readout * 100)) + "...")
                sys.stdout.flush()
                next_readout += 0.1
        sys.stdout.write('Percent completion...')
        sys.stdout.flush()
        sys.stdout.write('100')
        sys.stdout.flush()

class ChannelDownSlopeWithSmoothing(ChannelSlopeWithSmoothing):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_flow_direction'),
                                   )
    scale_factor = 2.0

    def points_along_path(self, points, i, j):
        position_of_center = list([ind[0] for ind in zip(range(len(points)),points) if (ind[1][0] == i and ind[1][1] == j)])[0]
        return points[0:position_of_center] if len(points[0:position_of_center]) > 0 else None

class ChannelUpSlopeWithSmoothing(ChannelSlopeWithSmoothing):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',), '_create'),
                                   (('ai_ascii_filename', 'EPSGprojectionCode'), '_read_ai'),
                                   (('gdal_filename',), '_read_gdal'),
                                   (('elevation', 'area', 'flow_direction',  'vertical_interval'),
                                    '_create_from_elevation_area_flow_direction'),
                                   (('elevation',  'area', 'flow_direction',  'horizontal_interval'),
                                    '_create_from_elevation_flow_direction'),
                                   )
    scale_factor = 2.0

    def points_along_path(self, points, i, j):
        position_of_center = list([ind[0] for ind in zip(range(len(points)),points) if (ind[1][0] == i and ind[1][1] == j)])[0]
        return points[position_of_center:] if len(points[position_of_center:]) > 0 else None

class GeographicKsFromChiWithSmoothing(GeographicGridMixin, KsFromChiWithSmoothing):
    pass

class GeographicThetaFromChiWithSmoothing(GeographicGridMixin, KsFromChiWithSmoothing):
    pass

class MultiscaleCurvatureValleyWidth(BaseSpatialGrid):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('elevation', 'area', 'area_cutoff', 'max_width', 'min_width'), '_create_from_inputs'),)

    class Utilities(object):

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

            from numpy.fft import fft2 as fft2
            from numpy.fft import ifft2 as ifft2
            from numpy import flipud as flipud
            from numpy import fliplr as fliplr
            from numpy.fft import ifftshift

            Xt = X
            Yt = Y

            FZ = fft2(fliplr(fliplr(Z._griddata)))
            FX1 = fft2(np.power(Xt, 2))
            FX2 = fft2(np.power(Yt, 2))
            FX3 = fft2(Xt * Yt)
            FX4 = fft2(Xt)
            FX5 = fft2(Yt)
            FX6 = fft2(K)

            g = np.real(ifftshift(ifft2(FZ * FX1))) - np.sum(np.power(Xt, 2)) * Z._griddata
            h = np.real(ifftshift(ifft2(FZ * FX2))) - np.sum(np.power(Yt, 2)) * Z._griddata
            i = np.real(ifftshift(ifft2(FZ * FX3))) - np.sum(Xt * Yt) * Z._griddata
            j = np.real(ifftshift(ifft2(FZ * FX4))) - np.sum(Xt) * Z._griddata
            k = np.real(ifftshift(ifft2(FZ * FX5))) - np.sum(Yt) * Z._griddata
            if not fix_center:
                l = np.real(ifftshift(ifft2(FZ * FX6))) - np.sum(K) * Z._griddata
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
        def _Cmin_for_scale(cls, Z, de, fix_center=False):

            X, Y, N, K = cls._create_kernel(Z, de)
            H = cls._calc_inv_G_for_kernel(X, Y, N, fix_center)
            g, h, i, j, k, l = cls._convolve(X, Y, Z, K, fix_center)
            Cmin = cls._Cmin(H, g, h, i, j, k, l, fix_center)
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

        [X,Y] = np.meshgrid(np.arange(xllcenter+dx/2, xllcenter + (nx-0.5)*dx, dx),np.arange(yllcenter+dx/2, yllcenter + (ny-0.5)*dx, dx))
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
        
        (area_cutoff, max_width, min_width, normalize, fix_center, use_dask) = (kwargs['area_cutoff'], kwargs['max_width'], kwargs['min_width'], kwargs.get('normalize', False), kwargs.get('fix_center', False), kwargs.get('use_dask', False))
        
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
            wrapper = partial(self.Utilities._Cmin_for_scale, Z, fix_center = fix_center)
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
                minC, _ = self.Utilities._Cmin_for_scale(Z, scale, fix_center)
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
        
        import copy
        unpadded = copy.deepcopy(self)
        unpadded._griddata =  unpadded._griddata[ypadding:-ypadding,xpadding:-xpadding]
        unpadded._minC = unpadded._minC[ypadding:-ypadding,xpadding:-xpadding]
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
        
        def get_band(gdal_dataset, band_number):
            band = gdal_dataset.GetRasterBand(band_number)
            nodata = band.GetNoDataValue()
            grid = band.ReadAsArray().astype(cls.dtype)
            if nodata is not None:
                nodata_elements = np.where(grid == nodata)
                from numpy import uint8
                if cls.dtype is not uint8:
                    grid[nodata_elements] = np.NAN
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
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('elevation', 'outlets'), '_create_from_elevation_outlets'),
                                   (('elevation', ), '_create_from_elevation'),)
    
    def _create_from_elevation(self, *args, **kwargs):
        elevation = Elevation()
        elevation._copy_info_from_grid(kwargs['elevation'])
        if kwargs.get('mask') is not None:
            i = np.where(kwargs['mask']._griddata != 1)
            elevation._griddata[i] = np.NaN
        kwargs['elevation'] = elevation
        i = np.where(~np.isnan(elevation._griddata))
        rcs = zip(i[0].tolist(),i[1].tolist())
        kwargs['outlets'] = elevation._rowscols_to_xy(rcs)
        self._create_from_elevation_outlets(*args, **kwargs)
        
    def _create_from_elevation_outlets(self, *args, **kwargs):
        
        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation, True)
        adjust = [(-1, -1), (0, -1), (1,-1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
        
        counter = 1
        
        for outlet in kwargs['outlets']:
            if kwargs.get('display_output') is True:
                print('Evaluating outlet ' + str(counter) + ' / ' + str(len(kwargs['outlets'])))
            counter += 1
            (ij, ) = elevation._xy_to_rowscols((outlet,))
            visited = (ij, )
            ij_a = [(ij[0] + x[0], ij[1] + x[1]) for x in adjust]
            e_a = np.array([elevation[i[0], i[1]] if elevation[i[0], i[1]] is not None else np.NaN for i in ij_a])
            new = True
            area = 0
            while None not in e_a and np.sum(np.isnan(e_a)) == 0 and new:
                area += self._area_per_pixel()[ij[0], ij[1]]
                if kwargs.get('terminations_only') is not True:
                    self._griddata[ij[0], ij[1]] = area
                sls = np.argsort(e_a)
                new = False
                for sl in sls:
                    this_ij = (ij[0] + adjust[sl][0], ij[1] + adjust[sl][1])
                    if this_ij not in visited:
                        visited += (this_ij, )
                        ij = this_ij
                        new = True
                        break
                if new:
                    ij_a = [(ij[0] + x[0], ij[1] + x[1]) for x in adjust]
                    e_a = np.array([elevation[i[0], i[1]] if elevation[i[0], i[1]] is not None else np.NaN for i in ij_a]) 
            if kwargs.get('terminations_only') is True:
                self._griddata[ij[0], ij[1]] = area  

    def _area_per_pixel(self, *args, **kwargs):
        return self._georef_info.dx**2 * np.ones((self._georef_info.ny, self._georef_info.nx))

    def _mean_pixel_dimension(self, *args, **kwargs):
        return self._georef_info.dx * np.ones_like(kwargs['elevation']._griddata, self.dtype)

class GeographicDiscreteFlowAccumulation(GeographicGridMixin, DiscreteFlowAccumulation):    
    pass
    
class FlowLength(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('flow_direction',), '_create_from_flow_direction_and_sorted_indexes'))
        
    def _create_from_flow_direction_and_sorted_indexes(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        self.__calculate_flow_length(*args, **kwargs)
        
    def __calculate_flow_length(self, *args, **kwargs):
        
        flow_dir = kwargs['flow_direction']
        if kwargs.get('sorted_indexes') is not None:
            idcs = kwargs.get('sorted_indexes')
        else:
            idcs = flow_dir.sort()
            
        self.__flow_directions = np.zeros_like(flow_dir._griddata, np.uint8)

        [ind_i, ind_j] = np.unravel_index(idcs, flow_dir._griddata.shape)            
        dx = self._mean_pixel_dimension(*args, **kwargs) * flow_dir.pixel_scale()

        for i, j in zip(ind_i, ind_j):  # Loop through all the data in sorted order  
              
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)
                
            if is_good:
                next_l = self._griddata[i,j] + dx[i,j]
                if next_l > self._griddata[i_next, j_next]:
                    self._griddata[i_next, j_next] = next_l
                    self.__flow_directions[i_next, j_next] = self.__flow_direction_for_length((i,j), (i_next,j_next))

    def __flow_direction_for_length(self, from_index, to_index):
        (i_from, j_from) = from_index
        (i_to, j_to) = to_index
        
        d_i = i_to - i_from
        d_j = j_to - j_from
                
        if d_i == 0:
            if d_j == 1:
                return 16
            if d_j == 0:
                return 0
            if d_j == -1:
                return 1
        elif d_i == 1:
            if d_j == 1:
                return 8
            if d_j == 0:
                return 4
            if d_j == -1:
                return 2
        elif d_i == -1:
            if d_j == 1:
                return 32
            if d_j == 0:
                return 64
            if d_j == -1:
                return 128
        
        return None
    
    def points_with_length(self, length, fd):
        
        tolerance = np.nanmin(self._mean_pixel_dimension())
        max_length = length
        min_length = length - tolerance*2.0
        
        indexes_of_locations = list()
        ind = np.where(np.logical_and((self._griddata >= min_length), (self._griddata <= max_length)))
        for (this_row,this_col) in zip(ind[0],ind[1]):
            (next_row, next_col, is_good) = fd.get_flow_to_cell(this_row, this_col)
            if next_row is not None and next_col is not None and is_good:
                this_value = self[this_row, this_col]
                next_value = self[next_row, next_col]
                if this_value <= length and next_value >= length:
                    indexes_of_locations.append((this_row, this_col))
        
        return self._rowscols_to_xy(indexes_of_locations)
    
    def indexes_along_flow_path_from_outlet(self, outlet):
        ij = self._xy_to_rowscols((outlet,))[0]
        return self.__get_upstream_indexes((ij,))
    
    def locations_along_flow_path_from_outlet(self, outlet):
        return self._rowscols_to_xy(self.indexes_along_flow_path_from_outlet(outlet))
    
    def locations_along_flow_path_from_outlets(self, outlets):
        points = tuple()
        for outlet in outlets:
            points += tuple(self.locations_along_flow_path_from_outlet(outlet))
        return points
        
    def __get_upstream_indexes(self, index):
        
        (i,j) = index[0]
        indexes = list(index)
        
        
        flow_code = self.__flow_directions[i,j]
        if flow_code == 1:
            indexes += self.__get_upstream_indexes(((i,j+1),))
        elif flow_code == 2:
            indexes += self.__get_upstream_indexes(((i-1,j+1),))
        elif flow_code == 4:
            indexes += self.__get_upstream_indexes(((i-1,j),))
        elif flow_code == 8:
            indexes += self.__get_upstream_indexes(((i-1,j-1),))
        elif flow_code == 16:
            indexes += self.__get_upstream_indexes(((i,j-1),))
        elif flow_code == 32:
            indexes += self.__get_upstream_indexes(((i+1,j-1),))
        elif flow_code == 64:
            indexes += self.__get_upstream_indexes(((i+1,j),))
        elif flow_code == 128:
            indexes += self.__get_upstream_indexes(((i+1, j+1),))
        
        return indexes
    
    def __map_flow_from_cell(self, index, **kwargs):
        
        (i,j) = index[0]
        
        return_dict = dict()
        flow_code = self.__flow_directions[i,j]
        
        return_dict['index'] = index[0]
        return_dict['distance'] = self._griddata[i,j]
        for arg in kwargs:
            return_dict[arg] = kwargs[arg][i,j]
        
        return_dict['next'] = []
        return_dict['distance_scale'] = 1.0
                
        if flow_code == 1:
            return_dict['distance_scale'] = 1.0
            child_dict = self.__map_flow_from_cell(((i,j+1),), **kwargs)
            return_dict['next'].append(child_dict)
        
        if flow_code == 2:
            return_dict['distance_scale'] = 1.4142135623730951
            child_dict = self.__map_flow_from_cell(((i-1,j+1),), **kwargs)
            return_dict['next'].append(child_dict)

        if flow_code == 4:
            return_dict['distance_scale'] = 1.0
            child_dict = self.__map_flow_from_cell(((i-1,j),), **kwargs)
            return_dict['next'].append(child_dict)
        
        if flow_code == 8:
            return_dict['distance_scale'] = 1.4142135623730951
            child_dict = self.__map_flow_from_cell(((i-1,j-1),), **kwargs)
            return_dict['next'].append(child_dict)
                
        if flow_code == 16:
            return_dict['distance_scale'] = 1.0
            child_dict = self.__map_flow_from_cell(((i,j-1),), **kwargs)
            return_dict['next'].append(child_dict)

        if flow_code == 32:
            return_dict['distance_scale'] = 1.4142135623730951
            child_dict = self.__map_flow_from_cell(((i+1,j-1),), **kwargs)
            return_dict['next'].append(child_dict)
            
        if flow_code == 64:
            return_dict['distance_scale'] = 1.0
            child_dict = self.__map_flow_from_cell(((i+1,j),), **kwargs)
            return_dict['next'].append(child_dict)
            
        if flow_code == 128:
            return_dict['distance_scale'] = 1.4142135623730951
            child_dict = self.__map_flow_from_cell(((i+1,j+1),), **kwargs)
            return_dict['next'].append(child_dict)
        
        if len(return_dict.get('next', 0)) == 0:
            return_dict.pop('next')
                   
        return return_dict

    def is_along_flow_length(self, from_index, to_index):
        
        (i_to, j_to) = to_index
        flow_code = self.__flow_directions[i_to, j_to]
        
        return flow_code == self.__flow_direction_for_length(from_index, to_index)
     
    def clip_to_extent(self, extent):
        
        return_grid = super(FlowLength, self).clip_to_extent(extent)
        
        import copy
        return_grid.__flow_directions = copy.deepcopy(self.__flow_directions)
        lower_left = (extent[0]-self._georef_info.dx, extent[2]-self._georef_info.dx)
        upper_right = (extent[1]+self._georef_info.dx, extent[3]+self._georef_info.dx)
        idx = self._xy_to_rowscols((lower_left, upper_right))
        return_grid.__flow_directions = return_grid._griddata[idx[1][0]+1:idx[0][0]+1,idx[0][1]:idx[1][1]]
                
        return return_grid

    def map_values_to_recursive_list(self, outlet, **kwargs):
        
        v = (outlet, )
        (ij_outlet, ) = self._xy_to_rowscols(v)
        return self.__map_flow_from_cell((ij_outlet,), **kwargs)
    
    def save(self, filename):
        
        super(FlowLength, self).save(filename)
        flow_dir_name = filename + "_directions"
        self._create_gdal_representation_from_array(self._georef_info, 'GTiff', self.__flow_directions, np.uint8, flow_dir_name, ['COMPRESS=LZW'])
        
    @classmethod
    def load(cls, filename):
        
        return_object_bsp = BaseSpatialGrid.load(filename)
        return_object = cls()
        return_object._georef_info = return_object_bsp._georef_info
        return_object._griddata = return_object_bsp._griddata
        flow_dir_filename = filename + "_directions"        
        gdal_dataset = gdal.Open(flow_dir_filename)
        band = gdal_dataset.GetRasterBand(1)
        return_object.__flow_directions = band.ReadAsArray().astype(np.uint8)
                
        gdal_file = None
        return return_object
        

class GeographicFlowLength(GeographicGridMixin, FlowLength):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('flow_direction',), '_create_from_flow_direction_and_sorted_indexes'))
    
class MaxFlowLengthTrackingMixin(object):
    
    def _calculate_by_tracking_down_max_flow_length(self, *args, **kwargs):
        flow_dir = kwargs['flow_direction']
        flow_length = kwargs['flow_length']
        if kwargs.get('sorted_indexes') is not None:
            idcs = kwargs.get('sorted_indexes')
        else:
            idcs = flow_dir.sort()
                
        [ind_i, ind_j] = np.unravel_index(idcs, flow_dir._griddata.shape)
        
        for i, j in zip(ind_i, ind_j):  # Loop through all the data in sorted order    
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)  
            if is_good and flow_length.is_along_flow_length((i,j), (i_next, j_next)):
                self._calculate_grid_value((i,j), (i_next, j_next), *args, **kwargs)        
        
class Relief(BaseSpatialGrid, MaxFlowLengthTrackingMixin):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('flow_direction', 'elevation', 'flow_length'), '_create_from_flow_direction_sorted_indexes_and_elevation'))   
    
    def _create_from_flow_direction_sorted_indexes_and_elevation(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['elevation'])
        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)
        self._griddata = self._griddata - kwargs['elevation']._griddata
    
    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        (i,j) = pos
        (i_next, j_next) = next_pos

        if kwargs.get('area') is not None and kwargs.get('Ao') is not None:
            if kwargs.get('area')[i,j] <= kwargs['Ao']:
                self._griddata[i_next, j_next] = self._griddata[i_next, j_next]
                return
            
        self._griddata[i_next, j_next] = self._griddata[i,j]
        
class ScaledRelief(Relief):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('flow_direction', 'flooded_dem', 'elevation', 'flow_length', 'Ao', 'theta'), '_create_scaled_from_flow_direction_flooded_dem_and_elevation'),
                           (('flow_direction', 'elevation', 'flow_length', 'Ao', 'theta'), '_create_scaled_from_flow_direction_sorted_indexes_and_elevation'))
    
    def _create_scaled_from_flow_direction_flooded_dem_and_elevation(self, *args, **kwargs):
        self._create_from_flow_direction_flooded_dem_and_elevation(*args, **kwargs)
        self._griddata = self._griddata * (np.power(kwargs['Ao'],kwargs['theta']))
    
    def _create_scaled_from_flow_direction_sorted_indexes_and_elevation(self, *args, **kwargs):
        self._create_from_flow_direction_sorted_indexes_and_elevation(*args, **kwargs)
        self._griddata = self._griddata * (np.power(kwargs['Ao'],kwargs['theta']))

class Ksi(BaseSpatialGrid, MaxFlowLengthTrackingMixin):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('area','flow_direction','theta', 'Ao', 'flow_length'), '_create_from_inputs'))
    
    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        area_grid = kwargs['area']._griddata - kwargs['Ao']
        area_grid[area_grid <= 0] = np.nan
        self._griddata = np.zeros_like(area_grid)
        i = np.where(area_grid > 0)
        self._griddata[i] = ( (kwargs['Ao'] / area_grid[i]) ** kwargs['theta']) * self._mean_pixel_dimension(*args, **kwargs)[i] * kwargs['flow_direction'].pixel_scale()[i]
        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)
        
    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        (i,j) = pos
        (i_next, j_next) = next_pos
        self._griddata[i_next, j_next] += self._griddata[i,j]
        

class GeographicKsi(GeographicGridMixin, Ksi):
    pass

class RestoredElevation(BaseSpatialGrid):
    
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
        v = self._xy_to_rowscols(outlets)
        pixel_dimension = self._mean_pixel_dimension(*args, **kwargs)
        
        
        mask = Mask(flow_direction = flow_direction, outlets = outlets)
        external_divides = copy.deepcopy(mask)
        
        if kwargs.get('fix_external_outlets') == True:
            from scipy.ndimage import binary_erosion as erosion
            external_divides._griddata = external_divides._griddata - erosion(external_divides._griddata, np.ones((3,3)), border_value = 0).astype(int)
        else:
            external_divides._griddata = np.zeros_like(external_divides._griddata, int)
        if randomize:
            filled = FilledElevation(elevation = self, mask = mask, randomize = randomize, outlets = outlets)
            flow_direction.update_flow_codes_in_mask(filled, mask)
            i = np.where(mask._griddata == 1)
            self._griddata[i] = filled._griddata[i]
            area = self.__recalculate_area(area, flow_direction, pixel_dimension, outlets)
        for i in range(iterations):
            last_grid = self._griddata.copy()
            print('Iteration {0}'.format(i))
            print('Filling outlets.')
            divides = self.__fill_outlets(area, flow_direction, pixel_dimension, v, ks, theta)
            print('Migrating divides')
            flow_direction = self.__migrate_divides(flow_direction, divides, external_divides)
            area = self.__recalculate_area(area, flow_direction, pixel_dimension, outlets)
            change_in_elevation = np.mean((self._griddata - last_grid)**2)
            print('Change: {0}'.format(change_in_elevation))
    
    def __is_unit_area(self, idx, flow_direction):
        
        (i,j) = idx
        
        if flow_direction[i,j+1] == 16 or flow_direction[i+1, j+1] == 32 or \
                flow_direction[i+1, j] == 64 or flow_direction[i+1, j-1] == 128 or \
                flow_direction[i, j-1] == 1 or flow_direction[i-1, j-1] == 2 or \
                flow_direction[i-1, j] == 4 or flow_direction[i-1, j+1] == 8:
            return False
        else:
            return True
        
    def __migrate_divides(self, *args, **kwargs):
        
        flow_direction = args[0]
        divides = args[1]
        external_divides = args[2]
        migrated = 0

        for (i, j) in divides:
            if external_divides[i,j] != 1:
                try:
                    if self[i,j+1] > self[i,j] and external_divides[i,j+1] != 1:
                        flow_direction[i,j+1] = 16
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i+1, j+1] > self[i,j] and external_divides[i+1,j+1] != 1:
                        flow_direction[i+1, j+1] = 32
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i+1, j] > self[i,j] and external_divides[i+1,j] != 1:
                        flow_direction[i+1, j] = 64
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i+1, j-1] > self[i,j] and external_divides[i+1,j-1] != 1:
                        flow_direction[i+1, j-1] = 128
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i, j-1] > self[i,j] and external_divides[i,j-1] != 1:
                        flow_direction[i, j-1] = 1
                        migrated += 1
                except:
                    pass
                
                try:        
                    if self[i-1,j-1] > self[i,j] and external_divides[i-1,j-1] != 1:
                        flow_direction[i-1, j-1] = 2
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i-1, j] > self[i,j] and external_divides[i-1,j] != 1:
                        flow_direction[i-1, j] = 4
                        migrated += 1
                except:
                    pass
                
                try:
                    if self[i-1, j+1] > self[i,j] and external_divides[i-1,j+1] != 1:
                        flow_direction[i-1, j+1] = 8
                        migrated += 1
                except:
                    pass
                
        print("Migrated " + str(migrated) + " divides.")
        return flow_direction
    
    def __fill_outlets(self, area, flow_direction, pixel_dimension, outlet_indexes, ks, theta):
        
        divides = list()
        for ind in outlet_indexes:
            divides += self.__fill_upstream_points(ind, ks, theta, area, flow_direction, pixel_dimension)
        return divides
        
    def __fill_upstream_points(self, ind, ks, theta, area, flow_direction, pixel_dimension):
        
        (i, j) = ind
        elevation_for_cell = self[i,j]
        is_divide = True
        divides = list()
        
        try:
            if flow_direction[i, j+1] == 16:
                self._griddata[i,j+1] = elevation_for_cell + pixel_dimension[i,j]*1.00*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i,j+1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:
            if flow_direction[i+1, j+1] == 32:
                self._griddata[i+1,j+1] = elevation_for_cell + pixel_dimension[i,j]*1.414*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i+1,j+1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:                
            if flow_direction[i+1, j] == 64:
                self._griddata[i+1,j] = elevation_for_cell + pixel_dimension[i,j]*1.00*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i+1,j), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:
            if flow_direction[i+1, j-1] == 128:
                self._griddata[i+1,j-1] = elevation_for_cell + pixel_dimension[i,j]*1.414*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i+1,j-1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:        
            if flow_direction[i, j-1] == 1:
                self._griddata[i,j-1] = elevation_for_cell + pixel_dimension[i,j]*1.00*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i,j-1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try: 
            if flow_direction[i-1, j-1] == 2:
                self._griddata[i-1,j-1] = elevation_for_cell + pixel_dimension[i,j]*1.414*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i-1,j-1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:    
            if flow_direction[i-1, j] == 4:
                self._griddata[i-1,j] = elevation_for_cell + pixel_dimension[i,j]*1.00*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i-1,j), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        try:    
            if flow_direction[i-1, j+1] == 8:
                self._griddata[i-1,j+1] = elevation_for_cell + pixel_dimension[i,j]*1.414*ks*area[i,j]**(-theta)
                divides += self.__fill_upstream_points((i-1,j+1), ks, theta, area, flow_direction, pixel_dimension)
                is_divide = False
        except:
            pass
        
        if is_divide:
            return [(i,j)]
        else:
            return divides
    
    def __recalculate_area(self, area, flow_direction, pixel_dimension, outlets):
        outlet_indexes = self._xy_to_rowscols(outlets)
        for outlet in outlet_indexes:
            self.__recurse_area(outlet, area, flow_direction, pixel_dimension)
        return area
        
    def __recurse_area(self, ind, area, flow_direction, pixel_dimension):
        
        (i, j) = ind
        
        this_area = pixel_dimension[i,j]**2
        
        try:
            if flow_direction[i, j+1] == 16:
                ind = (i, j+1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:
            if flow_direction[i+1, j+1] == 32:
                ind = (i+1,j+1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:                
            if flow_direction[i+1, j] == 64:
                ind = (i+1,j)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:
            if flow_direction[i+1, j-1] == 128:
                ind = (i+1, j-1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:        
            if flow_direction[i, j-1] == 1:
                ind = (i,j-1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try: 
            if flow_direction[i-1, j-1] == 2:
                ind = (i-1, j-1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:    
            if flow_direction[i-1, j] == 4:
                ind = (i-1, j)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass
        
        try:    
            if flow_direction[i-1, j+1] == 8:
                ind = (i-1, j+1)
                this_area += self.__recurse_area(ind, area, flow_direction, pixel_dimension)
        except:
            pass    

        area[i,j] = this_area
        return this_area
   
class GeographicRestoredElevation(GeographicGridMixin, RestoredElevation):
    pass

class ChiScaledRelief(BaseSpatialGrid):
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('elevation', 'flow_direction', 'theta', 'Ao', 'outlets'), '_create_from_inputs'),
                           (('elevation', 'flow_direction', 'flow_length', 'theta', 'Ao', 'basin_length'), '_create_from_basin_length'))
    
    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['elevation'], True)
        outlet_indexes = self._xy_to_rowscols(kwargs['outlets'])
        elevation = kwargs['elevation']
        scale = np.power(kwargs['Ao'],kwargs['theta'])
        outlet_number = 1
        for outlet_index in outlet_indexes:
            indexes = kwargs['flow_direction'].get_indexes_of_upstream_cells(outlet_index[0], outlet_index[1])
            elevation_of_outlet = kwargs['elevation'][outlet_index[0],outlet_index[1]]
            for index in indexes:
                self[index[0],index[1]] = (elevation[index[0], index[1]] - elevation_of_outlet) * scale
            if kwargs.get('output_flag', False):
                print('Outlet ' + str(outlet_number) + '/' + str(len(outlet_indexes)) + ' completed.')
            outlet_number = outlet_number + 1
            
    def _create_from_basin_length(self, *args, **kwargs):
        kwargs['outlets'] = kwargs['flow_length'].points_with_length(kwargs['basin_length'],kwargs['flow_direction'])
        kwargs['output_flag'] = True
        return self._create_from_inputs(*args, **kwargs)       

class Chi(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('area','flow_direction', 'theta', 'Ao', 'outlets'), '_create_from_inputs'),
                           (('area', 'flow_direction', 'flow_length', 'theta', 'Ao', 'basin_length'), '_create_from_basin_length'))
    
    
    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        self.__calculate_chi(*args, **kwargs)
    
    def _create_from_basin_length(self, *args, **kwargs):
        kwargs['outlets'] = kwargs['flow_length'].points_with_length(kwargs['basin_length'],kwargs['flow_direction'])
        kwargs['output_flag'] = True
        return self._create_from_inputs(*args, **kwargs)
            
    def __calculate_chi(self, *args, **kwargs):
        pixel_dimension = self._mean_pixel_dimension(*args, **kwargs)
        area = kwargs['area']
        flow_direction = kwargs['flow_direction']
        Ao = kwargs['Ao']
        theta = kwargs['theta']
        outlet_indexes = self._xy_to_rowscols(kwargs['outlets'])
        self.__occupied = np.zeros_like(area._griddata, dtype=int)
        outlet_number = 1
        lmax = kwargs.get('maximum_length')
        
        for outlet in outlet_indexes:
            self.__recurse_chi(outlet, area, flow_direction, pixel_dimension, Ao, theta, 0, 1, kwargs.get('mask'), 0.0, lmax = lmax)
            if kwargs.get('output_flag', False):
                print('Outlet ' + str(outlet_number) + '/' + str(len(outlet_indexes)) + ' completed.')
            outlet_number += 1
        self.__occupied = None
        try:
            self._griddata = self._griddata * kwargs['mask']._griddata
        except:
            pass
        
    def __recurse_chi(self, v, area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = None):
        
        (i,j) = v
        
        if mask is not None:
            if mask[i,j] == 0:
                if area[i,j] > 1E4:
                    return
           
        if self.__occupied[i,j] == 1:
            return
        
        dl = pixel_dimension[i,j] * scale
        l = l + dl
        
        if lmax is not None and l >= lmax:
            return
        
        self._griddata[i,j] = chi + (Ao / area[i,j])**theta * dl
        self.__occupied[i,j] = 1;
        chi = self._griddata[i,j]
        

        try:
            if flow_direction[i, j+1] == 16:
                scale = 1.0
                self.__recurse_chi((i, j+1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:
            if flow_direction[i+1, j+1] == 32:
                scale = 1.414
                self.__recurse_chi((i+1, j+1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:                
            if flow_direction[i+1, j] == 64:
                scale = 1.0
                self.__recurse_chi((i+1, j), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:
            if flow_direction[i+1, j-1] == 128:
                scale = 1.414
                self.__recurse_chi((i+1, j-1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:        
            if flow_direction[i, j-1] == 1:
                scale = 1.0
                self.__recurse_chi((i, j-1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try: 
            if flow_direction[i-1, j-1] == 2:
                scale = 1.414
                self.__recurse_chi((i-1, j-1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:    
            if flow_direction[i-1, j] == 4:
                scale = 1.0
                self.__recurse_chi((i-1, j), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass
        
        try:    
            if flow_direction[i-1, j+1] == 8:
                scale = 1.414
                self.__recurse_chi((i-1, j+1), area, flow_direction, pixel_dimension, Ao, theta, chi, scale, mask, l, lmax = lmax)
        except:
            pass  

class GeographicChi(GeographicGridMixin, Chi):
    pass


class CrossDivideDChi(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('chi','flow_direction'), '_create_from_inputs'))
    
    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        paired_divides = kwargs['flow_direction'].paired_divides()
        divide_columns = zip(*paired_divides)
        from_point = divide_columns[0]
        to_point = divide_columns[1]
        
        from_idx = self._xy_to_rowscols(from_point)
        to_idx = self._xy_to_rowscols(to_point)
        
        from_idx_column = zip(*from_idx)
        to_idx_column = zip(*to_idx)
        
        minchi = np.minimum(kwargs['chi']._griddata[to_idx_column[0],to_idx_column[1]], kwargs['chi']._griddata[from_idx_column[0],from_idx_column[1]])
        maxchi = np.maximum(kwargs['chi']._griddata[to_idx_column[0],to_idx_column[1]], kwargs['chi']._griddata[from_idx_column[0],from_idx_column[1]])
        
        self._griddata[from_idx_column[0],from_idx_column[1]] = (maxchi - minchi) * (minchi != 0).astype(float)
                
class NormalizedCrossDivideDChi(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('chi','flow_direction'), '_create_from_inputs'))
    
    def _create_from_inputs(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        paired_divides = kwargs['flow_direction'].paired_divides(mask = kwargs['chi'])
        divide_columns = zip(*paired_divides)
        from_point = divide_columns[0]
        to_point = divide_columns[1]
        
        from_idx = self._xy_to_rowscols(from_point)
        to_idx = self._xy_to_rowscols(to_point)
        
        from_idx_column = zip(*from_idx)
        to_idx_column = zip(*to_idx)
        
        minchi = np.minimum(kwargs['chi']._griddata[to_idx_column[0],to_idx_column[1]], kwargs['chi']._griddata[from_idx_column[0],from_idx_column[1]])
        maxchi = np.maximum(kwargs['chi']._griddata[to_idx_column[0],to_idx_column[1]], kwargs['chi']._griddata[from_idx_column[0],from_idx_column[1]])
        
        rangechi = (maxchi - minchi) * (minchi != 0).astype(float)
        meanchi = (maxchi + minchi) / 2
        
        self._griddata[from_idx_column[0],from_idx_column[1]] = rangechi / meanchi

         
class Deflection(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('elevation', 'D', 'rho_m', 'rho_c', 'g'), '_deflect'))
    
    def _deflect(self, *args, **kwargs):
        if kwargs.get('restored_elevation'):
            self._griddata = self.__deflect_elevation(kwargs['elevation'], **kwargs) - self.__deflect_elevation(kwargs['restored_elevation'], **kwargs)
        else:
            self._griddata = self.__deflect_elevation(kwargs['elevation'], **kwargs)
        
    def __deflect_elevation(self, *args, **kwargs):
        elevation = args[0]
        self._copy_info_from_grid(elevation, False)
        kwargs['flow_direction'] = elevation
        dx = np.mean(self._mean_pixel_dimension(*args, **kwargs))
        gr = -self._griddata * kwargs['rho_c'] * kwargs['g']
        grf = np.fft.fft2(gr)
        Rx = np.arange(((self._georef_info.nx-1)*dx)/2,0,-dx)
        Ry = np.arange(((self._georef_info.ny-1)*dx)/2,0,-dx)
        lx = len(Rx)
        ly = len(Ry)
        wn_x = np.array((1 / dx) * np.arange(0,lx,1) / lx)
        wn_y = np.array((1 / dx) * np.arange(0,ly,1) / ly)
        wn_x = np.concatenate((wn_x,np.flipud(wn_x)))
        wn_y = np.concatenate((wn_y,np.flipud(wn_y)))
        [WN_x, WN_y] = np.meshgrid(wn_x, wn_y)
        
        R = 1 / ((kwargs['rho_m']-kwargs['rho_c'])*kwargs['g'] + (2*3.141)**4 * kwargs['D']*(WN_x**2 + WN_y**2)**2)
        
        grfR = grf * R
        
        return np.real(np.fft.ifft2(grfR))
        
class GeographicDeflection(GeographicGridMixin, Deflection):
    pass

                
class ChannelSlope(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('flow_direction','elevation'), '_create_from_flow_direction_and_elevation'))
    
    def _create_from_flow_direction_and_elevation(self, *args, **kwargs):
        
        flow_dir = kwargs['flow_direction']
        elevation = kwargs['elevation']
        self._copy_info_from_grid(flow_dir, True)
        if flow_dir.__class__ == FlowDirectionD8:
            self.__calc_D8_slope(flow_dir, elevation)            

    def __calc_D8_slope(self, fd, dem):

        dx = self._georef_info.dx
        idcs = fd.sort() # Get the sorted indices of the array in reverse order (e.g. largest first)
    
        for idx in idcs:  # Loop through all the data in sorted order
            [i, j] = np.unravel_index(idx, fd._griddata.shape)  # Get the row/column indices
    
            i_next, j_next, is_good = fd.get_flow_to_cell(i,j)
            if is_good:
                dist = np.sqrt((i-i_next)**2 + (j-j_next)**2)*dx
                self._griddata[i,j] = (dem[i,j] - dem[i_next, j_next]) / dist

    
def mosaicFolder(folderPath, fileSuffix, outfile):
    #This runs the gdal utility gdal_merge on the command line, is mainly here so that I don't have to continue looking
    #up how to actually accomplish this
    #use os to get files, existing gdal functions to merge them
    files = glob.glob1(folderPath, '*'+fileSuffix) # get all the files in the path
    ## Need to use existing gdal_merge.py.... need to check on interpolation options.
    # could also use gdal_warp i think. Should consider whether I care about the possibility of using the c
    #utility...
    argument = ['python', 'gdal_merge.py', '-o', outfile]
    for file in files:
        argument.append(os.path.join(folderPath, file))
    # sys.argv = argv
    subprocess.call(argument)

def plot(*args, **kwargs):
    
    grid1 = args[0]._griddata
    grid2 = args[1]._griddata

    interactive = kwargs.pop('interactive', True)
    if kwargs.get('xlabel') is None:
        xlabel = 'Grid 1'
    
    if kwargs.get('ylabel') is None:
        ylabel = 'Grid 2'
    
    if kwargs.get('symbol') is None:
        symbol = 'k.'
    
    if kwargs.get('indexes') is None:    
        valid_indexes = np.where( ~np.isnan(grid1+grid2) )
    else:
        valid_indexes = kwargs.get('indexes')
        
    if kwargs.get('decimation_factor') is not None:
        valid_indexes = (valid_indexes[0][::kwargs.get('decimation_factor')],valid_indexes[1][::kwargs.get('decimation_factor')])
    
    plt.plot(grid1[valid_indexes], grid2[valid_indexes], symbol)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if interactive:
        plt.ion()
        plt.show(block=False)
    else:
        plt.ioff()
        plt.show(block=True)

    ax = plt.gca()
    
    return ax


