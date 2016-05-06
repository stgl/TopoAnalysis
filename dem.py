#Functions for demMethods calculations typical of geomorphology
#Developed by Sam Johnstone, January 2015, samuelj@stanford.edu , johnsamstone@gmail.com
#Not all functions are benchmarked
#If you somehow have this, and have never spoken to me, please reach out
#- I'd be curious to hear what you are doing (maybe I can even help with something!)

from osgeo import gdal # Used to load gis in files
import osr # Used with gis files
from osgeo import ogr
import os # Used to join file paths, iteract with os in other ways..
import glob  # Used for finding files that I want to mosaic (by allowing wildcard searches of the filesystem)
import heapq # Used for constructing priority queue, which is used for filling dems
import numpy as np # Used for tons o stuff, keeping most data stored as numpy arrays
import subprocess # Used to run gdal_merge.py from the command line
import Error
from numpy import uint8
from matplotlib.mlab import dist
from matplotlib import pyplot as plt
from operator import pos


        
class GDALMixin(object):
    
    def _get_projection_from_EPSG_projection_code(self, EPSGprojectionCode):
        # Get raster projection
        srs = osr.SpatialReference()
        return srs.ImportFromEPSG(EPSGprojectionCode).ExportToWkt()

    def _get_gdal_type_for_numpy_type(self, numpy_type):
    
        from numpy import float64, uint8, uint16, int16, uint32, int32, float32, complex64
        from gdal import GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32, GDT_Float32, GDT_Float64, GDT_CFloat64, GDT_Unknown
        
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
        from gdal import GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32, GDT_Float32, GDT_Float64, GDT_CFloat64
                    
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

    def _create_gdal_representation_from_array(self, georef_info, GDALDRIVERNAME, array_data, dtype, outfile_path='name'):
        #A function to write the data in the numpy array arrayData into a georeferenced dataset of type
        #  specified by GDALDRIVERNAME, a string, options here: http://www.gdal.org/formats_list.html
        #  This is accomplished by copying the georeferencing information from an existing GDAL dataset,
        #  provided by createDataSetFromArray
    
        #Initialize new data
        drvr = gdal.GetDriverByName(GDALDRIVERNAME)  #  Get the desired driver
        outRaster = drvr.Create(outfile_path, georef_info.nx, georef_info.ny, 1 , self._get_gdal_type_for_numpy_type(dtype))  # Open the file
    
        #Write geographic information
        outRaster.SetGeoTransform(georef_info.geoTransform)  # Steal the coordinate system from the old dataset
        if georef_info.projection != 0:
            outRaster.SetProjection(georef_info.projection)   # Steal the Projections from the old dataset
    
        #Write the array
        outRaster.GetRasterBand(1).WriteArray(array_data)   # Writes my array to the raster
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
            for _ in xrange(5):
                line = ascii_file.readline()
                (key, value) = (line.split()[1], int(line.split()[-1]))
                georef_data[key.lower()] = value

        required_values = ('ncols', 'nrows', 'cellsize','nodata_value')
        
        if len(set(required_values).subtract(set(georef_data.keys()))) != 0:
            raise Error.InputError('A/I ASCII grid error','The following properties are missing: ' + set(required_values).subtract(set(georef_data.keys())))
        
        if georef_data.get('xllcorner') is None and georef_data.get('xllcenter') is None:
            raise Error.InputError('A/I ASCII grid error','Neither XLLCorner nor XLLCenter is present.')
        
        if georef_data.get('yllcorner') is None and georef_data.get('yllcenter') is None:
            raise Error.InputError('A/I ASCII grid error','Neither YLLCorner nor YLLCenter is present.')
        
        dx = georef_data.get('cellsize')
        nx = georef_data.get('ncols')
        ny = georef_data.get('nrows')
        
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
    
        header = "ncols     %s\n" % georef_info.ncols
        header += "nrows    %s\n" % georef_info.nrows
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
    
        data = dataOut.ReadAsArray(dtype = self.dtype)
        dataOut = None
        
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

        #Returns a grid of the area of each pixel within the domain specified by the gdalDataset
    
        re = 6371.0 * 1000.0 #radius of the earth in meters
    
        dLong = self._georef_info.geoTransform[1]
        dLat = self._georef_info.geoTransform[5]
        
        lats = self._georef_info.yllcenter + np.float64(range(self._georef_info.ny))*dLat
        longs = self._georef_info.xllcenter + np.float64(range(self._georef_info.nx))*dLong
            
        #Get the size of pixels in the lat and long direction

        [LONG, LAT] = np.meshgrid(longs, lats)
        LAT1 = np.radians(LAT - dLat / 2.0)
        LAT2 = np.radians(LAT + dLat / 2.0)
        LONG1 = np.radians(LONG - dLong / 2.0)
        LONG2 = np.radians(LONG + dLong / 2.0)
        
        initialAreas = np.abs((re**2)*(LONG1 - LONG2)*(np.sin(LAT2) - np.sin(LAT1)))
        
        return initialAreas
    
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
        Zbc = self.assignBCs(grid)
    
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
    
    from numpy import float64
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('nx', 'ny', 'dx'), '_create_random_grid'),)
    dtype = float64
    

    _georef_info = Georef_info()
        
    def __init__(self, *args, **kwargs):
        
        from numpy import zeros
        super(BaseSpatialGrid,self).__init__()

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

        if i < 0 or j < 0 or i >= self._georef_info.ny or j >= self._georef_info.nx:
            return None
        
        return self._griddata[i, j]
        
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
            self._griddata = np.zeros(shape = (self._geref_info.ny,self._georef_info.nx), dtype = self.dtype)
        else:
            self._griddata = kwargs.get('grid')

    def _create_random_grid(self, *args, **kwargs):
        self._georef_info.dx = kwargs['dx']
        self._georef_info.nx = kwargs['nx']
        self._georef_info.ny = kwargs['ny']
        self._georef_info.xllcenter = 0
        self._georef_info.yllcenter = 0
        self._griddata = np.random.rand(self._georef_info.ny, self._georef_info.nx)
        
    def _read_ai(self, *args, **kwargs):
        
        self._georef_info = []
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
        outRows = rowKernel + row
        outCols = colKernel + col
    
        #Determine which indices are out of bounds
        inBounds = (outRows >= 0)*(outRows < self._georef_info.ny)*(outCols >= 0)*(outCols < self._georef_info.nx)
        return (outRows[inBounds], outCols[inBounds], dxMults[inBounds])
    
    def _xy_to_rowscols(self, v):
        l = list()
        for (x,y) in v:
            col = round((x-self._georef_info.xllcenter)/self._georef_info.dx)
            row = self._georef_info.ny - round((y-self._georef_info.yllcenter)/self._georef_info.dx)
            if col > self._georef_info.nx or row > self._georef_info.ny:
                l.append((None, None))
            else:
                l.append((row,col))
        return tuple(l)
    
    def _rowscols_to_xy(self, l):
        from numpy import float64
        v = list()
        for(row,col) in l:
            x = float64(col)*self._georef_info.dx + self._georef_info.xllcenter
            y = (float64(self._georef_info.ny) - float64(row))*self._georef_info.dx + self._georef_info.yllcenter
            v.append((x,y))
        return tuple(v)
    
    def sort(self, reverse=True):
        if not self._sorted:
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

    def plot(self, **kwargs):
        
        extent = [self._georef_info.xllcenter, self._georef_info.xllcenter+(self._georef_info.nx-1)*self._georef_info.dx, self._georef_info.yllcenter, self._georef_info.yllcenter+(self._georef_info.ny-1)*self._georef_info.dx]
        plt.imshow(self._griddata, extent = extent, **kwargs)
        plt.ion()
        plt.show()
    
    def find_nearest_cell_with_value(self, index, value, pixel_radius):

        (i, j) = index
        nearest = tuple()
        minimum_difference = np.nan
        for y in range(int(i-pixel_radius),int(i+pixel_radius)):
            for x in range(int(j-pixel_radius),int(j+pixel_radius)):
                if np.sqrt( (y-i)**2 + (x-j)**2) <= pixel_radius:
                    value_difference = np.abs(value - self._griddata[y,x])
                    if np.isnan(minimum_difference) or minimum_difference > value_difference:
                        minimum_difference = value_difference
                        nearest = (y,x)
        
        (y,x) = nearest
        return y, x
    
    def save(self, filename):
        
        self._create_gdal_representation_from_array(self._georef_info, 'ENVI', self._griddata, self.dtype, filename)
    
    @classmethod
    def load(cls, filename):
        
        return_object = cls()
        
        gdal_dataset = gdal.Open(filename)
        band = gdal_dataset.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        return_object._griddata = band.ReadAsArray().astype(cls.dtype)
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
        return_object._georef_info.yllcenter = return_object._georef_info.geoTransform[3]-(return_object._georef_info.dx*(return_object._georef_info.ny-0.5))
        return_object._georef_info.nx = nx
        return_object._georef_info.ny = ny
        
            
        gdal_file = None
        return return_object
    
        return cls._    
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
        for i in range(self._georef_info.ny-1):
            for j in range(self._georef_info.nx-1):
                flow_code = self.__flow_code_for_position(flooded_dem, i, j)
                self._griddata[i,j] = flow_code
    
        self._sort_indexes = flooded_dem.sort(reverse = False)
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

    def __get_flow_from_cell(self, i, j):
        
        i_source = [i]
        j_source = [j]
        
        if self[i, j+1] == 16:
            [i_append, j_append] = self.__get_flow_from_cell(i,j+1)
            i_source = i_source + i_append
            j_source = j_source + j_append
        
        if self[i+1, j+1] == 32:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j+1)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i+1, j] == 64:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j)
            i_source = i_source + i_append
            j_source = j_source + j_append
        
        if self[i+1, j-1] == 128:
            [i_append, j_append] = self.__get_flow_from_cell(i+1,j-1)
            i_source = i_source + i_append
            j_source = j_source + j_append
                
        if self[i, j-1] == 1:
            [i_append, j_append] = self.__get_flow_from_cell(i,j-1)
            i_source = i_source + i_append
            j_source = j_source + j_append
         
        if self[i-1, j-1] == 2:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j-1)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i-1, j] == 4:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j)
            i_source = i_source + i_append
            j_source = j_source + j_append
            
        if self[i-1, j+1] == 8:
            [i_append, j_append] = self.__get_flow_from_cell(i-1,j+1)
            i_source = i_source + i_append
            j_source = j_source + j_append
               
        return i_source, j_source
    
    
    def get_indexes_of_upstream_cells(self, i, j):
        
        i,j = self.__get_flow_from_cell(i,j)
        return zip(i,j)
    
    def get_indexes_of_upstream_cells_for_location(self, x, y):
        
        ((i, j)) = self._xy_to_rowscols((x,y),)
        return self.get_indexes_of_upstream_cells(i, j)
        
    def searchDownFlowDirection(self, start):
    
        l = list()
        (row, col) = self._xy_to_rowscols(start)
        l.append((row,col))
        #So long as we are not at the edge of the demMethods
        while not (row == 0 or row == self._georef_info.ny-1 or col == 0 or col == self._georef_info.nx - 1):
            # Get the neighbors in the flow direction corresponding order - note, I might be in danger the way I handle this...
            # Flow directions are indices of arrays which may be different sizes, this is not truly giving me flow directions
            # Because in the getNeighbor function I don't return values at edges.... might want to think about improving this,
            # although it should work in this application, and damn if it isn't clean
            row,col,inBounds = self.__getFlowToCell(row, col) # Find the indices of the cell we drain too, only need first two inputs
            if not inBounds:
                break    
            l.append((row, col))
            
        return tuple()

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
        
            
class Elevation(CalculationMixin, BaseSpatialGrid):
    
    def findDEMedge(self):
        # Function to find the cells at the edge of a dem. Dem is a ny x nx array, but may be largely padded
        # by nans. Determines where the edge of the real data is. Does this by finding the maximum value within a 3x3 kernel,
        # if that is not zero, but the original data at the corresponding location is, then that is an edge cell
            
        #Pad the data so that we can take a windowed max
        padded = np.zeros((self._georef_info.ny+2, self._georef_info.nx+2))
        padded[1:-1, 1:-1] = self._griddata
        padded[padded == 0] = np.nan
    
        # windowMax = np.zeros_like(padded)
        borderCells = np.zeros_like(padded)
    
        #Iterate through all the data, find the max in a 3 x 3 kernel
        for i in range(self._georef_info.ny):
            for j in range(self._georef_info.nx):
                # windowMax[i+1, j+1] = np.nanmax(padded[i:i+3, j:j+3])
                borderCells[i+1, j+1] = np.any(np.isnan(padded[i:i+3, j:j+3]))*(~np.isnan(padded[i+1,j+1]))
    
    
        return np.where(borderCells[1:-1, 1:-1]) # Return edge rows and columns as a tuple

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
        Sx, Sy = self._calcFiniteSlopes(self._griddata, self._georef_info.dx, self._georef_info.nx, self._georef_info.ny)  # Calculate slope in X and Y directions
    
        AspectRad = np.arctan2(Sy, Sx) # Angle of aspect
        SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))  # magnitude of slope in radians
    
        self._griddata = 255.0 * ((np.cos(elevRad) * np.cos(SmagRad)) + (np.sin(elevRad)* np.sin(SmagRad) * np.cos(azRad - AspectRad)))
            
class FilledElevation(Elevation):
    
    aggradation_slope = 0.000000001
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('elevation',), '_create_from_elevation'))
    
    def _create_from_elevation(self, *args, **kwargs):

        elevation = kwargs['elevation']
        self._copy_info_from_grid(elevation)
        self.__flood(self.aggradation_slope) 
            
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


    def __flood(self, aggSlope = 0.0):
        # dem is a numpy array of elevations to be flooded, aggInc is the minimum amount to increment elevations by moving upstream
        # use priority flood algorithm described in  Barnes et al., 2013
        # Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models
        # NOTE: They have another algorithm to make this more efficient, but to use that and a slope makes things more
        # complicated
                            
        priority_queue = FilledElevation.priorityQueue() # priority queue to sort filling operation
    
        #Create a grid to keep track of which cells have been filled        
        closed = np.zeros_like(self._griddata)
    
        #Add all the edge cells to the priority queue, mark those cells as draining (not closed)
        edgeRows, edgeCols = self.findDEMedge()
            
        for i in range(len(edgeCols)):
            row, col = edgeRows[i], edgeCols[i]
    
            closed[row, col] = True
            priority_queue.put(self._griddata[row, col], (row, col)) # store the indices as a vector of row column, in the priority queue prioritized by the dem value
    
        #While there is anything left in the priority queue, continue to fill holes
        while not priority_queue.isEmpty():
                        
            elevation, (row, col) = priority_queue.get()

            neighborRows, neighborCols, dxMults = self._getNeighborIndices(row, col)
            
            dxs = self._georef_info.dx * dxMults
    
            #Look through the upstream neighbors
            for i in range(len(neighborCols)):
                if not closed[neighborRows[i], neighborCols[i]]:
                    #Do I need to increment(ramp) things or can I leave things flat? I think I can increment b/c my priority queue is stabley sorted
    
                    #If this was a hole (lower than the cell downstream), fill it
                    if self._griddata[neighborRows[i], neighborCols[i]] <= elevation:                        
                        self._griddata[neighborRows[i], neighborCols[i]] = elevation + aggSlope*dxs[i]
    
                    closed[neighborRows[i], neighborCols[i]] = True
                    priority_queue.put(self._griddata[neighborRows[i], neighborCols[i]], [neighborRows[i], neighborCols[i]])
                    
class Area(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                                   (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                                   (('gdal_filename',), '_read_gdal'), 
                                   (('flow_direction',), '_create_from_flow_direction'))    

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
        
        import itertools
        
        for i, j in itertools.izip(ind_i, ind_j):  # Loop through all the data in sorted order    
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)
            if is_good:
                area[i_next, j_next] += area[i,j]
    
        self._griddata = area # Return non bc version of area

    def _area_per_pixel(self, *args, **kwargs):
        return self._georef_info.dx**2 * np.ones((self._georef_info.ny, self._georef_info.nx))

    def _mean_pixel_dimension(self, *args, **kwargs):
        return self._georef_info.dx * np.ones_like(kwargs['flow_direction']._griddata, self.dtype)
    
class GeographicArea(GeographicGridMixin, Area):
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
        dx = self._mean_pixel_dimension() * flow_dir.pixel_scale()

        import itertools
        for i, j in itertools.izip(ind_i, ind_j):  # Loop through all the data in sorted order  
              
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
    
    def is_along_flow_length(self, from_index, to_index):
        
        (i_to, j_to) = to_index
        flow_code = self.__flow_directions[i_to, j_to]
        
        return flow_code == self.__flow_direction_for_length(from_index, to_index)
        
    def _mean_pixel_dimension(self, *args, **kwargs):
        return self._georef_info.dx * np.ones_like(kwargs['flow_direction']._griddata, self.dtype)

class GeographicFlowLength(GeographicGridMixin, FlowLength):
    pass
    
class MaxFlowLengthTrackingMixin(object):
    
    def _calculate_by_tracking_down_max_flow_length(self, *args, **kwargs):
        flow_dir = kwargs['flow_direction']
        flow_length = kwargs['flow_length']
        if kwargs.get('sorted_indexes') is not None:
            idcs = kwargs.get('sorted_indexes')
        else:
            idcs = flow_dir.sort()
                
        [ind_i, ind_j] = np.unravel_index(idcs, flow_dir._griddata.shape)
        
        import itertools
        for i, j in itertools.izip(ind_i, ind_j):  # Loop through all the data in sorted order    
            i_next, j_next, is_good = flow_dir.get_flow_to_cell(i,j)  
            if is_good and flow_length.is_along_flow_length((i,j), (i_next, j_next)):
                self._calculate_grid_value((i,j), (i_next, j_next), *args, **kwargs)
                
        
class Relief(BaseSpatialGrid, MaxFlowLengthTrackingMixin):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                           (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                           (('gdal_filename',), '_read_gdal'), 
                           (('flow_direction', 'elevation', 'flow_length'), '_create_from_flow_direction_sorted_indexes_and_elevation'))   
    
    def _create_from_flow_direction_sorted_indexes_and_elevation(self, *args, **kwargs):
        self._copy_info_from_grid(kwargs['flow_direction'], True)
        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)
    
    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        (i,j) = pos
        (i_next, j_next) = next_pos
        drf = kwargs['elevation'][i,j] - kwargs['elevation'][i_next, j_next]
        if drf < 0:
            drf = 0
        self._griddata[i_next, j_next] = self._griddata[i,j] + drf
        
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
        self._griddata = np.power( (kwargs['Ao'] / kwargs['area']._griddata) ** kwargs['theta']) * self._mean_pixel_dimension() * kwargs['flow_direction'].pixel_scale()
        self._calculate_by_tracking_down_max_flow_length(*args, **kwargs)
        
    def _calculate_grid_value(self, pos, next_pos, *args, **kwargs):
        (i,j) = pos
        (i_next, j_next) = next_pos
        self._griddata[i_next, j_next] += self._griddata[i,j]
        

    def _area_per_pixel(self, *args, **kwargs):
        
        return np.power( (kwargs['Ao'] / kwargs['area']._griddata) ** kwargs['theta']) * self._mean_pixel_dimension() * kwargs['flow_direction'].pixel_scale()

class GeographicKsi(Ksi, GeographicGridMixin):
    pass
            
class ChannelSlope(BaseSpatialGrid):
    
    required_inputs_and_actions = ((('nx', 'ny', 'projection', 'geo_transform',),'_create'),
                               (('ai_ascii_filename','EPSGprojectionCode'),'_read_ai'),
                               (('gdal_filename',), '_read_gdal'), 
                               (('flow_direction','elevation'), '__create_from_flow_direction_and_elevation'))
    
    def __create_from_flow_direction_and_elevation(self, *args, **kwargs):
        
        flow_dir = kwargs['flow_direction']
        elevation = kwargs['elevation']
        self._copy_info_from_grid(flow_dir, True)
        if flow_dir.__class__ == FlowDirectionD8:
            self.__calc_D8_slope(flow_dir, elevation)            

    def calc_D8_slope(self, dem,fd):

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
    plt.ion()
    plt.show()
    ax = plt.gca()
    
    return ax

