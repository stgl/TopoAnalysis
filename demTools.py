#Functions for DEM calculations typical of geomorphology
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
import re #regex for checking strings
from matplotlib import pyplot as plt # Used to get user input on plots
from gdalconst import *


class gdalGeoRefInfo:
    def __init__(self, ncols, nrows, projection, geoTransform):
        #A class for storing the gdal information needed for creating georeferenced raster data
        #the geotransform is (xUL, dx, 0, yUL, 0, -dx)
        self.ncols = ncols
        self.nrows = nrows
        self.projection = projection
        self.geoTransform = geoTransform
        self.dx = geoTransform[1]
        self.xllcenter = geoTransform[0]+self.dx/2.0
        self.yllcenter = geoTransform[3]-(self.dx*(nrows-0.5))

class priorityQueue:
    #Implements a priority queue using heapq. Python has a priority queue module built in, but it
    # is not stabley sorted (meaning that two items who are tied in priority are treated arbitrarily, as opposed to being
    # returned on a first in first out basis). This circumvents that by keeping a count on the items inserted and using that
    # count as a secondary priority

    def __init__(self):
        # A counter and the number of items are stored seperately to ensure that items remain stabley sorted and to
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

def writeArcAsciiRaster(geoRefInfo, outfilePath, npArrayData, noDataValue, formatString):
    #A function to write the data stored in the numpy array npArrayData to a ArcInfo Text grid. Gdal doesn't
    #allow creation of these types of data for whatever reason

    header = "ncols     %s\n" % geoRefInfo.ncols
    header += "nrows    %s\n" % geoRefInfo.nrows
    header += "xllcenter %s\n" % geoRefInfo.xllcenter
    header += "yllcenter %s\n" % geoRefInfo.yllcenter
    header += "cellsize %s\n" % geoRefInfo.dx
    header += "NODATA_value %s" % noDataValue

    np.savetxt(outfilePath, npArrayData, header=header, fmt=formatString, comments='')

def createDataSetFromArray(geoRefInfo, outfilePath, GDALDRIVERNAME, arrayData):
    #A function to write the data in the numpy array arrayData into a georeferenced dataset of type
    #  specified by GDALDRIVERNAME, a string, options here: http://www.gdal.org/formats_list.html
    #  This is accomplished by copying the georeferencing information from an existing GDAL dataset,
    #  provided by createDataSetFromArray

    #Get info needed to initialize new dataset
    ncols = geoRefInfo.ncols
    nrows = geoRefInfo.nrows

    #Initialize new data
    drvr = gdal.GetDriverByName(GDALDRIVERNAME)  #  Get the desired driver
    outRaster = drvr.Create(outfilePath, ncols, nrows, 1 , gdal.GDT_Float32)  # Open the file

    #Write geographic information
    outRaster.SetGeoTransform(geoRefInfo.geoTransform)  # Steal the coordinate system from the old dataset
    outRaster.SetProjection(geoRefInfo.projection)   # Steal the Projections from the old dataset

    #Write the array
    outRaster.GetRasterBand(1).WriteArray(arrayData)   # Writes my array to the raster

    outRaster = None

def getGeoRefInfo(gdalDataset):
    #Get info needed to initialize new dataset
    ncols = gdalDataset.RasterXSize
    nrows = gdalDataset.RasterYSize

    #Write geographic information
    geoTransform = gdalDataset.GetGeoTransform()  # Steal the coordinate system from the old dataset
    projection = gdalDataset.GetProjection()  # Steal the Projections from the old dataset

    return gdalGeoRefInfo(ncols, nrows, projection, geoTransform)

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

def getRasterGeoTransformFromAsciiRaster(fileName):

    #Read in the components of the geotransform from the raster, BEWARE! This
    #is specific to how some matlab/ C scripts I have write these rasters. I believe
    #this is that standard arc ascii raster export format, but could be wrong
    with open(fileName, "r") as file:
        line = file.readline()
        nx = int(line.split()[-1])
        line = file.readline()
        ny = int(line.split()[-1])
        line = file.readline()
        xllcenter = float(line.split()[-1])
        line = file.readline()
        yllcenter = float(line.split()[-1])
        line = file.readline()
        dx = float(line.split()[-1])

    xUL = xllcenter - (dx/2.0)
    yUL = (yllcenter - (dx/2.0)) + dx*ny

    return (xUL, dx, 0, yUL, 0, -dx), nx, ny

def AsciiRasterToMemory(fileName, EPSGprojectionCode):

    # the geotransfrom structured as (xUL, dx, skewX, yUL, scewY, -dy)
    gt, nx, ny = getRasterGeoTransformFromAsciiRaster(fileName)

    #Open gdal dataset
    ds = gdal.Open(fileName)

    #Prep output
    memDrv = gdal.GetDriverByName('MEM')  # Create a gdal driver in memory
    dataOut = memDrv.CreateCopy('name', ds, 0) #Copy data to gdal driver

    #Set the location
    dataOut.SetGeoTransform(gt)  # Set the new geotransform

    # Get raster projection
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSGprojectionCode)
    wkt = srs.ExportToWkt()

    # Set projection
    dataOut.SetProjection(wkt)

    return dataOut

def searchDownFlowDirection(flowDir, area, dem, dx, start = None, doTrim = True):
    #NOTE: I am copying the output areas each iteration... should really pre-allocate
    #  i,j,a,z = searchDownFlowDirection(np.array() area,np.array() dem, tuple (rowStart,ColumnStart) start, boolean doTrim)
    #  Searches down the grid of flow Directions, flowDir, to construct a numpy array of the row/column indices, drainage areas,
    #  and elevations of the values along a steepest descent path starting from start. If start is not specified the user is prompted
    #  to select a point. If doTrim is set to true, the final path is displayed and the user is prompted to select a point to trim
    #  the profile at.

    hs = None

    nrows, ncols = dem.shape

    if start == None:
        #If no starting point is specified prompt the user for one
        hs = calcHillshade(dem, dx, 315, 45) #Create hillshade
        plt.title('Pick a point to route down from')
        plt.imshow(dem, interpolation = 'bilinear', cmap='coolwarm')
        plt.imshow(hs, interpolation = 'bilinear', cmap='gray', alpha = 0.5)
        plt.gca().invert_yaxis()
        start = plt.ginput(1)[0] #ginput returns a list of tuples... index out the first item

    print(start)
    rows, cols = np.array([int(round(start[1]))]), np.array([int(round(start[0]))]) # insure starting indices are ints
    a = np.array([area[rows[-1],cols[-1]]]) # add the first items to the list of long profile info
    z = np.array([dem[rows[-1],cols[-1]]])

    #So long as we are not at the edge of the DEM
    while not (rows[-1] == 0 or rows[-1] == nrows-1 or cols[-1] == 0 or cols[-1] == ncols - 1):
        # Get the neighbors in the flow direction corresponding order - note, I might be in danger the way I handle this...
        # Flow directions are indices of arrays which may be different sizes, this is not truly giving me flow directions
        # Because in the getNeighbor function I don't return values at edges.... might want to think about improving this,
        # although it should work in this application, and damn if it isn't clean
        rowNeighbs, colNeighbs = getNeighborIndices(rows[-1], cols[-1], nrows, ncols)[0:2] #find the neighbors, only need the first two returns


        thisRow,thisCol,inBounds = getFlowToCell(rows[-1],cols[-1],flowDir[rows[-1],cols[-1]],nrows,ncols) # Find the indices of the cell we drain too, only need first two inputs

        if not inBounds:
            break

        rows = np.append(rows, thisRow)
        cols = np.append(cols, thisCol)

        # add the latest items to the list of long profile info
        a = np.append(a, dem[rows[-1], cols[-1]])
        z = np.append(z, dem[rows[-1], cols[-1]])

    #Trim off the end of the data if requested
    if doTrim:

        if hs == None:
            hs = calcHillshade(dem, dx, 315, 45) #Create hillshade

        plt.title('Pick a point to trim at')
        plt.imshow(dem, interpolation = 'bilinear', cmap='coolwarm')
        plt.imshow(hs, interpolation = 'bilinear', cmap='gray', alpha = 0.5)
        plt.plot(cols,rows,lw = 2, color = 'black')
        plt.gca().invert_yaxis()
        end = plt.ginput(1)[0]
        dists = np.sqrt((rows - end[1])**2 + (cols - end[0])**2)
        end = np.argmin(dists)

        rows, cols, a, z = rows[:end], cols[:end], a[:end], z[:end]


    return rows, cols, a, z

def buildSearchKernel(dx, windowRadius, isCircular):
    # Function to build a search kernel that is either square or circular

    pxlRadius = int(round(windowRadius/dx)) #Is this right...

    relCoords = np.arange(1 + 2*pxlRadius)-pxlRadius #Relative coordinates of neighbors within this row, col distance

    searchKernelRow, serchKernelCol = np.meshgrid(relCoords,relCoords)

    if isCircular:
        dists = np.sqrt(searchKernelRow**2 + searchKernelRow**2)
        searchKernelRow = searchKernelRow[dists<pxlRadius]
        serchKernelCol = serchKernelCol[dists<pxlRadius]

    return searchKernelRow.flatten(), serchKernelCol.flatten()

def adjustKernel(row,col,grid,searchKernelRows,searchKernelCols):
    # Function to turn the rows and columns of the search kernel specified by searchKernelRows & ...Cols as just indices of the
    # data that is actually within the grid

    #Get the current absolute positions within the grid
    theseRows = searchKernelRows+row
    theseCols = searchKernelCols+col

    #Do any points in searchKernel extend to before the start of the grid? or off the end?
    gdIndices = (theseCols < grid.shape[1]) & (theseRows < grid.shape[0]) & (theseRows >= 0) & (theseCols >= 0)

    #Which rows are within the grid
    goodRows = theseRows[gdIndices]
    goodCols = theseCols[gdIndices]

    #Now that we are all within the grid, of those points within the grid, do any include
    #Do any points in searchKernel include nans ?
    gdIndices = np.logical_not(np.isnan(grid[goodRows, goodCols]))

    return goodRows[gdIndices], goodCols[gdIndices]

def movingWindow(grid, searchKernelRows, searchKernelCols,function):
    # Function to scan the moving window specified by Kernel across the dem, a numpy grid, and apply the specified function
    # to each location. function is a method that returns a single value given any number of inputs, e.g.
    # movingWindow(i) = function(dem[Kernel@i])

    outgrid = np.zeros_like(grid)

    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            theseRows,theseCols = adjustKernel(i,j,grid,searchKernelRows,searchKernelCols) #trim the searchKernel and return it as a list of indices
            outgrid[i,j] = function(grid[theseRows, theseCols]) #apply the specified function to the dem values for this search location

    return outgrid


def findDEMedge(dem):
    # Function to find the cells at the edge of a dem. Dem is a ny x nx array, but may be largely padded
    # by nans. Determines where the edge of the real data is. Does this by finding the maximum value within a 3x3 kernel,
    # if that is not zero, but the original data at the corresponding location is, then that is an edge cell
    ny, nx = dem.shape

    #Pad the data so that we can take a windowed max
    padded = np.zeros((ny+2, nx+2))
    padded[1:-1, 1:-1] = dem
    padded[padded == 0] = np.nan

    # windowMax = np.zeros_like(padded)
    borderCells = np.zeros_like(padded)

    #Iterate through all the data, find the max in a 3 x 3 kernel
    for i in range(ny):
        for j in range(nx):
            # windowMax[i+1, j+1] = np.nanmax(padded[i:i+3, j:j+3])
            borderCells[i+1, j+1] = np.any(np.isnan(padded[i:i+3, j:j+3]))*(~np.isnan(padded[i+1,j+1]))

    #Border cells are those which were not values in the original DEM, but did have values within the 3x3 search kernel
    # borderCells = (~np.isnan(windowMax)) * (np.isnan(padded))

    #This is slower...
    #     # Function to find the cells at the edge of a dem. Dem is a ny x nx array, but may be largely padded
    # # by nans. Determines where the edge of the real data is. Does this by finding the maximum value within a 3x3 kernel,
    # # if that is not zero, but the original data at the corresponding location is, then that is an edge cell
    # ny, nx = dem.shape
    # borderCells = []
    # append = borderCells.append
    #
    # #Iterate through all the data, find the max in a 3 x 3 kernel
    # for i in range(ny):
    #     for j in range(nx):
    #         neighRet = getNeighborIndices(i,j,ny,nx)
    #         rows, cols = neighRet[0], neighRet[1]
    #
    #         if ~i==0 and ~i==ny-1 and ~j==0 and ~j==ny-1:
    #             if np.any(np.isnan(dem[rows,cols]))*(~np.isnan(dem[i,j])):
    #                 append((i, j))
    #         else:
    #             if ~np.isnan(dem[i,j]):
    #                 append((i, j))


    return np.where(borderCells[1:-1, 1:-1]) # Return edge rows and columns as a tuple

def findRidgeTops(dem, area, minConnectedPixels):
    #function to find all the connected hilltops in a dem... not sure how well this will work
    #the idea is to find all the points that are both convex and have zero drainage area and are connected with
    #some number of other pixels that meet this criteria
    return None

def getNeighborIndices(row, col, ny, nx):
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
    inBounds = (outRows >= 0)*(outRows < ny)*(outCols >= 0)*(outCols < nx)
    return outRows[inBounds], outCols[inBounds], dxMults[inBounds]

def priorityFlood(dem, dx, aggSlope = 0.0):
    # dem is a numpy array of elevations to be flooded, aggInc is the minimum amount to increment elevations by moving upstream
    # use priority flood algorithm described in  Barnes et al., 2013
    # Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models
    # NOTE: They have another algorithm to make this more efficient, but to use that and a slope makes things more
    # complicated

    ny, nx = dem.shape  # Size of dem grid

    open = priorityQueue() # priority queue to sort filling operation

    #Create a grid to keep track of which cells have been filled
    closed = np.zeros_like(dem)

    #Add all the edge cells to the priority queue, mark those cells as draining (not closed)
    edgeRows, edgeCols = findDEMedge(dem)

    for i in range(len(edgeCols)):
        row, col = edgeRows[i], edgeCols[i]

        closed[row, col] = True
        open.put(dem[row, col], [row, col]) # store the indices as a vector of row column, in the priority queue prioritized by the dem value

    #While there is anything left in the priority queue, continue to fill holes
    while not open.isEmpty():
        elevation, rowCol = open.get()
        row, col = rowCol
        neighborRows, neighborCols, dxMults = getNeighborIndices(row, col, ny, nx)
        dxs = dx * dxMults

        #Look through the upstream neighbors
        for i in range(len(neighborCols)):
            if not closed[neighborRows[i], neighborCols[i]]:
                #Do I need to increment(ramp) things or can I leave things flat? I think I can increment b/c my priority queue is stabley sorted

                #If this was a hole (lower than the cell downstream), fill it
                if dem[neighborRows[i], neighborCols[i]] <= elevation:
                    dem[neighborRows[i], neighborCols[i]] = elevation + aggSlope*dxs[i]

                closed[neighborRows[i], neighborCols[i]] = True
                open.put(dem[neighborRows[i], neighborCols[i]], [neighborRows[i], neighborCols[i]])

    return dem

def getUTMZone(dataset):
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


def approximateDxFromGeographicData(dataset):
    #Function to return the approximate grid spacing in Meters. Will return the closest integer value
    tVect = dataset.GetGeoTransform()  # Get the coordinate transform vector, (ulx, dx, xRot, uly, yRot, -dx)
    dTheta = tVect[1] #Angular grid spacing
    metersPerDegree = 110000 #110 km per degree (approximation)
    return int(dTheta*metersPerDegree) #convert degrees to meters


def convertToUTM(dataset, dx, utmZone):

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

def clipRasterToRaster(srcFilename, dataExtentFilename, outputFilename):

    # Source
    src = gdal.Open(srcFilename, GA_ReadOnly)
    src_proj = src.GetProjection()
    src_geotrans = src.GetGeoTransform()

    # We want a section of source that matches this:
    match_ds = gdal.Open(dataExtentFilename, GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize

    # Output / destination
    dst = gdal.GetDriverByName('GTiff').Create(outputFilename, wide, high, 1, GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, GRA_Bilinear)
    # gdal.ReprojectImage(src, dst, None, None, GRA_Bilinear)
    # gdal.ReprojectImage(dst, src, None, None, GRA_Bilinear)

    #
    return dst

def clipRasterToShape(srcFilename,shpFilename,outputFilename):

    # drvr = gdal.GetDriverByName(GdalDriver)
    # drvr.Create(outputFilename,1,1,1)
#    warp= 'gdalwarp -cutline \'%s\' -crop_to_cutline -dstalpha \'%s\' \'%s\'' % (shpFilename, srcFilename, outputFilename)
    warp= 'gdalwarp -cutline \'%s\' -crop_to_cutline \'%s\' \'%s\'' % (shpFilename, srcFilename, outputFilename)

    os.system(warp)
    drvr = None

def createMaskFromShape(geoRefInfo,shpFilename,noDataValue = 0):

    #Open Shapefile
    source_ds = ogr.Open(shpFilename)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()


    maskSrc = gdal.GetDriverByName('MEM').Create('name',geoRefInfo.ncols,geoRefInfo.nrows, 1, gdal.GDT_Byte)
    maskSrc.SetGeoTransform(geoRefInfo.geoTransform)
    maskBand = maskSrc.GetRasterBand(1)
    maskBand.SetNoDataValue(noDataValue)


    # 5. Rasterize why is the burn value 0... isn't that the same as the background?
    gdal.RasterizeLayer(maskSrc, [1], source_layer, burn_values=[1])

    grid = maskSrc.ReadAsArray().astype(np.float)

    maskSrc = None
    source_ds = None

    return grid


def getDEMcoords(GdalData, dx):

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

def assignBCs(elevGrid):
    # Pads the boundaries of a grid
    # Boundary condition pads the boundaries with equivalent values
    # to the data margins, e.g. x[-1,1] = x[1,1]
    # This creates a grid 2 rows and 2 columns larger than the input

    ny, nx = elevGrid.shape  # Size of array
    Zbc = np.zeros((ny + 2, nx + 2))  # Create boundary condition array
    Zbc[1:-1,1:-1] = elevGrid  # Insert old grid in center

    #Assign boundary conditions - sides
    Zbc[0, 1:-1] = elevGrid[0, :]
    Zbc[-1, 1:-1] = elevGrid[-1, :]
    Zbc[1:-1, 0] = elevGrid[:, 0]
    Zbc[1:-1, -1] = elevGrid[:,-1]

    #Assign boundary conditions - corners
    Zbc[0, 0] = elevGrid[0, 0]
    Zbc[0, -1] = elevGrid[0, -1]
    Zbc[-1, 0] = elevGrid[-1, 0]
    Zbc[-1, -1] = elevGrid[-1, 0]

    return Zbc

def calcFiniteSlopes(elevGrid, dx):
    # sx,sy = calcFiniteDiffs(elevGrid,dx)
    # calculates finite differences in X and Y direction using the
    # 2nd order/centered difference method.
    # Applies a boundary condition such that the size and location
    # of the grids in is the same as that out.

    # Assign boundary conditions
    Zbc = assignBCs(elevGrid)

    #Compute finite differences
    Sx = (Zbc[1:-1, 2:] - Zbc[1:-1, :-2])/(2*dx)
    Sy = (Zbc[2:,1:-1] - Zbc[:-2, 1:-1])/(2*dx)

    return Sx, Sy

def calcAverageSlopeOfGridSubset(gridSubset,dx):
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

def calcFiniteSlopesOverWindow(elevGrid, dx,N):
    # sx,sy = calcFiniteDiffs(elevGrid,dx)
    # calculates finite differences in X and Y direction using a finite difference
    # kernel that extends N cells from the center. The width of the kernel is then
    # 2N + 1 by 2N + 1. Applies a boundary condition such that the size and location
    # of the grids in is the same as that out. However, the larger N is, the more NoData
    #Will be around the edges .

    #Compute finite differences
    Sx = (elevGrid[N:-N, (2*N):] - elevGrid[N:-N, :-(2*N)])/(((2*N)+1)*dx)
    Sy = (elevGrid[(2*N):,N:-N] - elevGrid[:-(2*N), N:-N])/(((2*N)+1)*dx)

    print(Sx.shape)
    print(Sy.shape)

    #Create two new arrays of the original DEMs size
    SxPadded = np.empty(elevGrid.shape)
    SxPadded[:] = np.NAN
    SyPadded = np.empty(elevGrid.shape)
    SyPadded[:] = np.NAN

    SyPadded[N:-N, N:-N] = Sy
    SxPadded[N:-N, N:-N] = Sx

    return SxPadded, SyPadded

def calcFiniteCurv(elevGrid, dx):
    #C = calcFiniteCurv(elevGrid, dx)
    #calculates finite differnces in X and Y direction using the centered difference method.
    #Applies a boundary condition such that the size and location of the grids in is the same as that out.

    #Assign boundary conditions
    Zbc = assignBCs(elevGrid)

    #Compute finite differences
    Cx = (Zbc[1:-1, 2:] - 2*Zbc[1:-1, 1:-1] + Zbc[1:-1, :-2])/dx**2
    Cy = (Zbc[2:, 1:-1] - 2*Zbc[1:-1, 1:-1] + Zbc[:-2, 1:-1])/dx**2

    return Cx+Cy


def calcFiniteCurvOverWindow(grid,dx, winRad):
    #C = calcFiniteCurv(elevGrid, dx)
    #calculates finite differnces in X and Y direction using the centered difference method.
    #Applies a boundary condition such that the size and location of the grids in is the same as that out.

    #Assign boundary conditions
    Curv = np.zeros_like(grid)*np.nan

    #Compute finite differences
    Cx = (grid[winRad:-winRad, (2*winRad):] - 2*grid[winRad:-winRad, winRad:-winRad] + grid[winRad:-winRad, :-(2*winRad)])/(2*dx*winRad)**2
    Cy = (grid[(2*winRad):, winRad:-winRad] - 2*grid[winRad:-winRad, winRad:-winRad] + grid[:-(2*winRad), winRad:-winRad])/(2*dx*winRad)**2

    Curv[winRad:-winRad,winRad:-winRad] = Cx+Cy
    return Curv

def calcContourCurvature(elevGrid,dx):
    # kt = (fxx*fy^2 - 2*fxyfxfy + fyy*fx^2)/((fx^2 + fy^2)*sqrt((fx^2 + fy^2)+1)

    #Preallocate
    Kt = np.zeros_like(elevGrid)*np.nan

    #First derivatives, 2nd order centered difference
    fx = (elevGrid[1:-1,2:] - elevGrid[1:-1,:-2])/(dx*2)
    fy = (elevGrid[2:,1:-1] - elevGrid[:-2,1:-1])/(dx*2)

    #Second derivatives, 2nd order centered differece
    fxx = (elevGrid[1:-1,2:] - 2*elevGrid[1:-1,1:-1] + elevGrid[1:-1,:-2])/(dx**2)
    fyy = (elevGrid[2:,1:-1] - 2*elevGrid[1:-1,1:-1] + elevGrid[:-2,1:-1])/(dx**2);

    #Partial derivative
    fxy = (elevGrid[2:,2:] - elevGrid[2:,1:-1] - elevGrid[1:-1,2:] + 2*elevGrid[1:-1,1:-1] - elevGrid[:-2,1:-1] - elevGrid[1:-1,:-2] + elevGrid[:-2,:-2])
    fxy = fxy/(4*dx**2)

    #Contour curvature
    Kt[1:-1, 1:-1] = (fxx*fy**2 - 2*fxy*fx*fy + fyy*fx**2)/((fx**2 + fy**2)*np.sqrt((fx**2 + fy**2)+1))

    return Kt

def calcHillshade(elevGrid,dx,az,elev):
    #Hillshade = calcHillshade(elevGrid,az,elev)
    #Esri calculation for generating a hillshade, elevGrid is expected to be a numpy array

    # Convert angular measurements to radians
    azRad, elevRad = (360 - az + 90)*np.pi/180, (90-elev)*np.pi/180
    Sx, Sy = calcFiniteSlopes(elevGrid, dx)  # Calculate slope in X and Y directions

    AspectRad = np.arctan2(Sy, Sx) # Angle of aspect
    SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))  # magnitude of slope in radians

    return 255.0 * ((np.cos(elevRad) * np.cos(SmagRad)) + (np.sin(elevRad)* np.sin(SmagRad) * np.cos(azRad - AspectRad)))

def calcD8Area(elevGrid,dx):

    # I am returning area and flowDirections but NOTE!I might be in danger the way I handle this...
    # Flow directions are indices of arrays which may be different sizes, this is not truly giving me flow directions
    # Because in the getNeighbor function I don't return values at edges.... might want to think about improving this,
    # although it should work in this application
    # Calculate the D8 drainage area for the numpy array representing a DEM in elevGrid, assumes that elevGrid has already been filled
    # Assigns BCs to deal with edges - these should be changed to handle flow differently, currently it will not force flow through the edges

    pxlArea = dx**2  # area of a pixel

    idcs = elevGrid.argsort(axis= None)[::-1] # Get the sorted indices of the array in reverse order (e.g. largest first)
    area = pxlArea*np.ones_like(elevGrid)  # All pixels have at least their own area

    [nrows, ncols] = elevGrid.shape  # How big is the BC array? We don't need to do calulate area on the boundarys...
    flowDir = np.zeros_like(elevGrid, dtype=int)

    for idx in idcs:  # Loop through all the data in sorted order
        [i, j] = np.unravel_index(idx, elevGrid.shape)  # Get the row/column indices

        if not np.isnan(elevGrid[i, j]):
            iNeighbs, jNeighbs, dxMults = getNeighborIndices(i, j, nrows, ncols) # Find the actual indices of the neighbors

            #Find the distance to each of the neighbors
            dxs = dx*dxMults

            #Assign the flow direction of the current point
            thisFD,iRel,jRel = assignD8FlowDir(iNeighbs-i, jNeighbs-j, (elevGrid[i, j] - elevGrid[iNeighbs, jNeighbs])/dxs)

            if not np.isnan(thisFD):
                # accumulate current area, downstream area
                flowDir[i,j] = thisFD
                area[i+iRel, j+jRel] += area[i, j]


    return area, flowDir # Return non bc version of area

def calcD8AreaSlope(filledDem,elevGrid,dx):

    # I am returning area and flowDirections but NOTE!I might be in danger the way I handle this...
    # Flow directions are indices of arrays which may be different sizes, this is not truly giving me flow directions
    # Because in the getNeighbor function I don't return values at edges.... might want to think about improving this,
    # although it should work in this application
    # Calculate the D8 drainage area for the numpy array representing a DEM in elevGrid, assumes that elevGrid has already been filled
    # Assigns BCs to deal with edges - these should be changed to handle flow differently, currently it will not force flow through the edges

    pxlArea = dx**2  # area of a pixel

    idcs = filledDem.argsort(axis= None)[::-1] # Get the sorted indices of the array in reverse order (e.g. largest first)
    area = pxlArea*np.ones_like(elevGrid)  # All pixels have at least their own area
    slope = np.ones_like(elevGrid)*np.nan

    [nrows, ncols] = elevGrid.shape  # How big is the BC array? We don't need to do calulate area on the boundarys...
    flowDir = np.zeros_like(elevGrid, dtype=int)

    for idx in idcs:  # Loop through all the data in sorted order
        [i, j] = np.unravel_index(idx, elevGrid.shape)  # Get the row/column indices

        if not np.isnan(elevGrid[i, j]):
            iNeighbs, jNeighbs, dxMults = getNeighborIndices(i, j, nrows, ncols) # Find the actual indices of the neighbors

            #Find the distance to each of the neighbors
            dxs = dx*dxMults

            #Assign the flow direction of the current point
            thisFD,iRel,jRel = assignD8FlowDir(iNeighbs-i, jNeighbs-j, (filledDem[i, j] - filledDem[iNeighbs, jNeighbs])/dxs)

            if not np.isnan(thisFD):
                # accumulate current area, downstream area
                flowDir[i,j] = thisFD
                area[i+iRel, j+jRel] += area[i, j]

                #Calculate slope to downstream cell
                thisDx = np.sqrt((i-iRel)**2 + (j-jRel)**2)*dx
                slope[i,j] = (elevGrid[i,j] - elevGrid[i+iRel, j+jRel])/thisDx


    return area, flowDir, slope # Return non bc version of area


def assignD8FlowDir(iRel,jRel,slopes):
    ## iRel and jRel are the relative indices from the current point to the surrounding points, the slopes of which are
    ## stored in 'slopes'

    #Search kernel for D8 flow routing, the relative indices of each of the 8 points surrounding a pixel, this is
    # ArcGIS convection
    # |i-1,j-1  i-1,j  i-1,j+1|  |32 64 128|
    # |i,j-1     i,j     i,j+1|  |16  X  1 |
    # |i+1,j-1  i+1,j  i+1,j+1|  |8   4  2 |

    idx = np.argmax(slopes)  # Find steepest surrounding slope
    iOut = iRel[idx]
    jOut = jRel[idx] # Find the index of the steepest surrounding slope

    fd = np.nan

    if iOut == 0 and jOut == 1:
        fd = 1
    elif iOut == 1 and jOut == 1:
        fd = 2
    elif iOut == 1 and jOut == 0:
        fd = 4
    elif iOut == 1 and jOut == -1:
        fd = 8
    elif iOut == 0 and jOut == -1:
        fd = 16
    elif iOut == -1 and jOut == -1:
        fd = 32
    elif iOut == -1 and jOut == 0:
        fd = 64
    elif iOut == -1 and jOut == 1:
        fd = 128

    return fd, iOut, jOut

def calcD8SlopeGrid(dem,fdGrid,dx):

    gridShape = dem.shape
    slopes = np.zeros_like(dem)*np.nan #Preallocate with nans

    for i in range(gridShape[0]):
        for j in range(gridShape[1]):
            slopes[i, j] = getD8slope(dem,fdGrid[i,j],dx,i,j,gridShape[0],gridShape[1])

    return slopes

def getD8slope(dem,fd,dx,i,j,nRows,nCols):

    #Function to get the slope cell in the down flow direction
    iOut,jOut,isGood = getFlowToCell(i,j,fd,nRows,nCols)

    if isGood:
        dist = np.sqrt((i-iOut)**2 + (j-jOut)**2)*dx
        return (dem[i, j]-dem[iOut, jOut])/dist
    else:
        return np.nan

def getFlowToCell(i,j,fd,nRows,nCols):
    #Function to get the indices of the cell that is drained to based on the flow direction specified in fd

    iOut = None
    jOut = None
    isGood = False

    if fd == 1 and j+1 < nCols:
        iOut = i
        jOut = j+1
    elif fd == 2 and i+1 < nRows and j+1 < nCols:
        iOut = i+1
        jOut = j+1
    elif fd == 4 and i+1 < nRows:
        iOut = i+1
        jOut = j
    elif fd == 8 and i+1 < nRows and j-1 >= 0:
        iOut = i+1
        jOut = j-1
    elif fd == 16 and j-1 >= 0:
        iOut = i
        jOut = j-1
    elif fd == 32 and i-1 >= 0 and j-1 >= 0:
        iOut = i-1
        jOut = j-1
    elif fd == 64 and i-1 >= 0:
        iOut = i-1
        jOut = j
    elif fd == 128 and i-1 >= 0 and j+1 < nCols:
        iOut = i-1
        jOut = j+1

    if not(iOut is None):
        isGood = True

    return iOut, jOut, isGood

def convertRiverToolsFlwDirToArc(flowDir,noData):
    # Function to convert river tools flow directions to arcGisFlowDirections
      # ArcGIS convection
    # |i-1,j-1  i-1,j  i-1,j+1|  |32 64 128|
    # |i,j-1     i,j     i,j+1|  |16  X  1 |
    # |i+1,j-1  i+1,j  i+1,j+1|  |8   4  2 |
    # In river tools convention the 1 is in the top right

    convertedFlowDir = int(np.log2(flowDir))
    convertedFlowDir -= 1
    convertedFlowDir[convertedFlowDir == -1] = 7
    convertedFlowDir = 2**convertedFlowDir
    convertedFlowDir[flowDir == noData] = noData

    return convertedFlowDir


