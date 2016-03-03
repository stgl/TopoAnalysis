import gdal
import numpy as np


def loadGrids(basePath):
    accExt = 'acc_15s.asc'
    fdExt = 'fdir_15s.asc'
    demExt = 'dem_15s.asc'
    # accExt = 'Acc.tif'
    # fdExt = 'Direc.tif'
    # demExt = 'DEM.tif'

    src = gdal.Open(basePath+accExt)
    acc = src.ReadAsArray().astype(np.float)
    src = None

    src = gdal.Open(basePath+fdExt)
    fd = src.ReadAsArray().astype(np.float)
    src = None

    src = gdal.Open(basePath+demExt)
    dem = src.ReadAsArray().astype(np.float)

    lats,longs = getLatsLongsFromGdalDataset(src)
    geoTransform = src.GetGeoTransform()

    src = None

    return dem, acc, fd, lats, longs, geoTransform

def getLatsLongsFromGdalDataset(gdalData):
    gt = gdalData.GetGeoTransform()
    dLong = gt[1]
    dLat = gt[5]

    #Determine the location of the center of the first pixel of the data
    xllcenter = gt[0]+dLong/2.0
    yllcenter = gt[3]+dLat/2.0

    #Assign the latitudinal (i.e. y direction) coordinates
    lats = np.zeros(gdalData.RasterYSize)
    for i in range(len(lats)):
        lats[i] = yllcenter + i*dLat

    #Assign the longitudinal (i.e. x direction) coordinates
    longs = np.zeros(gdalData.RasterXSize)
    for i in range(len(longs)):
        longs[i] = xllcenter + i*dLong

    return lats,longs

def calcPixelAreas(lats,longs,gt):

    #Returns a grid of the area of each pixel within the domain specified by the gdalDataset


    #Get the size of pixels in the lat and long direction
    dLong = gt[1]
    dLat = gt[5]

    #Preallocate space for the areas
    initialAreas = np.zeros((len(lats),len(longs)))

    #Iterate through the grid, calculating the pixel dimensions in length based on the center coordinate of that pixel
    for i in range(len(lats)):
        for j in range(len(longs)):

            #Get the x and y dimension of the pixels (in km), the first two returns of the function
            dx,dy,dist = getDistFromLatLong(lats[i]-dLat/2.0, longs[j]-dLong/2.0, lats[i]+dLat/2.0, longs[j]+dLong/2.0)

            #Calculate the area of the current pixel
            initialAreas[i,j] = np.abs(dx*dy)

    return initialAreas

def calculateArea(flowDirectionGrid,flowAcc,dx,areas,noDataValue,mask = None):

    # Get the sorted indices of the array in reverse order (e.g. largest first)
    idcs = flowAcc.argsort(axis= None)

    #Mask out the indexes outside of the mask
    if not mask is None:
        allRows,allCols = np.unravel_index(idcs,flowAcc.shape)
        gdIdcs = mask[allRows,allCols]
        idcs = idcs[gdIdcs]

    #Loop through all the indices in sorted order
    for idx in idcs:
        #Get the currents row column index
        [i, j] = np.unravel_index(idx, flowAcc.shape)  # Get the row/column indices

        #Get the index of the point downstream of this and the distance to that point
        [downI,downJ,inBounds] = getFlowToCell(i,j,flowDirectionGrid[i,j],flowAcc.shape[0],flowAcc.shape[1])

        #So long as we are not draining to the outskirts, proceed
        if inBounds and not flowAcc[downI,downJ] == noDataValue:

            #Accumulate area
            areas[downI,downJ] += areas[i,j]

    return areas

def calculateAll(flowDirectionGrid,flowAcc,dem,lats,longs,areas, theta,Ao,noDataValue,mask = None):

    #returns five grids of the same size as those input, maxTransportLength, meanDirs (the mean direction to that point),
    # and the maximum elevation along the longest flow path (relief along this flow path is difference in elev b/w that and dem) .
    #Input grids are  1) flowDirectionGrid, a grid of flow directions in ArcGIS convention
    # 2) flowAcc, a grid of the flow accumulations used to sort the grid from the top down, 3) dem, a grid of
    # the elevations 4) lats and 5) longs, the coordinates of the center of the pixels, 4) areas a preallocated grid of drainage
    #areas specifying the area of each pixel, 6) noDataValue, the value used to represent no data, mask a grid of which points to ignore

    # Get the sorted indices of the array in reverse order (e.g. largest first)
    idcs = flowAcc.argsort(axis= None)

    #Mask out the indexes outside of the mask
    if not mask is None:
        allRows,allCols = np.unravel_index(idcs,flowAcc.shape)
        gdIdcs = mask[allRows,allCols]
        idcs = idcs[gdIdcs]

    #Preallocate space for the final answers
    maxZalongMaxL = dem.copy() #The highest elevation point along the flow path initially (before routing flow), is the elevation of that point
    maxTransportLength = np.zeros_like(flowDirectionGrid,dtype=np.float)
    delXPath = np.zeros_like(flowDirectionGrid)
    delYPath = np.zeros_like(flowDirectionGrid)
    intA = np.zeros_like(flowDirectionGrid)

    #Loop through all the indices in sorted order
    for idx in idcs:

        #Get the currents row column index
        [i, j] = np.unravel_index(idx, flowAcc.shape)  # Get the row/column indices

        #Get the index of the point downstream of this and the distance to that point
        [downI,downJ,inBounds] = getFlowToCell(i,j,flowDirectionGrid[i,j],flowAcc.shape[0],flowAcc.shape[1])

        #So long as we are not draining to the outskirts, proceed
        if inBounds and not dem[downI,downJ] == noDataValue:

            #Accumulate area
            areas[downI,downJ] += areas[i,j]

            #How far do we move to the next cell?
            dx,dy,dist =  getDistFromLatLong(lats[i],longs[j],lats[downI],longs[downJ])
            newL = maxTransportLength[i,j] + dist

            #Keep track of the direction, length, and number of cells of the longest flow path
            if maxTransportLength[downI,downJ] < newL:
                maxTransportLength[downI,downJ] = newL
                maxZalongMaxL[downI,downJ] = maxZalongMaxL[i,j]
                delXPath[downI,downJ] = delXPath[i,j] + dx
                delYPath[downI,downJ] = delYPath[i,j] + dy
                intA[downI,downJ] = intA[i,j] + ((Ao/areas[i,j])**theta)*dist

    meanDirs = np.arctan2(delYPath,delXPath)

    #If this cell doesn't have flow into it, don't count it as having a direction
    meanDirs[maxTransportLength == 0] = np.nan

    #If the original dem was a no data value, make sure this is too
    badData = dem == noDataValue
    maxTransportLength[badData] = np.nan
    meanDirs[badData] = np.nan
    maxZalongMaxL[badData] = np.nan
    areas[badData] = np.nan
    intA[badData] = np.nan


    return maxTransportLength, meanDirs, maxZalongMaxL, areas, intA

def getDistFromLatLong(lat0,long0,lat1,long1):
        re = 1000*6371.0 #Earth radius in m

        #Convert to radians
        lat0,long0,lat1,long1 = np.radians(lat0), np.radians(long0),np.radians(lat1),np.radians(long1)

        dx = re*(long1 - long0)*np.cos((lat1+lat0)/2)
        dy = re*(lat1 - lat0)
        dist = np.sqrt(dx**2 + dy**2)

        return dx,dy, dist

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

def calcHillshade(elevGrid,dx,az,elev):
    #Hillshade = calcHillshade(elevGrid,az,elev)
    #Esri calculation for generating a hillshade, elevGrid is expected to be a numpy array

    # Convert angular measurements to radians
    azRad, elevRad = (360 - az + 90)*np.pi/180, (90-elev)*np.pi/180
    Sx, Sy = calcFiniteSlopes(elevGrid, dx)  # Calculate slope in X and Y directions

    AspectRad = np.arctan2(Sy, Sx) # Angle of aspect
    SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))  # magnitude of slope in radians

    return 255.0 * ((np.cos(elevRad) * np.cos(SmagRad)) + (np.sin(elevRad)* np.sin(SmagRad) * np.cos(azRad - AspectRad)))

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