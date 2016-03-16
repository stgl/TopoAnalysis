import gdal
import numpy as np


def loadGrids(basePath):
    accExt = 'acc_15s'
    fdExt = 'dir_15s'
    demExt = 'dem_15s'
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

    #Get the Geo transform from the gdal data, this provides information about pixel dimensions, locations, and orientations
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

def calcPixelAreas(lats,longs,gt, mask = None):

    #Returns a grid of the area of each pixel within the domain specified by the gdalDataset

    re = 6371.0 * 1000.0 #radius of the earth in meters

    if mask is None:
        rows = range(len(lats))
        cols = range(len(longs))
    else:
        idcs = np.where(mask.flatten())
        rows,cols = np.unravel_index(idcs,mask.shape)


    #Get the size of pixels in the lat and long direction
    dLong = gt[1]
    dLat = gt[5]

    #Preallocate space for the areas
    initialAreas = np.zeros((len(lats),len(longs)))

    #Iterate through the grid, calculating the pixel dimensions in length based on the center coordinate of that pixel
    for i in rows:
        for j in cols:

            #Get the bounding coordinates of these pixels
            lat1 = np.radians(lats[i]-dLat/2.0)
            lat2 = np.radians(lats[i]+dLat/2.0)
            long1 = np.radians(longs[j]-dLong/2.0)
            long2 = np.radians(longs[j]+dLong/2.0)

            #Calculate the area of that pixel, assuming a spherical earth
            initialAreas[i,j] = np.abs((re**2)*(long2 - long1)*(np.sin(lat2) - np.sin(lat1)))

    return initialAreas

def calculateArea(flowDirectionGrid,flowAcc,areas,noDataValue,mask = None):

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

def calculateAll(flowDirectionGrid,flowAcc,dem,lats,longs,areas,theta,Ao,noDataValue,mask = None):

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

def calculateAllMultipleThetas(flowDirectionGrid,flowAcc,dem,lats,longs,areas,thetas,Ao,noDataValue,mask = None):

    #The same as calculate all, but will return a m by n by l grid for intA, where l is the length of the input thetas.
    #In other words, this calculates multiple values of intA for each of the theta values present in theta

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
    intA = np.zeros((dem.shape[0],dem.shape[1],len(thetas)))

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

                #Calculate the downstream integrated drainage area for each of the theta values provided
                for k in range(len(thetas)):
                    intA[downI,downJ,k] = intA[i,j,k] + ((Ao/areas[i,j])**thetas[k])*dist

    #Calculate the running mean direction along each flow path (in radians) based on the integrated X and Y components of motion
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
        re = 1000.0*6371.0 #Earth radius in m

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