import hydroshedsProcessing as hp
import gdal as gd
import numpy as np
import os
from time import strftime
import demTools as dt
from matplotlib import pyplot as plt

#What is the base path for the file we need to run? This assumes that the extensions are specified as written
# in the hydroshesProccessing.py function 'loadGrids'
basePath = '/Users/johnstone/Box Sync/GlobalSteepness/Datasets/Grids'
# basePath = '/data/sesfs/scarp1/GlobalSteepness'

#What is the name of the grids we are running this on?
baseNames = ('af','as','au','ca','eu','na','sa')
# baseNames = ('af','au','ca')

#What is the no data value? Should improve this to automatically find it
ndv = -32768.0

#What parameters do we want to use?
thetas = np.array([0.3, 0.5, 0.7]) #Scaling power of drainage area
thetaStrings = ('0p3','0p5','0p7')#String representation of those thetas for file names
Ao = 500**2 #Approximate drainage area of a pixel, used to normalized drainage area when calculating the integrated area

#Get a gdal driver to use to save the grid
drvrName = 'GTIFF'

for name in baseNames:

    print('Working on continent: %s' %name)
    print('At the tone, current time will be: %s' %strftime("%m/%d/%Y %H:%M:%S"))

    #Where are these files being stored
    thisBasePath = os.path.join(basePath,name,name+'_')


    #Load in some source data (e.g. that has geo referencing info)
    src = gd.Open(thisBasePath+'dem_15s')
    dem = src.ReadAsArray().astype(np.float)


    #Get a structure that stores the georeferencing info
    geoRefInfo = dt.getGeoRefInfo(src)

    # ## These lines will convert the numpy grids to Tiffs and delete the numpy grids
    # maxTransportLength = np.load(thisBasePath+'maxL.npy')
    # dt.createDataSetFromArray(geoRefInfo, thisBasePath+'maxL.tif', drvrName, maxTransportLength)
    # os.remove(thisBasePath+'maxL.npy')
    #
    # meanDirs = np.load(thisBasePath+'meanDir.npy')
    # dt.createDataSetFromArray(geoRefInfo, thisBasePath+'meanDirs.tif', drvrName, meanDirs)
    # os.remove(thisBasePath+'meanDir.npy')
    #
    # maxZalongMaxL = np.load(thisBasePath+'maxZalongL.npy')
    # dt.createDataSetFromArray(geoRefInfo, thisBasePath+'maxZalongMaxL.tif', drvrName, maxZalongMaxL)
    # os.remove(thisBasePath+'maxZalongL.npy')
    #
    # areas = np.load(thisBasePath+'area.npy')
    # dt.createDataSetFromArray(geoRefInfo, thisBasePath+'areas.tif', drvrName, areas)
    # os.remove(thisBasePath+'area.npy')
    #
    # intA = np.load(thisBasePath+'intA.npy')
    # for i in range(len(thetas)):
    #     dt.createDataSetFromArray(geoRefInfo, thisBasePath+'intA_theta_%s.tif'%thetaStrings[i], drvrName, intA[:,:,i])
    # os.remove(thisBasePath+'intA.npy')

    #Load in the other grids
    #
    # src = gd.Open(thisBasePath+'maxL.tif')
    # maxTransportLength = src.ReadAsArray.astype(np.float)
    #
    # src = gd.Open(thisBasePath+'meanDir.tif')
    # meanDirs = src.ReadAsArray().astype(np.float)

    src = gd.Open(thisBasePath+'maxZalongMaxL.tif')
    maxZalongMaxL = src.ReadAsArray().astype(np.float)

    # src = gd.Open(thisBasePath+'area.tif')
    # area = src.ReadAsArray().astype(np.float)

    intA = np.zeros((dem.shape[0],dem.shape[1],len(thetas)))
    for i in range(thetas.shape[0]):
        src = gd.Open(thisBasePath+'intA_theta_%s.tif'%thetaStrings[i])
        intA[:,:,i] = src.ReadAsArray().astype(np.float)

    #Calculate the downstream integrated fluvial relief
    dsRelief = maxZalongMaxL - dem


    ## Create figures
    k = 1 #Which theta to plot?
    thisIntA = intA[:,:,k]
    theta = thetas[k]
    #Make evenly spaced pins of intA
    nBins = 100
    intABinWidth = (np.nanmax(thisIntA) - np.nanmin(thisIntA))/nBins
    binCents =(np.arange(nBins) + 1)*intABinWidth

    #For each of the bin centers, find the top p percent of reliefs at the intA value
    p = 0.01

    #Store a mask of where the steepest reliefs are
    mostSteepestMask = np.zeros_like(thisIntA) == 1

    for x in binCents:
        #What range of intAs are we looking at this iteration?
        intARange = (thisIntA>=(x-intABinWidth)) & (thisIntA<(x+intABinWidth))
        rangeRows,rangeCols = np.where(intARange)

        #what are the reliefs and intAs in that range?
        rlfInRange = dsRelief[intARange]

        #Obtain the sorted idcs of reliefs in order of largest to smallest
        topPIdcs = rlfInRange.argsort(axis=None)[::-1]

        #Clip out the largest p percent
        topPIdcs = topPIdcs[:np.round(len(rangeRows)*p)]

        mostSteepestMask[rangeRows[topPIdcs],rangeCols[topPIdcs]] = True


    #Create a random mask of a subset of the data that contains elveation
    subSetPerc = 0.05

    #Find the coordinates of the points that are valid
    gdRows,gdCols = np.where(dem>0)

    #Select a random subset of the valid points
    randIdcs = np.random.randint(0,len(gdRows),np.round(subSetPerc*len(gdRows)))

    #Define a range of Ksn values to plot, will plot twice and double this max
    maxKsn = 300
    maxChiSlope = maxKsn/(Ao**theta)

    #Prep the plot axis
    plt.rc('text',usetex=True)

    #Plot a random subset of the data
    plt.plot(thisIntA[gdRows[randIdcs],gdCols[randIdcs]],dsRelief[gdRows[randIdcs],gdCols[randIdcs]], '.k')

    #Plot the highest relief subset of the data
    plt.plot(thisIntA[mostSteepestMask], dsRelief[mostSteepestMask],'or',label = 'Top %.0f %%'%(p*100))

    #Plot reference values of ksn
    plt.plot(binCents,binCents*maxChiSlope/2.0,'-g',label='Ksn = %s' % (maxKsn/2.0))
    plt.plot(binCents,binCents*maxChiSlope,'-r',label='Ksn = %s' % maxKsn)
    plt.plot(binCents,binCents*maxChiSlope*2.0,'-b',label='Ksn = %s' % (maxKsn*2.0))

    #Label axes
    plt.xlabel(r'Downstream integrated $A^{m/n}, \xi, [m], \theta = %.1f $'%theta)
    plt.ylabel('Relief [m]')
    plt.legend(loc='upper left')
    plt.ylim((0,9000))
    plt.savefig(thisBasePath+'ReliefVsXi_theta%s.png'%thetaStrings[k])
    plt.clf()
    plt.close()


    maxTransportLength = None
    meanDirs = None
    maxZalongMaxL = None
    areas = None
    intA = None
    src = None

    print('Finished continent: %s' %name)
    print('At the tone, current time will be: %s' %strftime("%m/%d/%Y %H:%M:%S"))

