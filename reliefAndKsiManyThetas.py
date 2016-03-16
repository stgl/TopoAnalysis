import hydroshedsProcessing as hp
import numpy as np
import os
from time import strftime
from matplotlib import pyplot as plt

#What is the base path for the file we need to run? This assumes that the extensions are specified as written
# in the hydroshesProccessing.py function 'loadGrids'
basePath = '/Users/johnstone/Box Sync/GlobalSteepness/Datasets/Grids'
# basePath = '/data/sesfs/scarp1/GlobalSteepness'

#What is the name of the grids we are running this on?
baseNames = ('af','as','au','ca','eu','na','sa')
baseNames = ('na', 'sa')
# baseNames = ('ca','fail')

#What is the no data value? Should improve this to automatically find it
ndv = -32768.0

#What parameters do we want to use?
thetas = (0.3, 0.5, 0.7) #Scaling power of drainage area
Ao = 500**2 #Approximate drainage area of a pixel, used to normalized drainage area when calculating the integrated area

for name in baseNames:

    print('Working on continent: %s' %name)
    print('At the tone, current time will be: %s' %strftime("%m/%d/%Y %H:%M:%S"))

    #Where are these files being stored
    thisBasePath = os.path.join(basePath,name,name+'_')

    #Load in the different grids
    dem, flowAcc, flowDir, lats, longs, geoTransform = hp.loadGrids(thisBasePath)

    #Calculate the areas of each pixels
    areas = hp.calcPixelAreas(lats, longs, geoTransform, np.logical_not(dem==ndv))

    #Do all the things!
    maxTransportLength, meanDirs, maxZalongMaxL, areas, intA = hp.calculateAllMultipleThetas(flowDir,flowAcc,dem,lats,longs,areas, thetas,Ao,ndv,mask = None)

    #Save the grids
    np.save(thisBasePath+'maxL.npy',maxTransportLength)
    np.save(thisBasePath+'meanDir.npy',meanDirs)
    np.save(thisBasePath+'maxZalongL.npy',maxZalongMaxL)
    np.save(thisBasePath+'area.npy',areas)
    np.save(thisBasePath+'intA.npy',intA)

    maxTransportLength = None
    meanDirs = None
    maxZalongMaxL = None
    areas = None
    intA = None

    print('Finished continent: %s' %name)
    print('At the tone, current time will be: %s' %strftime("%m/%d/%Y %H:%M:%S"))

