import hydroshedsProcessing as hp
import numpy as np
from matplotlib import pyplot as plt

#What is the base path for the file we need to run? This assumes that the extensions are specified as written
# in the hydroshesProccessing.py function 'loadGrids'
# basePath = '/Users/johnstone/Box Sync/GlobalSteepness/Datasets/TestCases/mad_'
# outPath = '/Users/johnstone/Box Sync/GlobalSteepness/Datasets/TestCases/Derivatives/mad_'

basePath = 'D:\Box Sync\GlobalSteepness\Datasets\TestCases\cors_'
outPath = 'D:\Box Sync\GlobalSteepness\Datasets\TestCases\Derivatives\cors_'

#Load in the different grids
dem, flowAcc, flowDir, lats, longs,geoTransform = hp.loadGrids(basePath)

#What is the no data value? Should improve this to automatically find things
ndv = -9999.0

#What parameters do we want to use?
theta = 0.5 #Scaling power of drainage area
Ao = 500**2 #Approximate drainage area of a pixel

#Calculate the areas of each pixels
areas = hp.calcPixelAreas(lats,longs,geoTransform)

#Do all the things!
maxTransportLength, meanDirs, maxZalongMaxL, areas, intA = hp.calculateAll(flowDir,flowAcc,dem,lats,longs,areas, theta,Ao,ndv,mask = None)

#Save the grids
# np.save(outPath+'maxL.npy',maxTransportLength)
# np.save(outPath+'meanDir.npy',meanDirs)
# np.save(outPath+'maxZalongL.npy',maxZalongMaxL)
# np.save(outPath+'area.npy',areas)
# np.save(outPath+'intA.npy',intA)



#Makes some plots!
maxKsn = 300
maxChiSlope = maxKsn/(Ao**theta)

plt.rc('text',usetex=True)
plt.plot(intA,(maxZalongMaxL - dem), '.k')
chiRange = np.linspace(0,np.nanmax(intA))
plt.plot(chiRange,chiRange*maxChiSlope,'-r',label='Ksn = %s' % maxKsn)
plt.xlabel(r'Downstream integrated $A^{m/n}, \xi, [m]$')
plt.ylabel('Relief [m]')
plt.legend(loc='upper left')
plt.ylim((0,np.nanmax((maxZalongMaxL - dem))))
# plt.savefig(outPath+'ReliefVsXi.eps')



plt.figure()
plt.title('Relief along longest flow path, [m]')
plt.imshow((maxZalongMaxL - dem))
plt.colorbar()

plt.figure()
plt.title(r'Downstream integrated $\frac{A_0}{A}^{m/n}, \xi, [m]$')
plt.imshow(intA)
plt.colorbar()

plt.show()