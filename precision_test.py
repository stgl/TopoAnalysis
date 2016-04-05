# precision test

import gdal
import numpy as np

# Read in files
fnASCII = 'malta_area_15s.asc'
fnTIFF = 'malta_area_15s.tif'

#gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float32')

src = gdal.Open(fnASCII)
accASCII = src.ReadAsArray().astype(np.float64)
band = src.GetRasterBand(1)
typeASCII = gdal.GetDataTypeName(band.DataType)
src = None

src = gdal.Open(fnTIFF)
accTIFF = src.ReadAsArray().astype(np.float64)
band = src.GetRasterBand(1)
typeTIFF = gdal.GetDataTypeName(band.DataType)
src = None

# Find some nonzero cell
(m, n) = accASCII.shape

while True:
    r = np.random.randint(1, m) 
    c = np.random.randint(1, n) 
    if(accTIFF[r][c] != 0 and ~np.isnan(accTIFF[r][c])):
        break

# or perform quic kcheck with first element (Malta)
r = 0
c = 0
ref = 0.173059156637050

# Check precision
print('Element: %d, %d' % (r, c))
print('TIFF datatype: ' + typeTIFF)
print('ASCII datatype: ' + typeASCII) + '\n'

print('Actual value:\t%.20f ' % ref)
print('TIFF value:\t%.20f ' % accTIFF[r][c])
print('ASCII value:\t%.20f \n'% accASCII[r][c])

print('Error:\t\t%.20f ' % abs(ref - accASCII[r][c]))
