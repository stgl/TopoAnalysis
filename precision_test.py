# precision test

import gdal
import numpy as np

# Read in files
fnASCII = 'malta_acc_15s.asc'
fnTIFF = 'malta_acc_15s.tif'

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

# Check precision
print('Element: %d, %d' % (r, c))
print('TIFF datatype: ' + typeTIFF)
print('ASCII datatype: ' + typeASCII)

print('TIFF value: %.20f ' % accTIFF[r][c])
print('ASCII value: %.20f '% accASCII[r][c])

print('Error: %.20f ' % (accTIFF[r][c] - accASCII[r][c]))
