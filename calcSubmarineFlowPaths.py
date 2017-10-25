#! /usr/bin/env python2.7

import dem as d
import scipy.ndimage.morphology as morph
import numpy as np

prefix = 'ca'

threshold = 5

dem = d.Elevation.load(prefix + '_elevation')
mask = d.Mask.load(prefix + '_mask')

mask2 = d.Mask()
mask2._copy_info_from_grid(mask)
mask2._griddata = mask2._griddata = morph.binary_dilation(mask._griddata).astype(np.float) + mask._griddata.astype(np.float)
mask2._griddata[mask2._griddata != 1] = 0  
mask2._griddata[dem._griddata >= 0] = 0
mask2._griddata[np.isnan(dem._griddata)] = 0

i = np.where(mask2._griddata == 1)
rc = zip(i[0].tolist(),i[1].tolist())
outlets = dem._rowscols_to_xy(rc) 

discrete_flow_terminations = d.GeographicDiscreteFlowAccumulation(elevation = dem, outlets = outlets)

i = np.where(discrete_flow_terminations._griddata >= threshold)
rc = zip(i[0].tolist(),i[1].tolist())
outlets = dem._rowscols_to_xy(rc)
