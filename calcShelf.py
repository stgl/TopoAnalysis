import dem as d
import numpy as np
import scipy.ndimage.morphology as morph

grid_name = 'gebco_bath_near'
slope_threshold_degrees = 1.0
base_elevation = -300.0

dem = d.Elevation.load(grid_name)

mask1 = d.Mask()
s = d.GeographicMaxSlope(elevation = dem)
threshold_slope = np.tan(slope_threshold_degrees * np.pi / 180.0)
mask1._copy_info_from_grid(dem, True)
mask1.dtype = np.uint8
mask1._griddata = mask1._griddata.astype(np.uint8)
mask1._griddata[np.isnan(dem._griddata)] = 1
mask1._griddata = mask1._griddata + morph.binary_dilation(mask1._griddata)
coast_mask = d.BaseSpatialGrid()
coast_mask._copy_info_from_grid(dem, True)
coast_mask.dtype = np.uint8
coast_mask._griddata = coast_mask._griddata.astype(np.uint8)
coast_mask._griddata[(mask1._griddata == 1) & (dem._griddata >= base_elevation)] = 1
i = np.where(coast_mask._griddata == 1)
rc = zip(i[0].tolist(),i[1].tolist())
outlets = coast_mask._rowscols_to_xy(rc)
mask = d.Mask()
mask._copy_info_from_grid(dem, True)
i = np.where(s._griddata <= threshold_slope)
mask._griddata[i] = 1
shelf = d.PriorityFillGrid(mask = mask, outlets = outlets)
shelf._griddata = morph.binary_fill_holes(shelf._griddata)
        

shelf.save(grid_name + '_classified')

