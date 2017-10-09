import dem as d
import numpy as np
import scipy.ndimage.morphology as morph
from dem import BaseSpatialGrid
from numpy import int8

grid_name = 'gebco_bath_near'
slope_threshold_degrees = 1.0
base_elevation = -300.0

dem = d.Elevation.load(grid_name)
dem.dtype = np.float32
dem._griddata = dem._griddata.astype(np.float32)

shelf_mask = d.Mask()
shelf_mask._copy_info_from_grid(dem, True)

# Need to implement tiling scheme:

x_beg = range(-180, 160, 20)
y_beg = range(-90, 70, 20)

#x_beg = [-180]
#y_beg = [-90]

for x_s in x_beg:
    for y_s in y_beg:
        extent = [x_s, x_s+20, y_s, y_s+20]
        coords = ( (x_s, y_s), (x_s+20, y_s+20))
        rc_bounds = dem._xy_to_rowscols(coords)
        this_dem = dem.clip_to_extent(extent)
        mask1 = d.BaseSpatialGrid()
        mask2 = d.BaseSpatialGrid()
        s = d.GeographicMaxSlope(elevation = this_dem)
        s.dtype = np.float32
        s._griddata = s._griddata.astype(np.float32)
        threshold_slope = np.tan(slope_threshold_degrees * np.pi / 180.0)
        mask1._copy_info_from_grid(this_dem, True)
        mask1.dtype = np.uint8
        mask1._griddata = mask1._griddata.astype(np.uint8)
        mask1._griddata[np.isnan(this_dem._griddata)] = 1
        mask2._copy_info_from_grid(mask1, False)
        mask2._griddata = morph.binary_dilation(mask2._griddata)
        mask2._griddata = mask1._griddata + mask2._griddata
        coast_mask = d.BaseSpatialGrid()
        coast_mask._copy_info_from_grid(this_dem, True)
        coast_mask.dtype = np.uint8
        coast_mask._griddata = coast_mask._griddata.astype(np.uint8)
        coast_mask._griddata[(mask2._griddata == 1) & (this_dem._griddata >= base_elevation)] = 1
        i = np.where(coast_mask._griddata == 1)
        rc = zip(i[0].tolist(),i[1].tolist())
        outlets = coast_mask._rowscols_to_xy(rc)
        mask = d.Mask()
        mask._copy_info_from_grid(this_dem, True)
        i = np.where(s._griddata <= threshold_slope)
        mask._griddata[i] = 1
        shelf = d.PriorityFillGrid(mask = mask, outlets = outlets)
        shelf._griddata = morph.binary_fill_holes(shelf._griddata)
        shelf_mask._griddata[rc_bounds[1][0]:rc_bounds[0][0],rc_bounds[0][1]:rc_bounds[1][1]] = shelf._griddata
        

shelf_mask.save(grid_name + '_classified')

