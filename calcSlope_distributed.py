import dem as d
import numpy as np
import scipy.ndimage.morphology as morph

grid_name = 'gebco_bath_near'
shelf_grid_name = 'gebco_bath_near_classified'

slope_threshold_degrees = 1.0

dem = d.Elevation.load(grid_name)
dem.dtype = np.float32
dem._griddata = dem._griddata.astype(np.float32)

shelf_mask = d.Mask.load(shelf_grid_name)

slope_mask = d.Mask()
slope_mask._copy_info_from_grid(dem, True)

# Need to implement tiling scheme:

x_beg = range(-180, 160, 20)
y_beg = range(-90, 70, 20)

for x_s in x_beg:
    for y_s in y_beg:
        extent = [x_s, x_s+20, y_s, y_s+20]
        coords = ( (x_s, y_s), (x_s+20, y_s+20))
        rc_bounds = dem._xy_to_rowscols(coords)
        this_dem = dem.clip_to_extent(extent)
        this_shelf = shelf_mask.clip_to_extent(extent)
        shelf_edge = d.BaseSpatialGrid()
        shelf_edge._copy_info_from_grid(this_shelf, True)
        i = np.where(this_shelf._griddata == 1)
        shelf_edge._griddata[i] = 1
        shelf_edge._griddata = shelf_edge._griddata + morph.binary_dilation(this_shelf._griddata).astype(np.float64)
        shelf_edge._griddata = shelf_edge._griddata + morph.binary_erosion(this_shelf._griddata).astype(np.float64)
        i = np.where((shelf_edge._griddata == 2) & (this_dem._griddata <= -50))
        rc = zip(i[0].tolist(),i[1].tolist())
        outlets = shelf_edge._rowscols_to_xy(rc)
        
        s = d.GeographicMaxSlope(elevation = this_dem)
        s.dtype = np.float32
        s._griddata = s._griddata.astype(np.float32)
        threshold_slope = np.tan(slope_threshold_degrees * np.pi / 180.0)
        
        mask = d.Mask()
        mask._copy_info_from_grid(this_dem, True)
        mask._griddata = (s._griddata >= threshold_slope)
        mask._griddata = morph.binary_fill_holes(mask._griddata)

        slope = d.PriorityFillGrid(mask = mask, outlets = outlets)
        
        slope._griddata[(this_dem._griddata >-50) & (np.isnan(this_dem._griddata))] = 0
        slope._griddata = morph.binary_fill_holes(slope._griddata)
                
        slope_mask._griddata[rc_bounds[1][0]:rc_bounds[0][0],rc_bounds[0][1]:rc_bounds[1][1]] = slope._griddata
        

shelf_mask.save(grid_name + '_slope_classified')

