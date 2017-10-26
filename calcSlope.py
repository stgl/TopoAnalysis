import dem as d
import numpy as np
import scipy.ndimage.morphology as morph

grid_name = 'gebco_bath_near'
shelf_grid_name = 'gebco_bath_near_classified'

slope_threshold_degrees = 1.0

dem = d.Elevation.load(grid_name)
shelf_mask = d.Mask.load(shelf_grid_name)

shelf_edge = d.BaseSpatialGrid()
shelf_edge._copy_info_from_grid(shelf_mask, True)
i = np.where(shelf_mask._griddata == 1)
shelf_edge._griddata[i] = 1
shelf_edge._griddata = shelf_edge._griddata + morph.binary_dilation(shelf_mask._griddata).astype(np.float64)
shelf_edge._griddata = shelf_edge._griddata + morph.binary_erosion(shelf_mask._griddata).astype(np.float64)
i = np.where((shelf_edge._griddata == 2) & (dem._griddata <= -50))
rc = zip(i[0].tolist(),i[1].tolist())
outlets = shelf_edge._rowscols_to_xy(rc)
        
s = d.GeographicMaxSlope(elevation = dem)
threshold_slope = np.tan(slope_threshold_degrees * np.pi / 180.0)
        
mask = d.Mask()
mask._copy_info_from_grid(dem, True)
mask._griddata = (s._griddata >= threshold_slope)
mask._griddata = morph.binary_fill_holes(mask._griddata)

slope = d.PriorityFillGrid(mask = mask, outlets = outlets)
        
slope._griddata[(dem._griddata >-50) & (np.isnan(dem._griddata))] = 0
slope._griddata = morph.binary_fill_holes(slope._griddata)
                        

slope.save(grid_name + '_slope_classified')

