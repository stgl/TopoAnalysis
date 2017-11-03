import dem as d
import numpy as np
import scipy.ndimage.morphology as morph

grid_name = 'gebco'
slope_threshold_degrees = 1.0
base_elevation = -300.0
maximum_pit_depth = 200


dem = d.Elevation.load(grid_name + '_dem')

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
        
shelf.save(grid_name + '_shelf')



threshold = 5*1E6

mask = shelf

mask2 = d.Mask()
mask2._copy_info_from_grid(mask)
mask2._griddata = mask2._griddata = morph.binary_dilation(mask._griddata).astype(np.float) + mask._griddata.astype(np.float)
mask2._griddata[mask2._griddata != 1] = 0  
mask2._griddata[dem._griddata >= 0] = 0
mask2._griddata[np.isnan(dem._griddata)] = 0

i = np.where(mask2._griddata == 1)
rc = zip(i[0].tolist(),i[1].tolist())
outlets = dem._rowscols_to_xy(rc) 

discrete_flow_terminations = d.GeographicDiscreteFlowAccumulation(elevation = dem, outlets = outlets, terminations_only = True, display_output = True)

i = np.where(discrete_flow_terminations._griddata >= threshold)
rc = zip(i[0].tolist(),i[1].tolist())
outlets = dem._rowscols_to_xy(rc)

i = np.where(shelf._griddata == 1)
dem._griddata[i] = np.NaN

filled_dem = d.FilledElevation(elevation = dem, outlets = outlets, display_output = True, maximum_pit_depth = maximum_pit_depth)
fd = d.FlowDirectionD8(flooded_dem = filled_dem)
area = d.GeographicArea(flow_direction = fd)
logarea = d.LogArea(area = area)

filled_dem.save(grid_name + '_filled')
fd.save(grid_name + '_flow_direction')
area.save(grid_name + '_area')
logarea.save(grid_name + '_logarea')
