import dem as d
import numpy as np
import scipy.ndimage.morphology as morph

grid_name = 'ca_bath_near'
slope_grid_name = 'ca_bath_near_slope'

dem = d.Elevation.load(grid_name)
slope_mask = d.Mask.load(slope_grid_name)

base_of_slope_mask = d.Mask()
base_of_slope_mask._copy_info_from_grid(slope_mask, True)

base_of_slope_mask._griddata = slope_mask._griddata.astype(np.float64) + morph.binary_erosion(slope_mask._griddata).astype(np.float64)
i = np.where((base_of_slope_mask._griddata == 1) & (dem._griddata <= -600))

rc = zip(i[0].tolist(),i[1].tolist())
outlets = base_of_slope_mask._rowscols_to_xy(rc)

dem._griddata[slope_mask._griddata == 0] = np.NaN

filled_dem = d.FilledElevation(elevation = dem, outlets = outlets, mask = slope_mask)
fd = d.FlowDirectionD8(flooded_dem = filled_dem)
area = d.GeographicArea(flow_direction = fd)

filled_dem.save(grid_name + '_filled')
fd.save(grid_name + '_flow_direction')
area.save(grid_name + '_area')

