from dem import FilledElevation
def processAll(prefix_name, Ao, theta, base_name = '.'):
    
    from dem import Elevation, FlowDirectionD8, GeographicArea, Area, GeographicFlowLength, GeographicKsi, ScaledRelief
    
    elevation_name = base_name + "/" + prefix_name + "_dem_15s"
    area_name = base_name + "/" + prefix_name + "_acc_15s"
    d8_name = base_name + "/" + prefix_name + "_dir_15s"
    
    elevation = Elevation(gdal_filename = elevation_name)
    elevation.save(prefix_name + "_elevation")
    area = Area(gdal_filename = area_name)
    d8 = FlowDirectionD8(gdal_filename = d8_name)
    d8.save(prefix_name + "_flow_direction")
    
    idx = area.sort(reverse = False)
    area = GeographicArea(flow_direction = d8, sorted_indexes = idx)
    area.save(prefix_name + "_area")
    flow_length = GeographicFlowLength(flow_direction = d8, sorted_indexes = idx)
    flow_length.save(prefix_name + "_flow_length")
    ksi = GeographicKsi(area = area, flow_direction = d8, theta = theta, Ao = Ao, flow_length = flow_length, sorted_indexes = idx)
    ksi.save(prefix_name + "_ksi_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))
    relief = ScaledRelief(flow_direction = d8, elevation = elevation, flow_length = flow_length, Ao = Ao, theta = theta, sorted_indexes = idx)
    relief.save(prefix_name + "_relief_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))

def processAllUTM(prefix_name, EPSGprojectionCode, Ao, theta, base_name = '.'):
    
    from dem import Elevation, FlowDirectionD8, Area, FlowLength, Ksi, ScaledRelief
    
    full_path_without_suffix = base_name + "/" + prefix_name
    elevation_unfilled_ascii_filename = full_path_without_suffix + ".txt"
    elevation = Elevation(ai_ascii_filename = elevation_unfilled_ascii_filename, EPSGprojectionCode= EPSGprojectionCode)
    elevation.save(full_path_without_suffix + "_elevation")
    filled = FilledElevation(elevation = elevation)
    filled.save(full_path_without_suffix + "_filled")
    d8 = FlowDirectionD8(flooded_dem = filled)
    d8.save(full_path_without_suffix + "_flow_direction")
    area = Area(flow_direction = d8)
    area.save(full_path_without_suffix + "_area")
    
    flow_length = FlowLength(flow_direction = d8)
    flow_length.save(full_path_without_suffix + "_flow_length")
    ksi = Ksi(area = area, flow_direction = d8, theta = theta, Ao = Ao, flow_length = flow_length)
    ksi.save(full_path_without_suffix + "_ksi_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))
    relief = ScaledRelief(flow_direction = d8, elevation = elevation, flow_length = flow_length, Ao = Ao, theta = theta)
    relief.save(full_path_without_suffix + "_relief_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))
        
def processForTheta(prefix_name, Ao, theta, base_name = '.'):
    
    from dem import FlowDirectionD8, GeographicArea, GeographicFlowLength, GeographicKsi, ScaledRelief, Elevation
    
    area_name = base_name + "/" + prefix_name + "_area"
    d8_name = base_name + "/" + prefix_name + "_flow_direction"
    flow_length_name = base_name + "/" + prefix_name + "_flow_length"
    elevation_name = base_name + "/" + prefix_name + "_elevation"
 
    area = GeographicArea.load(area_name)
    d8 = FlowDirectionD8.load(d8_name)
    flow_length = GeographicFlowLength.load(flow_length_name)
    elevation = Elevation.load(elevation_name) 
    idx = area.sort(reverse = False)
    ksi = GeographicKsi(area = area, flow_direction = d8, theta = theta, Ao = Ao, flow_length = flow_length, sorted_indexes = idx)
    ksi.save(prefix_name + "_ksi_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))
    relief = ScaledRelief(flow_direction = d8, elevation = elevation, flow_length = flow_length, Ao = Ao, theta = theta, sorted_indexes = idx, area = area) 
    relief.save(prefix_name + "_relief_" + str(Ao).replace('.','_') + "_" + str(theta).replace('.','_'))
    
def plotGrids(x_grid, y_grid, plot_string, **kwargs):
    
    import numpy as np
    
    x_vec = np.ndarray.flatten(x_grid._griddata)
    y_vec = np.ndarray.flatten(y_grid._griddata)
    
    ind = np.where(y_vec >= 0)
    x_vec = x_vec[ind]
    y_vec = y_vec[ind]
    
    from matplotlib import pyplot as plt
    
    plt.plot(x_vec, y_vec, plot_string, **kwargs)

    return x_vec, y_vec
    
def create_density(x_grid, y_grid, x_boundaries, y_boundaries):
    
    import numpy as np
    x_vec = np.ndarray.flatten(x_grid._griddata)
    y_vec = np.ndarray.flatten(y_grid._griddata)
    
    #x_bin_boundaries = (x_centers[1:] + x_centers[0:-1]) / 2.0
    #print(x_bin_boundaries)
    #print((3.0*x_centers[-1]/2.0 - x_bin_boundaries[-1]/2.0))
    #x_bin_boundaries = np.concatenate((x_bin_boundaries[0] - x_centers[0], x_bin_boundaries, (3.0*x_centers[-1]/2.0 - x_bin_boundaries[-1]/2.0)))
    
    #y_bin_boundaries = (y_centers[1:] + y_centers[0:-1]) / 2.0
    #y_bin_boundaries = np.concatenate((y_bin_boundaries[0] - y_centers[0], y_bin_boundaries, (3.0*y_centers[-1]/2.0 - y_bin_boundaries[-1]/2.0)))
    
    
    
    H, xedges, yedges = np.histogram2d(x_vec, y_vec, bins = (x_boundaries, y_boundaries))
    
    return H, xedges, yedges
    
    
