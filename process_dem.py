import dem as d

def process_dem(dem_name, EPSGcode, folder_name = '.'):

    full_prefix = folder_name + '/' + dem_name
    ascii_name = full_prefix + '.txt'
    elevation = d.Elevation(ai_ascii_filename = ascii_name, EPSGprojectionCode = EPSGcode)
    elevation.save(full_prefix + '_elevation')

  
    filled_elevation = d.FilledElevation(elevation = elevation)
    filled_elevation.save(full_prefix + '_filled')
    d8 = d.FlowDirectionD8(flooded_dem = filled_elevation)
    d8.save(full_prefix + '_d8')
    area = d.Area(flow_direction = d8)
    area.save(full_prefix + '_area')
    flow_length = d.FlowLength(flow_direction = d8)
    flow_length.save(full_prefix + '_length')

def calc_ks_and_associated_grids(dem_name, Ao, theta, v, folder_name = '.', use_mask = False):

    full_prefix = folder_name + '/' + dem_name
    area = d.Area.load(full_prefix + '_area')
    d8 = d.FlowDirectionD8.load(full_prefix + '_d8')
    length = d.FlowLength.load(full_prefix + '_length')
    filled = d.FilledElevation.load(full_prefix + '_filled')
    d8.sort()
    d8._sort_indexes = filled.sort(reverse = False, force = True)
    elevation = d.Elevation.load(full_prefix + '_elevation')
    ksi = d.Ksi(area = area, flow_direction = d8, theta = theta, Ao = Ao, flow_length = length)
    relief = d.ScaledRelief(flow_direction = d8, elevation = elevation, flow_length = length, Ao = Ao, theta = theta, area=area)
    ksi.save(full_prefix + '_' + str(Ao) + '_' + str(theta) + '_ksi')
    relief.save(full_prefix + '_' + str(Ao) + '_' + str(theta) + '_relief')
    if use_mask:
        try:
            mask = d.BaseSpatialGrid.load(full_prefix + '_mask')
        except:
            mask = None
    else:
        mask = None
        
    from denudationRateAnalysis import calculate_ks_for_sample
    return calculate_ks_for_sample(v, d8, ksi, relief, area, Ao = Ao, mask = mask)


    
  

