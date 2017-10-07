def relocate_p_and_b(data, pixel_radius = 5):
       
    import dem as d
        
    prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']
    
    data_processed = dict()
        
    for sample_name, lat, lon, dr, dr_sig, a in data:
        
        data_processed[sample_name] = {'lat': float(lat),
                                       'lon': float(lon),
                                       'dr': float(dr),
                                       'dr_sig': float(dr_sig),
                                       'area_measured': float(a)*1.0E6}
        
    for prefix in prefixes:
        print('Loading prefix: ' + prefix)
        area = d.GeographicArea.load(prefix + '_area')
        print('Done loading prefix: ' + prefix)
        
        for sample in data_processed:
            xy = (data_processed[sample]['lon'], data_processed[sample]['lat'])
            rc = area._xy_to_rowscols((xy, ))[0]
            if rc[0] is not None and rc[1] is not None:
                area_measured = data_processed[sample]['area_measured']
                locations_snap = area.snap_locations_to_closest_value((xy,), (area_measured,), pixel_radius = pixel_radius)
                data_processed[sample]['lat_adjusted'] = locations_snap[0][1]
                data_processed[sample]['lon_adjusted'] = locations_snap[0][0]
                locations_snap_indexes = area._xy_to_rowscols(locations_snap)
                data_processed[sample]['area_dem'] = area[locations_snap_indexes[0][0], locations_snap_indexes[0][1]]
                data_processed[sample]['prefix'] = prefix
    return data_processed

import denudationRateAnalysis as dra

data = dra.read_csv('portengadata.csv')
del(data[0])

data_processed = relocate_p_and_b(data)

final_list = list()

for sample in data_processed:
    
    final_list.append([sample, data_processed[sample].get('lat'), data_processed[sample].get('lon'), data_processed[sample].get('dr'), data_processed[sample].get('dr_sig'), data_processed[sample].get('area_measured') / 1.0E6, data_processed[sample].get('lat_adjusted',-9999), data_processed[sample].get('lon_adjusted', -9999), data_processed[sample].get('area_dem', -9999) / 1.0E6, data_processed[sample].get('prefix','Not Found')])
    
import csv
with open('portenta_relocated.csv', 'wb') as f:
    writer = csv.writer(f, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for row in final_list:
        writer.writerow(row)
        
        
    