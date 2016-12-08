import dem as d
base_path = '/data/sesfs/scarp1/GlobalSteepness/Datasets/Grids/'

prefixes = ['au', 'ca', 'eu', 'na', 'sa']

for prefix in prefixes:
    area = d.GeographicArea(gdal_filename = base_path + prefix + '/' + prefix + '_acc_15s')
    print('Successfully read: ' + prefix)
    d8 = d.FlowDirectionD8.load(prefix + '_flow_direction')
    idx = area.sort(reverse = False)
    print('Sorting done.')
    area = d.GeographicArea(flow_direction = d8, sorted_indexes = idx)
    print('Area recomputed.')
    area.save(prefix + "_area")
    print('File saved.')
