import dem as d

prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']

for prefix in prefixes:
    area = d.GeographicArea.load(prefix + '_area')
    d8 = d.FlowDirectionD8.load(prefix + '_flow_direction')
    idx = area.sort(reverse = False)
    area = d.GeographicArea(flow_direction = d8, sorted_indexes = idx)
    area.save(prefix + "_area")
    
