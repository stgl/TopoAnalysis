import dem as d
from sys.stdout import flush 

prefixes = ['af', 'as', 'au', 'ca', 'eu', 'na', 'sa']

for prefix in prefixes:
    
    print('Started ' + prefix)
    flush()
    area = d.GeographicArea.load(prefix + '_area')
    idx = area.sort(reverse = False)
    area = None
    d8 = d.FlowDirectionD8.load(prefix + '_flow_direction')
    flow_length = d.GeographicFlowLength(flow_direction = d8, sorted_indexes = idx)
    flow_length.save(prefix + '_flow_length')
    print('Done with ' + prefix)
    flush()
    
