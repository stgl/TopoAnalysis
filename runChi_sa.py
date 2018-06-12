from denudationRateAnalysis import create_chi_grid_for_geographic_prefix as create_chi

prefix = 'sa' 

thetas = [0.4, 0.5, 0.6]
Ao = 1000000
basin_lengths = [50000, 100000, 200000, 400000]

create_chi(prefix, thetas, Ao, basin_lengths)

