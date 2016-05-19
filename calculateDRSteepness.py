#! /usr/bin/env python2.7

import denudationRateAnalysis as dra
data = dra.read_csv('portengadata.csv')
del(data[0])
ksn_vec, area_vec = dra.calculate_ksn_for_data(data)
