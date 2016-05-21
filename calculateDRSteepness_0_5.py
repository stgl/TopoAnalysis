#! /usr/bin/env python2.7

import denudationRateAnalysis as dra
import numpy as np

data = dra.read_csv('portengadata.csv')
del(data[0])
ksn_vec, area_vec = dra.calculate_ksn_for_data(data)

np.savez_compressed('ksn_area_data.npz', ksn_vec = ksn_vec, area_vec = area_vec)
