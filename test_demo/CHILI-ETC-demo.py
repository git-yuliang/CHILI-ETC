# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:43:05 2022

@author: DELL
"""

import matplotlib.pyplot as plt
from chili_etc.sp.chili_config import build_default_calc
from chili_etc.sp.chili_perform_calculation import perform_calculation

chili_config = build_default_calc()
chili_config['obst'] = 300
chili_config['repn'] = 3
chili_config['source']['normalization']['value'] = 18.0
chili_config['source']['spectrum']['name'] = 'SFgal_texp_FeH0_tau5_Ew10_AGN1.fits'
report = perform_calculation(chili_config)
print(report.snr)

plt.plot(report.snr)