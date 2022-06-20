# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 09:34:57 2022

@author: DELL
"""

from ifs_etc.etc1d.config import build_default_calc
from ifs_etc.etc1d.perform_calculation import perform_calculation

config = build_default_calc()
config['obst'] = 300
config['repn'] = 3
config['source']['normalization']['value'] = 18.0
config['source']['spectrum']['name'] = 'SFgal_texp_FeH0_tau5_Ew10_AGN1.fits'
report = perform_calculation(config)
print(report.snr)

#
config = build_default_calc()
config['targetsnr'] = 10
config['obst'] = 300
config['source']['normalization']['value'] = 18.0
config['source']['spectrum']['name'] = 'SFgal_texp_FeH0_tau1_Ewd.fits'
report = perform_calculation(config, calculation_mode='snr2exptime')
print(report.exptime)