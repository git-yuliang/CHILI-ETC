# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 10:33:40 2022

@author: Liang Yu
This work is modified on the basis of previous work(by Lin Lin@SHAO: https://ifs-etc.readthedocs.io/en/latest/quickstart.html) 

This code is used for test the chili-etc codes.
by YuLiang 
yuliang@shao.ac.cn
"""

import matplotlib.pyplot as plt
from astropy.io import fits
from chili_etc.sp.chili_config import build_default_calc
from chili_etc.sp.chili_perform_calculation import perform_calculation

chili_config = build_default_calc()
chili_config['obst'] = 300
chili_config['repn'] = 3
chili_config['source']['normalization']['value'] = 18.0
chili_config['source']['spectrum']['name'] = 'SFgal_texp_FeH0_tau5_Ew10_AGN1.fits'
report = perform_calculation(chili_config)
print(report.snr)


hdu = fits.open("../refdata/sed/SFgal_texp_FeH0_tau5_Ew10_AGN1.fits")
print(hdu[1].header)
imagn = hdu[1].data
print(imagn.shape)
template_wave = hdu[1].data['wavelength']
template_flux = hdu[1].data['flux'] * 1e-12

plt.figure(figsize=[16,7])
plt.subplot(121)
plt.plot(template_wave,template_flux)
plt.title('TEMP_FLUX')
plt.subplot(122)
plt.plot(report.snr)
plt.title('report.SNR')