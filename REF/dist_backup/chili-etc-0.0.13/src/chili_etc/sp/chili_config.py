# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 14:41:12 2022

@author: DELL
This work is modified on the basis of previous work(by Lin Lin@SHAO: https://ifs-etc.readthedocs.io/en/latest/quickstart.html) 

This code is used for setting the chili-etc parameters.
by YuLiang 
yuliang@shao.ac.cn
"""

import os
import pandas as pd
import json

path = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
refdata = os.path.join(path, 'refdata/')

# dict_t : means dict_telescope
def get_telescope_config():
    dict_t = {}
    dict_t['diameter'] = 240                          # unit: cm?
    dict_t['obscure'] = 0.                            # 
    dict_t['coll_area'] = 31415.926535897932          # Unit: 1E4 * Pi ?
    return dict_t

# dict_i : means dict_instrument, need be updated by CHILI parameters
def get_instrument_config(readout_xbin=1, readout_ybin=1, gain=1, dark=0.017):
    dict_i = dict()
    # user parameter
    dict_i['readout_xbin'] = readout_xbin
    dict_i['readout_ybin'] = readout_ybin
    dict_i['gain'] = gain                           # e/ADU
    dict_i['dark'] = dark                           # e/s/pix
    dict_i['readout_noise'] = 4.0                   # e/pix
    dict_i['QE'] = 1.0                              # e/pix
    dict_i['efficiency_file'] = 'chili_IFU_throughput.dat'
    # hidden parameters
    dict_i['ccd_xsize'] = 1.755555                              # wavelenght dircection, A, delta_lambda_per_pixel
    dict_i['ccd_ysize'] = 0.1                                   # spatial direction, arcsec, spatial axis
    dict_i['extract_aper_width'] = 2 * dict_i['readout_ybin']   # binning in spatial dirction, extract spectrum with aperture of 2 pixels
    dict_i['spaxel_xsize'] = 0.2                                # arcsec
    dict_i['spaxel_ysize'] = dict_i['extract_aper_width'] * 0.1   #pixel index to unit of arcsec
    dict_i['fov_xsize'] = 6                                     # arcsec
    dict_i['fov_ysize'] = 6                                     # arcsec
    dict_i['wave_start'] = 3500
    dict_i['wave_end'] = 10000
    dict_i['wave_delta'] = dict_i['ccd_xsize'] * dict_i['readout_xbin']
    return dict_i


def get_throughput():

    """
    read IFS throughput file

    Returns:
        throughputwave: unit in Angstrom
        throughputvalue: unit in energy fraction or photon fraction ?
        
    Chili`s throughput data can be got by the observed flat data.

    """

    # throughput_file = os.path.join(os.getenv('SNPP_refdata'), 'csst', 'ifs', 'IFU_throughput.dat')
    # throughput = pd.read_csv(throughput_file,
    #                          sep='\s+', skiprows=1, header=None, names=['nm', 'q'])
    # lambdaq = throughput.nm.values * 10  # nm to A
    # qifs = throughput.q.values  # ; throughput of the whole system,
    # # ;assuming the total throughput cannot reach the theory value, 0.3 is the upper limit.
    # qifs[qifs >= 0.3] = 0.3
    # qinput = 1.0                    # the throughput of the telescope
    # qtot = qifs * qinput            # *qe ;qtot of CSST already includes the CCD efficiency


    throughput_file = os.path.join(refdata, 'chili', 'ifs', 'IFU_throughput.dat')
    cat = pd.read_csv(throughput_file,
                      comment='#', sep='\s+', header=None, names=['nm', 'q'])
    throughputwave = cat['nm'].values * 10  # nm to A
    throughputvalue = cat['q'].values       # energe fraction or photon fraction??

    return throughputwave, throughputvalue

# ?
def build_default_source():
    dict_s = {}      # dict_s means dict of the default source
    dict_s['spectrum'] = {'name': 'SFgal_texp_FeH0_tau5_Ew10_AGN1.fits', 'redshift': 0.0, 'ebv': 0.0, 'lines':[]}
    dict_s['normalization'] = {'value': 17.7, 'unit':'mag/arcsec^2', 'band': 'sdss_g'}

    # out_file = open(refdata + "source/config.json", "w")
    # json.dump(dict, out_file, indent=2)
    # out_file.close()

    return dict_s


# What's the meaning?
def build_default_calc():
    calc = {}
    calc['configuration'] = get_instrument_config()
    calc['source'] = build_default_source()
    calc['background'] = 'from_msc'
    calc['background_level'] = 'from_msc'
    calc['repn'] = 20.
    calc['obst'] = 300.
    calc['targetsnr'] = 10      # only be used in calculate_type = 'limitmag',
                                # in this case, normalization value is invalid
    return calc