# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 10:33:40 2022

@author: Liang Yu
This work is modified on the basis of previous work(by Lin Lin@SHAO: https://ifs-etc.readthedocs.io/en/latest/quickstart.html) 

This code is used for setting the chili-etc parameters.
by YuLiang 
yuliang@shao.ac.cn
"""

from setuptools import setup, find_packages

setup(
    name="chili-etc",
    version="0.1.1",
    description='exposure time calculator of CHILI',
    long_description=('CHILI-ETC is used for setting the CHILI exposure time parameters. by YuLiang yuliang@shao.ac.cn .\
                      This work is based on the work of the predecessors(by Lin Lin@SHAO: https://ifs-etc.readthedocs.io/en/latest/quickstart.html), \
                      and has been modified and completed on the basis of it.'),
    author="Yu Liang",
    author_email="yuliang@shao.ac.cn",
    url="https://github.com/git-yuliang/CHILI-ETC/",
    license="MIT Licence",
    packages=find_packages(where="src"),
    install_requires=['pandas>=1.3.3',
                      'numpy>=1.21.6',
                      'h5py>=2.8.0',
                      'einops>=0.3.2',
                      'matplotlib>=3.2.2',
                      'astropy>=4.2.1',
                      'scipy>=1.1.0',
                      'extinction>=0.4.0'],
    #package_dir={'': 'src'},
    #include_package_data=True,
    package = {'chili_etc'},
    package_dir={'chili_etc': 'src/chili_etc'},
    package_data={'chili_etc': ['refdata/sed/*.fits',
                                'refdata/source/*.json',
                                'refdata/csst/background/*.csv',
                                'refdata/csst/background/*.dat',
                                'refdata/csst/ifs/*.dat',
                                'refdata/csst/ifs/*.json',
                                'refdata/csst/telescope/*.json',
                                'refdata/chili/ifs/*.dat',
                                'refdata/normalization/filters/*.par',
                                'demo/*.ipynb',
                                'demo/*.py',
                                ]},
    #include_package_data=True,

    python_requires='>=3.7',
)

