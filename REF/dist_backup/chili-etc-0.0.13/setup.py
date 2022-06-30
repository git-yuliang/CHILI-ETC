# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 16:02:45 2022

@author: DELL
"""

from setuptools import setup, find_packages

setup(
    name="chili-etc",
    version="0.0.13",
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
    package_dir={"": "src"},
# =============================================================================
#     # 安装过程中，需要安装的静态文件，如配置文件、service文件、图片等
#     data_files=[
#         ('src/referdata'),
#                ],
# =============================================================================
    include_package_data = True,
    include_package = {"": "src/refdata/"},
    python_requires='>=3.7',
)

