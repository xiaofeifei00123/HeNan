#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画时间高度廓线, 值是区域平均值(33-34N,111.5-113E)
变量是w
-----------------------------------------
Time             :2022/01/05 17:31:02
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/ERA5/era5.pl.20210720.nc'
# ds = xr.open_dataset(flnm)
# ds

