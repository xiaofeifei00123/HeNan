#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
寻找降水质心
-----------------------------------------
Time             :2021/11/09 13:19:10
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m/rain.nc'
ds = xr.open_dataset(flnm)
da = ds.max(dim=['south_north', 'west_east']).to_array()
da.plot()






