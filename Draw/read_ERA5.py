#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2021/09/09 11:41:05
Author           :Forxd
Version          :1.0
'''

import xarray as xr
# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/rain.nc'
ds = xr.open_dataset(flnm)
ds

# %%
ds['tp'].sel(expver=1)

