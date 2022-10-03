#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
测试github
-----------------------------------------
Author           :Forxd
Version          :1.0
Time：2022/05/11/ 15:49
'''
# %%
import xarray as xr

# %%
flnm1 = '/mnt/zfm_18T/fengxiang/aa/2014_07_11-12/radar/ob_radar/data/REF_mosaic.nc'
ds = xr.open_dataset(flnm1)
ds

cc = ds.expand_dims(dim='model')
dd = cc.assign_coords({'model':['obs']})
dd



