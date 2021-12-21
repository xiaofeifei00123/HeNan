#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取EC细网格数据
-----------------------------------------
Time             :2021/09/09 11:41:05
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import meteva.base as meb
import numpy as np
import os
import pandas as pd
# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/EC_thin_grid/21071508.003'
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/grapse_3km/21071608.013'
grd = meb.read_griddata_from_micaps4(flnm)
grd
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station/20210720200000.000'
# station = meb.read_stadata_from_micaps3(flnm)
# grd.plot()
# meb.tool.plot_tools.scatter_sta(station)
# station
# grd
# grid1 = meb.grid([108,117,0.25],[31,38,0.25])
# grd2 = meb.interp_sg_idw(station, grid1)
# grd2.astype('float32')

# %%
def get_rain_era():
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/EC_thin_grid'
    fl_list = os.popen('ls {}/21071920*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    dds_list = []
    for fl in fl_list:
        print(fl)
        da = meb.read_griddata_from_micaps4(fl)
        dds_list.append(da)
    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate


dd = get_rain_era()
dd
# %%
# dd.dtime
tt = pd.date_range(start='2021-07-19 23', end='2021-07-22 20', freq='3h')
# tt.shape
dd.coords['time'] = tt
dd.squeeze().to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_era.nc')
# dd.max()
# da = xr.open_dataarray(flnm)
# da.dtime
# pd.date_range(start='2016-07-20')

# %%



# %%
def read_one_station(flnm):
    station = meb.read_stadata_from_micaps3(flnm)
    grid1 = meb.grid([90,130,0.25],[28,50,0.25])
    grd2 = meb.interp_sg_idw(station, grid1)
    da_return = grd2.astype('float32')
    return da_return

def get_rain_station(path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station/'):
    path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station'
    fl_list = os.popen('ls {}/2021*.000'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    dds_list = []
    for fl in fl_list:
        da = read_one_station(fl)
        dds_list.append(da)
    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate

c = get_rain_station()

# %%
# dd.squeeze().to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_era.nc')
c.squeeze().to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
# c.to_netcdf('')








# %%
# grd2.coords
# aa = grd2.values.astype('float32')
# grd2.plot()
# grd2.squeeze().plot()
# meb.plot_tools.contourf_2d_grid(grd2)
# print(grid1)
# ds = xr.open_dataset(flnm)
# ds = xr.open_dataset(flnm, engine='cfgrib',)
#                     # backend_kwargs={'filter_by_keys':
#                         # {'typeOfLevel': 'isobaricInhPa'}})
# ds

# %%
# ds['tp'].sel(expver=1)

