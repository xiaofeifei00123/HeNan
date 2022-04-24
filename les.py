#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
3d理想试验的图
-----------------------------------------
Time             :2022/04/15 08:50:05
Author           :Forxd
Version          :1.0
'''
# %%
import os
import xarray as xr
from Draw.Draw_lunwen.draw_distance_height_cross import Draw
import matplotlib.pyplot as plt
import netCDF4 as nc
import wrf
import numpy as np
import pandas as pd

# %%
def get_les():
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/3dles/wrfout_d01_2021-07-20_00:00:00'
    wrfnc = nc.Dataset(path)
    x = 75
    y = 75
    p = wrf.getvar(wrfnc, 'pres', units='hpa')[:,x,y]
    u = wrf.getvar(wrfnc, 'ua', units='m/s')[:,x,y]
    v = wrf.getvar(wrfnc, 'va', units='m/s')[:,x,y]
    hagl = wrf.getvar(wrfnc, 'height', units='m')[:,x,y]  # 先纬度后经度
    ds = xr.merge([u, v, p, hagl])
    # ds = xr.merge([t,td, u, v, wind_speed, wind_angle, p, hagl])
    dds = ds.set_coords(['height', 'pressure'])
    ds_return = dds.swap_dims({'bottom_top':'pressure'})
    ds_return = ds_return.rename({'ua':'u', 'va':'v'})
    return ds_return
ds_les = get_les()
#%%
ds_les.pressure
# dds_list = []
# flnm2 = fl_list[12]
# wrfnc = nc.Dataset(flnm)
# wrfnc
# flnm2



# %%
path = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/3dles/'
tt = pd.date_range('2021-07-20 0000', '2021-07-20 1200', freq='1H')
# tt
fl_list = []
for t in tt:
    fl = 'wrfout_d01_'+t.strftime('%Y-%m-%d_%H:%M:%S')
    flnm = path+fl
    print(flnm)
    fl_list.append(flnm)




# %%
ds = xr.open_dataset(flnm)
# da = ds['W'].isel(Time=0).sel(south_north=75)
da = ds['W'].isel(Time=0).sel(south_north=75)
# x = da.west_east
# y = da.bottom_top_stag
# y
ds

# %%
def sounding_1station_1time(flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d03_2021-07-19_00:00:00'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d04_2021-07-19_00:00:00'
    print(flnm[-19:])
    wrfnc = nc.Dataset(flnm)
    # lat, lon = get_centroid()

    # lat = self.sta_dic['lat']
    # lon = self.sta_dic['lon']

    # x,y = wrf.ll_to_xy(wrfnc, 34.71, 113.66)
    # x,y = wrf.ll_to_xy(wrfnc, lat, lon)
    x = 75
    y = 75

    ## 高度， 各层海拔高度, 单位m, 和探空资料保持一致
    hagl = wrf.getvar(wrfnc, 'height', units='m')[:,x,y]  # 先纬度后经度
    # print(hagl)
    p = wrf.getvar(wrfnc, 'pres', units='hpa')[:,x,y]
    u = wrf.getvar(wrfnc, 'ua', units='m/s')[:,x,y]
    v = wrf.getvar(wrfnc, 'va', units='m/s')[:,x,y]
    w = wrf.getvar(wrfnc, 'wa', units='m/s')[:,x,y]
    
    # ## 根据u,v风计算风向和风速
    deg = 180.0/np.pi # 角度和弧度之间的转换
    rad = np.pi/180.0

    # wind_speed = xr.ufuncs.sqrt(u**2+v**2)
    wind_speed = np.sqrt(u**2+v**2)
    wind_speed.name = 'wind_speed'
    wind_angle = 180.0+np.arctan2(u, v)*deg
    wind_angle.name = 'wind_angle'
    u = u.rename('u')
    v = v.rename('v')
    w = v.rename('w')
    

    t = wrf.getvar(wrfnc, 'temp', units='degC')[:,x,y]
    td = wrf.getvar(wrfnc, 'td', units='degC')[:,x,y]
    ds = xr.merge([t,td, u, v, p, hagl])
    ds = xr.merge([t,td, u, v,w, wind_speed, wind_angle, p, hagl])
    dds = ds.set_coords(['height', 'pressure'])
    ds_return = dds.swap_dims({'bottom_top':'pressure'})
    # print(ds_return)
    return ds_return
# ds = sounding_1station_1time(flnm)
# ds
# %%
# ds['P_TOP']
# %%
path = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/3dles/'
tt = pd.date_range('2021-07-20 0000', '2021-07-20 1200', freq='1H')
fl_list = []
for t in tt:
    fl = 'wrfout_d01_'+t.strftime('%Y-%m-%d_%H:%M:%S')
    flnm = path+fl
    # print(flnm)
    fl_list.append(flnm)

ds_list = []
for fl in fl_list:
    # print(fl)
#     # ds = 
    ds = sounding_1station_1time(fl)
    ds_list.append(ds)
dds = xr.concat(ds_list, dim='Time')
dds

# %%
# dds['w']
# dds['w'].interp(pressure=[1000,900,800])

# flnm == flnm2
# print(flnm2)
# print(flnm)
# xr.open_dataset(flnm)
# flnm
# hagl = wrf.getvar(wrfnc, 'height', units='m')[:,75,75]  # 先纬度后经度
# hagl
# ds.Time.values.dt.strftime('%d%h%m')
# %%
# ds['W']
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/3dles/wrfout_d01_2021-07-19_12:00:00'
# ds = xr.open_dataset(flnm)
# da = ds['W'].isel(Time=0).sel(south_north=75)
da = ds['W'].isel(Time=0).sel(south_north=75)
x = da.west_east
y = da.bottom_top_stag
y
dr = Draw()
cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 7*cm), dpi=300)
ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])

cf = dr.draw_contourf_no_terrain(ax,da,x,y)

cb = fig.colorbar(
    cf,
    # cf_list[i],
    ax=ax,
    orientation='horizontal',
    ticks=dr.colorlevel[1:-1],
    fraction = 0.05,  # 色标大小,相对于原图的大小
    pad=0.1,  #  色标和子图间距离
)
cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小
# %%
# fig.savefig('0720')
