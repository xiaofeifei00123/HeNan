#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取高空填图数据
插值成格点
聚合成一个文件
和wrfout插值出的格点数据作差
-----------------------------------------
Time             :2021/10/06 16:59:52
Author           :Forxd
Version          :1.0
'''

# %%
from time import strftime
from cartopy.crs import Projection
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import xarray as xr
import meteva.base as meb
import metpy.interpolate as interp
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xesmf as xe


# %%
def get_plot(dic):
    """读取micaps 2类数据，高空填图数据"""
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/'
    # flnm = path+dic['level']+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/500/20210720080000.000'
    df = read_micaps_2(flnm)
    sta = meb.sta_data(df, columns = [
                'id', 'lon', 'lat', 'alt', 'grade', 'height', 'temperature', 'dewpoint',
                'wind_angle', 'wind_speed', 'time', 'level'])

    return sta


def interp_metpy(sta, ):
    """
    站点插值到格点
    反距离权重插值

    Args:
        sta (DataFrame): [lon,lat,height]

    Returns:
        [type]: [description]
    """
    h = sta['height']
    lon = sta['lon']
    lat = sta['lat']

    area = {
        'lon1':107,
        'lon2':135,
        'lat1':20,
        'lat2':40,
        'interval':0.1,
    }

    ds_regrid = xe.util.grid_2d(area['lon1'], area['lon2'], area['interval'], area['lat1'], area['lat2'], area['interval'])
    mx = ds_regrid['lon'].values
    my = ds_regrid['lat'].values

    z = interp.inverse_distance_to_grid(lon, lat, h, mx, my, r=5, min_neighbors=1)

    ## 重新设置coords属性
    lon1 = mx[0,:]
    lat1 = my[:,0]

    da = xr.DataArray(
        z,
        coords={
            'lon':lon1,
            'lat':lat1,
        },
        dims=['lat','lon']
    )
    return da


t = pd.Timestamp('2021-07-19 2000')
dic_t = {
    'var':'temp',
    'level':500,
    'time':t
}
sta = get_plot(dic_t)
# sta
c = interp_metpy(sta, )
# a



# ttt = pd.date_range(start='2021-07-18 08', end='2021-07-20 20', freq='12H')
# for level in ['200', '500', '850']:
#     for t in ttt:
#         draw_all(t, level)

# TODO  现在插值一个时次，一个层次， 一个变量已经写好了
# TODO  下一步是对多个时次, 多个层次，多个变量进行插值聚合













# %%
# c
# c.shape
# b[:,0]
# %%
# ds_regrid = xe.util.grid_2d(area['lon1'], area['lon2'], area['interval'], area['lat1'], area['lat2'], area['interval'])
# ds_regrid['lon'].values

# a[1]
# %%
# a
fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d02_latlon.nc'
ds = xr.open_dataset(fl_wrf)
# %%
da = ds.isel(time=0).isel(pressure=0)['geopt'].squeeze()
da
# %%
# c.dims
# da.dims
# dd = da-c
# dd.min()
# da
# xr.DataArray(
#     c,
#     coords={
#         'lat':my
#     }
# )