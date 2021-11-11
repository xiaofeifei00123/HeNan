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

def get_centroid(flnm= '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc', t_start='2021-07-20 00', t_end='2021-07-20 12'):
    """获得一段时间降水的质心

    Args:
        flnm (str, optional): 原始wrfout数据的聚合. Defaults to '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'.
        t_start (str, optional): 起始时间. Defaults to '2021-07-20 00'.
        t_end (str, optional): 结束时间. Defaults to '2021-07-20 12'.

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'
    da = xr.open_dataarray(flnm)
    # tt = slice('2021-07-20 00', '2021-07-20 12')
    tt = slice(t_start, t_end)
    rain = da.sel(time=tt).sum(dim='time')
    lat = sum(sum(rain*rain.lat))/sum(sum(rain))
    lon = sum(sum(rain*rain.lon))/sum(sum(rain))
    lon = lon.round(3)
    lat = lat.round(3)
    return lat, lon

lat, lon = get_centroid()
# %%
lat, lon