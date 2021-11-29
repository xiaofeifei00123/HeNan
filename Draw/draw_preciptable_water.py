#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
通过micaps探空资料计算大气可降水量，绘图
-----------------------------------------
Time             :2021/11/24 10:42:18
Author           :Forxd
Version          :1.0
'''

# %%
from scipy.ndimage.measurements import label
import xarray as xr
import numpy as np
import pandas as pd
import metpy.calc as ca
from metpy.units import units  # units是units模块中的一个变量名(函数名？类名？)
import matplotlib.pyplot as plt
import datetime
from cycler import cycler
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
# %%

def get_pw(da):
    """单个时次的数据，获得它的课降水量
    """
    td = da.sel(vars='dewpoint').dropna(dim='pressure')
    # pr = td.pressure.values*units.hPa
    # td = td.values*units.degC
    pr = td.pressure*units.hPa
    td = td.values*units.degC
    pw = ca.precipitable_water(pr, td)
    pw
    return pw

def get_one_station(flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_jiuquan.nc'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station2.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_lanzhou.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_jiuquan.nc'
    ds = xr.open_dataarray(flnm)
    ds.time
    da = ds.sel(time=slice('2021-07-01', '2021-07-31'))
    tt = da.time.values
    pw_list = []
    for t in tt:
        dda = da.sel(time=t)
        pw = get_pw(dda).magnitude
        pw_list.append(pw)
    na = np.array(pw_list)
    da = xr.DataArray(
        na,
        coords={'time':tt},
        dims=['time'],
    )
    return da



###### 画图  ######

# %%
def draw_one_station(da):
    t2 = da.time.sel(time=datetime.time(int(00)))
    dda = da.sel(time=t2)
    dda.to_netcdf('/mnt/Disk4T_4/fengxiang_file/jiuquan.nc')

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    x_label = dda.coords['time']+pd.Timedelta('+8H')
    # x_label = str(x_label.dt.strftime('%d%H').values).split()
    x_label = x_label.dt.strftime('%d')
    x_label
    ax.plot(x_label, dda.values, lw=2.5, color='black',label='郑州')
    ft = 11
    ax.set_xticks(x_label[::4])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_yticks(np.arange(10,76,5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_ylim(10,75)
    ax.xaxis.set_tick_params(labelsize=ft*2.0)
    ax.yaxis.set_tick_params(labelsize=ft*2.0)
    ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
    ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.set_xlabel("Date(BJT)", fontsize=ft*2.5)
    ax.set_ylabel("Preciptable Water (mm)", fontsize=ft*2.5)
    # ax.legend(fontsize=ft*3.0, edgecolor='white')
    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/pw_jiuquan.png')

if __name__ == '__main__':
    # main()
    
    da = get_one_station()
    draw_one_station(da)





