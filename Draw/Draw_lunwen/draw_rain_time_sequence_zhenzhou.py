#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
郑州站逐小时降水变化柱状图
-----------------------------------------
Time             :2021/06/04 14:32:20
Author          :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd

# import salem  # 插值
# import cartopy.crs as ccrs
# import cartopy.feature as cfeat
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# from cartopy.io.shapereader import Reader, natural_earth
# import matplotlib as mpl
# from matplotlib.path import Path
# import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt

## 显示中文
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# %%

def get_rain_zhenzhou():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
    ds = xr.open_dataset(flnm)
    da = ds.to_array().squeeze()
    da_return = da.sel(id='57083')
    return da_return

def draw_bar(dazz, damax):
    # da = da.sel(time=slice('2021-07-19 12', '2021-07-21 12'))
    da = dazz
    cm = 1/2.54
    fig = plt.figure(figsize=(9*cm, 8*cm), dpi=600)
    ax = fig.add_axes([0.16,0.25,0.78, 0.6])
    # ax.bar(da.time, da)
    # import datetime
    tt = da.time.dt.strftime('%d/%H')
    x_label = tt
    ax.bar(tt, da, label='郑州站降水')
    ax.plot(tt, damax, color='red',lw=1.5,  label='最大站点降水')
    # ax.hlines(y=30, xmin=tt[0], xmax=tt[-1])
    ax.axhline(y=50, color='black')
    ax.set_xlim(tt[0], tt[-1])
    # ax.hlines(y=30, xmin=tt[0], xmax=tt[-1])
    # ax.plot(tt, da)
    ax.set_xticks(x_label[::6],)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_yticks(np.arange(0, 220+1, 20))

    ax.xaxis.set_tick_params(labelsize=10, rotation=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.tick_params(which='major',length=4,width=1.0) # 控制标签大小 
    # ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
    # ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
    # ax.set_title('逐小时降水', fontsize=32)
    ax.set_xlabel('Date/Hour (UTC)', fontsize=10)
    ax.set_ylabel('Precipitation (mm)', fontsize=10)
    ax.legend(fontsize=10, edgecolor='white', loc='upper left')
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig_name = 'obs_rian_time_sequence.png'
    fig.savefig(fig_path+fig_name)


def get_max_station():
    """读取csv格式的站点数据

    Args:
        df_station (DataFrame): 输入csv格式站点数据

    Returns:
        [DataArray]: 小时雨强最大值
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
    ds = xr.open_dataset(flnm)
    da = ds.to_array().squeeze()
    rain_obs_max = da.max(dim='id')
    return rain_obs_max

def get_rain():
    """获得郑州站降水
    河南区域最大降水
    Returns:
        [type]: [description]
    """
    dazz = get_rain_zhenzhou()
    dazz = dazz.sel(time=slice('2021-07-19 12', '2021-07-21 00'))
    damax = get_max_station()
    damax = damax.sel(time=slice('2021-07-19 12', '2021-07-21 00'))
    ## 因为郑州站有部分时次的值是0，没有在列表中
    # t1 = dazz.time
    # t2 = damax.time
    rain = xr.Dataset()
    rain = xr.concat([damax, dazz], pd.Index(['max', 'zz'], name='model'))
    # da_rain = rain.fillna(0)
    # ds_rain = da_rain.to_dataset(dim='model')
    ds_rain = rain.to_dataset(dim='model')
    return ds_rain


if __name__ == '__main__':
    
    ds_rain = get_rain()
    draw_bar(ds_rain['zz'], ds_rain['max'])
