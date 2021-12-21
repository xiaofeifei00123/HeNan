#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
郑州站逐小时降水变化柱状图
画累积降水的柱状图
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
from cycler import cycler
import datetime
import sys,os
import xarray as xr
import numpy as np
import pandas as pd

import salem  # 插值
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import geopandas
import cmaps
from get_cmap import get_cmap_rain2
from multiprocessing import Pool

## 显示中文
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# %%

def get_rain_zhenzhou():
    df_station = pd.read_csv('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.csv')
    df = df_station
    df1 = df[(df['id']==57083)]
    df2 = pd.DataFrame(df1, columns=['time','lon', 'lat', 'data0'])
    df2['time'] = pd.to_datetime(df2['time'])
    da = xr.DataArray(
        df2['data0'].values,
        coords={
            'time':df2['time'].values,
            'lon':df2['lon'].values[0],
            'lat':df2['lat'].values[0],
        },
        dims=['time'],
    )
    return da

def draw_bar(dazz, damax):
    # da = da.sel(time=slice('2021-07-19 12', '2021-07-21 12'))
    da = dazz
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_axes([0.1,0.15,0.85, 0.79])
    # ax.bar(da.time, da)
    # import datetime
    tt = da.time.dt.strftime('%d/%H')
    x_label = tt
    ax.bar(tt, da, label='郑州')
    ax.plot(tt, damax, color='red',lw=1.5,  label='最大降水站点')
    # ax.hlines(y=30, xmin=tt[0], xmax=tt[-1])
    ax.axhline(y=30, color='black')
    ax.set_xlim(tt[0], tt[-1])
    # ax.hlines(y=30, xmin=tt[0], xmax=tt[-1])
    # ax.plot(tt, da)
    ax.set_xticks(x_label[::4],)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_yticks(np.arange(0, 220+1, 20))

    ax.xaxis.set_tick_params(labelsize=2.0*12, rotation=30)
    ax.yaxis.set_tick_params(labelsize=2.0*12)
    ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
    ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.set_title('逐小时降水', fontsize=32)
    ax.set_xlabel('Date/Hour (UTC)', fontsize=26)
    ax.set_ylabel('Precip (mm)', fontsize=26)
    ax.legend(fontsize=20, edgecolor='white')
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/'
    fig_name = 'obs_rian_time_sequence.png'
    fig.savefig(fig_path+fig_name)


def get_max_station():
    """读取csv格式的站点数据

    Args:
        df_station (DataFrame): 输入csv格式站点数据

    Returns:
        [DataArray]: 小时雨强最大值
    """
    ## 才读的数据，它的时间格式是object
    df_station = pd.read_csv('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.csv')
    # df_station = get_max_dataframe(aa)
    df_station['time']= pd.to_datetime(df_station['time'])
    t = pd.date_range(start='2021-07-18 00', end='2021-07-21 12', freq='1H')
    rain_list = []
    for tt in t:
        cc = df_station[df_station['time']==tt]
        rain_max = cc[(cc['lat']>32)&(cc['lat']<37)&(cc['lon']>110)&(cc['lon']<116)]['data0'].max()
        # rain_max = cc[(cc['lat']>32)&(cc['lat']<37)&(cc['lon']>110)&(cc['lon']<116)]['data0'].mean()
        rain_list.append(rain_max)
    rain_list
    ps = pd.Series(rain_list, index=t)
    da = xr.DataArray.from_series(ps)
    rain_obs_max = da.rename({'index':'time'})
    return rain_obs_max

def get_rain():
    """获得郑州站降水
    河南区域最大降水
    Returns:
        [type]: [description]
    """
    dazz = get_rain_zhenzhou()
    dazz = dazz.sel(time=slice('2021-07-19 00', '2021-07-21 08'))
    damax = get_max_station()
    damax = damax.sel(time=slice('2021-07-19 00', '2021-07-21 08'))
    ## 因为郑州站有部分时次的值是0，没有在列表中
    t1 = dazz.time
    t2 = damax.time
    rain = xr.Dataset()
    rain = xr.concat([damax, dazz], pd.Index(['max', 'zz'], name='model'))
    da_rain = rain.fillna(0)
    ds_rain = da_rain.to_dataset(dim='model')
    return ds_rain


if __name__ == '__main__':
    
    ds_rain = get_rain()
    draw_bar(ds_rain['zz'], ds_rain['max'])
