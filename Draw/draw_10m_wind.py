#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
插值后观测的10m风场
-----------------------------------------
Time             :2021/12/30 23:37:35
Author           :Forxd
Version          :1.0
'''

# %%
from time import strftime
from cartopy.crs import Projection
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from pandas._libs.tslibs.timestamps import Timestamp
import xarray as xr
import meteva.base as meb
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
# import matplotlib as mpl
import cartopy.crs as ccrs
from baobao.map import Map
import os


# %%
def get_wind_obs(flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
    ds = xr.open_dataset(flnm)
    t = '2021-07-20 00'
    ds1 = ds.sel(time=t)
    u = ds1['u']
    v = ds1['v']
    return u,v

def get_wind_wrf(flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3-test/10m_wind_station.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.swap_dims({'time':'Time'})
    ds = ds.drop_vars('time')
    ds = ds.rename({'Time':'time'})
    t = '2021-07-20 00'
    ds1 = ds.sel(time=t)
    da = ds1['uvmet10']
    u = da.sel(u_v='u')
    v = da.sel(u_v='v')
    return u,v
# %%
u1,v1 = get_wind_obs()
u2, v2 = get_wind_wrf()
u1
u2
# %%

def draw_quiver(u,v, ax):
    '''
    绘制风矢图
    '''
    y = u.lat.values
    x = u.lon.values
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=1.55, Y=0.05, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢

def draw(u,v,pic_dic={'model':'obs'}):

    # u,v = get_wind_obs()    
    # u,v = get_wind_wrf()

    proj = ccrs.PlateCarree()  # 创建坐标系
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地')
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)
    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    mb.set_extent([110, 116, 32, 36])

    
    
    
            
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        'extent':[110.5, 116, 32, 36.5],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    station = {
        'ZhengZhou': {
            'abbreviation':'ZZ',
            'lat': 34.76,
            'lon': 113.65
        },
        'NanYang': {
            'abbreviation':'NY',
            'lat': 33.1,
            'lon': 112.49,
        },
        'LuShi': {
            'abbreviation':'LS',
            'lat': 34.08,
            'lon': 111.07,
        },
    }
    mp.add_station(ax, station, justice=True, delx=0.1)
    ax.set_title(pic_dic['model'], fontsize=30,loc='left')
    
    

    draw_quiver(u,v,ax)
    fig_name = '10mwind_'+pic_dic['model']
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_10mwind/'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save)


def draw_wrf():
    # model_list = 
    model_list = ['gwd0', 'gwd1', 'gwd3','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    fpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    for model in model_list:
        fname = os.path.join(fpath, model)
        flnm = os.path.join(fname, '10m_wind_station.nc')
        u,v = get_wind_wrf(flnm)
        pic_dic = {'model':model}
        draw(u,v,pic_dic)

def draw_obs():
    # model_list = 
    # model_list = ['gwd0', 'gwd1', 'gwd3','gwd3-test','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    # fpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # for model in model_list:
    # fname = 
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
    # flnm = os.path.join(fname, '10m_wind_station.nc')
    u,v = get_wind_obs(flnm)
    pic_dic = {'model':'obs'}
    draw(u,v,pic_dic)
def main():
    draw_obs()
    draw_wrf()

# %%
if __name__ == '__main__':
    # u,v = get_wind_wrf()
    # draw(u,v)
    main()
    