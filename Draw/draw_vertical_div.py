#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
不同模式d04区域， 850hPa水汽通量的比较
各模式的图
各模式的差值的图
绘制的要素:
500hPa散度
使用的数据：
    wrfout的格点数据
    模式预报插值的格点数据
-----------------------------------------
Time             :2021/10/28 
Author           :Forxd
Version          :1.0
'''

# %%
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader
import xarray as xr
import meteva.base as meb
import numpy as np
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from wrf import constants, get_cartopy, smooth2d, getvar
# from baobao.caculate import Qv, QvDiv
from baobao.map import Map   # 一个类名
from baobao.caculate import caculate_div

# %%



# %%
def draw_contourf(ax, da):
    """在地图上绘制填色图
    """
    x = da.lon
    y = da.lat
    contour_levels = [-90,-2,-1.25,-0.75,-0.25,0.25,0.75,1.25,2,90]
    # contour_levels = [-10,-7,-5,-3,-1,1,3,5,7,10]
    # contour_levels = [-10,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,10]
    # contour_levels = [-10,-3.5,-2.5,-1.5,-0.25,0.25,1.5,2.5,3.5,10]
    # color_li = ['white', '#6CA6CD', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','#EE7600','red']
    color_li=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    crx = ax.contourf(x,
                        y,
                        da*10**3,
                        colors=color_li,
                        levels=contour_levels,
                        # levels = 10,
                        transform=ccrs.PlateCarree())
    return crx

def draw_contour(ax, da, **kw):
    """在地图上绘制等温线
    """
    x = da.lon
    y = da.lat
    level_contour = [6, 8, 10, 12, 14, 16] 
    
    da = smooth2d(field=da, passes=2)
    crx = ax.contour(x,
                        y,
                        da,
                        colors = 'red',
                        levels=level_contour,
                        linestyles = 'dashed',
                        transform=ccrs.PlateCarree())

    ax.clabel(crx,inline=1, fontsize=20, colors='red') # 等值线的标注
    return crx

def draw_quiver(u,v, ax):
    '''
    绘制风矢图
    '''
    u = u[::12,::12]
    v = v[::12,::12]
    y = u.lat.values
    x = u.lon.values
    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=30,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u.values,v.values,units='inches',linewidth=1.5, scale=30,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢

def draw(da):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.1,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    
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
    }
    mp.add_station(ax, station, justice=True)
    # w = ds['wa']
    cs = draw_contourf(ax,da)
    # contour_levels = [10, 12, 14, 16, 18, 20, 22, 24, 30]
    contour_levels = [-90,-2,-1.25,-0.75,-0.25,0.25,0.75,1.25,2,90]
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.05,  # 色标大小
        pad=0.12,  # colorbar和图之间的距离
        ticks=contour_levels[1:-1],
    )
    cb.ax.tick_params(labelsize=20)  # 设置色标标注的大小

    # u = ds['ua']
    # v = ds['va']
    # draw_quiver(u,v,ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/'
    fig_name = fig_path + 'div.png'
    fig.savefig(fig_name)

# %%
if __name__ == '__main__':
    # main()
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/upar.nc'
    ds = xr.open_dataset(flnm)
    ds1 = ds.sel(time='2021-07-20 00').sel(pressure=500)
    ds2 = ds1.rename({'ua':'u', 'va':'v'})
    div = caculate_div(ds2)
    # aa
    draw(div)
    
    # draw(ds1)


