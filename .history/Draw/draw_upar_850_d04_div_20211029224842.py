#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
不同模式d04区域， 850hPa水汽通量散度的比较
各模式的图
各模式的差值的图
绘制的要素:
    850hPa
    水汽通量(填色)
    风场(矢量)
使用的数据：
    模式预报插值的格点数据
-----------------------------------------
Time             :2021/10/28 
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
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cmaps
from wrf import constants, get_cartopy, smooth2d, getvar
import metpy.calc as ca
import metpy
from metpy.units import units
from pyproj import CRS
from metpy import constants
# from metdig.io.cassandra import get_obs_station
from caculate_diag import QvDiv



### 测试
# %%
def test():
    flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
    ds = xr.open_dataset(flnm_wrf)
    ds = ds.isel(time=0).isel(pressure=2)
    lon = ds.lon
    lat = ds.lat
    u = ds.u*units('m/s')
    v = ds.v*units('m/s')
    q = ds.q*10**3*units('g/kg')
    qv_u = q*u/constants.g
    qv_v = q*u/constants.g

    dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
    qv = ca.divergence(u=qv_u, v=qv_v, dx=dx, dy=dy)




# %%
def caculate_data(dds):
    """ 抽取并计算出需要的数据
    u, v, qu, qv

    Args:
        dds ([type]): [description]

    Returns:
        [type]: [description]
    """
    dict = {}

    lon = dds.lon
    lat = dds.lat
    u = dds.u*units('m/s')
    v = dds.v*units('m/s')
    q = dds.q*10**3*units('g/kg')
    qv_u = q*u/constants.g
    qv_v = q*u/constants.g

    dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
    qv = ca.divergence(u=qv_u, v=qv_v, dx=dx, dy=dy)
    qf = xr.ufuncs.sqrt(qv_u**2+qv_v**2)
    
    

    dict['q'] = q
    dict['qu'] = qv_u
    dict['qv'] = qv_v
    dict['qf'] = qf
    dict['u'] = u
    dict['v'] = v
    dict['qv_div'] = qv*10**3
    return dict

def get_data(dic):
    """获得模式和观测
    水汽和风场数据

    Args:
        dic ([dict]): 时间和层次信息

    Returns:
        [type]: [description]
    """
    flnm_wrf = dic['flnm']
    ds_wrf = xr.open_dataset(flnm_wrf)
    t = dic['time']
    level = dic['level']
    ds2 = ds_wrf.sel(time=t, pressure=level)
    qvd = QvDiv()
    # dic2 = caculate_data(ds2)
    dic2 = {}
    dic2 = caculate_data(ds2)

    dic_return = {
        'u_model':dic2['u'],
        'v_model':dic2['v'],
        'qf_model':dic2['qf'],
        'qv_model':dic2['qv_div']
    }

    return dic_return 

def draw_station(ax):
    pass
    # station = station_dic
    station = {
        'ZhengZhou': {
            'abbreviation':'郑州',
            'lat': 34.75,
            'lon': 113.62
        },
    }
    values = station.values()
    station_name = list(station.keys())
    station_name = []
    x = []
    y = []
    for i in values:
        y.append(float(i['lat']))
        x.append(float(i['lon']))
        station_name.append(i['abbreviation'])

    # 标记出站点
    ax.scatter(x,
                y,
               color='black',
                # color='green',
                transform=ccrs.PlateCarree(),
                alpha=1.,
                linewidth=5,
                s=20,
                marker='o',
                )
    # 给站点加注释
    # for i in range(len(x)):
    #     # print(x[i])
    #     ax.text(x[i]-0.2,
    #             y[i] + 0.2,
    #             station_name[i],
    #             transform=ccrs.PlateCarree(),
    #             alpha=1.,
    #             fontdict={
    #             'size': 28,
    #     })

def add_ticks(ax,):
    """添加坐标标签"""
    # mb.set_extent([110, 116, 32, 36])
    ax.set_yticks(np.arange(32, 36 + 1, 1))
    ax.set_xticks(np.arange(110, 116 + 1, 1))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.tick_params(which='major',length=10,width=2.0) # 控制标签大小 
    ax.tick_params(which='minor',length=5,width=1.0)  #,colors='b')
    ax.tick_params(axis='both', labelsize=25, direction='out')

def draw_south_sea(fig,):
    pass
    ax2 = fig.add_axes([0.798, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2.add_geometries(Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='black',linewidth=0.8)

def draw_contourf(ax, da):
    """在地图上绘制填色图
    """
    x = da.lon
    y = da.lat
    # contour_levels = [-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3,  -0.1, 0.1,  0.3,  0.5, 0.7, 0.9,  1.1, 1.3, 1.5]
    contour_levels = np.arange(-1.9, 1.9+0.1, 0.2)
    colormap = cmaps.ViBlGrWhYeOrRe
    color_li = ['white', '#6CA6CD', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','#EE7600','red']
    crx = ax.contourf(x,
                        y,
                        da,
                        cmap=colormap,
                        # colors=color_li,
                        levels=contour_levels,
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
    u = u[::3,::3]
    v = v[::3,::3]
    y = u.lat.values
    x = u.lon.values

    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=20,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢

def draw(qdif, qu,qv, dic):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地')
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    tt = (dic['time']).strftime('%Y-%m-%d %H%M')
    # print(tt)
    ax.set_title(tt, loc='center', fontsize=25)
    ax.set_title(dic['model'], loc='left', fontsize=25)
    # ax.set_title('OBS', loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=25)

    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    mb.set_extent([110, 116, 32, 36])
    cs = draw_contourf(ax,qdif)
    contour_levels = [-1.5, -1.1,   -0.7,  -0.3, 0.3,  0.7,   1.1, 1.5]
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.06,  # 色标大小
        pad=0.12,  # colorbar和图之间的距离
        ticks=contour_levels,

    )
    cb.ax.tick_params(labelsize=24)  # 设置色标标注的大小
    add_ticks(ax)
    draw_station(ax)

    draw_quiver(qu,qv,ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/850/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'div'
    fig.savefig(fig_name)


def draw_all(dic_model):
    """
    在一张图上画所有需要的要素
    获取数据，调用画图程序
    高度场
    温度场
    风场差

    Args:
        t ([type]): [description]
        level (str, optional): [description]. Defaults to '500'.
    """
    # t = pd.Timestamp('2021-07-19 2000')
    level = dic_model['level']
    t = dic_model['time']
    dic_t = {
        'var':'temp',
        'level':level,
        'time':t
    }
    dic_h = {
        'var':'height',
        'level':level,
        'time':t
    }
    dic_v = {
        'var':'height',
        'level':level,
        'time':t-pd.Timedelta('8H')
    }
    dic_model['time'] = dic_model['time']-pd.Timedelta('8H')
    dic_data = get_data(dic_model)
    qdif = dic_data['qv_model']
    qu = dic_data['u_model']
    qv = dic_data['v_model']
    draw(qdif, qu, qv, dic_model)



def draw_model_once():
    """画所有模式的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    # flnm = 'high_resolution_high_hgt_upar_d04_latlon.nc'
    flnm = 'high_resolution_upar_d04_latlon.nc'
    path_out = path_main+flnm
    dic_model = {'model': '1km', 'flnm':path_out}
    # dic_model = {'model': '1km_hr', 'flnm':path_out}

    print("画 [%s] 的图"%(dic_model['model']))
    ### 北京时
    dic_model['time'] = pd.Timestamp('2021-07-20 08')
    dic_model['level'] = 850
    draw_all(dic_model)
    
def draw_model_dual():
    """画所有模式的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    flnm1 = 'high_resolution_high_hgt_upar_d04_latlon.nc'
    flnm2 = 'high_resolution_upar_d04_latlon.nc'
    fl_list = [flnm1, flnm2]
    model_list = ['1km_hr', '1km']
    ### 北京时
    t_list = pd.DatetimeIndex(['2021-07-20 08','2021-07-20 16'])
    level_list = [850, 925]

    i = 0
    for fl in fl_list:

        for t in t_list:
            for level in level_list:
                path_out = path_main+fl
                if i == 0:
                    model = '1km_hr'
                else: 
                    model = '1km'
                dic_model = {'model': model, 'flnm':path_out}
                print("画北京时[%s], [%s]高度的图"%(t, model))
                dic_model['time'] = t
                dic_model['level'] = level
                draw_all(dic_model)
        i += 1


### 测试结束
# %%
if __name__ == '__main__':
    pass
    # draw_model_once()
    draw_model_dual()



# %%