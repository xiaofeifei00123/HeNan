#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
在进行相关的绘图或者数据处理前， 注意数值的单位，和有效位数
不同模式d03区域， 850hPa水汽通量散度和风场的比
各模式的图
各模式的差值的图
绘制的要素:
    850hPa
    水汽通量(填色)
    风场(矢量)
使用的数据：
    模式预报插值的格点数据
-----------------------------------------
Time             :2022/02/17
Author           :Forxd
Version          :1.1
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
# from caculate_diag import QvDiv
from baobao.caculate import QvDiv
from baobao.caculate import caculate_q_rh_thetav
from baobao.map import Map   # 一个类名
# from baobao.caculate.



# %%

# path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
# f1 = path_main+'gwd0'+'/upar.nc'
# f2 = path_main+'gwd3'+'/upar.nc'
# ds1 = xr.open_dataset(f1)
# ds2 = xr.open_dataset(f2)
# ds2
# # %%
# # ds2['ua']-ds1['ua']
# ddd = ds2-ds1
# ddd
# # %%
# ddd['ua']







# flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
# ds = xr.open_dataset(flnm_wrf)
# ds['q'].max()*1000
### 测试
# %%
def test():
    # flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
    # flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar_latlon.nc'
    flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
    ds = xr.open_dataset(flnm_wrf)

    ds = ds.isel(time=0).isel(pressure=0)
    lon = ds.lon
    lat = ds.lat
    u = ds.ua
    u.max()
    #.round(1)# *units('m/s')
    u = u*units('m/s')
    u.max()
    u.max()
    # u.metpy.dequantify(), 变成xarray

    q = ds.q*10**3*units('g/kg')
    qv_u = q*u/constants.g
    qv_v = q*u/constants.g

    dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
    qv = ca.divergence(u=qv_u, v=qv_v, dx=dx, dy=dy)
    qv





# %%

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
                # linewidth=5,
                s=10,
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
    ax.tick_params(which='major',length=6,width=1.0) # 控制标签大小 
    ax.tick_params(which='minor',length=3,width=0.5)  #,colors='b')
    ax.tick_params(axis='both', labelsize=10, direction='out')

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
    # contour_levels = [0.2, 0.6, 1.0, 1.4, 1.8, 2.2]
    # contour_levels = np.arange(-2.6, 2.6+0.1, 0.4)
    # contour_levels = np.arange(-5.6, 5.6+0.1, 0.4)
    # contour_levels = np.arange(-2.6, 2.6+0.1, 0.4)*10
    # print(da.max().values, da.min().values)
    print('最大值是{}, 最小值是{}'.format(da.max().values, da.min().values))

    # colormap = cmaps.ViBlGrWhYeOrRe
    # color_li = ['white', '#6CA6CD', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','#EE7600','red']
    # crx = ax.contourf(x,
    #                     y,
    #                     da,
    #                     cmap=colormap,
    #                     # colors=color_li,
    #                     levels=contour_levels,
    #                     # levels=10,
    #                     transform=ccrs.PlateCarree())

    # colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    # colorlevel= [-90,-2,-1.25,-0.75,-0.25,0.25,0.75,1.25,2,90]
    # colorlevel= [-90,-3,-2,-1,-0.5,0.5,1,2,3,90]
    colorlevel= [-90,-4,-3,-2,-1,1,2,3,4,90]
    # colorlevel=[-80, -3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3, 80]#雨量等级
    # colorticks=[-3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3,]#雨量等级
    # colorlevel=[-80, -8, -6, -3, -2, -1, 1, 2, 3, 6, 8, 80]#雨量等级
    # colorlevel=[-80, -5, -4, -3, -2, -1, 1, 2, 3,4, 5, 80]#雨量等级
    # colorlevel=[-80, -5, -4, -3, -2, -1, 1, 2, 3,4, 5, 80]#雨量等级
    # colorlevel=[-80, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 80]#雨量等级
    # colorticks=[-3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3,]#雨量等级
    # colorticks=colorlevel[1:-2]
    crx = ax.contourf(x,
                        y,
                        da.values,
                        colors=colordict,
                        levels=colorlevel,
                        transform=ccrs.PlateCarree()
    )




    return crx, colorlevel

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

    ax.clabel(crx,inline=1, fontsize=10, colors='red') # 等值线的标注
    return crx

def draw_quiver(u,v, ax):
    '''
    绘制风矢图
    '''
    # u = u[::3,::3]
    # v = v[::3,::3]
    u = u[::10,::10]
    v = v[::10,::10]
    # u = u[::20,::20]
    # v = v[::20,::20]
    y = u.lat.values
    x = u.lon.values

    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢

def draw(qdif, qu,qv, dic):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    # fig = plt.figure(figsize=[10,8])
    # ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    # mb = mapview.BaseMap()
    # # mb.set_extent('中国陆地')
    # mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    
    fig = plt.figure(figsize=[3.8,3.6], dpi=300)
    ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
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
    
    

    tt = (dic['time']).strftime('%Y-%m-%d %H%M')
    # print(tt)
    ax.set_title(tt, loc='center', fontsize=10)
    ax.set_title(dic['model'], loc='left', fontsize=10)
    # ax.set_title('OBS', loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=10)

    # mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    # mb.set_extent([110, 116, 32, 36])
    cs, colorlevel = draw_contourf(ax,qdif)
    colorticks = colorlevel[1:-1]
    # print(colorticks)
    # colorticks=[-3, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3,]#雨量等级
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.06,  # 色标大小
        pad=0.12,  # colorbar和图之间的距离
        ticks=colorticks,

    )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    # add_ticks(ax)
    draw_station(ax)

    draw_quiver(qu,qv,ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/850/div/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'div'
    fig.savefig(fig_name)


def get_data(dic):
    """获得模式和观测
    水汽和风场数据

    Args:
        dic ([dict]): 时间和层次信息

    Returns:
        [type]: [description]
    """
    ## 筛选数据
    flnm_wrf = dic['flnm']
    # print(flnm_wrf)
    ds_wrf = xr.open_dataset(flnm_wrf)
    ds_wrf = ds_wrf.rename({'ua':'u', 'va':'v'})
    t = dic['time']
    level = dic['level']
    ds2 = ds_wrf.sel(time=t, pressure=level)

    ## 计算水汽通量和散度
    qvd = QvDiv()
    ds3 = qvd.caculate_qfdiv(ds2)
    ds_return = xr.merge([ds2,ds3])
    return ds_return

def draw_model_once():
    """画所有模式的数据
    筛选数据
    时间
    高度
    模式
    """
    path_out = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/upar.nc'

    dic_model = {
        'model':'gwd0',
        'level':925,
        'time':pd.Timestamp('2021-07-20 00'),
        'flnm':path_out,
    }    
    print("画 [%s] 的图"%(dic_model['model']))
    ds = get_data(dic_model)
    draw(ds['qf_div']*10, ds['u'], ds['v'], dic_model)
    
    
def draw_model_dual():
    """画所有模式的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    model_list = ['gwd0', 'gwd3']
    fl_list = []
    for model in model_list:
        fl = path_main+model+'/upar.nc'
        fl_list.append(fl)
    t_list = pd.DatetimeIndex(['2021-07-20 00','2021-07-20 06', '2021-07-20 12'])
    t_list = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='1H')
    level_list = [925]

    i = 0
    for fl in fl_list:
        for t in t_list:
            for level in level_list:
                # path_out = path_main+fl
                model = fl.split('/')[-2]
                # dic_model = {'model': model, 'flnm':fl}
                # dic_model['time'] = t
                # dic_model['level'] = level
                dic_model = {
                    'model':model,
                    'level':level,
                    'time':t,
                    'flnm':fl,
                }    
                print("画北京时[%s],[%s]hPa高度, [%s]分辨率的图"%(t,level, model))
                # draw_all(dic_model)
                ds = get_data(dic_model)
                draw(ds['qf_div']*10, ds['u'], ds['v'], dic_model)
        i += 1


def draw_minus():
    """画所有模式的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # model_list = ['gwd0', 'gwd3']
    # fl_list = []
    # for model in model_list:
    f1 = path_main+'gwd0'+'/upar.nc'
    f2 = path_main+'gwd3'+'/upar.nc'
        # fl_list.append(fl)
    # t_list = pd.DatetimeIndex(['2021-07-20 00','2021-07-20 06', '2021-07-20 12'])
    t_list = pd.DatetimeIndex(['2021-07-20 00','2021-07-20 03','2021-07-20 06', '2021-07-20 09','2021-07-20 12'])
    # t_list = pd.DatetimeIndex(['2021-07-20 00'])
    level_list = [700, 850, 900]
    # level_list = [925]

    # i = 0
    # for fl in fl_list:
    for t in t_list:
        for level in level_list:
            # path_out = path_main+fl
            model1 = f1.split('/')[-2]
            model2 = f2.split('/')[-2]
            # dic_model = {'model': model, 'flnm':fl}
            # dic_model['time'] = t
            # dic_model['level'] = level
            dic_model1 = {
                'model':model1,
                'level':level,
                'time':t,
                'flnm':f1,
            }    
            dic_model2 = {
                'model':model2,
                'level':level,
                'time':t,
                'flnm':f2,
            }    
            print("画北京时[%s],[%s]hPa高度, [%s]分辨率的图"%(t,level, 'minus'))
            # draw_all(dic_model)
            ds1 = get_data(dic_model1)
            # print(ds1)
            ds2 = get_data(dic_model2)
            ds = ds2-ds1
            # print(ds)
            dic_model = {
                'model':'minus',
                'level':level,
                'time':t,
                'flnm':f1,
            }    
            draw(ds['qf_div']*10, ds['u'], ds['v'], dic_model)
    # print(ds1)
    # return ds1,ds2
        # i += 1
### 测试结束
# %%
if __name__ == '__main__':
    pass
    # draw_model_once()
    draw_model_dual()
    # draw_minus()
    # aa
