#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画850hPa,
wrf, 观测，ERA5的高空风对比
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
# import meteva.base as meb
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
# import meteva.base as meb
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
# from Draw.draw_upar_850_d03 import draw_one_model
# from metdig.io.cassandra import get_obs_station
# from caculate_diag import QvDiv
# from baobao.caculate import QvDiv
# from baobao.caculate import caculate_q_rh_thetav
# from baobao.caculate import caculate_vo_div, caculate_div
from baobao.map import Map   # 一个类名
# from baobao.caculate.
from dask import delayed, compute
from baobao.timedur import timeit



# %%
# flnmec = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/era5.pl.20210720.nc'
# ds = xr.open_dataset(flnmec)
# ds


# %%
def get_uvst(pre=900, tt='2021-07-20 00'):
    def get_uv(ds, staname='zhengzhou', pre=900, tt='2021-07-20 00'):
        ds1 = ds.sel(station=staname)
        ds2 = ds1.sel(time=tt)
        u = ds2['u']
        u1 = u.dropna(dim='pressure')
        u2 = u1.sel(pressure=pre, method='nearest')
        v = ds2['v']
        v1 = v.dropna(dim='pressure')
        v2 = v1.sel(pressure=pre, method='nearest')
        return u2.values, v2.values

    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_all.nc'
    ds = xr.open_dataset(flnm)
    u_zz, vzz = get_uv(ds, 'zhengzhou', pre, tt)
    u_ny, vny = get_uv(ds, 'nanyang', pre, tt)
    u_ls, vls = get_uv(ds, 'lushi', pre, tt)
    ust = {'ZhengZhou':u_zz, 'NanYang':u_ny, 'LuShi':u_ls}
    vst = {'ZhengZhou':vzz, 'NanYang':vny, 'LuShi':vls}
    return ust, vst

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

def draw_station(ax, ust, vst):
    pass
    # station = station_dic
    # station = {
    #     'ZhengZhou': {
    #         'abbreviation':'郑州',
    #         'lat': 34.75,
    #         'lon': 113.62
    #     },
    # }

    station = {
            'ZhengZhou': {
                'abbreviation':'郑州',
                'lat': 34.76,
                'lon': 113.65
            },
            'NanYang': {
                'abbreviation':'南阳',
                'lat': 33.1,
                'lon': 112.49,
            },
            'LuShi': {
                'abbreviation':'卢氏',
                'lat': 34.08,
                'lon': 111.07,
            },
        }
    
    

    values = station.values()
    keys = list(station.keys())
    station_name = []
    x = []
    y = []
    for i,j in zip(values,keys):
        y.append(float(i['lat']))
        x.append(float(i['lon']))
        station_name.append(i['abbreviation'])
        # u = u
    # for j in station.keys():
        # Q = ax.quiver(i['lon'],i['lat'],ust[j],vst[j], units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree(), color='red')
        Q = ax.quiver(i['lon'], i['lat'], ust[j],vst[j],units='width',scale=200,pivot='middle',  width=0.003, transform=ccrs.PlateCarree(), color='red')  # 绘制风矢
        # Q = ax.quiver(i['lon'],i['lat'],ust[j],vst[j], scale=80,pivot='middle', transform=ccrs.PlateCarree(), color='red')
        # qk = ax.quiverkey(Q,
        #                 X=0.85, Y=0.1, 
        #                 U=10 ,
        #                 label=r'$10 m/s$', 
        #                 labelpos='S',  # label在参考箭头的哪个方向
        #                 labelsep=0.05, # 箭头和标签之间的距离
        #                 coordinates='figure',  
        #                 fontproperties={'size':8}
        #                 )   # 设置参考风矢
        # ax.quiver(x, y, u.values,v.values,units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢

    # 标记出站点
    # ax.scatter(x,
    #             y,
    #            color='black',
    #             # color='green',
    #             transform=ccrs.PlateCarree(),
    #             alpha=1.,
    #             # linewidth=5,
    #             s=5,
    #             marker='o',
    #             )
    # 给站点加注释
    # ax.quiver(x,y,u2,v2)
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


def draw_south_sea(fig,):
    pass
    ax2 = fig.add_axes([0.798, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2.add_geometries(Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='black',linewidth=0.8)

def draw_contourf(ax, da):
    """在地图上绘制填色图
    """
    x = da.lon
    y = da.lat
    # print(da.max().values, da.min().values)
    print('最大值是{}, 最小值是{}, 平均值是{}'.format(da.max().values, da.min().values, da.mean().values))


    # colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    colordict=['#0000fb','#3232fd','#6464fd','white','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    # colordict=['#0000fb','#5a5af9','#a9a9fd','#dfdffd','white','white','#ffdbdb', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    # colorlevel= [-90,-2.5,-1.5,-0.6,-0.3,0, 0.3,0.6,1.5,2.5,90]# 散度的色标
    # colorlevel= [-90,-12,-9,-6,-3, 0, 3,6,9,12,90]  # 垂直速度的色标
    colorlevel= [-90,-20,-10,-3, 0, 3,10,20,90]  # 垂直速度的色标
    # colorlevel= [-90,-4,-3,-2,-1,0,1,2,3,4,90]# 水汽通量散度的色标
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
    # u = u[::10,::10]
    # v = v[::10,::10]
    u = u[::20,::20]
    v = v[::20,::20]
    y = u.lat.values
    x = u.lon.values

    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u,v,units='width',scale=200,pivot='middle',  width=0.003, transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q,
                      X=0.05, Y=0.15, 
                      U=10 ,
                      label=r'$10 m/s$', 
                    #   label='10 m/s', 
                      labelpos='S',  # label在参考箭头的哪个方向
                      labelsep=0.05, # 箭头和标签之间的距离
                      coordinates='figure',  
                      fontproperties={'size':8}
                      )   # 设置参考风矢

def draw_quiver_ec(u,v, ax):
    '''
    绘制风矢图
    '''
    u = u[::3,::3]
    v = v[::3,::3]
    # u = u[::10,::10]
    # v = v[::10,::10]
    # u = u[::20,::20]
    # v = v[::20,::20]
    y = u.lat.values
    x = u.lon.values

    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u,v,units='width',scale=200,pivot='middle',  width=0.003, transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q,
                      X=0.05, Y=0.15, 
                      U=10 ,
                      label=r'$10 m/s$', 
                    #   label='10 m/s', 
                      labelpos='S',  # label在参考箭头的哪个方向
                      labelsep=0.05, # 箭头和标签之间的距离
                      coordinates='figure',  
                      fontproperties={'size':8}
                      )   # 设置参考风矢

def draw(u,v, dic):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    
    ax = fig.add_axes([0.12,0.18,0.8,0.7], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        # 'extent':[110.5, 116, 32, 36.5],
        'extent':[107, 118, 28, 38],
        'extent_interval_lat':1,
        'extent_interval_lon':2,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])
    
    tt = (dic['time']).strftime('%Y-%m-%d %H')
    ax.set_title(tt, loc='center', fontsize=10)
    if dic['model'] == 'gwd0':
        dic['model'] = 'CTRL'
    ax.set_title(dic['model'], loc='left', fontsize=10)
    # ax.set_title('OBS', loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=10)
    ust, vst = get_uvst(dic['level'], dic['time'])
    if dic['model'] == 'EC':
        draw_quiver_ec(u, v, ax)
    else:
        draw_quiver(u,v,ax)

    draw_station(ax, ust, vst)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/upar_wind/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'w_speed'
    fig.savefig(fig_name, bbox_inches = 'tight')
    # fig.savefig(fig_name,)

# def draw_ec():
    # pass
    
# %%
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

    return ds2

def draw_model_once(tt=pd.Timestamp('2021-07-20 00'), level=500, model='GWD3'):
    """画所有模式的数据
    筛选数据
    时间
    高度
    模式
    """
    if model == 'GWD3':
        path_out = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
    else:
        path_out = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/upar.nc'

    dic_model = {
        'model':model,
        # 'model':'gwd0',
        'level':level,
        'time':tt,
        'flnm':path_out,
    }    
    print("画 [%s] [%s]的图"%(dic_model['model'], dic_model['time'].strftime('%Y%m%d-%H')))
    ds = get_data(dic_model)
    draw(ds['u'], ds['v'], dic_model)

# tt=pd.Timestamp('2021-07-20 00')
@timeit
def draw_ec():
    def draw_ec1(tt, level):
        # tt=pd.Timestamp('2021-07-19 12')
        dic_model = {
            'model':'EC',
            # 'model':'gwd0',
            'level':level,
            'time':tt,
        }    
        print("画 [%s] [%s] [%s]的图"%(dic_model['model'], dic_model['time'].strftime('%Y%m%d-%H'), level))
        flnmec = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/era5.nc'
        ds = xr.open_dataset(flnmec)
        dss = ds.rename({'latitude':'lat', 'longitude':'lon'}).sel(time=tt).sel(level=level)
        u = dss['u']
        v = dss['v']
        # u
        draw(u, v, dic_model)
    
    t_list = pd.DatetimeIndex(['2021-07-19 12','2021-07-20 00','2021-07-20 06','2021-07-20 12'])
    level_list = [200,500, 700, 850, 900, 925]
    for t in t_list:
        for level in level_list:
            draw_ec1(t, level)
            # c_list.append(aa)
    # compute(c_list)
    



# %%
    
def draw_model_dual():
    """画所有模式的数据
    """
    t_list = pd.DatetimeIndex(['2021-07-19 12','2021-07-20 00','2021-07-20 06','2021-07-20 12','2021-07-21 00'])
    # t_list = pd.DatetimeIndex(['2021-07-19 20','2021-07-20 00','2021-07-20 06','2021-07-20 12','2021-07-20 20','2021-07-21 00'])
    level_list = [200,500, 700, 850, 900, 925]

    for t in t_list:
        for level in level_list:
            for model in ['GWD3', 'CTRL']:
                draw_model_once(t, level, model)
    
### 测试结束
if __name__ == '__main__':
    # draw_model_once()
    draw_model_dual()
    draw_ec()

# %%
