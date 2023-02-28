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
# from time import strftime
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
# from metdig.io.cassandra import get_obs_station
# from caculate_diag import QvDiv
from baobao.caculate import QvDiv
from baobao.caculate import caculate_q_rh_thetav
from baobao.caculate import caculate_vo_div, caculate_div
from baobao.map import Map   # 一个类名
from baobao.map import get_rgb
# from draw_rain_distribution import Draw
from common import Common



# %%

path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
f1 = path_main+'gwd0'+'/upar.nc'
f2 = path_main+'gwd3'+'/upar.nc'
ds1 = xr.open_dataset(f1)
ds1
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
    # print('最大值是{}, 最小值是{}, 平均值是{}'.format(da.max().values, da.min().values, da.mean().values))
    print('最大值是{}, 最小值是{}, 平均值是{}'.format(np.nanmax(da.values), np.nanmin(da.values), np.nanmean(da.values)))


    # colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    # colordict=['#0000fb','#5a5af9','#a9a9fd','#dfdffd','white','white','#ffdbdb', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
    # colorlevel= [-90,-2.5,-1.5,-0.6,-0.3,0, 0.3,0.6,1.5,2.5,90]# 散度的色标
    # colorlevel= [-90,-12,-9,-6,-3, 0, 3,6,9,12,90]  # 垂直速度的色标
    colordict=['#0000fb','#3232fd','#6464fd','white','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    # fnrgb = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/11colors_blue_red.rgb'
    fnrgb = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/7colors.rgb'
    # fnrgb = ['#0000fb','#3232fd','#6464fd','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    rgb = get_rgb(fnrgb)
    colordict = rgb

    colordict=['#0000fb','#3232fd','#6464fd','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    # colorlevel= [-2000, -40,-30,-20,-10,-5, 5,10,20,30, 40, 2000]  # 垂直速度的色标
    # colorlevel= [-2000, -120,-80,-50,-30,-10,10,30, 50, 80,120,2000]  # 垂直速度的色标
    # colorlevel= [-3000, -70,-40,-5,5,40,70,3000]  # 垂直速度的色标
    # colorlevel= [-3000, -50,-20,-5,5,20,50,3000]  # 垂直速度的色标
    colorlevel= [-3000, -80,-40,-10,10,40,80,3000]  # 垂直速度的色标

    # colorlevel= [-90,-20,-10,-3, 0, 3,10,20,90]  # 垂直速度的色标
    # colorlevel= [-90,-20,-10,-1, 0, 1,10,20,90]  # 垂直速度的色标
    # colorlevel= [-90,-20,-10,-5, 0, 5,10,20,90]  # 垂直速度的色标
    # colorlevel= [-90,-4,-3,-2,-1,0,1,2,3,4,90]# 水汽通量散度的色标
    # colorticks=colorlevel[1:-2]
    # da = smooth2d(field=da, passes=16)
    da = smooth2d(field=da, passes=6)
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
    qk = ax.quiverkey(Q,
                      X=0.05, Y=0.15, 
                      U=10 ,
                      label=r'$10 m/s$', 
                      labelpos='S',  # label在参考箭头的哪个方向
                      labelsep=0.05, # 箭头和标签之间的距离
                      coordinates='figure',  
                      fontproperties={'size':8}
                      )   # 设置参考风矢

def draw(qdif, qu,qv, dic):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    cm = round(1/2.54, 2)
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
    
    ax = fig.add_axes([0.12,0.1,0.8,0.8], projection=ccrs.PlateCarree())
    # mb = mapview.BaseMap()
    # mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    
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
    ax.set_title(tt, loc='center', fontsize=10)
    if dic['model'] == 'gwd0':
        dic['model'] = 'no-gwd'
    ax.set_title(dic['model'], loc='left', fontsize=10)
    # ax.set_title('OBS', loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=10)

    cs, colorlevel = draw_contourf(ax,qdif)
    colorticks = colorlevel[1:-1]
    # ax.set_xlabel('vertical velocity $(w, 10^{-1}m/s)$', labelpad=0.01, fontsize=8)  # 和图片的距离
    # ax.set_xlabel('散度 $10^{-5}s^{-1}$', labelpad=0.01, fontsize=8)  # 水汽通量散度
    ax.set_xlabel('垂直速度 $10^{-2}m \cdot s^{-1}$', labelpad=0.01, fontsize=8)  # 水汽通量散度
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.05,  # 色标大小
        pad=0.15,  # colorbar和图之间的距离
        ticks=colorticks,
    )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小
    # draw_station(ax)

    draw_quiver(qu,qv,ax)

    
    # dr = Draw(fig, ax)
    dr = Common()
    ax.plot(np.linspace(dr.cross_start[0], dr.cross_end[0], 10), np.linspace(dr.cross_start[1], dr.cross_end[1], 10), color='black')
    # mp.add_station(dr.station)
    mp.add_station(ax, dr.station_zz, justice=True, ssize=30)
    

    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/850/alltime/'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_upar/vs/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'div'
    # fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'div'
    # fig.savefig(fig_name, bbox_inches = 'tight')
    fig.savefig(fig_name)
    # fig.savefig(fig_name,)


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
    # qvd = QvDiv()
    # ds3 = qvd.caculate_qfdiv(ds2)
    ## 计算涡度和散度
    ds3 = caculate_vo_div(ds2)
    ds_return = xr.merge([ds2,ds3])
    return ds_return

def draw_model_once():
    """画所有模式的数据
    筛选数据
    时间
    高度
    模式
    """
    # path_out = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
    # path_out = '/home/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/upar.nc'
    path_out = '/home/fengxiang/HeNan/Data/GWD/d03/DA/GWD3/upar.nc'


    # for t in pd.date_range('2021-07-17 00', '2021-07-23 00', freq='12H'):
    for t in pd.date_range('2021-07-20 01', '2021-07-20 01', freq='12H'):
        lev_list = [900, 925, 850, 700, 600, 500, 400,300,200]
        for lev in lev_list:
            dic_model = {
                'model':'GWD3',
                # 'level':850,
                'level':lev,
                'time':t,
                'flnm':path_out,
            }    
            print("画 [%s] 的图"%(dic_model['model']))
            ds = get_data(dic_model)
            # print(ds)
            # daa = ds['div']
            # da = daa*10**5

            daa = ds['wa']
            da = daa*10**2

            da1 = xr.DataArray(da,
                            coords=daa.coords,
                            dims=daa.dims,
                            )

            draw(da1, ds['u'], ds['v'], dic_model)

    
    
def draw_model_dual():
    """画所有模式的数据
    """
    # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    path_main = '/home/fengxiang/HeNan/Data/GWD/d03/newall/'
    model_list = ['GWD3', 'CTRL', 'FD', 'SS']
    fl_list = []
    for model in model_list:
        fl = path_main+model+'/upar.nc'
        fl_list.append(fl)
    # t_list = pd.DatetimeIndex(['2021-07-20 06','2021-07-20 06'])
    # t_list = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='6H')
    t_list = pd.date_range('2021-07-17 00', '2021-07-23 00', freq='12H')
    # level_list = [500, 700, 850, 900, 925]
    # level_list = [500, 700, 850, 900, 925]
    # level_list = [500, 850,]
    # level_list = [700,]
    level_list = [850,]

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
                # dic_model2 = {
                #     'model':model,
                #     'level':900,
                #     'time':t,
                #     'flnm':fl,
                # }    
                print("画北京时[%s],[%s]hPa高度, [%s]分辨率的图"%(t,level, model))
                # draw_all(dic_model)
                ds = get_data(dic_model)

                
                # print(ds)
                daa = ds['div']
                da = daa*10**5

                da1 = xr.DataArray(da,
                                coords=daa.coords,
                                dims=daa.dims,
                                )

                draw(da1, ds['u'], ds['v'], dic_model)
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
    # t_list = pd.DatetimeIndex(['2021-07-20 00','2021-07-20 03','2021-07-20 06', '2021-07-20 09','2021-07-20 12'])
    # t_list = pd.date_range('2021-07-20 01', '2021-07-21 00', freq='1H')
    t_list = pd.DatetimeIndex(['2021-07-20 00'])
    # level_list = [700, 850, 900]
    level_list = [700]

    # i = 0
    # for fl in fl_list:

    ds_list = []
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
            ds_list.append(ds)
            # print(ds)
            dic_model = {
                'model':'minus',
                'level':level,
                'time':t,
                'flnm':f1,
            }    
            draw(ds['qf_div']*10, ds['u'], ds['v'], dic_model)

    # dic_model = {
    #     'model':'mean',
    #     'level':level,
    #     'time':t,
    #     'flnm':f1,
    # }    
    # aa = xr.concat(ds_list, dim='time')
    # bb = aa.mean(dim='time')
    # draw(bb['qf_div']*10, bb['u'], bb['v'], dic_model)
    
    # return ds_list
    # print(ds1)
    # return ds1,ds2
        # i += 1
### 测试结束
if __name__ == '__main__':
    pass
    draw_model_once()
    # draw_model_dual()
    # aa = draw_minus()

# %%
