#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图
实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''

# %%
from cmath import pi
import sys,os
import xarray as xr
import numpy as np
import pandas as pd

# import salem  # 插值
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
from baobao.map import Map


# from matplotlib import rcParams

# config = {
#     "font.family": 'serif', # 衬线字体
#     "font.size": 12, # 相当于小四大小
#     "font.serif": ['stix'], # 宋体
#     "mathtext.fontset": 'stix', # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
#     'axes.unicode_minus': False # 处理负号，即-号
# }
# rcParams.update(config)


# %%
class Draw(object):

    def __init__(self) -> None:
        super().__init__()
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'


    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
        colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        x = data.lon
        y = data.lat
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                          levels=colorlevel,
                          colors = colordict,
                          transform=ccrs.PlateCarree())
        return crx

    def draw_single(self, da, picture_dic):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        cm = round(1/2.54, 2)
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.08,0.85,0.85], projection=proj)
        date = picture_dic['date']
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
               'time':str(date)}
               
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
            # 'ZhengZhou': {
            #     'abbreviation':'郑州',
            #     'lat': 34.76,
            #     'lon': 113.65
            # },
            # 'NanYang': {
            #     'abbreviation':'南阳',
            #     'lat': 33.1,
            #     'lon': 112.49,
            # },
            # 'LuShi': {
            #     'abbreviation':'卢氏',
            #     'lat': 34.08,
            #     'lon': 111.07,
            # },
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
        mp.add_station(ax, station, justice=True, delx=-0.1)
        # mp.add_station(ax, station)
        # ax.set_title(date, fontsize=35,)

        
        if 'south_north' in da.dims:
            rain_max = da.max(dim=['south_north', 'west_east'])        
        elif 'lat' in da.dims:
            rain_max = da.max(dim=['lat', 'lon'])        
        else:
            print("出错啦")
            
        # rain_mean = da.mean(dim=['south_north', 'west_east'])        
        # # ax.set_title('Max = %s, Avg = %s'%(rain_max.values.round(2),rain_mean.values.round(2)), fontsize=35,loc='left')
        ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='left')
        
        # ax.set_title(picture_dic['initial_time'], fontsize=30,loc='left')
        if picture_dic['type'] == 'gwd0':
            picture_dic['type'] = 'no-gwd'

        ax.set_title(picture_dic['type'], fontsize=10,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
        colorticks = colorlevel[1:-1]
        
        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.1,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
        fig_name = picture_dic['type']+'_'+picture_dic['initial_time']+'typhoon'
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_24h_gwd/'
        fig.savefig(fig_path+fig_name, bbox_inches = 'tight')


def draw_tricontourf(rain):

    """rain[lon, lat, data],离散格点的DataArray数据

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    # from nmc_met_graphics.plot import mapview
    # mb = mapview.BaseMap()
    # fig = plt.figure(figsize=[3.8, 3.7], dpi=300)
    cm = round(1/2.54,2)
    fig = plt.figure(figsize=[8*cm, 8*cm], dpi=300)
    ax = fig.add_axes([0.1,0.08,0.85,0.85], projection=ccrs.PlateCarree())
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        'extent':[110.5, 116, 32, 36.5],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
    colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000',]#颜色列表
    # cs = ax.tricontour(rain.lon, rain.lat, rain, levels=colorlevel, transform=ccrs.PlateCarree())
    cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=colorlevel,colors=colordict, transform=ccrs.PlateCarree())
    colorticks = colorlevel[1:-1]
    
    cb = fig.colorbar(
        cs,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小

    mp = Map()
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
    mp.add_station(ax, station, justice=True)


    rain_max = rain.max()        
    ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='left')
    # ax.set_title('2021-07 20/00--21/00', fontsize=35,)
    ax.set_title('OBS', fontsize=10,loc='right')


    # mp.add_station(ax)
    fig_name = 'obs_distribution' 
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_24h_gwd/'
    fig.savefig(fig_path+fig_name, bbox_inches = 'tight')


def draw_obs():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
    da = xr.open_dataarray(flnm)
    da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da = da.sum(dim='time') 
    draw_tricontourf(da)
    

def draw_onemodel(model='gwd3'):
    pass

    dr = Draw()
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Typhoon/'+model+'/'+'rain.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
    da = xr.open_dataarray(flnm)
    da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da = da.sum(dim='time') 
    # picture_dic = {'date':'2021-07 20/00--21/00', 'type':model, 'initial_time':''}
    picture_dic = {'date':'2021-07 20/00--21/00', 'type':model, 'initial_time':''}
    dr.draw_single(da, picture_dic)

def draw_dual_model():    
    """对于重力波试验
    """
    # model_list = ['gwd3-LS']
    # model_list = ['gwd3-BL', 'gwd3-FD', 'gwd3-LS', 'gwd3-SS']
    # model_list = ['gwd0','gwd1', 'gwd3','gwd3-BL', 'gwd3-FD', 'gwd3-LS', 'gwd3-SS']
    model_list = ['strengthen_typhoon', 'weak_typhoon']
    # model_list = ['gwd3-test']
    # model_list = ['gwd0', 'gwd1','gwd3']
    for model in model_list:
        print("画 %s 模式"%model)
        draw_onemodel(model)



def draw_EC(model='gwd3'):
    """绘制EC细网格降水

    Args:
        model (str, optional): _description_. Defaults to 'gwd3'.
    """
    pass

    dr = Draw()
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_ec.nc'
    da = xr.open_dataarray(flnm)
    da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
    da = da.sum(dim='time').squeeze()
    picture_dic = {'date':'2021-07 20/01--21/00', 'type':'EC', 'initial_time':''}
    dr.draw_single(da, picture_dic)



if __name__ == '__main__':

    draw_obs()
    draw_dual_model()
    # draw_EC()
