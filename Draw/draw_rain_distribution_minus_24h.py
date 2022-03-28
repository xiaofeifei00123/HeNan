#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图
# 实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
模式降水的差值
-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''

# %%
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
from baobao.interp import rain_station2grid   # 站点插值成格点，这里插到和EC网格点一样



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

        # colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
        # colorlevel=[-700, -200, -100, -50, -10, 10, 50 , 100, 200,700 ]#雨量等级
        colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
        # colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红

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
        # fig = plt.figure(figsize=(3.8, 3.7), dpi=300)
        cm = round(1/2.54, 2)
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
        # fig = plt.figure()
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
        mp.add_station(ax, station, justice=True, delx=-0.2)
        ax.set_title(picture_dic['type'], fontsize=10,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
        colorticks = colorlevel[1:-1]
        
        # cb = fig.colorbar(
        #     cf,
        #     orientation='vertical',
        #     ticks=colorticks,
        #     fraction = 0.038,  # 色标大小,相对于原图的大小
        #     pad=0.02,  #  色标和子图间距离
        # )

        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.1,  #  色标和子图间距离
        )
        
        
        
        

        cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
        fig_name = picture_dic['type']+'_'+picture_dic['initial_time']
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_24h_gwd/'
        fig.savefig(fig_path+fig_name)




def draw_minus(model='gwd3'):
    """原始投影格式的数据，直接相减

    Args:
        model (str, optional): _description_. Defaults to 'gwd3'.
    """
    pass

    dr = Draw()
    flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain.nc'
    flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'

    da3 = xr.open_dataarray(flnm3)
    da3 = da3.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da3 = da3.sum(dim='time') 

    da0 = xr.open_dataarray(flnm0)
    da0 = da0.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da0 = da0.sum(dim='time') 
    
    da = da3-da0
    # picture_dic = {'date':'2021-07 20/00--21/00', 'type':model, 'initial_time':''}
    picture_dic = {'date':'2021-07 20/00--21/00', 'type':'gwd3-nogwd', 'initial_time':''}
    dr.draw_single(da, picture_dic)

    
def draw_minus_latlon(model='gwd3'):
    """插值过后的数据相减

    Args:
        model (str, optional): _description_. Defaults to 'gwd3'.
    """
    pass

    dr = Draw()
    flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
    flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
    da3 = xr.open_dataarray(flnm3)
    da3 = da3.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da3 = da3.sum(dim='time') 

    da0 = xr.open_dataarray(flnm0)
    da0 = da0.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    da0 = da0.sum(dim='time') 
    
    da = da3.values-da0.values

    da = xr.DataArray(
        da,
        coords={
            'lat':da3.lat.values,
            'lon':da3.lon.values
        },
        dims=['lat', 'lon']
        )
    picture_dic = {'date':'2021-07 20/00--21/00', 'type':'gwd3-obs', 'initial_time':''}
    dr.draw_single(da, picture_dic)
    

def draw_minus_EC(model='gwd3'):
    """原始投影格式的数据，直接相减

    Args:
        model (str, optional): _description_. Defaults to 'gwd3'.
    """
    pass

    dr = Draw()

    flnm = '/home/fengxiang/HeNan/Data/OBS/rain_ec.nc'
    da2 = xr.open_dataarray(flnm)
    da2 = da2.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
    da2 = da2.sel(lat=slice(32,36.5)).sel(lon=slice(110.5,116))
    da2 = da2.sum(dim='time')


    area = {
        'lon1':110.5,
        'lon2':116,
        'lat1':32,
        'lat2':36.5,
        'interval':0.125,
    }
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc')
    da1 = rain_station2grid(da, area)
    da1 = da1.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
    da1 = da1.sum(dim='time')
    
    da = da2-da1
    # picture_dic = {'date':'2021-07 20/00--21/00', 'type':model, 'initial_time':''}
    picture_dic = {'date':'2021-07 20/00--21/00', 'type':'EC-OBS', 'initial_time':''}
    dr.draw_single(da, picture_dic)



if __name__ == '__main__':

    draw_minus()
    # draw_minus_latlon()
    # draw_minus_EC()
