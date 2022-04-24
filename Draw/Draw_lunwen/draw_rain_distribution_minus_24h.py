#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图
# 实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
模式降水的差值
原始的图片
原始的数据
封装色标， 封装子图

基本图像
基本数据

主函数， 对它们进行拼接
保证可以对单幅图有效

-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''

# %%
import sys,os





from matplotlib.ticker import Locator
import xarray as xr
import numpy as np
import pandas as pd
import wrf
import netCDF4 as nc

# import salem  # 插值
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
import matplotlib.patches as patches
import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import geopandas
import cmaps
from multiprocessing import Pool
from baobao.map import Map
from baobao.interp import rain_station2grid   # 站点插值成格点，这里插到和EC网格点一样

# %%

station = {
        'ZhengZhou': {
            # 'abbreviation':'郑州',
            'abbreviation':'ZZ',
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

# station[list(station.keys())[0]]
# station['ZhengZhou']
# list(station)



# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/gwd0/wrfout_d04_2021-07-19_12:00:00'
def get_hgt(
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/gwd0/wrfout_d04_2021-07-19_12:00:00'
    flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d02_2021-07-19_12:00:00'
    ):
    """从wrfout文件中获取海拔高度数据
    """
    pass
    wrfnc = nc.Dataset(flnm)
    h = wrf.getvar(wrfnc, 'ter')
    h = h.rename({'XLAT':'lat', 'XLONG':'lon'})
    return h


# %%
class Draw(object):

    def __init__(self, fig, ax) -> None:
        super().__init__()
        self.fig = fig
        self.ax = ax
        self.colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
        self.colorticks = self.colorlevel[1:-1]
        self.colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'
        self.map_dic = {
                'proj':ccrs.PlateCarree(),
                'extent':[110.5, 116, 32, 36.5],
                'extent_interval_lat':1,
                'extent_interval_lon':1,
            }
        self.station = {
                'ZhengZhou': {
                    'abbreviation':'郑州',
                    # 'abbreviation':'ZZ',
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



    def draw_single(self, da):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        ax = self.ax
        mp = Map()

        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])

        # for sta in self.station.keys():
        mp.add_station(ax, {'ZhengZhou':self.station['ZhengZhou']}, justice=False, delx=-0.2, marker='o', ssize=15)
        mp.add_station(ax, {'NanYang':self.station['NanYang']}, justice=False, delx=-0.2, marker='s', ssize=15)
        mp.add_station(ax, {'LuShi':self.station['LuShi']}, justice=False, delx=-0.2, marker='x', ssize=15)
        # mp.add_station(ax, self.station['NanYang'], justice=False, delx=-0.2, marker='^')
        # mp.add_station(ax, self.station['LuShi'], justice=False, delx=-0.2, marker='h')

        x = da.lon
        y = da.lat
        
        crx = ax.contourf(x,
                          y,
                          da,
                          corner_mask=False,
                          levels=self.colorlevel,
                          colors = self.colordict,
                          transform=ccrs.PlateCarree()
                          )
        
        return crx
        

class GetData():
    def model_model(self,):
        """原始投影格式的数据，直接相减
        gwd3-nogwd

        Args:
            model (str, optional): _description_. Defaults to 'gwd3'.
        """
        pass

        # dr = Draw()
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
        # picture_dic = {'date':'2021-07 20/00--21/00', 'type':'gwd3-nogwd', 'initial_time':''}
        return da
        # cf = dr.draw_single(da, picture_dic)
        # return cf

        
    def model_obs(self,model):
        """插值过后的数据相减
        gwd3-obs
        gwd0-obs


        Args:
            model (str, optional): _description_. Defaults to 'gwd3'.
        """
        pass

        # dr = Draw()
        if model == 'gwd3':
            flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
        elif model == 'gwd0':
            flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain_latlon.nc'
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
        # picture_dic = {'date':'2021-07 20/00--21/00', 'type':'gwd3-obs', 'initial_time':''}
        # cf = dr.draw_single(da, picture_dic)
        return da
        

    def EC_OBS(self,):
        """EC和观测插值过后的数据相减
        EC-OBS

        Args:
            model (str, optional): _description_. Defaults to 'gwd3'.
        """
        pass

        # dr = Draw()

        flnm = '/home/fengxiang/HeNan/Data/OBS/rain_ec.nc'
        da2 = xr.open_dataarray(flnm)
        da2 = da2.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
        da2 = da2.sel(lat=slice(32,36.5)).sel(lon=slice(110.5,116))
        da2 = da2.sum(dim='time')

        ## 将站点数据重新插值到EC的网格点上
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
        return da




if __name__ == '__main__':

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    # %%
    def get_dr():
        cm = 1/2.54
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.08,0.85,0.85], projection=proj)
        dr = Draw(fig, ax)
        return dr

    ## EC和观测降水差值
    # dr = get_dr()
    # gd = GetData()
    # da = gd.EC_OBS()
    # cf = dr.draw_single(da)
    # cb = dr.fig.colorbar(
    #     cf,
    #     # cax=ax6,
    #     orientation='horizontal',
    #     ticks=dr.colorticks,
    #     fraction = 0.05,  # 色标大小,相对于原图的大小
    #     pad=0.1,  #  色标和子图间距离
    #     )
    # fig_name = 'EC-OBS'
    # dr.fig.savefig(fig_path+fig_name+'.png', rasterized=True)
    
    ## 模式和模式降水差值
    dr = get_dr()
    gd = GetData()
    da = gd.model_model()
    cf = dr.draw_single(da)
    cb = dr.fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小
    title_name = '(a)'
    dr.ax.set_title(title_name, loc='left', fontsize=10,y=0.98)
    fig_name = 'gwd3-no_gwd'

    ## 画降水变化最大的框
    area_left = {
        'lon1':111.6,
        'lon2':112.4,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    ## 伏牛山右边
    area_right = {
        'lon1':112.4,
        'lon2':113.3,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    def ad_patch(area):
        xy = (area['lon1'], area['lat1'])
        width = area['lon2']-area['lon1']
        height = area['lat2']-area['lat1']
        rect = patches.Rectangle(xy=xy, width=width, height=height, edgecolor='red', fill=False, lw=1.5, ) # 左下角的点的位置
        dr.ax.add_patch(rect)
    ad_patch(area_left)
    ad_patch(area_right)

    # mpl.font_manager._rebuild()
    # font1 = {
    #     # 'weight':'black',
    #     'fontweight':'normal',
    #     'size': 10,
    #     'fontstyle':'oblique',
    # }
    dr.ax.text(111.7,33.5, '$\\bf{A}$', color='black', fontweight='bold')
    # dr.ax.text(111.7,33.5, '$\\bf{A}$', color='black',fontdict=font1, fontweight='bold')
    # dr.ax.text(110.7,33.5, '$\\bf{this}$', color='black',)
    aa = dr.ax.text(112.9,33.5, '$\\bf{B}$', color='black')
    # print(type(aa))
    # aa.color = 'blue'
    # aa.font_wight = 200
    ## 画这个300m的等高线
    h_level = [80,300,]
    h = get_hgt()
    CS = dr.ax.contour(h.lon,h.lat, h, levels=h_level, colors='black')    
    manual_locaton = [(111,32.5),(112.5,35), (114,33)]  # 指定线上值的标注位置
    # fmt = {300:'300m',80:'80m'} # 将这个值变为相关的值
    fmt = {300:'300',80:'80'} # 将这个值变为相关的值
    dr.ax.clabel(CS, inline=True, fontsize=10, manual=manual_locaton, fmt=fmt)

    dr.fig.savefig(fig_path+fig_name+'.png')


    ## 模式和观测降水差值 
    # for model in ['gwd3', 'gwd0']:
    #     dr = get_dr()
    #     gd = GetData()
    #     da = gd.model_obs(model)
    #     cf = dr.draw_single(da)
    #     cb = dr.fig.colorbar(
    #         cf,
    #         # cax=ax6,
    #         orientation='horizontal',
    #         ticks=dr.colorticks,
    #         fraction = 0.05,  # 色标大小,相对于原图的大小
    #         pad=0.1,  #  色标和子图间距离
    #         )
    #     cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    #     fig_name = model+'-obs'
    #     dr.fig.savefig(fig_path+fig_name+'.png', rasterized=True)
# %%

# %%
