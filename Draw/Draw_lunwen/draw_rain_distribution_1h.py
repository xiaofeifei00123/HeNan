#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图, 小时降水
实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
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
# import cartopy.feature as cfeat
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# from cartopy.io.shapereader import Reader, natural_earth
# import matplotlib as mpl
# from matplotlib.path import Path
# import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import cmaps
# from multiprocessing import Pool

from baobao.map import Map



# %%
class Draw(object):

    def __init__(self) -> None:
        super().__init__()
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_micaps/City.shp'
        # self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        # self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'


    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140,  700]#雨量等级
        colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
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
        cm = 1/2.54
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        # ax.set_extent([])
        print("画{}时刻的图".format(str(picture_dic['date'])))
        date = picture_dic['date']
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
               'time':str(date)}
               
        # ax = self.create_map(ax)
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
        }
        mp.add_station(ax, station, justice=True)
        # ax.set_title(date, fontsize=10,)
        # ax.set_title(picture_dic['initial_time'], fontsize=10,loc='left')
        # ax.set_title(picture_dic['type'], fontsize=10,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        # colorticks=[0.1,5,15,30.0,70,140]#雨量等级
        colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
        colorticks = colorlevel[1:-1]
        
        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
        fig_name = 'rain1h'+picture_dic['type']
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
        fig.savefig(fig_path+fig_name)


def draw_one(model='gwd0'):
    pass

    dr = Draw()
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'
    da = xr.open_dataarray(flnm)
    da = da.sel(time=slice('2021-07-20 12', '2021-07-20 12'))
    # da = da.sel(time=slice('2021-07-19 16', '2021-07-20 04 '))
    # da = da.sel(time=slice('2021-07-19 17', '2021-07-20 05 '))
    # da = da.sum(dim='time') 
    tt = da.time
    for t in tt:
        rain1h = da.sel(time=t) 
        picture_dic = {'date':t.dt.strftime("%Y%m%d-%H").values, 'type':model, 'initial_time':''}
        dr.draw_single(rain1h, picture_dic)

def draw_dual():
    model_list = ['gwd0', 'gwd3']
    # model_list = ['gwd0','gwd1', 'gwd3', 'gwd3-LS', 'gwd3-BL', 'gwd3-SS', 'gwd3-FD']
    # model_list = ['1912_90m_gwd3']
    for model in model_list:
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/'
        draw_one(model)


if __name__ == '__main__':

    # draw_one()
    draw_dual()
    # draw_obs()
    # draw_forecast()
    # gd = GetData()
    # dd = gd.get_rain_obs()
    # gd = GetData()
    # da = gd.get_rain_wrf()
    # dr = Draw()
    # dr.draw_single(da, '2000_2012')
# %%
