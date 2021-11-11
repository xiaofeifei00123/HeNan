# %%
import sys,os
import xarray as xr
import numpy as np
import pandas as pd

import salem  # 插值
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import geopandas
import cmaps
from get_cmap import get_cmap_rain2
from multiprocessing import Pool

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_era.nc'
da_rain_all = xr.open_dataarray(flnm)
# da_rain = da_rain_all.sel(time=slice('2021-07-20 08', '2021-07-20 20'))  # 24小时逐小时降水
da_rain = da_rain_all.sel(time=slice('2021-07-19 08', '2021-07-21 08'))  # 24小时逐小时降水
# dd
# da_rain = da_rain_all.isel(dtime=0)
dd = da_rain.sum(dim='time')


# %%


# %%


class Draw(object):

    
    def __init__(self) -> None:
        super().__init__()
        # self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        # self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'
        # self.levels = np.arange(0,301,20) 
        # self.levels = np.arange(0,50,1) 
        # self.levels = np.arange(0,301,20) 
        # self.levels = np.arange(0,531,35) 
        # self.levels = np.arange(0,151,10)
        # self.levels = [0, 10, 25, 50, 100, 250,500,]
        self.levels = [0, 10, 25, 50, 100, 250]
        # self.levels = np.linspace(0,660,16)

    def draw_station(self, ax):
        pass
        # station = station_dic
        station = {
            'ZhengZhou': {
                'abbreviation':'郑州',
                'lat': 34.76,
                'lon': 113.65
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
                #    color='black',
                   color='black',
                   transform=ccrs.PlateCarree(),
                   alpha=1.,
                   linewidth=5,
                   s=30,
                   )
        # 给站点加注释
        for i in range(len(x)):
            # print(x[i])
            ax.text(x[i]-0.2,
                    y[i] + 0.2,
                    station_name[i],
                    transform=ccrs.PlateCarree(),
                    alpha=1.,
                    fontdict={
                    'size': 28,
            })


    def create_map(self, ax):
        """创建地图对象
        ax 需要添加底图的画图对象

        Returns:
            ax: 添加完底图信息的坐标子图对象
        """
        # proj = ccrs.PlateCarree()
        proj = ccrs.PlateCarree()
        # --设置地图属性
        # 画省界
        provinces = cfeat.ShapelyFeature(
            Reader(self.path_province).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')
        
        city = cfeat.ShapelyFeature(
            Reader(self.path_city).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')


        # ax.set_extent([108, 118, 30, 40], crs=proj)
        ax.add_feature(provinces, linewidth=1, zorder=2)  # 添加青藏高原区域
        # ax.add_feature(city, linewidth=0.6, zorder=2)  # 添加青藏高原区域

        # --设置图像刻度
        # ax.set_xticks(np.arange(105, 120 + 2, 1))
        
        # ax.add_feature(provinces, linewidth=1, zorder=2)
        # ax.gridlines(draw_labels=True, dms=True)
        # ax.add_feature(provinces, linewidth=2, zorder=10)
        # ax.add_feature(city, linewidth=1, zorder=2)  # 添加青藏高原区域
        ax.set_yticks(np.arange(30, 40 + 2, 2))
        ax.set_xticks(np.arange(110, 120 + 2, 2))
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.tick_params(axis='both', labelsize=20, direction='out')
        # -- 设置图像范围
        # ax.set_extent([78, 98, 26, 38], crs=ccrs.PlateCarree())
        return ax

    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        # norm = mpl.colors.BoundaryNorm(self.levels, dic['cmap'].N, extend='both')
        # ax = self.create_map(ax)  # 创建坐标图像
        # ax.set_title("jjjj")
        x = data.lon
        y = data.lat
        # contour_levels = [0, 10, 25, 50, 100, 250,500,]
        colorlevel=[0, 0.1,10.0,25.0,50.0,100.0,250.0,400.0, 700]#雨量等级
        colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        # print(x.max())
        # crx = ax.contourf(x,
        #                   y,
        #                   data,
        #                   cmap=dic['cmap'],
        #                 #   norm=norm,
        #                 #   extend='both',
        #                 #   extend='max',
        #                 #   extendfrac='scalar',
        #                 #   levels=levels,
        #                   levels=self.levels,
        #                   transform=ccrs.PlateCarree())

        # print(x.max())
        # sns.kdeplot(data, x=x,y=y)
        
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                        #   levels=levels_rain,
                        #   levels=self.levels,
                            levels = colorlevel,
                            colors = colordict,
                          transform=ccrs.PlateCarree())
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          

        # ax.set_title(title[2], loc='left',  fontsize=12)
        # ax.set_title(dic['name'], loc='left',  fontsize=12)
        # ax.set_extent([70, 105, 25, 41], crs=ccrs.LambertConformal())
        # ax.xaxis.set_tick_params(labelsize=10)
        # ax.yaxis.set_tick_params(labelsize=10)
        # ax.text(99.5, 38.5, title[1])
        # ax.set_title(title[1], loc='right',  fontsize=12)
        # ax.set_tilte(title[1], loc='right')
        # ax.text(78,26.2,title[0], fontsize=12)
        self.draw_station(ax)
        return crx

    def draw_single(self, da, date):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        fig = plt.figure(figsize=(12, 12), dpi=600)
        # proj = ccrs.LambertConformal()  # 创建坐标系
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        # ax.set_extent([])
        ax.set_extent([108, 118, 30, 40], crs=ccrs.PlateCarree())
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
            #    'cmap':get_cmap_rain2(),
               'time':str(date)}
        ax = self.create_map(ax)
        # ax.set_title(date, fontsize=30)
        ax.set_title('OBS', fontsize=20,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        # ax.contourf(da.lon, da.lat, da)
        color_ticks=[0.1,10.0,25.0,50.0,100.0,250.0,400.0]#雨量等级
        
        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            # ticks=self.levels[:-1],
            ticks = color_ticks,
            # label='Precipitation mm',
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和图片间距离
        )
        cb.ax.tick_params(labelsize=18)  # 设置色标标注的大小
        cb.set_label('Precipitation (mm)', fontdict={'size':20})
        # # fig_name = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_12/'+date+'_YSU_1912_ERAI.png'
        # fig_name = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_12/'+date+'ERAsmall.png'
        fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/test.png')
        # fig.savefig(fig_name)
        # pass

# dd
dr = Draw()
# dr.draw_single(dd, '20_08-20_20')
dr.draw_single(dd, '19_08-21_00')