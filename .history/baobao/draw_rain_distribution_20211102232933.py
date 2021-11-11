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
# from get_cmap import get_cmap_rain2
from multiprocessing import Pool

from baobao.map import Map


# %%



# %%
class GetData():

    def get_rain_obs(self,):
        """获得站点插值后的降水

        Returns:
            [type]: [description]
        """
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_grid.nc'
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/baobao/rain_grid.nc'
        da_rain_all = xr.open_dataarray(flnm)
        # da_rain = da_rain_all.sel(time=slice('2021-07-20 00', '2021-07-20 12'))  # 24小时逐小时降水
        # dd = da_rain.sum(dim='time')
        dd = da_rain_all
        return dd
    
    def get_rain_ec(self,):
        """获得EC细网格的降水

        Returns:
            [type]: [description]
        """
        pass
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_ec.nc'
        da = xr.open_dataarray(flnm)
        da_rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12'))  # 24小时逐小时降水
        dd = da_rain.sum(dim='time')
        return dd

    def get_sum_rain(self, da_rain):
        """计算wrfout降水数据，某时间段内的降水和
        涉及到重新处理经纬度坐标的问题

        Args:
            da_rain ([type]): [description]

        Returns:
            [type]: [description]
        """
        dd = da_rain.sel(time=slice('2021-07-20 00', '2021-07-20 12'))  # 24小时逐小时降水
        da_sum = dd.sum(dim='time')  # 24小时总降水量
        lat = dd.isel(time=0).lat
        lon = dd.isel(time=0).lon

        ds_sum = xr.Dataset(
            {'RAINNC':(['south_north', 'west_east'], da_sum.values)},
            coords={
                'lon':(['south_north','west_east'], lon),
                'lat':(['south_north', 'west_east'], lat),
                'time':pd.Timestamp('2021-07-20 0008')
            },
        )
        return ds_sum['RAINNC']

    def get_rain_wrf_gfs(self,):
        """使用gfs的预报数据作为初始场，使用wrf转出的数据

        Returns:
            [type]: [description]
        """
        pass
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/YSU_GFS.nc'
        ds = xr.open_dataset(flnm)
        da_rain = ds['RAINNC']
        dda = da_rain.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
        da = dda.swap_dims({'Time':'time'})
        # da_rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12'))  # 24小时逐小时降水
        # dd = da_rain.sum(dim='time')
        dd = self.get_sum_rain(da)
        return dd 

    def get_rain_wrf(self, flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/YSU_rain_1km.nc'):
        """
        特定试验的降水数据
        使用ERA5和GDAS数据， 不同起报时间，转出的wrf降水数据

        Returns:
            [type]: [description]
        """
        # pass
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912_rain.nc'
        ds = xr.open_dataset(flnm)
        da_rain = ds['RAINNC']
        dd = self.get_sum_rain(da_rain)
        return dd 


    def get_rain_yyx(self,):
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/YYX/wrfout_d02_2021-07-19_12:00:00'
        ds = xr.open_dataset(flnm)
        da_rain = ds['RAINNC']
        dda = da_rain.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
        dc = dda.swap_dims({'Time':'time'})
        rain_list = []
        tt = dc.time
        # i = 0
        r = 0
        for t in tt:
            rr = dc.sel(time=t)  # 时间t时刻的累计降水
            rain = rr-r
            r = rr.values
            rain_list.append(rain.squeeze())
        rain_all = xr.concat(rain_list, pd.Index(tt.values, name='time'))

        da_rain = rain_all.sel(time=slice('2021-07-20 00', '2021-07-20 12'))  # 24小时逐小时降水
        dd = self.get_sum_rain(da_rain)
        return dd

class Draw(object):

    def __init__(self) -> None:
        super().__init__()
        # self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        # self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'
        # self.levels = np.arange(0,400,40)
        # self.levels=[0.1,10.0,25.0,50.0,100.0,250.0,500.0]#雨量等级
        # self.levels=[0, 0.1,10.0,25.0,50.0,100.0,250.0,400.0, 600]#雨量等级


    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        # norm = mpl.colors.BoundaryNorm(self.levels, dic['cmap'].N, extend='both')
        # ax = self.create_map(ax)  # 创建坐标图像
        # ax.set_title("jjjj")
        # levels_rain = [0, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400]
        # levels_rain = np.arange(0,400,40)
        # colorlevel=[0.1,10.0,25.0,50.0,100.0,250.0,500.0]#雨量等级
        # colorlevel=[0, 0.1,10.0,25.0,50.0,100.0,250.0,400.0, 700]#雨量等级
        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 400.0, 700]#雨量等级
        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 250.0, 700]#雨量等级
        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140,  700]#雨量等级
        colorlevel=[0, 0.1, 10, 20, 30, 40, 140,  700]#雨量等级
        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
        colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        # colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040']#颜色列表
        x = data.lon
        y = data.lat
        # print(x.max())
        # sns.kdeplot(data, x=x,y=y)
        
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                        #   levels=levels_rain,
                          levels=colorlevel,
                            colors = colordict,
                        #   colors=['33ff66', '33ff66', '33ff66', '339900', 'ffff00', 'cc9900', 'ff0000'],
                        #   colors=['white','#00ffff','#00ff99', '#009933', '#ffff00','#cc9900', '#ff0000','#ff00ff'],#,'#EE0000'],
                        #   cmap = cmaps.precip3_16lev,

                        #   cmap = 'Greys',
                        #   extend='both',
                        #   vmin=10,
                          transform=ccrs.PlateCarree())
        # # ax.set_title(title[2], loc='left',  fontsize=12)
        # # ax.set_title(dic['name'], loc='left',  fontsize=12)
        # # ax.set_extent([70, 105, 25, 41], crs=ccrs.LambertConformal())
        # # ax.xaxis.set_tick_params(labelsize=10)
        # # ax.yaxis.set_tick_params(labelsize=10)
        # # ax.text(99.5, 38.5, title[1])
        # # ax.set_title(title[1], loc='right',  fontsize=12)
        # # ax.set_tilte(title[1], loc='right')
        # # ax.text(78,26.2,title[0], fontsize=12)
        # self.draw_station(ax)
        return crx

    def draw_single(self, da, picture_dic):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        fig = plt.figure(figsize=(12, 12), dpi=600)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        # ax.set_extent([])
        date = picture_dic['date']
        # ax.set_extent([110, 116, 32, 37], crs=ccrs.PlateCarree())
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
            #    'cmap':get_cmap_rain2(),
               'time':str(date)}

               

        # ax = self.create_map(ax)
        mp = Map()
        map_dic = {
            'proj':ccrs.PlateCarree(),
            # 'extent':[110, 116, 32, 37],
            'extent':[116, 123, 30, 35],
            'extent_interval_lat':1,
            'extent_interval_lon':1,
        }
        ax = mp.create_map(ax, map_dic)

        station = {
            'ZhengZhou': {
                'abbreviation':'ZZ',
                'lat': 34.76,
                'lon': 113.65
            },
        }
        mp.add_station(ax, station, justice=True)
        
        

        ax.set_title(date, fontsize=35,)
        ax.set_title(picture_dic['initial_time'], fontsize=30,loc='left')
        ax.set_title(picture_dic['type'], fontsize=30,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        # ax.contourf(da.lon, da.lat, da)
        # colorticks = [0.1,10.0,25.0,50.0,100.0,250.0,400.0]#雨量等级
        # colorticks=[0.1, 5, 15.0, 30, 70, 140, 250.0]#雨量等级
        # colorticks=[0.1, 5, 15.0, 30, 70, 140]#雨量等级
        # colorticks=[0.1,5,15,30.0,70,140,250]#雨量等级
        colorticks=[0.1,5,15,30.0,70,140]#雨量等级
        
        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            # ticks=self.levels,
            ticks=colorticks,
            # label='Precipitation mm',
            # shrink=0.5, # 按比例缩放colorbar的大小
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.05,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=30)  # 设置色标标注的大小
        # cb.set_label('Precipitation (mm)', fontdict={'size':20})
        # # fig_name = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_12/'+date+'_YSU_1912_ERAI.png'
        # fig_name = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_12/'+date+'ERAsmall.png'
        fig_name = picture_dic['type']+'_'+picture_dic['initial_time']
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/'
        fig.savefig(fig_path+fig_name)
        # fig.savefig(fig_name)
        # pass

def draw_obs():
    pass

    dr = Draw()
    gd = GetData()
    ## obs
    da = gd.get_rain_obs()
    da = da.squeeze()
    picture_dic = {'date':'2000_2012', 'type':'OBS_baobao', 'initial_time':''}
    dr.draw_single(da, picture_dic)



def draw_one():
    pass

    dr = Draw()
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1km_1912_hgt/'
    ds = xr.open_dataset(flnm)
    da = ds['RAINNC'].sum(dim='time')
    ### latlon数据结束
    

    t = '1900'
    picture_dic = {'date':'2021-07-20 00-12', 'type':'90m', 'initial_time':t}
    dr.draw_single(da, picture_dic)



if __name__ == '__main__':

    # draw_yyx()
    # draw_one()
    draw_obs()
    # draw_forecast()
    # gd = GetData()
    # dd = gd.get_rain_obs()
    # gd = GetData()
    # da = gd.get_rain_wrf()
    # dr = Draw()
    # dr.draw_single(da, '2000_2012')
# %%
