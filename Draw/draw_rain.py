#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画郑州暴雨降水分布图
-----------------------------------------
Time             :2021/09/09 16:19:22
Author           :Forxd
Version          :1.0
'''

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


def nearest_position(stn_lat, stn_lon, xlat, xlon ):
    """获取最临近格点坐标索引, 距离目标站点最近的格点坐标
    根据二维的经纬度数据，获取某个点的索引值
    parameters:
    ----------------------
    stn_lon  : 单个站点经度
    stn_lat  : 单个站点纬度
    xlon : numpy.ndarray网格二维经度坐标
    xlat : numpy.ndarray网格二维纬度坐标 
    Return:
    -----------------------
    (xindx,yindx), 可以直接使用
    """
    difflat = stn_lat - xlat;
    difflon = stn_lon - xlon;
    rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)#difflat * difflat + difflon * difflon
    aa=np.where(rad==np.min(rad))
    # print(aa)
    #print(stn_lat, stn_lon)
    ind = np.squeeze(np.array(aa))
    print('距离该站点最近的格点经纬度索引为：',ind)  
    return tuple(ind)



# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912_rain_ERAI.nc'
ds_all = xr.open_dataset(flnm)
## 获取降水数据
da_rain = ds_all['RAINNC']
# da_rain = da_rain.rename({'XLAT':'lat','XLONG':'lon', 'XTIME':'time'})
da_rain = da_rain.rename({'XLAT':'lat','XLONG':'lon',})
# da_rain = da_rain.swap_dims({'Time':'time'})
## 改成北京时间
tt = da_rain.time+pd.Timedelta(hours=+8)
da_rain = da_rain.assign_coords({'time':('time',tt.values)})


## 筛选单个时次的数据
# da = da_rain.sel(time='2021-07-20 1600').squeeze()

# # da.lat
#                 # 'lat': 34.75,
#                 # 'lon': 113.62
# station = {
#     'ZhengZhou': {
#         'abbreviation':'郑州',
#         'lat': 34.75,
#         'lon': 113.62
#     },
# }
# stn_lat = station['ZhengZhou']['lat']
# stn_lon = station['ZhengZhou']['lon']

# cc = nearest_position(stn_lon=stn_lon, stn_lat=stn_lat,xlat=da.lat, xlon=da.lon)
# da[cc]






# %%
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
class Draw(object):

    
    def __init__(self) -> None:
        super().__init__()
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'
        self.levels = np.arange(0,301,20) 
        # self.levels = np.arange(0,531,35) 
        # self.levels = np.arange(0,151,10)
        # self.levels = [0, 10, 25, 50, 100, 250,500,]
        # self.levels = np.linspace(0,660,16)

    def draw_station(self, ax):
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

    def draw_patch(self, ax):
        """画矩形框
        """
        area = [None] * 3  # 设置一个维度为8的空列表
        area[0] = {"lat1": 33.5, "lat2": 40, "lon1": 80, "lon2": 90}  # north
        area[1] = {"lat1": 28, "lat2": 33,
                   "lon1": 83, "lon2": 94}  # south left
        area[2] = {"lat1": 26, "lat2": 33,
                   "lon1": 95, "lon2": 103}  # south right
        for i in range(3):
            lon = np.empty(4)
            lat = np.empty(4)
            lon[0], lat[0] = area[i]['lon1'], area[i]['lat1']
            lon[1], lat[1] = area[i]['lon2'], area[i]['lat1']
            lon[2], lat[2] = area[i]['lon2'], area[i]['lat2']
            lon[3], lat[3] = area[i]['lon1'], area[i]['lat2']
            x, y = lon, lat
            xy = list(zip(x, y))
            poly = plt.Polygon(xy, edgecolor="red", fc="none", lw=.9, alpha=1)
            ax.add_patch(poly)

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
        # ax.add_feature(provinces, linewidth=1, zorder=2)  # 添加青藏高原区域
        # ax.add_feature(city, linewidth=0.6, zorder=2)  # 添加青藏高原区域

        # --设置图像刻度
        ax.set_xticks(np.arange(105, 120 + 2, 1))
        
        # ax.add_feature(provinces, linewidth=1, zorder=2)
        # ax.gridlines(draw_labels=True, dms=True)
        ax.add_feature(provinces, linewidth=2, zorder=10)
        # ax.add_feature(city, linewidth=1, zorder=2)  # 添加青藏高原区域
        ax.set_yticks(np.arange(30, 40 + 2, 1))
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
        crx = ax.contourf(x,
                          y,
                          data,
                          cmap=dic['cmap'],
                        #   norm=norm,
                        #   extend='both',
                          extend='max',
                        #   extendfrac='scalar',
                        #   levels=levels,
                          levels=self.levels,
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
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
            #    'cmap':get_cmap_rain2(),
               'time':str(date)}
        ax = self.create_map(ax)
        ax.set_title(date, fontsize=30)
        ax.set_title('ERA5_1912', fontsize=20,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        # ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

        
        # bounds = [0,10, 25, 50, 100, 250, 300, 500]
        

        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=self.levels[:-1],
            # label='Precipitation mm',
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和图片间距离
        )
        cb.ax.tick_params(labelsize=18)  # 设置色标标注的大小
        cb.set_label('Precipitation (mm)', fontdict={'size':20})
        fig_name = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_12/'+date+'_YSU_1912_ERAI.png'
        # fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/test.png')
        fig.savefig(fig_name)
        pass

    def draw_dual(self, rain):
        shp = geopandas.read_file(self.path_tibet)
        model_list = ['QNSE', 'QNSE_EDMF', 'TEMF', 'ACM2', 'YSU']
        # --->画图
        proj = ccrs.LambertConformal()  # 创建坐标系
        fig = plt.figure(figsize=(9, 9), dpi=400)  # 创建页面
        grid = plt.GridSpec(3,
                            2,
                            figure=fig,
                            left=0.07,
                            right=0.96,
                            bottom=0.12,
                            top=0.96,
                            wspace=0.2,
                            hspace=0.1)

        axes = [None] * 6  # 设置一个维度为8的空列表
        axes[0] = fig.add_subplot(grid[0, 0:1], projection=proj)
        axes[1] = fig.add_subplot(grid[0, 1:2], projection=proj)
        axes[2] = fig.add_subplot(grid[1, 0:1], projection=proj)
        axes[3] = fig.add_subplot(grid[1, 1:2], projection=proj)
        axes[4] = fig.add_subplot(grid[2, 0:1], projection=proj)
        axes[5] = fig.add_subplot(grid[2, 1:2], projection=proj)

        # cmap = get_cmap_rain3()
        cmap = cmaps.precip_16v

        time_2005 = '2800_2906 Jul 2005'
        time_2014 = '1900_2000 Aug 2014'

        self.draw_contourf(rain['ACM2'].salem.roi(
            shape=shp), axes[0], cmap, ['ACM2', '(a)', time_2005])
        self.draw_contourf(rain['ACM2'].salem.roi(
            shape=shp), axes[1], cmap, ['YSU', '(b)', time_2005])
        self.draw_contourf(rain['QNSE'].salem.roi(
            shape=shp), axes[2], cmap, ['QNSE', '(c)', time_2005])
        self.draw_contourf(rain['QNSE_EDMF'].salem.roi(
            shape=shp), axes[3], cmap, ['QNSE_EDMF', '(d)', time_2005])
        self.draw_contourf(rain['TEMF'].salem.roi(
            shape=shp), axes[4], cmap, ['TEMF', '(e)', time_2005])
        cf = self.draw_contourf(rain['obs'].salem.roi(
            shape=shp), axes[5], cmap, ['OBS', '(f)', time_2005])
        ax6 = fig.add_axes([0.18, 0.06, 0.7, 0.02])  # 重新生成一个新的坐标图

        # cb = fig.colorbar(
        #     cf,
        #     cax=ax6,
        #     orientation='horizontal',
        #     ticks=self.levels,
        #     # fraction = 0.1,  # 色标大小
        #     pad=0.1,  #
        #     fontsize=30,
        #     label='preciptation'
        # )
        path = os.path.join(self.picture_path, self.fig_name)
        fig.savefig(path)

def get_sum_rain(da_rain):
    dd = da_rain.sel(time=slice('2021-07-20 08', '2021-07-20 20'))  # 24小时逐小时降水
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


    
if __name__ == '__main__':
    pass
    cc = get_sum_rain(da_rain)
    dr = Draw()
    dr.draw_single(cc, '20_08-20_20')
    # main()
    # dr = Draw()
    # aa = da_rain.time.to_series()
    # for i in aa:
    # dd = da_rain.sel(time=slice('2021-07-19 2000', '2021-07-20 2000'))
    # dd = dd.sum(dim='time')

    # print(i.strftime('%d_%H'))
    # dr.draw_single(dd, '1920-2020')
    # print(dd)
    # %%
    ### 单进程结束

    ## 多进程
    # dr = Draw()
    # time_series = da_rain.time.to_series()
    # pool = Pool(6)
    # result = []
    # # for i in range(50):
    # for i in time_series:
    #     dd = da_rain.sel(time=i)
    #     tr = pool.apply_async(dr.draw_single, args=(dd, i.strftime('%d_%H'),))
    #     # print(i.strftime('%d_%H'))
    #     # print("计算%d"%i)
    #     # result.append(tr)
    # pool.close()
    # pool.join()

    # c = []
    # for j in result:
    #     c.append(j.get())
    # print(c)
    ### 多进程结束

