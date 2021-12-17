#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取观测站点降水，并插值
观测站降水为0时，该时次观测站数据没有显示
需要对这些数据进行添加
保存：
    插值的数据，nc搁谁
    未插值的数据, csv
-----------------------------------------
Time             :2021/09/28 20:39:34
Author           :Forxd
Version          :1.0
'''

# %%

# from HeNan.Draw.read_ERA5 import get_rain_station
import xarray as xr
import meteva.base as meb
import numpy as np
import os
import pandas as pd
# from metpy.interpolate import inverse_distance_to_grid
# from metpy.interpolate import interpolate_to_grid
# from baobao.quick_draw import quick_contourf,quick_contourf_station
from baobao.interp import rain_station2grid
# %%


path='/mnt/zfm_18T/fengxiang/HeNan/bb/2017081500ref/level_dataref1.txt'
data_file = path
# df = pd.read_csv(path,sep='\s', col=['lon', 'lat', 'height', 'val'])
col_names=['lon', 'lat', 'height', 'val']
df = pd.read_table(
    data_file,  # 由字符串虚拟的文件
    sep='\\s+',
    skiprows=1,
    usecols=[0,1,2,3,],
    names=col_names,
    na_values= -9999.000
)
ddf = df.dropna(axis=0, how='any')
d1 = ddf['val']
# %%
# ddf
rain = ddf[['lon', 'lat','val']]
rain
# d1
# %%
# len(rain.lon)
# len(rain.lon.dropna())
# rain.lat.values
# rain
# %%
# d1.values
import matplotlib.pyplot as plt
from baobao.map import Map
import cartopy.crs as ccrs
import cmaps

def draw_tricontourf(rain):

    """rain[lon, lat, data],离散格点的DataArray数据

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    from nmc_met_graphics.plot import mapview
    mb = mapview.BaseMap()
    fig = plt.figure(figsize=[12, 12], dpi=600)
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.LambertConformal(central_latitude=34, central_longitude=113))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
    ax = fig.add_axes([0.1, 0.1, 0.85, 0.85], projection=ccrs.PlateCarree())
    # ax.plot(rain.lon.values, rain.lat.values, 'ko',ms=3,zorder=2, transform=ccrs.PlateCarree())
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        # 'extent':[110.5, 116, 32, 36.5],
        'extent':[119, 123,30 , 33],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
    # colorlevel=[-20, -15, -10, -5.0, 0,  5, 10, 15]#雨量等级
    colorlevel=np.arange(-16,64,8)
    # colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
    colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000',]#颜色列表
    cs = ax.tricontour(rain.lon.values, rain.lat.values, rain.val.values, levels=colorlevel, transform=ccrs.PlateCarree())
    # cs = ax.tricontourf(rain.lon.values, rain.lat.values, rain.val.values, levels=colorlevel,colors=colordict, transform=ccrs.PlateCarree())
    cs = ax.tricontourf(rain.lon.values, rain.lat.values, rain.val.values, levels=15,cmap=cmaps.radar, transform=ccrs.PlateCarree())
    # cs = ax.tricontour(rain.lon, rain.lat, rain, levels=colorlevel, transform=ccrs.PlateCarree())
    # cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=colorlevel,colors=colordict, transform=ccrs.PlateCarree())
    # # colorticks=[0.1,5,15,30.0,70,140]#雨量等级
    # colorticks = colorlevel[1:-1]
    # cb = fig.colorbar(
    #     cs,
    #     # cax=ax6,
    #     orientation='horizontal',
    #     ticks=colorticks,
    #     fraction = 0.05,  # 色标大小,相对于原图的大小
    #     pad=0.05,  #  色标和子图间距离
    # )
    # # mb.cities(ax, city_type='base_station', color_style='black', 
    # #             marker_size=5, font_size=16)
    # # mb.gridlines()
    # cb.ax.tick_params(labelsize=30)  # 设置色标标注的大小



    # mp = Map()
    # station = {
    #     'ZhengZhou': {
    #         'abbreviation':'ZZ',
    #         'lat': 34.76,
    #         'lon': 113.65
    #     },
    # }
    # mp.add_station(ax, station, justice=True)


    # ax.set_title('2021-07-20 00-12', fontsize=35,)
    # # ax.set_title('2021-07 19/16--20/04', fontsize=35,)
    # ax.set_title('OBS', fontsize=30,loc='left')


    # # mp.add_station(ax)
    # fig_name = 'obs_distribution' 
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_12h_gwd/'
    # fig.savefig(fig_path+fig_name)


draw_tricontourf(rain)



# %%
class RainStation():
    """站点数据聚合到一起"""

    def read_one_station(self, flnm):
        station = meb.read_stadata_from_micaps3(flnm)
        ## 转换为世界时
        station['time'] = station['time']+pd.Timedelta('-8H')

        df = station

        ## 将一个站点的数据，变为DataArray
        da = xr.DataArray(
            df['data0'].values,
            coords={
                'id':df['id'],
                'lat':('id',df['lat']),
                'lon':('id',df['lon']),
                'time':df['time'].values[0],
            },
            dims=['id',]
        )
        return da

    def get_rain_station(self,):
        # path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station'
        path='/mnt/zfm_18T/fengxiang/HeNan/bb/2017081500ref'
        fl_list = os.popen('ls {}/level*.txt'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        dds_list = []
        for fl in fl_list:
            da = self.read_one_station(fl)
            # print(da)
            dds_list.append(da)

        ## 针对micaps数据的各个站点数据进行聚合
        da_concat = xr.concat(dds_list, dim='time')
        lat = da_concat['lat'].mean(dim='time')  # 将多列数据，变成一列
        lon = da_concat['lon'].mean(dim='time')
        dda = da_concat.drop_vars(['lat', 'lon'])
        daa = dda.assign_coords({'lat':('id',lat.values), 'lon':('id',lon.values)})
        # dc = daa.fillna(0)
        # dc = daa.dropna(dim='id')
        dc = daa
        # print()
        
        return dc
        
        
    def save_rain_station(self,):
        ## 保存成micaps3格式, 保存站点数据
        # rs = rain_station()
        rain_st = self.get_rain_station()
        # print(rain_st)
        ## 保存所有站点的数据
        rain_st.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc') # 所有站点

        ## 筛选出河南范围内的站点数据
        area = {
            'lon1':110.5,
            'lon2':116,
            'lat1':32,
            'lat2':36.5,
            'interval':0.125,
        }
        # area = {
        #     'lon1':111.5,
        #     'lon2':113.5,
        #     'lat1':33.5,
        #     'lat2':35,
        #     'interval':0.125,
        # }
        da = rain_st
        index = ((da.lat<=area['lat2']) & (da.lat>=area['lat1']) & (da.lon>=area['lon1']) & (da.lon<=area['lon2']))
        da_obs = da.loc[:,index]  # 这里时间维度在前
        print(da_obs)
        da_obs.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc') # 所有站点
        
class RainGrid():
    """格点数据

    Returns:
        [type]: [description]
    """
    pass
    def rain_station2grid(self, da):
        """将站点数据，插值为格点数据

        Args:
            da ([type]): [description]

        Returns:
            [type]: [description]
        """
        area = {
            'lon1':110,
            'lon2':116,
            'lat1':31,
            'lat2':37,
            'interval':0.125,
        }
        ddc = rain_station2grid(da, area)
        return ddc

    def save_rain_grid(self,):
        ## 读取存储好的站点数据
        da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc')
        ## 插值为格点数据
        da_grid = self.rain_station2grid(da)
        da_grid.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc')

def save_rain():
    rs = RainStation()
    rs.save_rain_station()
    # rg = RainGrid()
    # rg.save_rain_grid()
            
        

# %%
if __name__ == '__main__':
    pass
    save_rain()
    
