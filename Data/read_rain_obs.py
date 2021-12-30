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
        path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station'
        fl_list = os.popen('ls {}/2021*.000'.format(path))  # 打开一个管道
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
            'interval':0.05,
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
            'lon1':110.5,
            'lon2':116,
            'lat1':32,
            'lat2':36.5,
            'interval':0.01,
        }
        ddc = rain_station2grid(da, area)
        ## 这个插值出来的数据太不靠谱了
        return ddc

    def save_rain_grid(self,):
        ## 读取存储好的站点数据
        da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc')
        ## 只保留24小时累积降水的插值
        da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
        # da = da.assign_coords({'time':[pd.Timestamp('2021-07-20 00'),]})
        da = da.expand_dims(dim='time')
        tt = pd.date_range('2021-07-20 00', '2021-07-20 00', freq='1H')
        da = da.assign_coords({'time':tt}) 
        print(da.time)
        # ## 插值为格点数据
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
    
