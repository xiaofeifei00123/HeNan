#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取史文茹她们的数据
-----------------------------------------
Time             :2021/11/12 09:14:06
Author           :Forxd
Version          :1.0
'''

# %%
import pandas as pd
import xarray as xr
import numpy as np
import time, datetime
import os
from baobao.interp import rain_station2grid


# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain-henan/rain2/rain-2021071900CST.txt'



class RainStation():

    def rain1h(self,flnm):
        """获得这一个小时的所有站点降水数据
        Returns:
            DataArray : 
        """
        ## 获得数据时间
        df = pd.read_csv(flnm, delim_whitespace=True, nrows=0)
        time1 = df.columns[1]
        time2 = datetime.datetime.strptime(time1,'%Y%m%d%H')
        # time3 = pd.DatetimeIndex([time2]).values[0]+pd.Timedelta('6H')
        # time3 = pd.DatetimeIndex([time2]).values[0]
        time3 = pd.DatetimeIndex([time2]).values[0]+pd.Timedelta('-8H')

        df = pd.read_csv(flnm, delim_whitespace=True, header=1, na_values=99999, usecols=['station','lon', 'lat', '1hrain'])
        da = xr.DataArray(
            df['1hrain'].values,
            coords={
                'station':df['station'],
                'lat':('station',df['lat']),
                'lon':('station',df['lon']),
                'time':time3,
            },
            dims=['station',]
        )
        return da


    def rainall(self,):
        path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain-henan/rain2'
        fl_list = os.popen('ls {}/rain*.txt'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        dds_list = []
        for fl in fl_list:
            print(fl)
            da = self.rain1h(fl)
            dds_list.append(da)

        ## 针对micaps数据的各个站点数据进行聚合
        da_concat = xr.concat(dds_list, dim='time')
        lat = da_concat['lat'].mean(dim='time')  # 将多列数据，变成一列
        lon = da_concat['lon'].mean(dim='time')
        dda = da_concat.drop_vars(['lat', 'lon'])
        daa = dda.assign_coords({'lat':('station',lat.values), 'lon':('station',lon.values)})
        # dc = daa.fillna(0)
        dc = daa.rename({'station':'id'}).dropna(dim='id')
        
        return dc

    def save_rain(self,):
            rain_st = self.rainall()
            ## 保存所有站点的数据
            rain_st.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc') # 所有站点
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
    rs.save_rain()
    # rg = RainGrid()
    # rg.save_rain_grid()
# %%
if __name__ == '__main__':

    save_rain()