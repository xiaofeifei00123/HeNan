#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取观测站点降水，并插值
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
# %%
class rain_station_grid():
    """格点数据

    Returns:
        [type]: [description]
    """
    pass
    def read_one_station(self,flnm):
        station = meb.read_stadata_from_micaps3(flnm)
        grid1 = meb.grid([110-1,116+1,0.125],[32-1,37+1,0.125],)
        grd2 = meb.interp_sg_idw(station, grid1, effectR=60, nearNum=3)
        da_return = grd2.astype('float32')
        ## 转换为世界时
        da_return['time'] = da_return.time+pd.Timedelta('-8H')
        return da_return

    def get_rain_station(self,):
        path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station'
        fl_list = os.popen('ls {}/2021*.000'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        dds_list = []
        for fl in fl_list:
            da = self.read_one_station(fl)
            dds_list.append(da)
        dds_concate = xr.concat(dds_list, dim='time')
        dds_concate
        return dds_concate.squeeze()



class rain_station():
    """站点数据"""
    pass
    def read_one_station(self,flnm):
        station = meb.read_stadata_from_micaps3(flnm)
        ## 转换为世界时
        station['time'] = station['time']+pd.Timedelta('-8H')
        return station

    def get_rain_station(self,):
        path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station'
        fl_list = os.popen('ls {}/2021*.000'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        dds_list = []
        for fl in fl_list:
            da = self.read_one_station(fl)
            dds_list.append(da)
        df_concat = pd.concat(dds_list, axis=0)
        # dds_concate = xr.concat(dds_list, dim='time')
        # dds_concate
        return df_concat


def save_rain_grid():
    ## 读取并插值, 保存格点数据
    grs = rain_station_grid()
    rain_grid = grs.get_rain_station()
    rain_grid.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    
def save_rain_station():
    ## 保存成micaps3格式, 保存站点数据
    rs = rain_station()
    rain_st = rs.get_rain_station()
    rain_st.to_csv('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.csv')

def get_max_dataframe(df_station):
    """读取csv格式的站点数据

    Args:
        df_station (DataFrame): 输入csv格式站点数据

    Returns:
        [DataArray]: 小时雨强最大值
    """
    ## 才读的数据，它的时间格式是object
    df_station['time']= pd.to_datetime(df_station['time'])
    t = pd.date_range(start='2021-07-20 00', end='2021-07-21 00', freq='1H')
    rain_list = []
    for tt in t:
        cc = df_station[df_station['time']==tt]
        rain_max = cc[(cc['lat']>32)&(cc['lat']<37)&(cc['lon']>110)&(cc['lon']<116)]['data0'].max()
        rain_list.append(rain_max)
    rain_list
    ps = pd.Series(rain_list, index=t)
    da = xr.DataArray.from_series(ps)
    rain_obs_max = da.rename({'index':'time'})
    return rain_obs_max

# %%
if __name__ == '__main__':
    pass
    save_rain_grid()
    # save_rain_station()
    
