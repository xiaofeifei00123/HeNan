#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取wrfout数据中的降水
将wrfout数据中的降水集合成一个文件
# TODO 将wrfout的格点数据插值为站点数据
-----------------------------------------
Time             :2021/09/09 10:53:08
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import os
import xesmf as xe
import numpy as np
from baobao.interp import regrid_xesmf
from baobao.coord_transform import xy_ll
# import wrf



# %%
def combine_rain(path):
    """
    由于wrfout数据的降水是累计降水，
    这里将它变为逐小时降水,同时进行合并
    又由于wrfout数据的坐标是x,y格点上的，
    通过wrf-python 库将其转为不规则的latlon格点坐标

    Args:
        path ([type]): 包含有wrfout数据的文件夹路径

    Returns:
        rain[DataArray] : 多时次聚合后的降水 
    """    

    fl_list = os.popen('ls {}/wrfout_d04*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    dds_list = []
    r = 0
    for fl in fl_list:
        print(fl[-18:])
        ds = xr.open_dataset(fl)
        da = ds['RAINNC']-r
        r = ds['RAINNC'].values

        dda = da.squeeze()  # 该是几维的就是几维的
        dc = dda.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
        dds_list.append(dc)
    da_concate = xr.concat(dds_list, dim='time') # 原始的， 未经坐标变化的降水

    ## 改变坐标, add latlon coords
    latlon = xy_ll(fl)  # 获得所有点的经纬度坐标
    rain = xr.DataArray(
                da_concate.values.round(1), # 降水精确到1位小数
                coords={
                    'time':da_concate.time.values,
                    'lat':latlon['lat'].round(3),  # 经纬度精确到3位
                    'lon':latlon['lon'].round(3),
                },
                dims=['time', 'lat', 'lon']
            )
    # ds.to_netcdf(path_dic['path_rain'])
    return rain

# %%
def grid2station(flnm_obs, flnm_wrf, area):
    """将格点数据插值为站点数据, 
    都是多时间尺度的数据

    Returns:
        [type]: [description]
    """

    # 特定区域范围内观测的站点数据
    da_obs = xr.open_dataarray(flnm_obs)
    # wrf输出后转化的latlon格点数据
    da_wrf = xr.open_dataarray(flnm_wrf)

    ## 构建插值区域
    cc = da_obs.isel(time=0)
    sta = cc.id.values
    lon = xr.DataArray(cc.lon.values, coords=[sta], dims=['sta'])
    lat = xr.DataArray(cc.lat.values, coords=[sta], dims=['sta'])
    ## 插值
    # rr = da_wrf.interp(lon=lon, lat=lat, method='nearest').round(1)
    rr = da_wrf.interp(lon=lon, lat=lat, method='linear').round(1)
    xr.DataArray.interp
    return rr


# %%
def regrid_latlon(flnm_rain, area):
    """
    将combine得到的数据，
    插值到格距较粗的latlon格点上, 
    也包含有投影转换的需求
    插值到latlon格点上
    将二维的latlon坐标水平插值到一维的latlon坐标上
    """
    ds = xr.open_dataset(flnm_rain)
    ds_out = regrid_xesmf(ds, area, rd=1)
    return ds_out
    

    
def save_one(path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'):
    """处理一个模式的数据

    Args:
        path_main (str, optional): [description]. Defaults to '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'.
    """

    path_dic = {
        'path_main':path_main,  # 模式数据文件夹
        'path_rain_wrf_grid':path_main+'rain.nc', # 原始降水数据存储路径+文件名
        'path_rain_wrf_latlon':path_main+'rain_latlon.nc',  # 插值到latlon之后的文件名
        'path_rain_wrf_station':path_main+'rain_station.nc',  # 插值到站点之后的文件名
        'path_rain_obs_station':'/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc', # 站点降水
    }
    area = {
        'lon1':110.5,
        'lon2':116,
        'lat1':32,
        'lat2':36.5,
        'interval':0.125,
    }

    ## 合并数据
    da = combine_rain(path_main)
    da.to_netcdf(path_dic['path_rain_wrf_grid'])

    ## 降低分辨率
    # regrid_rain_model(path_dic, area)
    da1 = regrid_latlon(path_dic['path_rain_wrf_grid'], area)
    da1.to_netcdf(path_dic['path_rain_wrf_latlon'])

    ## 插值到站点
    da2 = grid2station(path_dic['path_rain_obs_station'], path_dic['path_rain_wrf_grid'],area)
    da2.to_netcdf(path_dic['path_rain_wrf_station'])
    pass

def dual():
    """处理多个模式的数据
    """
    pass
    model_list = ['1900_90m', '1912_900m', '1912_90m', 'YJF']
    for model in model_list:
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/'
        save_one(path_main)
    
if __name__ == '__main__':

    pass
    # main()
    dual()




