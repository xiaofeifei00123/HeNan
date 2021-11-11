#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
将不同模式和观测的latlon降水聚合成一个文件
station降水聚合成一个文件
这样聚合成的文件有部分时次，部分模式是没有降水的
-----------------------------------------
Time             :2021/10/11 22:16:30
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr




# %%
def combine_latlon():
    ds = xr.Dataset()
    model_list = ['1900_90m','1900_900m', '1912_900m', '1912_90m', ]
    for model in model_list:
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/rain_latlon.nc'
        da = xr.open_dataarray(flnm)
        ds[model] = da

    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
    da_obs = xr.open_dataarray(flnm_obs)
    ds['OBS'] = da_obs
    ds.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_latlon.nc')
    # return ds

def combine_station():
    ds = xr.Dataset()
    model_list = ['1900_90m', '1900_900m','1912_900m', '1912_90m']
    for model in model_list:
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/rain_station.nc'
        da = xr.open_dataarray(flnm)
        ds[model] = da
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
    da_obs = xr.open_dataarray(flnm_obs)
    ds['OBS'] = da_obs.rename({'id':'sta'})
    # ds['OBS'] = xr.open_dataarray(flnm_obs)
    ds.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc')
    # return ds


if __name__ == '__main__':
    
    combine_latlon()
    combine_station()
