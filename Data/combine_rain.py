#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
将不同模式和观测的latlon降水聚合成一个文件
这样聚合成的文件有部分时次，部分模式是没有降水的
-----------------------------------------
Time             :2021/10/11 22:16:30
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import os

# %%
def get_rain_obs():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
    da = xr.open_dataarray(flnm)
    da = da.rename({'id':'sta'})
    return da

def get_rain_wrf(flnm):
    da = xr.open_dataarray(flnm)
    return da

def main():
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'
    # type_list = ['gwd3-BL', 'gwd3-FD', 'gwd3-LS', 'gwd3-SS']
    # type_list = ['gwd3', 'gwd1', 'gwd0', 'gwd3-test']
    type_list = ['gwd3-BL', 'gwd3-FD', 'gwd3-LS', 'gwd3-SS', 'gwd3', 'gwd1', 'gwd0', 'gwd3-test']

    ds = xr.Dataset()
    for type in type_list:
        fpath = os.path.join(path, type)
        flnm = os.path.join(fpath, 'rain_station.nc')
        # print(flnm)
        da = get_rain_wrf(flnm)
        ds[type] = da
    ds['OBS'] = get_rain_obs()
    print(ds)
    return ds
if __name__ == '__main__':
    
    dds = main()
    dds.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station1.nc')
