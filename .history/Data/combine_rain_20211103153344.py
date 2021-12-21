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


# %%
# def get_data_mean():
def get_rain_obs():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_grid.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.drop_vars(['member', 'level', 'dtime'])
    lat=slice(32,37)
    lon=slice(110, 116)
    da = ds.sel(lon=lon, lat=lat)
    # da = dds['data0']
    return da


def get_rain_wrf(flnm):

    # model_list = ['1900_90m', '1912_900m', '1912_90m', 'YJF']
    # for model in model_list:
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/rain_latlon.nc'
    ds = xr.open_dataset(flnm)
    lat=slice(32,37)
    lon=slice(110, 116)
    dds = ds.sel(lon=lon, lat=lat)
    da = dds['RAINNC']
    return da


def main():
    type_list = ['ERA5', 'GDAS']
    time_list = ['1800', '1812', '1900', '1912']
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    ds = xr.Dataset()
    model_list = ['1900_90m', '1912_900m', '1912_90m', 'YJF']
    for model in model_list:
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/rain_latlon.nc'
        da = get_rain_wrf(flnm)
        ds[model] = da
        # print(type+t)
    ds['OBS'] = get_rain_obs()
    return ds
if __name__ == '__main__':
    
    dds = main()
    dds.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc')
