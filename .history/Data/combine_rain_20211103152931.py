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
    dds = ds.sel(lon=lon, lat=lat)
    da = dds['data0']
    return da


def get_rain_wrf(flnm):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_rain_latlon.nc'

    model_list = ['1900_90m', '1912_900m', '1912_90m', 'YJF']
    for model in model_list:
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/'
    

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
    for type in type_list:
        for t in time_list:
            flnm = path+type+'/'+'YSU_'+t+'_rain_latlon.nc'
            # print(flnm)
            da = get_rain_wrf(flnm)
            ds[type+t] = da
            # print(type+t)
    ds['OBS'] = get_rain_obs()
    # ds['EC'] = get_rain_ec()
    # ds['GFS'] = get_rain_GFS()
    return ds
if __name__ == '__main__':
    
    dds = main()
    dds.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc')
