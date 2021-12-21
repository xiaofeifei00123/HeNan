#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取wrfout数据中的降水
将wrfout数据中的降水集合成一个文件
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
from read_global import caculate_diagnostic, regrid_xesmf

# %%
def get_rain(path):
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1900/'
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

        # dda = da.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
        # dc = dda.swap_dims({'Time':'time'})
        # da = da.swap_dims({'Time':'XTIME'})
        # da = da.rename({'XTIME':'time'})
        dds_list.append(dc)

    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate

def get_ERAI():
    pass
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/'
    # time_list = ['1800', '1812', '1900', '1912']
    # time_list = ['1900', '1912']
    # for t in time_list:
    # path_wrfout = path_main+'YSU_'+t+'_ERAI/'
    path_wrfout = path_main
    ds = get_rain(path_wrfout)
    print(ds.max())
    flnm = 'YSU'
    path_save = path_main+flnm+'_rain_1km.nc'
    ds.to_netcdf(path_save)


def regrid():
    """
    将combine得到的数据，插值到latlon格点上
    将二维的latlon坐标水平插值到一维的latlon坐标上
    """
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    interval = 0.125
    area = {
        'lon1':110-1-interval/2,
        'lon2':116+1,
        'lat1':32-1-interval/2,
        'lat2':35+1,
        'interval':interval,
    }
    # ## wrf的插值
    # for f in initial_file_list:
    #     path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
    #     for t in time_list:
    #         # path_wrfout = path_main+'YSU_'+t+'/'
    #         # ds = gu.get_upar(path_wrfout)
    #         flnm = 'YSU_'+t
    #         path_in = path_main+flnm+'_rain.nc'
    #         ds = xr.open_dataset(path_in)
    #         ds_out = regrid_xesmf(ds, area)
    #         path_out = path_main+flnm+'_rain_latlon.nc'
    #         # ds_out = ds_out.rename({'ua':'u', 'va':'v', 'geopt':'height'})
    #         ds_out.to_netcdf(path_out)
    ## GFS的插值
    path_high_in = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/YSU_rain_1km.nc'
    ds = xr.open_dataset(path_high_in)
    ds_out = regrid_xesmf(ds, area)
    path_GFS_out = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/YSU_rain_1km_latlon.nc'
    ds_out.to_netcdf(path_GFS_out)
    
# %%
### 测试开始
def test():
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912'
    ds = get_rain(path_main)
    ds.max()
### 测试结束
# %%
if __name__ == '__main__':

    pass
    # get_ERAI()
    regrid()


