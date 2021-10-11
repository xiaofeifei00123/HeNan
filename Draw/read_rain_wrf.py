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

# %%
def get_rain(path):
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1900/'
    fl_list = os.popen('ls {}/wrfout_d03*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    dds_list = []
    r = 0
    for fl in fl_list:
        print(fl[-18:])
        ds = xr.open_dataset(fl)
        da = ds['RAINNC']-r
        r = ds['RAINNC'].values
        dda = da.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
        dc = dda.swap_dims({'Time':'time'})
        # da = da.swap_dims({'Time':'XTIME'})
        # da = da.rename({'XTIME':'time'})
        dds_list.append(dc)

    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate

def get_ERAI():
    pass
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    # time_list = ['1800', '1812', '1900', '1912']
    time_list = ['1900', '1912']
    for t in time_list:
        path_wrfout = path_main+'YSU_'+t+'_ERAI/'
        ds = get_rain(path_wrfout)
        print(ds.max())
        flnm = 'YSU_'+t
        path_save = path_main+flnm+'_rain_ERAI.nc'
        ds.to_netcdf(path_save)

def get_ERA5():
    pass
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    time_list = ['1800', '1812', '1900', '1912']
    # time_list = ['1900', '1912']
    for t in time_list:
        print(t)
        path_wrfout = path_main+'YSU_'+t
        ds = get_rain(path_wrfout)
        # print(ds.max())
        print("小时最大降水是%s"%str(ds.max().values))
        flnm = 'YSU_'+t
        path_save = path_main+flnm+'_rain.nc'
        ds.to_netcdf(path_save)

def get_GDAS():
    pass
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GDAS/'
    time_list = ['1800', '1812', '1900', '1912']
    # time_list = ['1900', '1912']
    for t in time_list:
        # print(t)
        path_wrfout = path_main+'YSU_'+t
        ds = get_rain(path_wrfout)
        print("小时最大降水是%s"%str(ds.max().values))
        flnm = 'YSU_'+t
        path_save = path_main+flnm+'_rain.nc'
        ds.to_netcdf(path_save)


def get_GFS():
    pass
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/'
    # time_list = ['1800', '1812', '1900', '1912']
    # time_list = ['1900', '1912']
    # for t in time_list:
        # print(t)
    path_wrfout = path_main+'YSU_GFS'
    ds = get_rain(path_wrfout)
    print("小时最大降水是%s"%str(ds.max().values))
    flnm = 'YSU_GFS'
    path_save = path_main+flnm+'_rain.nc'
    ds.to_netcdf(path_save)

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
    # get_ERA5()
    get_GFS()
    # get_GDAS()


