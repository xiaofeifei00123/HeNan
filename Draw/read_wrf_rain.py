#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取wrfout数据中的降水
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
import netCDF4 as nc
# import wrf
# from multiprocessing import Pool

# %%
def get_rain(path):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800'
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'
    fl_list = os.popen('ls {}/wrfout*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()

    dds_list = []
    # fl_list = fl_list[0:2]

    r = 0
    for fl in fl_list:
        # fl = fl_list[0]
        ds = xr.open_dataset(fl)
        da = ds['RAINNC']-r
        dds_list.append(da)
        r = ds['RAINNC']
    # dds_concate = xr.concat(dds_list, dim='Time')

    # return dds_concate
    return dds_list


cc = get_rain('1')
# %%
# cc[1]



# %%
# fl = fl_list[0]
# ds = xr.open_dataset(fl)
# data_nc = nc.Dataset(fl)






# %%



# %%
if __name__ == '__main__':
    # main()
    # initial_file_list = ['ERA5', 'GDAS']
    # time_list = ['1800', '1812', '1900', '1912']
    time_list = ['1800']
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    # gu = GetUpar()
    for t in time_list:
        path_wrfout = path_main+'YSU_'+t+'/'
        # ds = gu.get_upar(path_wrfout)
        ds = get_rain_dual(path_wrfout)
        # flnm = 'YSU_'+t
        # path_save = path_main+flnm+'_upar.nc'
        # print(path_save)
        # ds.to_netcdf(path_save)
    # %%
    ds


