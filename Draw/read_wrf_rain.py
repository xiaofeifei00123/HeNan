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

# %%
def get_rain(path):
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1900/'
    fl_list = os.popen('ls {}/wrfout*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    dds_list = []
    r = 0
    for fl in fl_list:
        ds = xr.open_dataset(fl)
        da = ds['RAINNC']-r
        r = ds['RAINNC'].values
        da = da.swap_dims({'Time':'XTIME'})
        da = da.rename({'XTIME':'time'})
        dds_list.append(da)

    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate



# %%
if __name__ == '__main__':
    # main()
    # initial_file_list = ['ERA5', 'GDAS']
    # time_list = ['1800', '1812', '1900', '1912']
    time_list = ['1900', '1912']
    # time_list = ['1800']
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    # gu = GetUpar()
    for t in time_list:
        path_wrfout = path_main+'YSU_'+t+'_ERAI/'
        # ds = gu.get_upar(path_wrfout)
        ds = get_rain(path_wrfout)
        flnm = 'YSU_'+t
        path_save = path_main+flnm+'_rain_ERAI.nc'
        # print(path_save)
        ds.to_netcdf(path_save)
    # %%
    # ds.isel(time=10).max()


