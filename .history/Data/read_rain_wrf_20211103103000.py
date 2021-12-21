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
    """由于wrfout数据的降水是累计降水，
    这里将它变为逐小时降水,
    同时进行合并"""

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
    dds_concate = xr.concat(dds_list, dim='time')
    dds_concate
    return dds_concate

def get_rain_model(path_dic):
    pass
    # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'
    # path_wrfout = path_main
    ds = get_rain(path_wrfout)
    print(ds.max())
    path_save = path_main+'rain_90m.nc'
    ds.to_netcdf(path_save)

def regrid(path_dic):
    """
    将combine得到的数据，插值到latlon格点上
    将二维的latlon坐标水平插值到一维的latlon坐标上
    """
    interval = 0.125
    area = {
        'lon1':111-1-interval/2,
        'lon2':115+1,
        'lat1':32-1-interval/2,
        'lat2':36+1,
        'interval':interval,
    }
    ds = xr.open_dataset(path_dic['path_rain'])
    ds_out = regrid_xesmf(ds, area)
    ds_out.to_netcdf(path_dic['path_rain_latlon'])
    
def main():
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'
    path_dic = {
        'path_main':path_main,  # 模式数据文件夹
        'path_rain':path_main+'rain.nc', # 原始降水存储路径+文件名
        'path_rain_latlon':path_main+'rain_latlon.nc'  # 插值之后的文件名
    }
    pass
    

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
    # get_high_hgt()
    # get_ERAI()
    regrid()


