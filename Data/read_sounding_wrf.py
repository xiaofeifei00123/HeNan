#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
其实很简单，先插值到一个站点
然后再插值到垂直层
读取单个和多个站的wrf逐层数据
目的是为了多插值几层
参考：
https://wrf-python.readthedocs.io/en/latest/user_api/generated/wrf.interp1d.html#wrf.interp1d

## 对于读取的文件可以这样设置不同的纵坐标
dds = ds.set_coords(['height_agl', 'pressure'])
dds.swap_dims({'bottom_top':'pressure'})
-----------------------------------------
Time             :2021/10/05 22:53:08
Author           :fengxiang
Version          :1.1
'''

# %%
import xarray as xr
import os
import xesmf as xe
import numpy as np
import netCDF4 as nc
import wrf
from multiprocessing import Pool
# from multiprocessing import Pool
# from read_global import caculate_diagnostic, regrid_xesmf
# from baobao.caculate import caculate_q_rh_thetav
# from baobao.interp import regrid_xesmf
from baobao.caculate import caculate_q_rh_thetav
# from baobao.coord_transform import xy_ll
# from netCDF4 import MFDataset

# %%

def get_centroid(flnm= '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc', t_start='2021-07-20 00', t_end='2021-07-20 12'):
    """获得一段时间降水的质心
    根据降水聚合文件获取

    Args:
        flnm (str, optional): 原始wrfout降水数据的聚合. Defaults to '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'.
        t_start (str, optional): 起始时间. Defaults to '2021-07-20 00'.
        t_end (str, optional): 结束时间. Defaults to '2021-07-20 12'.

    Returns:
        [type]: [description]
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'
    da = xr.open_dataarray(flnm)
    # tt = slice('2021-07-20 00', '2021-07-20 12')
    tt = slice(t_start, t_end)
    rain = da.sel(time=tt).sum(dim='time')
    lat = sum(sum(rain*rain.lat))/sum(sum(rain))
    lon = sum(sum(rain*rain.lon))/sum(sum(rain))
    lon = lon.round(3)
    lat = lat.round(3)
    return lat, lon
# 我其实想要获得的是某一个点的多个变量, 但是变量只能一个一个获取

def sounding_1station_1time(flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d04_2021-07-19_00:00:00'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d04_2021-07-19_00:00:00'
    print(flnm[-19:])
    wrfnc = nc.Dataset(flnm)
    lat, lon = get_centroid()
    # x,y = wrf.ll_to_xy(wrfnc, 34.71, 113.66)
    x,y = wrf.ll_to_xy(wrfnc, lat, lon)

    hagl = wrf.getvar(wrfnc, 'height_agl', units='m')[:,x,y]  # 先纬度后经度
    # pj = str(hagl.attrs['projection'])
    pj = hagl.attrs['projection'].proj4()
    hagl = hagl.assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])
    p = wrf.getvar(wrfnc, 'pres', units='hpa')[:,x,y].assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])
    u = wrf.getvar(wrfnc, 'ua', units='m/s')[:,x,y].assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])

    v = wrf.getvar(wrfnc, 'va', units='m/s')[:,x,y].assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])
    t = wrf.getvar(wrfnc, 'temp', units='degC')[:,x,y].assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])
    td = wrf.getvar(wrfnc, 'td', units='degC')[:,x,y].assign_attrs({'projection':pj}).drop_vars(['latlon_coord'])
    ds = xr.merge([t,td, u, v, p, hagl])
    dds = ds.set_coords(['height_agl', 'pressure'])
    ds_return = dds.swap_dims({'bottom_top':'pressure'})
    return ds_return

def sounding_1station(fl_list):
    """单进程循环读取文件
    单个站点多个时次
    """
    pass
    dds_list = []
    for fl in fl_list:
        dds = sounding_1station_1time(fl)
    dds_list.append(dds)
    dds_concate = xr.concat(dds_list, dim='Time')
    dds_return = dds_concate.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time'}).drop_vars('XTIME')
    return dds_return

def sounding_1station_mp(fl_list):
    """多进程读取文件
    单个站点多个时次
    """
    pass
    pool = Pool(12)
    result = []
    for fl in fl_list:
        tr = pool.apply_async(sounding_1station_1time, args=(fl,))
        result.append(tr)
    pool.close()
    pool.join()

    dds_list = []
    for j in result:
        dds_list.append(j.get())

    dds_concate = xr.concat(dds_list, dim='Time')
    # ds_upar = dds_concate.rename({'level':'pressure', 'XLAT':'lat', 'XLONG':'lon', 'Time':'time'})
    dds_return = dds_concate.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time'}).drop_vars('XTIME')
    return dds_return


def get_upar(path):
    pass
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'
    fl_list = os.popen('ls {}/wrfout_d04*'.format(path))  # 打开一个管道
    fl_list = fl_list.read().split()
    ## 临时测试
    # fl_list = fl_list[0:2]
    dds = sounding_1station_mp(fl_list)
    print("开始计算诊断变量")
    # cc = caculate_diagnostic(dds)
    cc = caculate_q_rh_thetav(dds)
    print("合并保存数据")
    ds_upar = xr.merge([dds, cc])
    return ds_upar



def combine_one(model='1912_90m'):
    """
    将wrfout数据中需要的变量聚合成一个文件，并进行相关的垂直插值, 和诊断量的计算
    处理两种模式，不同时次的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    # gu = GetUpar()
    # path_wrfout = path_main+'1912_90m'
    path_wrfout = path_main+model
    ds = get_upar(path_wrfout)
    # ds = gu.get_upar_multi(path_wrfout)
    flnm = model+'/sounding.nc'
    path_save = path_main+flnm
    print(path_save)
    ds.to_netcdf(path_save)
    return ds

def combine():
    """
    将wrfout数据中需要的变量聚合成一个文件，并进行相关的垂直插值, 和诊断量的计算
    处理两种模式，不同时次的数据
    """
    model_list = ['1900_90m', '1900_900m','1912_90m', '1912_900m']
    for model in model_list:
        combine_one(model)


# %%


# %%
if __name__ == '__main__':
    ### combine和regrid一般不同时进行
    combine()
    # combine_one()
    # regrid_one()
    # combine() 
    # regrid()