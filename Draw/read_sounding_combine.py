#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
目的：
    合并郑州站观测和模式数据
    观测:micaps资料TlogP聚合成的文件
    模式:不同边界和起报时间
    
参考: 
    https://zhuanlan.zhihu.com/p/102032513
需要注意:
    这个是画单个探空站用的，所以画的要素比较全,

-----------------------------------------
Time             :2021/03/22 14:03:00
Author           :Forxd
Version          :1.0
'''

# %%
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # 这个还没用过
import numpy as np
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import xarray as xr
from read_global import caculate_diagnostic
import datetime

# %%
def get_wrf(flnm):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d03_latlon.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    ds = ds.drop_vars(['theta_e', 'XTIME'])
    ds = ds.rename({'ua':'u', 'va':'v'})
    t1 = ds.time.sel(time=datetime.time(int(0)))
    t2 = ds.time.sel(time=datetime.time(int(6)))
    t3 = ds.time.sel(time=datetime.time(int(12)))
    t4 = xr.concat([t1, t2, t3], dim='time')
    ds_model = ds.sel(time=t4, lat=34.71, lon=113.66, method='nearest')
    # ds_model
    return ds_model

    
def get_wrf_all():
    pass
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    ds2 = xr.Dataset()
    wrf_list = []
    name_list = []
    for f in initial_file_list:
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
        for t in time_list:
            dic_model = {'initial_time':t, 'file_type':f}
            model = f+t
            print(model)
            name_list.append(model)
            flnm = 'YSU_'+t
            path_in = path_main+flnm+'_upar_d03_latlon.nc'
            ds = get_wrf(path_in)
            wrf_list.append(ds)
    
    ds_wrf = xr.concat(wrf_list, pd.Index(name_list, name='model'))
    return ds_wrf

def get_micaps():
    """读取micaps数据
    计算诊断变量

    Returns:
        [type]: [description]
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    ds = ds.drop_vars('variable')
    
    ### 这里的时间处理不具有通用性
    t1 = pd.date_range('2021-07-18 00', '2021-07-21 12', freq='6H')
    t2 = pd.date_range('2021-07-18 18', '2021-07-21 12', freq='24H')
    t3 = t1.drop(t2)
    ds = ds.sel(time=t3)
    
    ds_diagnostic = caculate_diagnostic(ds)
    ds_obs = xr.merge([ds, ds_diagnostic])
    ## 增加一个model维度
    cc = ds_obs.expand_dims(dim='model')
    ds_return = cc.assign_coords({'model':['micaps']})
    return ds_return
    
if __name__ == '__main__':
    ds_wrf = get_wrf_all()
    ds_obs = get_micaps()    
    ds_all = xr.concat([ds_wrf, ds_obs], dim='model')
    ds_all.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou_all.nc')

