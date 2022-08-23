#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
时间高度剖面数据, wrfout数据
单个站点，或者区域平均的
气压层，高度层

参考：
https://wrf-python.readthedocs.io/en/latest/user_api/generated/wrf.interp1d.html#wrf.interp1d

## 对于读取的文件可以这样设置不同的纵坐标
dds = ds.set_coords(['height', 'pressure'])
dds.swap_dims({'bottom_top':'pressure'})
# TODO 气压或者高度坐标会不会随着时间变化
# TODO 只负责处理数据 
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
import pandas as pd
import netCDF4 as nc
import wrf
from multiprocessing import Pool
from baobao.caculate import caculate_q_rh_thetav, caculate_average_wrf, get_div_wrfout, get_vor_wrfout

# %%


class Sounding():
    """单个站点, 单个模式的
    诊断量计算
    数据聚合
    """
    def __init__(self,
                sta={'sta_num':'57083','sta_name':'zhengzhou','lon':113.66,'lat':34.71 },
                area = {
                    'lat1':34.2,
                    'lat2':35.2,
                    'lon1':113,
                    'lon2':114,
                },
                flag = 'area',
                 ):
        """ 存放一些公共使用的变量, 也就是每个函数都一样的哪些 将一些不同模式不一样的东西，全部定义到这里
        """
        self.flag = flag
        self.area = area
        self.sta = sta
        self.var_list = ['p', 'ua', 'va', 'wa', 'temp', 'td', 'theta', 'theta_e', 'height']
        pass

    def sounding_1station_1time(self, flnm, *args, **kwargs):
        """_summary_

        Args:
            flnm (_type_): 一个wrfout数据路径

        Returns:
            _type_: _description_
        """
        print(flnm[-19:])
        wrfnc = nc.Dataset(flnm)

        if self.flag == 'sta':
            lat, lon = self.sta['lon'], self.sta['lat']
            x,y = wrf.ll_to_xy(wrfnc, lat, lon)

        elif self.flag == 'area':
            area = self.area
            hagl = wrf.getvar(wrfnc, 'height', units='m')  # 先纬度后经度
            pj = hagl.attrs['projection'].proj4()
            va_list = []
            for var in self.var_list:
                if var in ['pres', 'pressure', 'p']:
                    va = wrf.getvar(wrfnc, var, units='hpa').assign_attrs({'projection':pj})
                elif var in ['ua', 'va', 'wa']:
                    va = wrf.getvar(wrfnc, var, units='m/s').assign_attrs({'projection':pj})
                elif var in ['temp', 'td', 'theta', 'theta_e']:
                    va = wrf.getvar(wrfnc, var, units='degC').assign_attrs({'projection':pj})
                elif var in ['height', 'z']:
                    va = wrf.getvar(wrfnc, var, units='m')  # 先纬度后经度
                else:
                    va = wrf.getvar(wrfnc, var).assign_attrs({'projection':pj})
                va2 = caculate_average_wrf(va, area=self.area)
                va_list.append(va2)
            div = get_div_wrfout(flnm)
            div1 = caculate_average_wrf(div, area=self.area)
            vor = get_vor_wrfout(flnm)
            vor1 = caculate_average_wrf(vor, area=self.area)
            va_list.append(div1)
            va_list.append(vor1)

        ds = xr.merge(va_list)
        dds = ds.set_coords(['height', 'pressure'])
        ds_return = dds.swap_dims({'bottom_top':'pressure'})
        return ds_return
        # return dds
        
    def sounding_1station(self, fl_list):
        """单进程循环读取文件
        单个站点多个时次
        """
        pass
        dds_list = []
        for fl in fl_list:
            dds = self.sounding_1station_1time(fl,)
        dds_list.append(dds)
        dds_concate = xr.concat(dds_list, dim='Time')
        dds_return = dds_concate.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time'}).drop_vars('XTIME')
        return dds_return

    def sounding_1station_mp(self, fl_list):
        """多进程读取文件
        单个站点多个时次
        """
        pass
        pool = Pool(13)
        result = []
        for fl in fl_list:
            tr = pool.apply_async(self.sounding_1station_1time, args=(fl,))
            result.append(tr)
        pool.close()
        pool.join()

        dds_list = []
        for j in result:
            dds_list.append(j.get())
        # print(dds_list)
        dds_concate = xr.concat(dds_list, dim='Time')
        # ds_upar = dds_concate.rename({'level':'pressure', 'XLAT':'lat', 'XLONG':'lon', 'Time':'time'})
        # dds_return = dds_concate.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time'}).drop_vars('XTIME')
        dds_return = dds_concate.rename({'Time':'time'}).drop_vars('XTIME')
        return dds_return

    
    def sounding_main(self, path_wrfout):
        """处理流程的主控制函数

        Args:
            path ([type]): [description]
            path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'

        Returns:
            [type]: [description]
        """
        # path_wrfout = os.path.join(self.path_main, self.model)
        fl_list = os.popen('ls {}/wrfout_d03*'.format(path_wrfout))  # wrfout文件的path
        fl_list = fl_list.read().split()
        # print(fl_list)

        print("1.合并不同时次数据")
        dds = self.sounding_1station_mp(fl_list)
        print("2. 开始计算诊断变量")
        cc = caculate_q_rh_thetav(dds)
        print("3. 合并诊断变量")
        ds_upar = xr.merge([dds, cc])
        return ds_upar


# sd = Sounding()
# sd.sounding_main()


# %%

def sounding_dual():
    """
    将wrfout数据中需要的变量聚合成一个文件，并进行相关的垂直插值, 和诊断量的计算
    处理两种模式，不同时次的数据
    多模式数据的合并
    """
    sd = Sounding()
    ds = sd.sounding_main()
        # ds_list.append(ds)

# %%
def combine():
    """将不同模式，不同站点的数据聚合到一起
    """

    model_list = ['gwd3']
    sta_dic_list = [
        {'sta_num':'57083','sta_name':'zhengzhou','lon':113.66,'lat':34.71 },
        {'sta_num':'57178','sta_name':'nanyang','lon':112.4,'lat':33.1 }, 
        {'sta_num':'57067','sta_name':'lushi','lon':111.04,'lat':34.05 }, 
    ]
    sta_list = ['zhengzhou', 'nanyang', 'lushi']
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    ds_sta_list = []
    for stadic in sta_dic_list:
        ds_model_list = []
        for model in model_list:
            flnm = 'sounding_'+stadic['sta_name']+'_'+model+'.nc'
            # print(path_main)
            path_save1 = os.path.join(path_main, model)
            path_save = os.path.join(path_save1,flnm)
            print(path_save)
            ds = xr.open_dataset(path_save)
            ds_model_list.append(ds)
        dds = xr.concat(ds_model_list, dim=pd.Index(model_list, name='model'))
        ds_sta_list.append(dds)
    ddds  = xr.concat(ds_sta_list, dim=pd.Index(sta_list, name='station'))
    flnm_all = os.path.join(path_main,'sounding_all.nc')
    ddds.to_netcdf(flnm_all)
    return ddds
    # dds


if __name__ == '__main__':
    pass
    # sounding_dual() # 分别存储
    sd = Sounding()
    ds = sd.sounding_main()
    # aa = combine()  # 合并为一个文件
    # print(aa)
    # sd = Sounding()
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-19_17:00:00'
    # sd.sounding_1station_1time(flnm, flag='area')
    

# %%
# import pyproj
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d01_2021-07-20_04:00:00'
# wrfnc = nc.Dataset(flnm)
# h = wrf.getvar(wrfnc, 'height', units='m')  # 先纬度后经度
# # h.to_netcdf('aa.nc', engine='h5netcdf')
# # h.attrs['projection']
# pj = h.attrs['projection'].proj4()
# # pj
# h.attrs['projection'] = pj
# h
# # %%
# h.to_netcdf('aa.nc')
# type(h.attrs['projection'])
# type(pj)
# h.encoding
# ds = xr.open_dataset(flnm, engine='netcdf4')
# ds.attrs
# %%
# pj