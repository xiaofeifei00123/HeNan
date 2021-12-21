#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取聚合wrf高空数据
保存
原始网格数据
插值成latlon网格数据
combine:
    从原始的wrfout数据中，
    读取需要的变量
    插值到需要的气压坐标上
    计算需要的诊断变量，例如q
    聚合成一个文件
regrid:
    将combine的数据，水平插值到需要的等经纬网格点上
    或者说是转换为latlon坐标
此程序读取和最终生成的wrfout变量有
u, v, q, temp, height_agl, geopt
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
from read_global import caculate_diagnostic, regrid_xesmf


# %%
class GetUpar():
    """获得wrfout高空数据，原始投影
    """
    def get_upar_one(self, fl):
        # dds_list = []
        pre_level = [1000, 850, 700, 500, 200]
        # pre_level = np.arange(1000, 90, -20)  # 对d03使用垂直方向10度的分辨率
        # pre_level = np.arange(1000, 90, -50)  # 对d02使用垂直方向100度的分辨率
        dds = xr.Dataset()
        data_nc = nc.Dataset(fl)
        print(fl[-19:])
        p = wrf.getvar(data_nc, 'pressure', squeeze=False)

        for var in ['ua', 'va', 'td', 'temp', 'theta_e', 'height_agl', 'geopt']:
            da = wrf.getvar(data_nc, var, squeeze=False)
            # dds[var] = da.expand_dims(dim='Time')
            dds[var] = wrf.interplevel(da, p, pre_level, squeeze=False)
            attr_proj = str(dds[var].projection)
            dds[var]=dds[var].assign_attrs({'projection':attr_proj})
        # dds['height_agl']=dds['height_agl'].assign_attrs({'projection':'lambert'})
        return dds

    def get_upar_dual(self, fl_list):
        """单进程循环读取文件
        """
        pass
        dds_list = []
        for fl in fl_list:
            # dds = xr.Dataset()
            dds = self.get_upar_one(fl)
        dds_list.append(dds)
        dds_concate = xr.concat(dds_list, dim='Time')
        dds_return = dds_concate.rename({'level':'pressure', 'XLAT':'lat', 'XLONG':'lon', 'Time':'time'})
        return dds_return

    def get_upar_multi(self, fl_list):
        """多进程读取文件
        """
        pass
        pool = Pool(2)
        result = []
        for fl in fl_list:
            tr = pool.apply_async(self.get_upar_one, args=(fl,))
            result.append(tr)
        pool.close()
        pool.join()

        dds_list = []
        for j in result:
            dds_list.append(j.get())

        dds_concate = xr.concat(dds_list, dim='Time')
        dds_return = dds_concate.rename({'level':'pressure', 'XLAT':'lat', 'XLONG':'lon', 'Time':'time'})
        dds_return = dds_return.drop_vars(['XTIME'])
        return dds_return

    def get_upar(self, path):
        pass
        # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'
        fl_list = os.popen('ls {}/wrfout_d02*'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        # fl_list = fl_list[0:2]
        dds = self.get_upar_multi(fl_list)
        print("开始计算诊断变量")
        cc = caculate_diagnostic(dds)
        print("合并保存数据")
        ds_upar = xr.merge([dds, cc])
        return ds_upar

def combine():
    """
    将wrfout数据中需要的变量聚合成一个文件，并进行相关的垂直插值, 和诊断量的计算
    处理两种模式，不同时次的数据
    """
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    for f in initial_file_list:
        # time_list = ['1800']
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
        gu = GetUpar()
        for t in time_list:
            path_wrfout = path_main+'YSU_'+t+'/'
            ds = gu.get_upar(path_wrfout)
            flnm = 'YSU_'+t
            path_save = path_main+flnm+'_upar_d02.nc'
            print(path_save)
            ds.to_netcdf(path_save)

def combine_one():
    """
    将wrfout数据中需要的变量聚合成一个文件，并进行相关的垂直插值, 和诊断量的计算
    处理两种模式，不同时次的数据
    """
    # time_list = ['1800', '1812', '1900', '1912']
    # initial_file_list = ['ERA5', 'GDAS']
    # for f in initial_file_list:
        # time_list = ['1800']
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    gu = GetUpar()
    # for t in time_list:
    path_wrfout = path_main+'High'
    ds = gu.get_upar(path_wrfout)
    flnm = 'YSU_'+t
    path_save = path_main+flnm+'_upar_d02.nc'
    print(path_save)
    ds.to_netcdf(path_save)

        

def regrid():
    """
    将combine得到的数据，插值到latlon格点上
    将二维的latlon坐标水平插值到一维的latlon坐标上
    """
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    for f in initial_file_list:
        # time_list = ['1800']
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
        # gu = GetUpar()
        area = {
            'lon1':107-1,
            'lon2':135+1,
            'lat1':20-1,
            'lat2':40+1,
            'interval':0.5,
        }
        for t in time_list:
            # path_wrfout = path_main+'YSU_'+t+'/'
            # ds = gu.get_upar(path_wrfout)
            flnm = 'YSU_'+t
            path_in = path_main+flnm+'_upar_d02.nc'
            ds = xr.open_dataset(path_in)
            ds_out = regrid_xesmf(ds, area)
            path_out = path_main+flnm+'_upar_d02_latlon.nc'
            ds_out = ds_out.rename({'ua':'u', 'va':'v', 'geopt':'height'})
            ds_out.to_netcdf(path_out)

            


# %%
if __name__ == '__main__':
    ### combine和regrid一般不同时进行
    combine() 
    # regrid()