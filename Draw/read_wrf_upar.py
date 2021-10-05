#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读高空数据，插值到气压坐标上
u, v, q, temp, height_agl
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
import wrf
from multiprocessing import Pool
from read_global import caculate_diagnostic




# %%
class GetUpar():
    """获得wrfout高空数据，原始投影
    """
    def get_upar_one(self, fl):
        # dds_list = []
        pre_level = [1000, 850, 700, 500, 250]
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
        pool = Pool(1)
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
        return dds_return

    def get_upar(self, path):
        pass
        # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'
        fl_list = os.popen('ls {}/wrfout_d02*'.format(path))  # 打开一个管道
        fl_list = fl_list.read().split()
        # fl_list = fl_list[0:2]
        dds = self.get_upar_multi(fl_list)
        cc = caculate_diagnostic(dds)
        ds_upar = xr.merge([dds, cc])
        return ds_upar


# %%
if __name__ == '__main__':
    # main()
    # initial_file_list = ['ERA5', 'GDAS']
    time_list = ['1800', '1812', '1900', '1912']
    # time_list = ['1800']
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
    gu = GetUpar()
    for t in time_list:
        path_wrfout = path_main+'YSU_'+t+'/'
        ds = gu.get_upar(path_wrfout)
        flnm = 'YSU_'+t
        path_save = path_main+flnm+'_upar_d02.nc'
        print(path_save)
        ds.to_netcdf(path_save)
