#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取wrfout的垂直剖面图的数据
-----------------------------------------
Time             :2021/11/09 20:51:01
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import numpy as np
import pandas as pd
from wrf import getvar, CoordPair, vertcross, get_cartopy
import wrf
from netCDF4 import Dataset
from multiprocessing import Pool
# import matplotlib.pyplot as plt
# from matplotlib.cm import get_cmap
# import cartopy.crs as crs
# from cartopy.feature import NaturalEarthFeature
from baobao.caculate import caculate_div3d


# %%
class CrossData():
    """获得垂直方向切向的数据
    提供剖面数据
    地形的剖面数据
    """
    def __init__(self, wrf_file) -> None:
        pass
        ## Create the start point and end point for the cross section
        ## 12时降水
        self.cross_start = CoordPair(lat=33.8, lon=111.7)

        self.cross_end = CoordPair(lat=33.5, lon=113.2)

        # 伏牛山
        # self.cross_start = CoordPair(lat=34, lon=110.5)
        # self.cross_end = CoordPair(lat=33.5, lon=113)


        # 嵩山
        # self.cross_start = CoordPair(lat=34.5, lon=112.5)
        # self.cross_end = CoordPair(lat=34, lon=115)



        # self.cross_start = CoordPair(lat=36, lon=111)
        # self.cross_end = CoordPair(lat=32.5, lon=113)

        # self.cross_start = CoordPair(lat=35.5, lon=113)
        # self.cross_end = CoordPair(lat=33.5, lon=113.5)
        ## read the ncfile
        self.ncfile = Dataset(wrf_file)
        ## 计算垂直坐标, 可以是离地高度、气压等
        # self.vert = getvar(self.ncfile, "height_agl")  # 离地高度坐标
        self.vert = getvar(self.ncfile, "z")  # 离地高度坐标
        # self.vert = getvar(self.ncfile, "pres")/100  # 气压坐标

    def get_vcross(self, var):
        """获得单个变量的切向数据, 竖着切

        Args:
            var ([type]): 变量名, 需要是wrf-python支持的

        Returns:
            [type]: [description]
        """

        # var =  getvar(self.ncfile, var)
        var_vcross = vertcross(var, self.vert, wrfin=self.ncfile,
                                     start_point=self.cross_start,
                                        end_point=self.cross_end, 
                                        latlon=True, )
        ## 改变投影的attrs的格式
        pj = var_vcross.attrs['projection'].proj4()
        var_vcross = var_vcross.assign_attrs({'projection':pj})


        ## 改变xy_loc的coords的存储格式
        coord_pairs = var_vcross.coords["xy_loc"].values
        x_labels = [pair.latlon_str(fmt="{:.3f}, {:.3f}")
                    for pair in coord_pairs]
        var_vcross = var_vcross.assign_coords({'xy_loc':('cross_line_idx',x_labels)})
        return var_vcross

    def get_ter(self,):
        """获得地形高度
        """
        ter = wrf.getvar(self.ncfile, "ter", timeidx=-1)
        ter_line = wrf.interpline(ter, wrfin=self.ncfile, 
                            start_point=self.cross_start,
                            end_point=self.cross_end)
        ter_line = ter_line.assign_attrs({'projection':'lambert'})
        return ter_line

    def get_cross_data(self, var_list=['ua', 'va', 'wa', 'theta_e']):
        """获得垂直切一刀的数据

        Returns:
            [type]: [description]
        """
        da_cross_list = []
        for var in var_list:
            var =  getvar(self.ncfile, var)  # 直接传变量, 不要传文件
            da = self.get_vcross(var)
            da_cross_list.append(da)

        ## 增加计算散度, 三维的散度计算
        u =  getvar(self.ncfile, 'ua')
        v =  getvar(self.ncfile, 'va')
        lon = u.XLONG
        lat = u.XLAT
        div = caculate_div3d(u,v, lon, lat) # 计算得到div
        div.attrs = u.attrs  # 因为get_vcross,对project做了统一处理，这里需要是一样的
        div_cross = self.get_vcross(div)  # 插值到剖面上
        da_cross_list.append(div_cross)
        
        ds = xr.merge(da_cross_list)
        return ds

    def get_proj(self):
        z = wrf.getvar(self.ncfile, "z")
        pj = get_cartopy(z)
        return pj

# %%
def save_one_model():
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'
    # tt = pd.date_range('2021-07-20 0000', '2021-07-20 1200', freq='3H')
    tt = pd.date_range('2021-07-20 0000', '2021-07-20 1200', freq='12H')
    # tt
    fl_list = []
    for t in tt:
        fl = 'wrfout_d03_'+t.strftime('%Y-%m-%d_%H:%M:%S')
        flnm = path+fl
        fl_list.append(flnm)

    ds_list = []
    for fl in fl_list:
        print(fl[-19:])
        cd = CrossData(fl)
        var_list = []
        ds = cd.get_cross_data()
        ds_list.append(ds)
    ds = xr.concat(ds_list, dim='Time')    

    ## 再把地形高度读出来, 利用了上面的fl和cd, 所以位置不能变
    ter = cd.get_ter()
    ds['ter'] = ter

    ds = ds.rename({'Time':'time'})
    save_name = path+'cross_rain.nc'
    ds.to_netcdf(save_name)


# %%
def __cross_1model_1time(flnm):
    """子函数，多进程中的每一个进程处理的过程
    """
    print(flnm[-19:])
    cd = CrossData(flnm)
    ds = cd.get_cross_data()
    return ds

def save_one_model_mp(path='/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'):


    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'
    # tt = pd.date_range('2021-07-20 0000', '2021-07-20 1200', freq='1H')
    tt = pd.date_range('2021-07-20 0000', '2021-07-21 0000', freq='1H')
    # tt
    fl_list = []
    for t in tt:
        fl = 'wrfout_d03_'+t.strftime('%Y-%m-%d_%H:%M:%S')
        flnm = path+fl
        fl_list.append(flnm)


    pool = Pool(13)
    result = []
    for fl in fl_list:
        tr = pool.apply_async(__cross_1model_1time, args=(fl,))
        result.append(tr)
    pool.close()
    pool.join()

    dds_list = []
    for j in result:
        dds_list.append(j.get())
    ds = xr.concat(dds_list, dim='Time')


    ## 将地形高度加到数据里面去
    fl = fl_list[0]
    cd = CrossData(fl)
    ter = cd.get_ter()
    ds['ter'] = ter

    ds = ds.rename({'Time':'time'})
    save_name = path+'cross_rain.nc'
    ds.to_netcdf(save_name)

def save_all_model():
    # model_list = ['1900_90m', '1900_900m','1912_90m', '1912_900m']
    # model_list = ['gwd3-NO','gwd3-CTL','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    model_list = ['gwd0', 'gwd3']
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    for model in model_list:
        path = path_main+model+'/'
        save_one_model_mp(path)
        # print(path)


if __name__ == '__main__':
    pass
    save_all_model()
