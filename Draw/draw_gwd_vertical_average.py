#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制重力波拖曳应力
区域平均的廓线
-----------------------------------------
Time             :2021/11/23 10:29:43
Author           :Forxd
Version          :1.0
'''
# %%
# from meteva.base.fun.statisticing import var_of_grd
import os
import xarray as xr
import pandas as pd
import numpy as np
import wrf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
# %%
# flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd0/wrfout_d01_2021-07-20_00:00:00'
# flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/wrfout_d01_2021-07-20_00:00:00'
# ds = xr.open_dataset(flnm)

# # ds
# da = ds['DTAUX3D']
# da

# %%
class Vp():
    """vertical profile of gravity dray
    """
    # def __init__(self):
    #     pass
    #     self.var_list = ['DTAUY3D_LS', 'DTAUY3D_SS', 'DTAUY3D_FD', 'DTAUY3D_BL']
    #     # self.var_list = ['DTAUX3D_LS', 'DTAUX3D_SS', 'DTAUX3D_FD', 'DTAUX3D_BL']
    #     self.label_list = ['LS-GWD', 'SS-GWD', 'TOFD', 'BL']
    #     self.color_list = ['green','orange',  'red', 'blue']

    def __init__(self):
        """gwd1
        """
        pass
        self.var_list = ['DTAUY3D']
        # self.var_list = ['DTAUX3D_LS', 'DTAUX3D_SS', 'DTAUX3D_FD', 'DTAUX3D_BL']
        self.label_list = ['gwd1']
        self.color_list = ['green','orange',  'red', 'blue']

class Data(Vp):
    def __init__(self, flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3/wrfout_d01_2021-07-20_00:00:00'):
        super().__init__()  # 几层父类的属性，没有这个好像继承不起来
        self.flnm  = flnm
        self.wrfin = Dataset(flnm)

    def get_da_vertical(self, da):
        """求区域平均值

        Args:
            da ([type]): 某一个变量的所有空间网格三维数据

        Returns:
            [type]: [description]
        """

        # area = {
        #     'lon1':111.5,
        #     'lon2':113,
        #     'lat1':33,
        #     'lat2':34,
        #     'interval':0.125,
        # }
        area = {
            'lon1':113,
            'lon2':114,
            'lat1':34.5,
            'lat2':35.5,
            'interval':0.125,
        }
        # wrfin = Dataset(self.flnm)
        wrfin = self.wrfin
        x1, y1 = wrf.ll_to_xy(wrfin,area['lat1'], area['lon1']).values
        x2, y2 = wrf.ll_to_xy(wrfin,area['lat2'], area['lon2']).values
        dda = da[:,y1:y2, x1:x2]
        ddda = dda.mean(dim={'south_north', 'west_east'})
        return ddda

    def get_da(self,):
        ds = xr.open_dataset(self.flnm)
        dds = xr.Dataset()
        for var in self.var_list:
            da = ds[var].squeeze()
            dds[var] = self.get_da_vertical(da)
        return dds

    def get_da_dual(self,):
        ds = xr.open_dataset(self.flnm)
        dds = xr.Dataset()
        # var_list1 = ['DTAUY3D_LS', 'DTAUY3D_SS', 'DTAUY3D_FD', 'DTAUY3D_BL']
        # var_list2 = ['DTAUX3D_LS', 'DTAUX3D_SS', 'DTAUX3D_FD', 'DTAUX3D_BL']
        var_list1 = ['DTAUY3D']
        var_list2 = ['DTAUX3D']
        for var1,var2,label in zip(var_list1,var_list2, self.label_list):
            da1 = ds[var1].squeeze()
            # da1 = self.get_da_vertical(da1)
            da2 = ds[var2].squeeze()
            # da2 = self.get_da_vertical(da2)
            da = np.sqrt(da1**2+da2**2)
            da = self.get_da_vertical(da)  # 可以先求合力，再插值, 比较大小是可以的
            dds[label] = da
        return dds

def get_drag_gwd0():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/wrfout_d03_2021-07-20_00:00:00'
    ds = xr.open_dataset(flnm)
    dax = ds['DTAUX3D']
    day = ds['DTAUY3D']
    da = np.sqrt(dax**2+day**2)
    da = da.squeeze()
    return da

# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_00:00:00'
# gd = Data(flnm)
# dds = gd.get_da_dual()
# # dds

# da1 = get_drag_gwd0()
# da1 = gd.get_da_vertical(da1)
# da1
# # %%
# dds['gwd0_all'] = da1
# dds

# %%
class Draw(Vp):
    def __init__(self,pic_dic) -> None:
        super().__init__()
        self.pic_dic = pic_dic
        self.pic_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/av_upar/'
        self.fig_name = os.path.join(self.pic_path, pic_dic['time']+pic_dic['domain']+'gwd1')

    def draw(self,ds):
        plt.rcParams['font.family'] = 'Nimbus Roman No9 L' # 作用：解决图上汉字和负号显示为方框的问题 
        # plt.rcParams['font.sans-serif'] = 'Nimbus Roman No9 L'
        fig = plt.figure(figsize=[10,10])
        ax = fig.add_axes([0.1,0.1, 0.8,0.8])
        # color_list = ['green','orange',  'red', 'blue']
        for var,color,label in zip(self.var_list,self.color_list, self.label_list):
            dav = ds[var]
            ax.plot(dav.values, dav.bottom_top, label=label, lw=5, color=color)

        ax.set_xlim(10**(-10),10**(-2))
        ax.set_ylim(0,50)
        ax.set_xscale('log')

        ax.xaxis.set_tick_params(labelsize=22)
        ax.yaxis.set_tick_params(labelsize=22)
        ax.legend(fontsize=24, edgecolor='white')
        ax.set_xlabel('Drag ($m/s^2$)', fontsize=26)
        ax.set_ylabel('Model Level', fontsize=26)
        fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_gwd/ogwd_average_distribution.png')

    def draw_dual(self,ds):
        """水平反向合力的图

        Args:
            ds ([type]): [description]
        """
        plt.rcParams['font.family'] = 'Nimbus Roman No9 L' # 作用：解决图上汉字和负号显示为方框的问题 
        # plt.rcParams['font.sans-serif'] = 'Nimbus Roman No9 L'
        fig = plt.figure(figsize=[10,10])
        ax = fig.add_axes([0.1,0.1, 0.8,0.8])
        # color_list = ['green','orange',  'red', 'blue']
        # for var,color,label in zip(self.var_list,self.color_list, self.label_list):
        for label,color in zip(self.label_list,self.color_list):
            dav = ds[label]
            ax.plot(dav.values, dav.bottom_top, label=label, lw=5, color=color)
        ax.set_xlim(10**(-10),10**(-2))
        ax.set_ylim(0,50)
        ax.set_xscale('log')

        ax.xaxis.set_tick_params(labelsize=22)
        ax.yaxis.set_tick_params(labelsize=22)
        ax.legend(fontsize=24, edgecolor='white')
        ax.set_xlabel('Drag ($m/s^2$)', fontsize=26)
        ax.set_ylabel('Model Level', fontsize=26)
        ax.set_title(self.pic_dic['domain'], loc='left', fontsize=22)
        ax.set_title(self.pic_dic['time'], loc='right', fontsize=22)
        
        # fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/aa.png')
        fig.savefig(self.fig_name+'right')
        plt.close()

def main():
    # f_path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/'
    f_path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/'
    domain_list = ['d01', 'd02', 'd03']
    time_range = pd.date_range('2021-07-20 1600', '2021-07-20 1600', freq='3H')
    for t in time_range:
        for d in domain_list:
            f_name = 'wrfout_'+d+'_'+t.strftime('%Y-%m-%d_%H:%M:%S')
            flnm = f_path+f_name

            pic_dic = {'time':t.strftime('%d-%H'), 'domain':d}            
            print(pic_dic['domain']+pic_dic['time']+'gwd1')
            

    
            dt = Data(flnm)
            ds = dt.get_da_dual()
            dr = Draw(pic_dic)
            dr.draw_dual(ds)
    
    

# def draw(flnm):
#     dt = Data(flnm)
#     ds = dt.get_da_dual()
#     dr = Draw()
#     dr.draw_dual(ds)


# %%
if __name__ == '__main__':
    main()
    # pass
    