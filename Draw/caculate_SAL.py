#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算并绘制24h总降水的SAL评分
-----------------------------------------
Time             :2022/01/06 
Author           :Forxd
Version          :1.0
'''

# %%
# from numpy.lib.shape_base import column_stack
# from pandas.core.frame import DataFrame
import xarray as xr
import meteva.method as mem
import meteva.base as meb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from baobao.select import select_area_latlon
# from read_data import TransferData
# from read_data import GetData

# from global_variable import station_dic

# %%

class SAL():

    def __init__(self, ):
        # self.rain = rain
        pass

    def caculate_threshold(self, rain):
        """"计算降水阈值"""
        rain_dim1 = rain.values.flatten()  # 二维降水数据变为1维
        # tt = rain_dim1[~np.isin(rain_dim1, 0.)]  # 删除0降水的点
        # rain_dim1 = np.sort(tt)  # 排序
        rain_dim1 = np.sort(rain_dim1)  # 这里有没有计算0降水量影响在2mm内
        max_num = len(rain_dim1) - 1  # 取最大降水序号
        # 最大降水序号*0.95为降水阈值
        rain_threshold = rain_dim1[int(max_num * 0.95)]
        # rain_dim1
        # rain_threshold/15
        return rain_threshold

    def caculate_space_scale(self, rain_obs, rain_model, rain_threshold):
        """计算一个时次降水的SAL评分
        """
        grid_obs = meb.xarray_to_griddata(rain_obs)
        grid_model = meb.xarray_to_griddata(rain_model)

        look_ff = mem.mode.feature_finder(grid_obs, grid_model, smooth=0, threshold=rain_threshold, minsize=5)
        look_match = mem.mode.centmatch(look_ff)
        look_merge = mem.mode.merge_force(look_match)
        sal = mem.sal(look_ff)
        ssal = pd.Series(sal)
        # return ssal, look_merge, look_ff
        return ssal
        

# %%
def get_rain_obs():
    pass



def get_sal_one_model(flnm_model):
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
    da = xr.open_dataarray(flnm_obs)
    # area = {
    #     'lon1':111.5,
    #     'lon2':113,
    #     'lat1':33,
    #     'lat2':34,
    #     'interval':0.05,
    # }
    area = {
        'lon1':112.5,
        'lon2':114.5,
        'lat1':34,
        'lat2':35,
        'interval':0.05,
    }

    rain_obs = da.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
    rain_obs = select_area_latlon(rain_obs, area)
    rain_obs['time'] = pd.Timestamp('2021-07-20 00')
    sal = SAL()
    # rain_threshold = ca.caculate_threshold(rain_obs)
    rain_threshold = 100
    # flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
    da = xr.open_dataarray(flnm_model)
    

    rain_model = da.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
    rain_model = select_area_latlon(rain_model, area) # 筛选指定范围内的数据
    # print("1111"*10)
    # print(rain_obs)
    # print("2222"*10)
    # print(rain_model)
    # print("3333"*10)
    rain_model['time'] = pd.Timestamp('2021-07-20 00')
    sal = sal.caculate_space_scale(rain_obs, rain_model, rain_threshold)
    sal = xr.DataArray.from_series(sal)
    # sal
    return sal


def get_sal_dual_model():
    pass
    model_list = ['gwd0', 'gwd1', 'gwd3']
    # flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
    flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'

    sal_list = []
    for model in model_list:
        pass
        fname = model+'/rain_latlon.nc'
        flnm = os.path.join(flpath, fname)
        # print(flnm)
        sal = get_sal_one_model(flnm)
        sal_list.append(sal)
    sal_all = xr.concat(sal_list, dim=pd.Index(model_list, name='model'))
    return sal_all
    

        # sal = get_sal_one_model(model)




class Draw():
    def __init__(self):
        pass
        self.fontsize=12

    def set_ticks(self, ax, ):    
        pass
        ax.xaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.yaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        # ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
        # ax.set_xlabel("Date/Hour (UTC)", fontsize=self.fontsize*2.5)
        # ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.5)
        # ax.legend(fontsize=self.fontsize*2.0, edgecolor='white')



    def draw_SAL(self,df, ax, title='SAL'):

        labels = ["Structure", "Amplitude", "Location"]
        x = np.arange(len(labels)) * 2
        width = 0.3

        colors = ['red', 'green', 'blue', 'orange', 'cyan']
        # 画柱状图

        rects1 = ax.bar(x - width * 2, df['gwd0'][0:3], width, label='gwd0', color=colors[0])
        rects2 = ax.bar(x - width, df['gwd1'][0:3], width, label='gwd1', color=colors[1])
        rects3 = ax.bar(x, df['gwd3'][0:3], width, label='gwd3', color=colors[2])
        # rects4 = ax.bar(x + width, df['TEMF'][0:3], width, label='TEMF', color=colors[3])
        # rects5 = ax.bar(x + width*2, df['MYJ'][0:3], width, label='MYJ', color=colors[4])

        # ax.set_ylabel('SAL')
        ax.set_xticks(x)
        # ax.xaxis.set_tick_params(labelsize=22)
        ax.set_yticks(np.arange(-2, 2.1, 0.5))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.set_ylim(-0.4, 0.2)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_ylim(-1, 2)  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.yaxis.set_tick_params(labelsize=20)
        ax.set_xticklabels(labels)
        ax.set_title(title, fontsize=26)
        ax.axhline(y=0, color='black') # 画0线
        ax.legend(fontsize=24, loc='upper left', edgecolor='white')
        self.set_ticks(ax)


aa = get_sal_dual_model()
aa.rename('sal')
bb = aa.to_dataset(dim='model')
bb
dr = Draw()
fig = plt.figure(figsize=(10,8), dpi=600)
ax = fig.add_axes([0.1,0.1,0.8,0.8])
dr.draw_SAL(bb,ax)
fig_name = 'sal_area2.png'
fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/'
fig_save = os.path.join(fig_path, fig_name)
fig.savefig(fig_save)




# sal
# xr.DataArray.from_series(sal)

# %%
# area = {
#     'lon1':112,
#     'lon2':113,
#     'lat1':33,
#     'lat2':34,
#     'interval':0.05,
# }
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
# da = xr.open_dataarray(flnm_obs)
# da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
# # da
# da1 = xr.where((da.lat.round(2)<=area['lat2'])&(da.lat.round(2)>=area['lat1']) , da, np.nan)
# da2 = da1.dropna(dim='lat')
# da3 = xr.where((da2.lon.round(2)<=area['lon2'])&(da2.lon.round(2)>=area['lon1']) , da2, np.nan)
# da4 = da3.dropna(dim='lon')
# da4

# # %%
# da4.lon
