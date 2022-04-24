#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
逐日降水的时间变化曲线
1. 这里的时间是从当天的08时到第二天的08时(BJT), 也就是世界时(00-23)
2. 24小时降水变化曲线(平均), 各站点分开
3. 多站点求平均，上面两张图应该都要有
站点降水的时间序列
-----------------------------------------
Time             :2021/06/04 14:32:20
Author          :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
# import datetime
# import sys,os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from baobao.caculate import caculate_average_wrf


# %%

# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'
# da = xr.open_dataarray(flnm)
# da  = da.rename({'lon':'XLONG', 'lat':'XLAT'})
# caculate_average_wrf(da)
# %%
def gd_mean():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station1.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
    ## 伏牛山
    area_all = {
        'lon1':111.6,
        'lon2':113.3,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    ## 伏牛山左边
    area_left = {
        'lon1':111.6,
        'lon2':112.4,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    ## 伏牛山右边
    area_right = {
        'lon1':112.4,
        'lon2':113.3,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }

    area_list = [area_all, area_left, area_right]
    # name_list = ['all', 'left', 'right']
    name_list = ['A+B', 'A', 'B']
    # name_list 
    rain_dict = {}
    for area,j in zip(area_list, name_list):
        index = ((ds.lat<=area['lat2']) & (ds.lat>=area['lat1']) & (ds.lon>=area['lon1']) & (ds.lon<=area['lon2']))
        rain = ds.sel(sta=index)
        rain_dict[j] = rain.mean(dim='sta')
    return rain_dict

def gd_mean_wrf():
    """wrf模拟区域，格点平均降水

    Returns:
        _type_: _description_
    """
    pass
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station1.nc'
    # ds = xr.open_dataset(flnm)
    # ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))

    ## gwd0的区域平均降水
    flnm1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'
    da1 = xr.open_dataarray(flnm1)
    da1  = da1.rename({'lon':'XLONG', 'lat':'XLAT'})
    # da_gwd0 = caculate_average_wrf(da1)
    da1 = da1.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
    
    ## gwd3的区域平均降水
    flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain.nc'
    da2 = xr.open_dataarray(flnm2)
    da2  = da2.rename({'lon':'XLONG', 'lat':'XLAT'})
    # da_gwd3 = caculate_average_wrf(da2)
    da2 = da2.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
    # ds = {'gwd0':da_gwd0,'gwd3':da_gwd3}
    

    ## 伏牛山
    area_all = {
        'lon1':111.6,
        'lon2':113,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    ## 伏牛山左边
    area_left = {
        'lon1':111.6,
        'lon2':112.3,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }
    ## 伏牛山右边
    area_right = {
        'lon1':112.3,
        'lon2':113.0,
        'lat1':33.4,
        'lat2':33.8,
        'interval':0.125,
    }

    area_list = [area_all, area_left, area_right]
    name_list = ['all', 'left', 'right']
    # name_list 
    rain_dict = {}
    for area,j in zip(area_list, name_list):
        da_gwd0 = caculate_average_wrf(da1, area)
        da_gwd3 = caculate_average_wrf(da2, area)
        ds = {'gwd0':da_gwd0,'gwd3':da_gwd3}
        # rain = 
        # index = ((ds.lat<=area['lat2']) & (ds.lat>=area['lat1']) & (ds.lon>=area['lon1']) & (ds.lon<=area['lon2']))
        # rain = ds.sel(sta=index)
        rain_dict[j] = ds
    return rain_dict

def set_ticks(ax):
    fontsize = 10 
    ax.xaxis.set_tick_params(labelsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=fontsize)
    # ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
    # ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlabel("Date/Hour (UTC)", fontsize=fontsize*1.2)
    ax.set_ylabel("Precipitation (mm)", fontsize=fontsize*1.2)
    # ax.legend(fontsize=fontsize, edgecolor='white')
        
def draw_time_sequence(ax, da, label = 'obs', color='red', line_style='solid', mmarker='O'):
    """[summary]
    Args:
        ax ([type]): [description]
        dic ([type]): 数据
        pic_dic ([type]): 图片元素控制
    """
    # print(dic)

    x_label = da.coords['time']
    x_label = x_label.dt.strftime('%d/%H')
    y = da.values.round(2)
    ax.plot(x_label, y, label=label, lw=1, markersize=5, color=color, linestyle=line_style)
    ax.set_xticks(x_label[::6])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    set_ticks(ax)

        

# rain = gd_mean_wrf()
rain = gd_mean()
cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)  # 创建页面
ax = fig.add_axes([0.15, 0.2, 0.82, 0.73])

# for i in 
# area_list = []
# area_list = ['all', 'left', 'right']
# style_list = ['-', '--', ':']
# # model_list = ['gwd0', 'gwd3', 'OBS']
# model_list = ['gwd0', 'gwd3', ]
# color_list = ['blue', 'red', 'black']
# for area,lnstyle in zip(area_list, style_list):
#     # line_style = j
#     for model,color in zip(model_list, color_list):
#         draw_time_sequence(ax,rain[area][model],label=area+model, color=color, line_style=lnstyle)

# area_list = ['all', 'left', 'right']
area_list = ['A+B', 'A', 'B']
style_list = ['-', '--', ':']
# model_list = ['gwd0', 'gwd3', 'OBS']
model_list = ['gwd0', 'gwd3', ]
color_list = ['green','blue', 'red', ]
for area,color in zip(area_list, color_list):
    # line_style = j
    for model,lnstyle in zip(model_list, style_list):
        if model == 'gwd0':
            label = area+'_'+'CTRL'
        elif model == 'gwd3':
            label = area+'_'+'GWD3'
        else:
            label = area+'_'+model
        draw_time_sequence(ax,rain[area][model],label=label, color=color, line_style=lnstyle)

ax.set_ylim(0,25)
ax.legend(ncol=2, fontsize=10, edgecolor='white', loc='upper left')
ax.set_title('(b)', loc='left', y=0.97, fontsize=10)
fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
fig_name = 'time_sequence'
fig.savefig(fig_path+fig_name)


