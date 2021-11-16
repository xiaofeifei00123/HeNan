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

## 显示中文
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# %%


flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
ds = xr.open_dataset(flnm)
# %%
ds.lat.max()



# %%
class Draw():
    
    def __init__(self,):
        self.fontsize = 10
        self.path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/'   # 这里要改

    def draw_time_sequence(self,ax, dic, pic_dic):
        """[summary]

        Args:
            ax ([type]): [description]
            dic ([type]): 数据
            pic_dic ([type]): 图片元素控制
        """
        print(dic)

        QNSE = dic['OBS']
        x_label = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        x_label = x_label.dt.strftime('%d/%H')


        ccolor = ['green', 'blue', 'orange', 'red', 'black', 'blue','green','green','green','green'  ]
        lline_style = ['-', '-', '-', '-', '-', '-','-.', '-.', '-.', '-.']
        mmarker = ['o', 'o', '^', '^', '^', '^', '*', '*', '*', '*']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style)           
            # cycler(label=module_list)
            # cycler(marker=mmarker))
                        )
        
        j = 0
        ax.set_prop_cycle(custom_cycler)
        for i in dic.data_vars:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            # ax.plot(x_label, y, label=i, color=ccolor[j])
            # ax.plot(x_label, y, label=i, lw=4)
            # if i == 'obs':
                # i = 'OBS'
            ax.plot(x_label, y, label=i, lw=2.5, markersize=5)
            j +=1 

        # ax.set_xticks(x_label[::24])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_xticks(x_label[::4])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.xaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.yaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.set_xlabel("Date/Hour (UTC)", fontsize=self.fontsize*2.5)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.5)
        ax.legend(fontsize=self.fontsize*2.0, edgecolor='white')
    

    def draw_single(self, rain, pic_dic):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        # fig = plt.figure(figsize=(12, 8), dpi=200)  # 创建页面
        fig = plt.figure(figsize=(12, 8), dpi=600)  # 创建页面
        # ax = fig.add_axes([0.12, 0.25, 0.83, 0.7])
        ax = fig.add_axes([0.12, 0.12, 0.83, 0.8])
        y_max = pic_dic['yticks'].max()
        y_min = pic_dic['yticks'].min()
        ax.set_ylim(y_min, y_max)
        self.draw_time_sequence(ax, rain, pic_dic)
        # ax.legend(ncol=5 ,bbox_to_anchor=(0.5,-0.3) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        title = str(pic_dic['title'])
        print(pic_dic['title'])
        ax.set_title(title, fontsize=30)
        flnm = self.path+'time_sequence'+"_"+title+'.png'   # 这里要改
        fig.savefig(flnm)

def draw_area_max():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station_cst.nc'
    ds = xr.open_dataset(flnm)
    ds_max = ds.max(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-21 00'))
    dr = Draw()
    # ## 区域最大降水
    pic_dic_max = {
        # 'title':'max_rain',
        'title':'最大小时降水',
        'yticks':np.arange(0,225+25,25),
    }
    dr.draw_single(ds_max, pic_dic_max)

def draw_area_mean():
    # ### 区域平均降水
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station_cst.nc'
    ds = xr.open_dataset(flnm)
    dr = Draw()
    pic_dic_mean = {
        # 'title':'mean_rain',
        'title':'平均小时降水',
        'yticks':np.arange(0,3.6,0.5),
        # 'yticks':np.arange(0,20,1),
    }
    # ds_mean = get_data_mean()
    # ds_mean = ds.mean(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-20 12'))
    ds_mean = ds.mean(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-21 00'))
    dr.draw_single(ds_mean, pic_dic_mean)

def draw_zhenzhou():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    ds = xr.open_dataset(flnm)
    rain  = ds.sel(sta='57083')
    pic_dic_rain= {
        # 'title':'mean_rain',
        'title':'郑州站逐小时降水',
        'yticks':np.arange(0,210,5),
        # 'yticks':np.arange(0,20,1),
    }
    dr = Draw()
    dr.draw_single(rain, pic_dic_rain)

def draw_zhenzhou_max():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    ds = xr.open_dataset(flnm)
    area = {
        'lon1':110,
        'lon2':116,
        'lat1':32,
        'lat2':37,
        'interval':0.125,
    }
    index = ((ds.lat<=area['lat2']) & (ds.lat>=area['lat1']) & (ds.lon>=area['lon1']) & (ds.lon<=area['lon2']))
    rain = ds.sel(sta=index)
    rain = rain.max(dim='sta')
    pic_dic_rain= {
        'title':'郑州站及周边最大小时降水',
        'yticks':np.arange(0,210,5),
    }
    dr = Draw()
    dr.draw_single(rain, pic_dic_rain)

def draw_zhenzhou_mean():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    ds = xr.open_dataset(flnm)
    area = {
        'lon1':112.7,
        'lon2':114.2,
        'lat1':34.3,
        'lat2':34.9,
        'interval':0.125,
    }
    index = ((ds.lat<=area['lat2']) & (ds.lat>=area['lat1']) & (ds.lon>=area['lon1']) & (ds.lon<=area['lon2']))
    rain = ds.sel(sta=index)
    rain = rain.mean(dim='sta')
    pic_dic_rain= {
        'title':'郑州站及其周边站点平均降水逐小时变化',
        'yticks':np.arange(0,26,5),
    }
    dr = Draw()
    dr.draw_single(rain, pic_dic_rain)
if __name__ == '__main__':
    # main()
    # draw_zhenzhou()
    # draw_zhenzhou_max()
    # draw_zhenzhou_mean()
    pass
    # %%
    # import xarray as xr
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station_cst.nc'
    # ds = xr.open_dataset(flnm)
    # ds_max = ds.max(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-21 00'))
    # dr = Draw()
    # # ## 区域最大降水
    # pic_dic_max = {
    #     # 'title':'max_rain',
    #     'title':'最大小时降水',
    #     'yticks':np.arange(0,225+25,25),
    # }
    # dr.draw_single(ds_max, pic_dic_max)

    # ### 区域平均降水
    # pic_dic_mean = {
    #     # 'title':'mean_rain',
    #     'title':'平均小时降水',
    #     'yticks':np.arange(0,3.6,0.5),
    #     # 'yticks':np.arange(0,20,1),
    # }
    # # ds_mean = get_data_mean()
    # # ds_mean = ds.mean(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-20 12'))
    # ds_mean = ds.mean(dim='sta').sel(time=slice('2021-07-19 12', '2021-07-21 00'))

    # dr.draw_single(ds_mean, pic_dic_mean)


