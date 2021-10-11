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
import datetime
import sys,os
import xarray as xr
import numpy as np
import pandas as pd

import salem  # 插值
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import geopandas
import cmaps
from get_cmap import get_cmap_rain2
from multiprocessing import Pool

## 显示中文
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# %%
def get_data_mean():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=slice('2021-07-20 00', '2021-07-20 12'))
    ds_all = ds.mean(dim=['lat', 'lon'])
    return ds_all


# %%
class GetDataMax():
    def get_rain_ec_max(self, ):
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_ec.nc'
        ds = xr.open_dataarray(flnm)

        tt = ds.time
        r_list = []
        for t in tt:
            daa = ds.sel(time=t)
            r_max = daa.max().values
            r_list.append(r_max)

        ps = pd.Series(r_list, index=tt.values)
        da = xr.DataArray.from_series(ps)
        rain_max = da.rename({'index':'time'}).astype('float32')
        return rain_max
    # print(aa)

    def get_max_dataframe(self,):
        """读取csv格式的站点数据

        Args:
            df_station (DataFrame): 输入csv格式站点数据

        Returns:
            [DataArray]: 小时雨强最大值
        """
        ## 才读的数据，它的时间格式是object
        df_station = pd.read_csv('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.csv')
        # df_station = get_max_dataframe(aa)
        df_station['time']= pd.to_datetime(df_station['time'])
        t = pd.date_range(start='2021-07-20 00', end='2021-07-21 00', freq='1H')
        rain_list = []
        for tt in t:
            cc = df_station[df_station['time']==tt]
            rain_max = cc[(cc['lat']>32)&(cc['lat']<37)&(cc['lon']>110)&(cc['lon']<116)]['data0'].max()
            # rain_max = cc[(cc['lat']>32)&(cc['lat']<37)&(cc['lon']>110)&(cc['lon']<116)]['data0'].mean()
            rain_list.append(rain_max)
        rain_list
        ps = pd.Series(rain_list, index=t)
        da = xr.DataArray.from_series(ps)
        rain_obs_max = da.rename({'index':'time'})
        return rain_obs_max


    def get_rain_gfs_max(self, ):
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/YSU_GFS_rain.nc'
        ds = xr.open_dataset(flnm)
        # print(ds)
        ds = ds['RAINNC']
        tt = ds.time
        r_list = []
        for t in tt:
            daa = ds.sel(time=t)
            r_max = daa.max().values
            # r_max = daa.mean().values
            r_list.append(r_max)

        ps = pd.Series(r_list, index=tt.values)
        da = xr.DataArray.from_series(ps)
        rain_max = da.rename({'index':'time'}).astype('float32')
        return rain_max

    def get_rain_wrf_max(self, flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/YSU_GFS_rain.nc'):
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/YSU_GFS_rain.nc'
        ds = xr.open_dataset(flnm)
        # print(ds)
        ds = ds['RAINNC']
        tt = ds.time
        r_list = []
        for t in tt:
            daa = ds.sel(time=t)
            print(daa)
            r_max = daa.max().values
            # r_max = daa.mean().values
            r_list.append(r_max)

        ps = pd.Series(r_list, index=tt.values)
        da = xr.DataArray.from_series(ps)
        rain_max = da.rename({'index':'time'}).astype('float32')
        return rain_max
# r1 = get_rain_gfs_max()
# r2 = get_rain_ec_max()
# aa = pd.read_csv('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.csv')
# r3 = get_max_dataframe(aa)
# print(aa)
# %%
# ds = 
# ds = xr.Dataset()
# ds['gfs'] = r1.sel(time=slice('2021-07-20 00', '2021-07-20 12'))
# # ds['ec'] = r2.sel(time=slice('2021-07-20 00', '2021-07-20 12'))
# ds['obs'] = r3.sel(time=slice('2021-07-20 00', '2021-07-20 12'))

# %%

def draw_forecast():
    pass
    type_list = ['ERA5', 'GDAS']
    time_list = ['1800', '1812', '1900', '1912']
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/'

    dr = Draw()
    gd = GetDataMax()
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_rain.nc'
    ds = xr.Dataset()
    for type in type_list:
        for t in time_list:
            flnm = path+type+'/'+'YSU_'+t+'_rain.nc'
            # print(flnm)
            # da = gd.get_rain_wrf(flnm)
            da_max = gd.get_rain_wrf_max(flnm)
            ds[type+t] = da_max
            # picture_dic = {'date':'2000_2012', 'type':type, 'initial_time':t}
            # dr.draw_single(da, picture_dic)
    return ds



# %%
class Draw():
    
    def __init__(self,):
        self.fontsize = 10
        self.path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/'   # 这里要改
        # self.module_list = ['obs', 'gfs']
        self.module_list = ['OBS', 'GFS', 
                            'ERA51800',
                            'ERA51812',
                            'ERA51900',
                            'ERA51912',
                            'GDAS1800',
                            'GDAS1812',
                            'GDAS1900',
                            'GDAS1912',
                            ]

    def draw_time_sequence(self,ax, dic, pic_dic):
        """[summary]

        Args:
            ax ([type]): [description]
            dic ([type]): 数据
            pic_dic ([type]): 图片元素控制
        """

        QNSE = dic['OBS']
        x_label = QNSE.coords['time']
        # x_label = str(x_label.dt.strftime('%d%H').values).split()
        x_label = x_label.dt.strftime('%d_%H')


        # x_label = x_label.dt.strftime("%H")  # 转换时间维字符串格式
        # y = QNSE.values
        # module_list = ['obs', 'ACM2','YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        # x = np.arange(len(x_label))
        # print(x)


        ccolor = ['black','m', 'green', 'cyan', 'blue', 'red', 'green','cyan','blue','red'  ]
        # ccolor = ['black','red', 'blue', 'blue', 'blue', 'blue', 'green','green','green','green'  ]
        # lline_style = ['-', '-', '-', '--', '-.', ':','-', '--', '-.', ':']
        lline_style = ['-', '-', '--', '--', '--', '--',':', ':', ':', ':']
        mmarker = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
        custom_cycler = (
            cycler(color=ccolor) +
            cycler(linestyle=lline_style)           
            # cycler(label=module_list)
            # cycler(marker=mmarker))
                        )
        
        j = 0
        ax.set_prop_cycle(custom_cycler)
        for i in self.module_list:
            # y = dr.loc[i,:].values
            y = dic[i].values
            y = np.around(y,2)
            # ax.plot(x_label, y, label=i, color=ccolor[j])
            # ax.plot(x_label, y, label=i, lw=4)
            if i == 'obs':
                i = 'OBS'
            ax.plot(x_label, y, label=i, lw=2.5, markersize=5)
            j +=1 

        # ax.set_xticks(x_label[::24])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        ax.set_xticks(x_label[::2])  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.set_yticks(np.arange(0, 230, 20))
        # ax.set_yticks(pic_dic['yticks'])
        
        # ax.set_yticks(np.arange(0, 20.1, 1))
        # ax.xaxis.set_tick_params(labelsize=15)
        # ax.xaxis.set_tick_params(labelsize=self.fontsize*1.8, rotation=45)
        ax.xaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.yaxis.set_tick_params(labelsize=self.fontsize*2.0)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
        # ax.set_yticks(np.arange(0, 5.01, 0.1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
        # ax.set_ylim(0,5.1)
        # ax.set_ylim(0.0,0.6)
        ax.set_xlabel("Time(Jul2021, UTC)", fontsize=self.fontsize*2.5)
        ax.set_ylabel("Precipitation (mm)", fontsize=self.fontsize*2.5)
        # ax.legend()
        # ax.set_title("201607", fontsize=18)
        # fig.savefig('/home/fengxiang/Project/Asses_pbl_July/Draw/Rain/time_sequecnce.png')
    

    def draw_single(self, rain, pic_dic):
        """[summary]

        Args:
            rain ([DataArray]): 一个模式的降水
        """
        fig = plt.figure(figsize=(12, 8), dpi=200)  # 创建页面
        ax = fig.add_axes([0.12, 0.2, 0.83, 0.7])
        self.draw_time_sequence(ax, rain, pic_dic)
        ax.legend(ncol=5 ,bbox_to_anchor=(0.5,-0.3) ,loc='lower center',fontsize=self.fontsize*2.0, edgecolor='white')
        title = pic_dic['title']
        ax.set_title(title, fontsize=30)
        # # flnm = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/rain_staion1_'+self.month+'.png'   # 这里要改
        flnm = self.path+'time_sequence'+"_"+title+'.png'   # 这里要改
        # fig.suptitle('The mean time sequence', fontsize=self.fontsize*2.0)
        fig.savefig(flnm)

def get_data_max():
    aa = draw_forecast()
    gx = GetDataMax()
    r1 = gx.get_max_dataframe()
    r2 = gx.get_rain_wrf_max()
    aa['OBS'] = r1
    aa['GFS'] = r2
    ds_all = aa.sel(time=slice('2021-07-20 00', '2021-07-20 12'))
    return ds_all
# # %%
# ds_all['ERA51900']
if __name__ == '__main__':
    # main()
    dr = Draw()

    ## 区域最大降水
    pic_dic_max = {
        'title':'max_rain',
        'yticks':np.arange(0,230,20),
    }
    ds_max = get_data_max()
    dr.draw_single(ds_max, 'max_rain')

    ## 区域平均降水
    pic_dic_mean = {
        'title':'mean_rain',
        'yticks':np.arange(0,4,0.5),
    }
    ds_mean = get_data_mean()
    dr.draw_single(ds_mean, pic_dic_mean)


