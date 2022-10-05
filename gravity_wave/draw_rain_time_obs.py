#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
格点观测降水的区域平均值随时间变化曲线
-----------------------------------------
Time             :2022/08/23 20:07:20
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from read_rain_wrf import GetData
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
from common import Common
from draw_rain_distribution_minus import Rain

# %%
class RainTime(Common):
    def __init__(self):
        super().__init__()  # 继承父类
        pass

class Data(RainTime):
    def __init__(self):
        super().__init__()

    def caculate_area_mean_obs(self, da,area):
        mask = (
            (da.coords['lat']>area['lat1'])
            &(da.coords['lat']<area['lat2'])
            &(da.coords['lon']<area['lon2'])
            &(da.coords['lon']>area['lon1'])
        )
        aa = xr.where(mask, 1, np.nan)
        db = da*aa
        dsr = db.mean(dim=['lat', 'lon'])
        return dsr

    def get_data_1area_obs(self, area):
        """
        area = {
            'lat1':32,
            'lat2':36.5,
            'lon1':110.5,
            'lon2':116,
            }        
        """
        flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
        ds_obs = xr.open_dataset(flnm_obs)
        ds_obs_mean  = self.caculate_area_mean_obs(ds_obs, area)
        ds = ds_obs_mean
        ds = ds.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
        return ds

    def get_data_dual_area_obs(self, ):
        """合并三个时次的降水数据
        """
        # com = Common()
        dsA = self.get_data_1area_obs(self.areaA)
        dsB = self.get_data_1area_obs(self.areaB)
        dsC = self.get_data_1area_obs(self.areaC)

        da1 = dsA['PRCP']
        da = dsA.resample(time='12H', closed='right', label='right').sum()
        da12 = da['PRCP']

        db1 = dsB['PRCP']
        db = dsB.resample(time='12H', closed='right', label='right').sum()
        db12 = db['PRCP']

        dc1 = dsC['PRCP']
        dc = dsC.resample(time='12H', closed='right', label='right').sum()
        dc12 = dc['PRCP']

        da1.name = 'A1'
        da12.name = 'A12'
        db1.name = 'B1'
        db12.name = 'B12'
        dc1.name = 'C1'
        dc12.name = 'C12'
        # ds = xr.merge([da1, da12, db1, db12])
        ds = xr.merge([da1, da12, db1, db12, dc1, dc12])
        return ds

class Draw(RainTime):
    def __init__(self):
        super().__init__()

    # def draw(self, ds, fig, ax, *args, **kw):
    #     # color_list = ['black', 'green', 'blue', 'red', 'orange']
    #     # color_list = ['black', 'blue', 'red','green', 'orange']
    #     color_list = ['black', 'red', 'blue','orange', 'green']
    #     color_list = ds.color_list
    #     linestyle_list = ds.line_list
    #     var_list = list(ds.data_vars)
    #     i = 0
    #     for var in var_list:
    #         da = ds[var]
    #         x = da.time.dt.strftime('%d/%H')
    #         y = da.values
    #         if var == 'PRCP':
    #             var = 'OBS'
    #         ax.plot(x,y, label=var, color='black',linestyle=linestyle_list[i], **kw)
    #         i+=1
    #     ax.legend(edgecolor='white')

    #     ax.set_xticks(x[::12])
    #     ax.set_xticklabels(x[::12].values, rotation=30, fontsize=10)
    #     ax.xaxis.set_minor_locator(AutoMinorLocator())
    #     ax.yaxis.set_minor_locator(AutoMinorLocator())
    #     ax.set_ylim(0, 25)

    def draw_rain_time(self, ds, ax2):

        # da = dsA['PRCP']
        ### 获得数据
        x = ds['A1'].time.dt.strftime('%d/%H')
        ya1 = ds['A1'].values
        ya12 = ds['A12'].values
        yb1 = ds['B1'].values
        yb12 = ds['B12'].values
        yc1 = ds['C1'].values
        yc12 = ds['C12'].values


        ## 小时降水
        # cm = 1/2.54
        # # fig = plt.figure(figsize=(12*cm, 5*cm), dpi=300)
        # # ax2  = fig.add_axes([0.14, 0.2, 0.73, 0.7])

        ax2.plot(x,yb1,  color='black', label='B_1h')
        ax2.plot(x,yc1,  color='#5e5e5e', linestyle='--', label='C_1h')


        ax2.set_ylabel('Hourly rainfall (mm)')
        ax2.set_xlabel('Time (Date/Hour)')
        ax2.set_xticks(x[::24])
        ax2.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.set_ylim(0, 30)
        ax2.yaxis.set_minor_locator(AutoMinorLocator())

        ## 12小时降水
        ax1 = ax2.twinx()
        width = 2
        x1 = np.arange(len(x))
        # ax1.bar(x1-width, yb12, width=2, color='red', label='B_12h', alpha=1.)
        # ax1.bar(x1, yc12, width=2, color='blue', label='C_12h')
        ax1.bar(x1-width/2, yb12, width=2, color='red', label='B_12h', alpha=1.)
        ax1.bar(x1+width/2, yc12, width=2, color='blue', label='C_12h')
        ax1.set_ylabel('12-h rainfall (mm)')
        ax1.set_ylim(0, 250)
        ax1.yaxis.set_minor_locator(AutoMinorLocator())

def draw_obs(ax2, ):
    gd = Data()
    ds = gd.get_data_dual_area_obs()
    dr = Draw()
    dr.draw_rain_time(ds, ax2)

def draw_obs_main():
    cm = 1/2.54
    fig = plt.figure(figsize=(9*cm, 5*cm), dpi=300)
    ax2  = fig.add_axes([0.12, 0.2, 0.73, 0.7])
    draw_obs(ax2,)

    fig.legend(loc='upper left', edgecolor='white', bbox_to_anchor=(0.13, 0.7, 0.2, 0.2))
    figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
    fig.savefig(figpath+'rain_time')
    pass

if __name__ == '__main__':
    draw_obs_main()
    
