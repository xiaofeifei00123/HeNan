#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
区域平均的降水变化
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
from common import Common
from draw_rain_distribution_minus import Rain
# %%
# flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
# ds_model = xr.open_dataset(flnm_model)
# ds_obs = xr.open_dataset(flnm_obs)
# flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/rain_model.nc'
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
# ds_model = xr.open_dataset(flnm_model)
# ds_obs = xr.open_dataset(flnm_obs)
# gd = GetData()
# area = {
#     'lat1':32,
#     'lat2':36.5,
#     'lon1':110.5,
#     'lon2':116,
#     }        
# ds_model_mean = gd.caculate_area_mean(ds_model, area)
# ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)
# # %%
# # ds_model_mean
# # ds_obs_mean
# xr.merge([ds_obs_mean, ds_model_mean])

# %%
# ds_obs
# ds2 = xr.merge([ds_model['GWD3'], ds_model['CTRL'], ds_])

# %%
def caculate_area_mean_obs(da,area):
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

def get_data(area):
    """
    area = {
        'lat1':32,
        'lat2':36.5,
        'lon1':110.5,
        'lon2':116,
        }        
    """
    # flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    # flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/rain_model.nc'
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model_da.nc'
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_model = xr.open_dataset(flnm_model)
    ds_obs = xr.open_dataset(flnm_obs)
    gd = GetData()
    ds_model_mean = gd.caculate_area_mean(ds_model, area)
    ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)

    # ds = xr.merge([ds_obs_mean, ds_model_mean])
    # dsall = dsall.sel(model=['GWD3', 'CTRL'])
    # ds = dsall.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    ds = xr.merge([ds_model_mean['GWD3'], ds_model_mean['CTRL'],ds_obs_mean['PRCP']])
    # ds = ds_model_mean

    ds = ds.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    # ds = ds.resample(time='3H').sum()
    return ds

def draw(ds, fig, ax, *args, **kw):
    # color_list = ['black', 'green', 'blue', 'red', 'orange']
    # color_list = ['black', 'blue', 'red','green', 'orange']
    # color_list = ['black', 'red', 'blue','orange', 'green']
    color_list = ds.color_list
    linestyle_list = ds.line_list
    var_list = list(ds.data_vars)
    i = 0
    for var in var_list:
        da = ds[var]
        x = da.time.dt.strftime('%d/%H')
        y = da.values
        if var == 'PRCP':
            var = 'OBS'
        # ax.plot(x,y, label=var, color=color_list[i], **kw)
        ax.plot(x,y, label=var, color=color_list[i],linestyle=linestyle_list[i], **kw)
        i+=1
    ax.legend(edgecolor='white')

    # ax.set_xticks(x[::24])
    # ax.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(6))

    ax.set_xticks(x[12::12])
    ax.set_xticklabels(x[12::12].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(47, 124)
    # ax.set_ylim(0, 20)
    # ax.set_xticks(x[::2])
    # ax.set_xticklabels(x[::2].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
# %%
# ds['FD'].resample(time='3H').sum()
def main():
    # area = {
    #     'lat1':33.4,
    #     'lat2':33.8,
    #     'lon1':112.6,
    #     'lon2':113.0,
    # }        
    # ra = Rain()

    areaA = {
        'lat1':34.4,
        'lat2':34.9,
        'lon1':113.0,
        'lon2':113.7,
    }        
    areaB = {
        'lat1':35.3,
        'lat2':35.8,
        'lon1':113.5,
        'lon2':114.2,
    }        
    dsA = get_data(areaA)
    dsB = get_data(areaB)
    dsA = dsA.rename({'GWD3':'GWD3A', 'CTRL':'CTRLA', 'PRCP':'OBSA'})
    dsA = xr.merge([dsA['GWD3A'], dsA['OBSA']])
    dsB = dsB.rename({'GWD3':'GWD3B', 'CTRL':'CTRLB', 'PRCP':'OBSB'})
    ds = xr.merge([dsA, dsB])
    # ds = dsA
    cm = 1/2.54
    # fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
    fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
    ax  = fig.add_axes([0.1, 0.15, 0.85, 0.8])
    ax.set_ylabel('Precipitation (mm)')
    ax.set_xlabel('Time (Date/Hour)')
    # ds.attrs['color_list'] = ['red', 'black', 'red', 'black']
    # ds.attrs['line_list'] = ['-', '-', '--', '--']
    # ds.attrs['color_list'] = ['red', 'green', 'black', 'red', 'green', 'black']
    ds.attrs['color_list'] = ['red', 'black', 'black', 'red', 'green', 'black']
    ds.attrs['line_list'] = ['-', '-', '-', '--', '--', '--']
    draw(ds, fig, ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
    # fig.savefig(fig_path+'time_sAB')
    fig.savefig(fig_path+'model_A')
# cm.areaA
if __name__ == "__main__":
    main()

# %%

# import xarray as xr
# # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/2021-07-19-12__2021-07-21-00/all.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/2021-07-20-12__2021-07-22-00/all.nc'
# da = xr.open_dataarray(flnm)
# da
# # %%
# gd = GetData()
# cm = Common()
# area = {
#     'lat1':33.5,
#     'lat2':36.0,
#     'lon1':112.5,
#     'lon2':114.5,
#     }        
# daa = gd.caculate_area_mean(da, area)
# daa
# # %%
# daa.plot()
# # da.mean(dim=['south'])
