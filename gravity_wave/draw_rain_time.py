#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

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
from read_rain_wrf import GetData
from common import Common

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

# %%
def get_data(area):
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_model = xr.open_dataset(flnm_model)
    ds_obs = xr.open_dataset(flnm_obs)
    # ds_obs.time.values = ds_obs.time.values+pd.Timedelta('12H')

    # tt = ds_obs.time.values+pd.Timedelta('8H')
    # ds_obs = ds_obs.assign_coords({'time':tt})

    cm = Common()
    gd = GetData()
    ds_model_mean = gd.caculate_area_mean(ds_model, area)
    # ds_model_mean.sel(model='gwd3')
    ds_model_mean = ds_model_mean['GWD3']
    # print(ds_model_mean.dims)
    ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)
    dsall = ds_obs_mean
    dsall = xr.merge([ds_obs_mean, ds_model_mean])
    ds = dsall.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    # ds = dsall.sel(time=slice('2021-07-17 00', '2021-07-19 00'))
    # ds = ds.resample(time='12H').sum()
    return ds

def draw(ds, fig, ax):
    color_list = ['black', 'green', 'blue', 'red', 'orange']
    var_list = list(ds.data_vars)
    i = 0
    for var in var_list:
        da = ds[var]
        x = da.time.dt.strftime('%d/%H')
        y = da.values
        if var == 'PRCP':
            var = 'OBS'
        ax.plot(x,y, label=var, color=color_list[i])
        i+=1
    ax.legend(edgecolor='white')
    ax.set_xticks(x[::24])
    ax.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
    ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
    # ax.set_xticks(x[::2])
    # ax.set_xticklabels(x[::2].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(1))

# ds['FD'].resample(time='3H').sum()
def main():
    area = {
        'lat1':33.5,
        'lat2':36.0,
        'lon1':112,
        'lon2':115,
        }        
    # area = {
    #     'lat1':32,
    #     'lat2':36.5,
    #     'lon1':110.5,
    #     'lon2':116,
    #     }        

    cm = Common()
    # area = cm.areaA
    ds = get_data(area)
    # ds = get_data(cm.areaB)
    cm = 1/2.54
    fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
    ax  = fig.add_axes([0.1, 0.15, 0.8, 0.8])
    ax.set_ylabel('Precipitation (mm)')
    ax.set_xlabel('Time (Date/Hour)')
    draw(ds, fig, ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/'
    fig.savefig(fig_path+'core')
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