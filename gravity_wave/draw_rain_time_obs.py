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
from turtle import width
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from read_rain_wrf import GetData
from common import Common
from draw_rain_distribution_minus import Rain

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
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_obs = xr.open_dataset(flnm_obs)
    ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)
    ds = ds_obs_mean
    ds = ds.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    return ds

def draw(ds, fig, ax, *args, **kw):
    # color_list = ['black', 'green', 'blue', 'red', 'orange']
    # color_list = ['black', 'blue', 'red','green', 'orange']
    color_list = ['black', 'red', 'blue','orange', 'green']
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
        ax.plot(x,y, label=var, color='black',linestyle=linestyle_list[i], **kw)
        i+=1
    ax.legend(edgecolor='white')

    ax.set_xticks(x[::12])
    ax.set_xticklabels(x[::12].values, rotation=30, fontsize=10)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylim(0, 25)
# %%
# def main():
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
dsB = get_data(areaB)
dsA = get_data(areaA)

db1 = dsB['PRCP']
db = dsB.resample(time='12H', closed='right', label='right').sum()
db12 = db['PRCP']

da1 = dsA['PRCP']
da = dsA.resample(time='12H', closed='right', label='right').sum()
da12 = da['PRCP']

da1.name = 'A1'
da12.name = 'A12'
db1.name = 'B1'
db12.name = 'B12'
ds = xr.merge([da1, da12, db1, db12])

# %%
cm = 1/2.54
fig = plt.figure(figsize=(12*cm, 5*cm), dpi=300)
ax2  = fig.add_axes([0.14, 0.2, 0.73, 0.7])
# da = dsA['PRCP']
x = ds['A1'].time.dt.strftime('%d/%H')
ya1 = ds['A1'].values
ya12 = ds['A12'].values
yb1 = ds['B1'].values
yb12 = ds['B12'].values

ax2.plot(x,ya1,  color='black', label='A_1h')
ax2.plot(x,yb1,  color='#5e5e5e', linestyle='--', label='B_1h')


# ax2.set_ylabel('Precip Hourly(mm)')
ax2.set_ylabel('Preciptation (mm)')
ax2.set_xlabel('Time (Date/Hour)')
ax2.set_xticks(x[::24])
ax2.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
# ax2.set_ylim(0, 230)
# draw(ds, fig, ax)
# x1 = ds['B'].time.dt.strftime('%d/%H')
# y1 = ds['B'].values
ax1 = ax2.twinx()
# ax1.bar(x,ya12, width=2, color='red', label='12h')
# ax1.bar(x,yb12, width=2, color='red', label='12h')
width = 2


x1 = np.arange(len(x))
x1
ax1.bar(x1-width/2, ya12, width=2, color='red', label='A_12h')
ax1.bar(x1+width/2, yb12, width=2, color='blue', label='B_12h')

# ax1.bar(x+width/2, yb12, width=2, color='red', label='12h')
# ax1.legend()
# ax2.legend()
# fig.legend(loc='upper left', edgecolor='white')#, bbox_to_anchor=(0.65, 0.7, 0.2, 0.2))
fig.legend(loc='upper left', edgecolor='white', bbox_to_anchor=(0.13, 0.7, 0.2, 0.2))
# ax1.set_ylim(0, 250)
# ax1.set_ylabel('Precip (mm)')
figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
# %%
fig.savefig(figpath+'rain_time')

# if __name__ == '__main__':
#     main()
    