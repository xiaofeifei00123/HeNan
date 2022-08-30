# %%
import imp
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from common import Common

from draw_rain_time import get_data as rain_gd
from draw_rain_time import draw as rain_dr

# %%

# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_A.nc'
def get_data(flnm):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_A.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
    hh = ds['height'].mean(dim='time').values
    ds2 = ds.assign_coords({'z':('pressure',hh)})
    ds3 = ds2.swap_dims({'pressure':'z'})
    da = ds3['div']#.sel(z=1000, method='nearest')
    db = da.sel(z = np.sort(da.z))
    dc = db.sel(z=1000, method='nearest')
    dc = dc*10**5
    return dc
# dc = get_data()

# dc
# %%
cm = 1/2.54
fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
ax  = fig.add_axes([0.1, 0.15, 0.8, 0.8])
area = {
    'lat1':33.5,
    'lat2':36.0,
    'lon1':112,
    'lon2':115,
    }        

cm = Common()
# area = cm.areaA
ds = rain_gd(area)
# ds = get_data(cm.areaB)
cm = 1/2.54
fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
ax  = fig.add_axes([0.1, 0.15, 0.8, 0.8])
ax.set_ylabel('Precipitation (mm)')
ax.set_xlabel('Time (Date/Hour)')
rain_dr(ds, fig, ax)



ax.set_ylabel('Precip (mm)/ Div($10^{-5}\cdot s^{-1}$)')
ax.set_xlabel('Time (Date/Hour)')
# color_list = ['black', 'green', 'blue', 'red', 'orange']
color_list = ['green', 'blue', 'red', 'orange']
model_list = ['CTRL', 'SS', 'FD', 'GWD3']

i = 0
for model in model_list:
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'+model+'/wrfout/time_cross_D.nc'
    # ds = xr.open_dataset(flnm)
    da = get_data(flnm)
    x = da.time.dt.strftime('%d/%H')
    y = da.values
    # ax.plot(x,y, label=model, color=color_list[i], ls='--')
    ax.plot(x,y, color=color_list[i], ls='--')
    i+=1
ax.legend(edgecolor='white')
ax.set_xticks(x[::24])
ax.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
ax.axhline(y=0, color='black')

fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/'
fig.savefig(fig_path+'div')
# #     # cm.areaA
# # if __name__ == "__main__":
# #     main()
# %%