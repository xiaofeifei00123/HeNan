# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_d03.nc'
flnm_all = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_d02_all.nc'
flnm_con = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_d02_con.nc'
flnm_grd = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_d02_grd.nc'
da_all = xr.open_dataarray(flnm_all)
da_con = xr.open_dataarray(flnm_con)
da_grd = xr.open_dataarray(flnm_grd)
# %%

# %%
# def get_xy(da):
#     rain = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
#     x = rain.time.dt.strftime('%H').values
#     y = rain.values
#     daa  = xr.where((rain.lat<35)&(rain.lat>33.5)&(rain.lon>112)&(rain.lon<113.5), rain, np.nan)
#     # daa = rain
#     y = daa.mean(dim=['south_north', 'west_east']).values
#     return x, y

def get_xy(da):    
    rain = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    rain  = xr.where((rain.lat<35)&(rain.lat>33.5)&(rain.lon>112)&(rain.lon<113.5), rain, np.nan)
    r3 = rain.mean(dim=['south_north', 'west_east'])
    df = r3.to_series()
    rr = df.resample('3H', label='right', closed='right').sum()  # 右端是闭合的
    x = rr.index.strftime('%H').values
    return x,rr

x1,y_all = get_xy(da_all)
x1,y_con = get_xy(da_con)
x1,y_grd = get_xy(da_grd)
y_list = [y_all, y_grd, y_con]
# %%
labels = x1
labels
x = np.arange(len(labels))

cm = 1/2.54
fig = plt.figure(figsize=(10*cm, 7*cm), dpi=300)
ax = fig.add_axes([0.2,0.2,0.7,0.7])
width = 0.3
rects = ax.bar(x-0.5*width, y_grd, width, label='grid', color='blue')
rects = ax.bar(x+0.5*width, y_con, width, label='convective', color='red')
ax.legend(edgecolor='white')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_yticks(np.arange(0,20,2))
ax.set_ylabel('precip (mm)')
ax.set_xlabel('Hour 2021-07-20~21 UTC')
fig.savefig('rain_time.png')
# %%
# rain = da_all.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
# r3 = rain.mean(dim=['south_north', 'west_east'])
# # %%
# df = r3.to_series()
# rr = df.resample('3H', label='right', closed='right').sum()  # 右端是闭合的
# # %%
# rr.index
# rr.plot(kind='barh')

# import pandas as pd
# import seaborn as sns
# sns.countplot(y='', data=rr)
# df.resample('3H', convention='start', label='right').sum()
# rr = rain.
# rr =  rain.resample(time='3H', convention='end').sum()
# rr
# %%
# rr[0,:,:]
# rr.sum(dim='time').max()
# rain.sum(dim='time').max()
# r3 = rr.mean(dim=['south_north', 'west_east'])
# r3.plot()





# %%

# %%


# color_list = ['red', 'blue', 'black']
# model_list = ['all', 'grid', 'con']

# x = np.arange(len(model_list))
# width = 0.1
# j = 0 
# for i in range(len(model_list)):
#     rects = ax.bar(x+width*i, y_list[i], width, label=model_list[i], color=color_list[i])
