# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
ds = xr.open_dataset(flnm)

# %%
t0 = slice('2021-07-20 00', '2021-07-20 03')
t1 = slice('2021-07-20 00', '2021-07-20 06')
t2 = slice('2021-07-20 00', '2021-07-20 09')
t3 = slice('2021-07-20 00', '2021-07-20 12')
t_list = [t0, t1, t2, t3]

## 最大降水站点的累计降水
da_list = []
for t in t_list:
    ds1 = ds.sel(time=t)
    rain = ds1.max(dim='sta').sum(dim='time')
    rain.data_vars
    da = rain.to_array()
    da_list.append(da)

daa = xr.concat(da_list, pd.Index(['3', '6', '9', '12'], name='sum_time'))
daa
# %%
# da_list[2]
# rain
# y = rain.to_array()
# y
# %%
# daa
# da_list[0]
# daa
# %%


# x_label = daa.sum_time.values.astype('str')
x_label = ['3h', '6h', '9h', '12h']
x = np.arange(len(x_label))
rain_ds = daa.to_dataset(dim='variable')
rain_ds
y1 = rain_ds['1900_90m']
y2 = rain_ds['1912_90m']
y3 = rain_ds['1912_900m']
y4 = rain_ds['OBS']
# y1
# daa


# %%
width = 0.2
fig = plt.figure(figsize=(12, 12))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.bar(x-width*2, y1.values, width, label='1900_90m')
ax.bar(x-width*1, y2.values, width, label='1912_90m')
ax.bar(x, y3.values, width, label='1912_900m')
ax.bar(x+width*1, y4.values, width, label='OBS')
ax.set_xticks(x)
ax.set_xticklabels(x_label)
ax.legend()
