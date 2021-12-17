# %%
import xarray as xr
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt

# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd0/wrfout_d03_2021-07-19_21:00:00'
# wrfnc = nc.Dataset(flnm)
# lat = 34.72
# lon = 113.65
# # x,y = wrf.ll_to_xy(wrfnc, 34.71, 113.66)
# x,y = wrf.ll_to_xy(wrfnc, lat, lon)

# ## 以位势高度
# # %%
# hagl = wrf.getvar(wrfnc, 'z', units='m')[:,x,y]  # 先纬度后经度
# hagl
# # %%
# p = wrf.getvar(wrfnc, 'p', units='hPa')[:,x,y]  # 先纬度后经度
# p
# # print(hagl)
# # print(p)


# %%
flnm1 = '/home/fengxiang/HeNan/Data/GWD/sounding_all.nc'
flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_all.nc'
ds1 = xr.open_dataset(flnm1)
ds2 = xr.open_dataset(flnm2)
# ds

# %%
# ds2['u'].sel(station='zhengzhou').isel(time=0).dropna(dim='pressure', how='all')
ds1['height'].sel(station='zhengzhou', model='gwd1').isel(time=0).dropna(dim='pressure', how='all')
# %%
# ds1['q'].sel(station='zhengzhou', model='gwd0').isel(time=0).dropna(dim='pressure', how='all')
# ds3 = ds2.assgin_coords({'model':'OBS'})
# ds3 = ds2.expand_dims(dim='model')
# xr.concat([ds1, ds3], dim='model')
# ds2['height']
# ds2['pressure']
ds2['q'].sel(station='zhengzhou').isel(time=0).dropna(dim='pressure', how='all')
# ds2['pressure']

# %%
fig = plt.figure(figsize=[12,8])
ax = fig.add_axes([0.1,0.1, 0.1, 0.4])
# theta_v = ds2['wind_angle'].sel(station='zhengzhou').sel(time='2021-07-20 00').dropna(dim='pressure', how='all')
# theta_v1 = ds1['wind_angle'].sel(station='zhengzhou', model='gwd1').sel(time='2021-07-20 00').dropna(dim='pressure', how='all')
theta_v = ds2['wind_speed'].sel(station='zhengzhou').sel(time='2021-07-20 00').dropna(dim='pressure', how='all')
theta_v1 = ds1['wind_speed'].sel(station='zhengzhou', model='gwd1').sel(time='2021-07-20 00').dropna(dim='pressure', how='all')
# theta_v.pressure.values[::-1]
ax.plot(theta_v.values, theta_v.pressure)
ax.plot(theta_v1.values, theta_v1.pressure)
ax.invert_yaxis()

da = theta_v1
ax2 = ax.twinx()
ax2.plot(da.values, da.pressure, alpha=0)
ax2.invert_yaxis()

dda = da.swap_dims({'pressure':'height'})
ddda = dda.interp(height=[500, 1000, 1500, 2000, 3000, 5000, 10000, 20000], method='linear', kwargs={'fill_value':'extrapolate'})
# ddda = dda
ax2.set_yticks(ddda.pressure.values)
ax2.set_yticklabels(ddda.height.values.round(1))
# ax2.invert_yaxis()



# ax.set_ylim(1000,700)
# %%
# theta_v1['height']
# ds1
# ds2
# theta_v1


