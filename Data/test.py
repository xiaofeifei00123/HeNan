# %%
import xarray as xr
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt

# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
ds = xr.open_dataarray(flnm)
# ds
da = ds.isel(time=0)
da
# da = ds.sel(time='2021-07-20 00')
# da = ds.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
# # da
# ds.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
# %%
south_north = da.lat.values
west_east = da.lon.values
# da = da.a
da = da.assign_coords({'south_north':('lat',south_north), 'west_east':('lon',west_east)})
da = da.swap_dims({'lat':'south_north', 'lon':'west_east'})
import sys
sys.path.append('/mnt/zfm_18T/fengxiang/HeNan/Draw')
import draw_rain_distribution_24h as d2
dr = d2.Draw()
picture_dic = {'date':'2021-07 20/00--21/00', 'type':'aa', 'initial_time':''}
dr.draw_single(da, picture_dic)

# import baobao.quick_draw as bq
# bq.quick_contourf(da)


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


