# %%
# import metdig
# import datetime
# from metdig.io.cassandra import get_obs_stations
import pandas as pd
import numpy as np
import xesmf as xe
from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units
import xarray as xr
import datetime
import meteva.base as meb
# %%
# flnm = ''
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GFS/YSU_GFS_rain.nc'
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
# flnm_ec = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_ec.nc'
# ds = xr.open_dataset(flnm_ec)
ds1 = xr.open_dataset(flnm)
ds1
# ds2 = xr.open_dataset(flnm_ec)
# ds1.lat
# ds2.lat
# ds1 = ds1.sel(lat=slice(32,37), lon=slice(110, 116))
# ds2 = ds2.sel(lat=slice(32,37), lon=slice(110, 116))
# ds2.lat
# ds1.lat
# ll = ds.lat.values
# ds1 = ds.sel(lat=slice(32,37), lon=slice(110, 116))
# ds.lat
# pressure = units.Quantity(200, "hPa")
# dew_point = units.Quantity(20.2, "degC")
# aa = specific_humidity_from_dewpoint(pressure, dew_point)
# aa
# %%
# ds.lat
path='/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_station/20210720170000.000'
station = meb.read_stadata_from_micaps3(path)
meb.tool.plot_tools.scatter_sta(station, value_column=0,map_extend = [110,116,32,37],
                             cmap = meb.cmaps.rain_1h)
# station

# %%
grid1 = meb.grid([110-1,116+1,0.125],[32-1,37+1,0.125],)
grd2 = meb.interp_sg_idw(station, grid1)
# ds2 = grd2.sel()
ds2 = grd2.sel(lat=slice(32,37), lon=slice(110, 116))
# grd2.lat
# %%
# ds1
# ds2
# grd2.sel(lat=32)
# ds2.lat
# ds1.lat
# ds1
# ds2.squeeze()
dds1 = ds1.sel(time='2021-07-20 1700')
dds2 = ds2.squeeze()

# dds2-dds1
# dds2.var
# dds2
da1 = dds1.to_array().squeeze()
# da2 = 
ds2
# dds2
# ds2.name = 'data'
# ds2
# ds1
# ds1
# dds2-dds1
# ds1.time
# ds2.lat
# ds1.lat
# %%
# grd2.lat
# t1 = pd.date_range('2021-07-18 00', '2021-07-21 00', freq='24H')
# t2 = pd.date_range('2021-07-18 00', '2021-07-21 00', freq='6H')
# t3 = pd.date_range('2021-07-18 18', '2021-07-21 00', freq='24H')
# t3
# t2.drop(t3)
# t2.sel()
# %%
t2 = pd.date_range('2021-07-18 00', '2021-07-21 00', freq='6H')
ds = xr.Dataset({'time':t2})
d1 = ds.time.sel(time=datetime.time(int(0)))
d2 = ds.time.sel(time=datetime.time(int(6)))
d3 = ds.time.sel(time=datetime.time(int(12)))
xr.concat([d1, d2, d3], dim='time')
# d4 = ds.time.sel(time=datetime.time(int(18)))

# ttt = pd.Timestamp('2021-07-20 18', )
# d5 = ds.drop_sel(time=ttt)


t = pd.DatetimeIndex(['2021-07-18 00',
                      '2021-07-18 06',
                      '2021-07-18 12',
                      '2021-07-19 00',
                      '2021-07-19 06',
                      '2021-07-19 12',
                      '2021-07-20 00',
                      '2021-07-20 06',
                      '2021-07-20 12',
                      ])

# t.append(ttt)

# bh = pd.offsets.BusinessHour(start='00:00', end='12:00', offset=ttt)
# bh = pd.offsets.BusinessHour(start='00:00', end='12:00')
# bh[0]
# bh
# ttt
# t1 = pd.Timestamp('2021-07-20 00')
# t1+bh+bh
# t[0]
# t1+t2
# pd.concat([t1, t2])
# t3 = t2.append(t1)
# t3
# %%
# t3
# pd.merge(t1, t2)

# tt.sel()
# id_a = tt.apply(lambda x:  x.hour==5)
# id_a
# area = {
#     'lon1':107,
#     'lon2':135,
#     'lat1':20,
#     'lat2':40,
#     'interval':0.5,
# }

# ds_regrid = xe.util.grid_2d(area['lon1'], area['lon2'], area['interval'], area['lat1'], area['lat2'], area['interval'])
# ds_regrid
# # lon = ds_regrid.lon_b
# # lat = ds_regrid.lat_b
# # # for i in lon:
# # #     print(i[0].values)
# # # aa = lon[0,:]
# # # for i in aa:
# #     # print(i.values)
# # ds_out = ds_regrid
# # lat = ds_out.lat.sel(x=0).values.round(2)
# # lon = ds_out.lon.sel(y=0).values.round(2)
# # lat

# %%
# for i in lat:
    # print(i)