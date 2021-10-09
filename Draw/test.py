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
# %%
# pressure = units.Quantity(200, "hPa")
# dew_point = units.Quantity(20.2, "degC")
# aa = specific_humidity_from_dewpoint(pressure, dew_point)
# aa

# %%
# t1 = pd.date_range('2021-07-18 00', '2021-07-21 00', freq='24H')
# t2 = pd.date_range('2021-07-18 00', '2021-07-21 00', freq='6H')
# t3 = pd.date_range('2021-07-18 18', '2021-07-21 00', freq='24H')
# t3
# t2.drop(t3)
# t2.sel()
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