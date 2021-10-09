# %%
# import metdig
# import datetime
# from metdig.io.cassandra import get_obs_stations
import pandas as pd
import numpy as np
import xesmf as xe
from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units
# %%
pressure = units.Quantity(200, "hPa")
dew_point = units.Quantity(20.2, "degC")
aa = specific_humidity_from_dewpoint(pressure, dew_point)
aa

# %%
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
for i in lat:
    print(i)