# %%
import xarray as xr
import os
from draw_rain_distribution_24h import Draw

# %%


flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
da = xr.open_dataarray(flnm)
rain = da.sel(time=slice('2021-07-20 01', '2021-07-21 00')).sum(dim='time')
rain
# %%

dr = Draw()
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
picture_dic = {'date':'2021-07 20/00--21/00', 'type':'aa', 'initial_time':''}
dr.draw_single(rain, picture_dic)







# %%
