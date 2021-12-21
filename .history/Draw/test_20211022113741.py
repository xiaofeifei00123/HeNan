# %%
import xarray as xr

# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/YSU_rain_1km_latlon.nc'
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/HighResolution/YSU_rain_1km.nc'
ds = xr.open_dataset(flnm)
ds
# %%
ds.isel(time=4).max()