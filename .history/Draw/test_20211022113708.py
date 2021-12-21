# %%
import xarray as xr

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_rain_latlon.nc'
ds = xr.open_dataset(flnm)
ds
# %%
ds.isel(time=4).max()