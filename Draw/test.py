# %%
import xarray as xr
# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
ds = xr.open_dataset(flnm)
ds
# %%
ds['1912_90m']