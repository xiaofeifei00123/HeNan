# %%
import xarray as xr
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Rain/rain_all_station.nc'
ds = xr.open_dataset(flnm)
ds
# %%
ds['OBS']