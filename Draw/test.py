# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_900m/rain.nc'
da = xr.open_dataarray(flnm)
da
tt = slice('2021-07-20 00', '2021-07-20 12')
rain = da.sel(time=tt).sum(dim='time')
rain
# %%
# rain.max()
# rain.lat
# rain.lon
latlon = np.meshgrid(rain.lon, rain.lat)
# latlon[0,:]
# latlon.shape()
# %%
# latlon[0]
# rain.lon.shape
# lon = latlon[0][0,:]
# lon
latlon[0]
# np.meshgrid([0,1,2], [1, 2, 3])


# %%
rain_shape = rain_select.shape
## 构建二维的经纬度
lat_one = np.ones(rain_shape[1]).reshape(1, rain_shape[1])
lat = rain_select['lat'].values.reshape(rain_shape[0], 1)
lat2d = lat * lat_one
lon_one = np.ones(rain_shape[0]).reshape(rain_shape[0], 1)
lon = rain_select['lon'].values.reshape(1, rain_shape[1])
lon2d = lon_one * lon
## 这里用两个sum,是因为每次求sum只可以求一行的
lat_return = sum(sum(rain * lat2d)) / sum(sum(rain))
lon_return = sum(sum(rain * lon2d)) / sum(sum(rain))

dict_return = {
    'rain': rain_select,
    'rain_threshold': rain_threshold,
    'lat': lat_return,
    'lon': lon_return
}







