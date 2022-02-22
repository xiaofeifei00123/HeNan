#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算水平散度的剖面数据
-----------------------------------------
Time             :2022/01/12 17:01:24
Author           :Forxd
Version          :1.0
'''
# %%

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cnmaps import get_map, draw_map
from baobao.interp import regrid_xesmf
from baobao.interp import rain_station2grid

# %%
flnm = '/home/fengxiang/HeNan/Data/OBS/rain_ec.nc'
da2 = xr.open_dataarray(flnm)
da2 = da2.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
da2 = da2.sel(lat=slice(32,36.5)).sel(lon=slice(110.5,116))
da2 = da2.sum(dim='time')


# %%
# da1
area = {
    'lon1':110.5,
    'lon2':116,
    'lat1':32,
    'lat2':36.5,
    'interval':0.125,
}
da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc')
da1 = rain_station2grid(da, area)
da1 = da1.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
da1 = da1.sum(dim='time')
# %%
# ddc.lat
# da2.lat
# dd
da2-da1



# da1-da_rain
# type(da_rain.dims)
# da_rain.plot()
# dd = da_rain.sum(dim='time').plot()
# dd
# da
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
# flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain.nc'
# flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'

# %%
# flnm3 = 
flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
da3 = xr.open_dataarray(flnm3)
da3 = da3.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
da3 = da3.sum(dim='time') 
da0 = xr.open_dataarray(flnm0)
da0 = da0.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
da0 = da0.sum(dim='time') 
da0
# %%
# da0
# da3-da0
# da3.lat.values
da = da3.values-da0.values
# da = xr.DataArray
da
# %%
dda = xr.DataArray(
    da,
    coords={
        'lat':da3.lat.values,
        'lon':da3.lon.values,
    },
    dims=['lat', 'lon']
    )
dda




# da0
# k
# da0










# %%
import xarray as xr
import numpy as np
import pandas as pd
from wrf import getvar, CoordPair, vertcross, get_cartopy
import wrf
from netCDF4 import Dataset
from multiprocessing import Pool
from baobao.caculate import caculate_div
from metpy.units import units  # units是units模块中的一个变量名(函数名？类名？)
from metpy import calc as ca  # calc是一个文件夹
import matplotlib.pyplot as plt

# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross.nc'
ds = xr.open_dataset(flnm)
da = ds['div_cross']
da = da.isel(time=0)
# %%

# self.cross_start = CoordPair(lat=33, lon=111)
# self.cross_end = CoordPair(lat=35.5, lon=114.5)

ldic = {
    'lat1':33,
    'lon1':111,
    'lat2':35.5,
    'lon2':114.5,
}
dy = (ldic['lat2']-ldic['lat1'])
dx = (ldic['lon2']-ldic['lon1'])
angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度

# np.cos(angle)
# deg = 180.0/np.pi # 角度和弧度之间的转换
# rad = np.pi/180.0


# angle
# np.arctan(dy/dx)




# %%
def draw_contour2(ax,da):
    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    cs = ax.contour(xs, ys, da.values*10**4, colors='red', zorder=3)
    ax.clabel(cs, inline=True, fontsize=10)

    

fig = plt.figure(figsize=[10,8])
ax = fig.add_axes([0.1,0.1,0.8,0.8])



# %%

def concat_dxy(var, index):
    """将二维数据在垂直方向上累加
    变成三维的

    Args:
        var ([type]): 二维变量
        index ([type]): 垂直方向的索引

    Returns:
        var_concat: 累加到一块的数据
    """
    var_list = []
    # index = u.bottom_top.values
    for i in index:
        # print(i)
        aa = var.magnitude  # 转为numpy
        bb = xr.DataArray(aa, dims=['south_north', 'west_east']) # 转为DataArray
        var_list.append(bb)
    var_concat = xr.concat(var_list, dim=pd.Index(index, name='bottom_top'))
    return var_concat

def caculate_div3d(u, v, lon, lat):
    """求wrfout数据中三维的u,v数据对应的散度
    就先原始的wrfout数据吧

    Args:
        u ([type]): 三维
        v ([type]): 三维
        lon ([type]): 二维
        lat ([type]): 二维
    """
    pass
    # u =  getvar(ncfile, 'ua')
    # v =  getvar(ncfile, 'va')
    # lon = u.XLONG
    # lat = u.XLAT
    u = u*units('m/s')
    v = v*units('m/s')
    dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
    ## 重组dx和dy, 其实就是把dx和dy的垂直维度加上，虽然每个垂直层上数据一样, 这是由于metpy计算时的问题导致的
    index = u.bottom_top.values
    ddx = concat_dxy(dx, index)
    ddy = concat_dxy(dy, index)
    dddx = ddx.values*units('m')
    dddy = ddy.values*units('m')
    ### 因为这个函数的问题，所以dx必须是和u维度相对应的
    div = ca.divergence(u=u, v=v, dx=dddx, dy=dddy)
    div
    return div
    # div.min()*10**3*10

# if __name__ == '__main__':
#     # main()
    
wrf_file = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d01_2021-07-19_18:00:00'
ncfile = Dataset(wrf_file)
u =  getvar(ncfile, 'ua')
v =  getvar(ncfile, 'va')
lon = u.XLONG
lat = u.XLAT
div = caculate_div3d(u,v, lon, lat)
div = div.rename('div')
div
# %%
# xr.merge([u,div])
u.projection
div.attrs = u.attrs
# %%
div.to_netcdf('aa.nc')
