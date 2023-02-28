
# %%
import netCDF4 as nc
from netCDF4 import Dataset
import wrf
from wrf import getvar
# from baobao.caculate import caculate_div3d
from baobao.caculate import get_div_wrfout
import xarray as xr
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
import pycwt
from netCDF4 import Dataset
from wrf import getvar, interpline, CoordPair, vertcross
from geopy.distance import distance  # 根据经纬度计算两点距离
# %%
# coord_pairs
def latlon2distance(da2):
    """将剖面数据的经纬度横坐标变为距离坐标

    Args:
        da2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
    # ds = xr.open_dataset(flnm)
    # ds1 = ds.sel(time='2021-07-20 08')
    # da = ds1['wa_cross']
    # da1 = da.interpolate_na(dim='vertical', method='linear',  fill_value="extrapolate")
    # da2 = da1.sel(vertical=2000, method='nearest')
    dd = da2.xy_loc
    def str_latlon(string):
        # d1 = dd.values[0]
        lat = float(string.split(',')[0])
        lon = float(string.split(',')[1])
        return lat, lon

    d2 = dd.values
    lat_list = []
    lon_list = []
    for i in d2:
        # print(i)
        lat, lon = str_latlon(i)
        lat_list.append(lat)
        lon_list.append(lon)

    dis_list = [0]
    di = 0
    for i in range(len(lat_list)-1):
        # print(i)
        lat1 = lat_list[i]
        lon1 = lon_list[i]
        loc1 = (lat1, lon1)
        lat2 = lat_list[i+1]
        lon2 = lon_list[i+1]
        loc2 = (lat2, lon2)
        dist = distance(loc1,loc2).km
        di = di+dist
        dis_list.append(di)
    dis_list
    dis_array = (np.array(dis_list)).round(1)
    dis_array
    da2 = da2.assign_coords({'distance':('line_idx',dis_array)})
    da3 = da2.swap_dims({'line_idx':'distance'})
    da3
    return da3

# %%

def get_rain_line(model):
    cross_start = CoordPair(lat=33, lon=111.5)
    cross_end = CoordPair(lat=36, lon=114.3)
    # flnm1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/wrfout_d03_2021-07-20_00:00:00'
    # flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/wrfout_d03_2021-07-20_01:00:00'
    flnm1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'+model+'/wrfout/wrfout_d03_2021-07-20_00:00:00'
    flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'+model+'/wrfout/wrfout_d03_2021-07-20_01:00:00'
    wrfnc1 = nc.Dataset(flnm1)
    wrfnc2 = nc.Dataset(flnm2)

    da1 = getvar(wrfnc1, 'RAINNC')+getvar(wrfnc1, 'RAINSH')+getvar(wrfnc1, 'RAINC')
    da2 = getvar(wrfnc2, 'RAINNC')+getvar(wrfnc2, 'RAINSH')+getvar(wrfnc2, 'RAINC')
    da = da2-da1
    var_vcross = interpline(da,wrfnc1, start_point=cross_start, end_point=cross_end, latlon=True)
    coord_pairs = var_vcross.coords["xy_loc"].values
    x_labels = [pair.latlon_str(fmt="{:.5f}, {:.5f}")
                for pair in coord_pairs]
    x_labels
    var = var_vcross.assign_coords({'xy_loc':('line_idx',x_labels)})
    var2 = latlon2distance(var)
    return var2



var1 = get_rain_line('GWD3')
var2 = get_rain_line('CTRL')
var3 = get_rain_line('FD')
var4 = get_rain_line('SS')


# %%
cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 4*cm), dpi=300)
ax = fig.add_axes([0.1,0.1, 0.8, 0.8])
# var2
var_list = [var1, var2, var3, var4]
# var_list = [var1]# var3, var4]
label_list = ['GWD3', 'CTRL', 'FD', 'SS']
color_list = ['red', 'black', 'blue', 'orange']
i = 0
for var in var_list:
    x = var.distance.values
    y = var.values
    ax.plot(x,y, label=label_list[i], color=color_list[i], )
    i+=1
ax.legend(edgecolor='white', loc='upper left', fontsize=5)
# ax.set_xlim(237, 238)