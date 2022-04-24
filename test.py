# %%
import xarray as xr
# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_all.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/wrfout_d06_2021-07-19_12:00:00'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/wrfout_d05_2021-07-20_00:00:00'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/LES/wrfout_d01_2021-07-19_12:00:00'
ds = xr.open_dataset(flnm)
da = ds.sel(station='zhengzhou').sel(time='2021-07-20 00').dropna(dim='pressure')

# %%
# da['q'].max()
df = da.to_dataframe()
ddf = df[['height', 'theta', 'q','u', 'v']]
ddf['q'] = ddf['q']*10**3
ddf
ddf = ddf.round(2)
ddf = ddf.sort_values(by='height', ascending=True)
ddf.to_csv('input_sounding',index=False, sep=' ')



# %%
import pandas as pd
flnm = './20200101000000.txt'
df = pd.read_csv(flnm, skiprows=1, sep='\s')
name = ['Station_Id_C', 'WIN_D_Avg_2mi']
df1 = df[name]
df1['val'] = df['WIN_D_Avg_2mi'].apply(lambda x: 1 if x==999999 else 0)
df1


flnm = './20200101000000.txt'
df = pd.read_csv(flnm, skiprows=1, sep='\s')
name = ['Station_Id_C', 'WIN_D_Avg_2mi']
df2 = df[name]
df2['val'] = df2['WIN_D_Avg_2mi'].apply(lambda x: 1 if x==999999 else 0)
df2
df1+df2










# %%
import xarray as xr
flnm = '/home/fengxiang/SOIL/GLDAS/2.0/GLDAS_CLSM025_D.A19480102.020.nc4'
ds = xr.open_dataset(flnm)
ds



# %%
import pandas as pd
flnm = './subset_GLDAS_NOAH025_3H_2.0_20220312_065359.txt'
aa = pd.read_table(flnm, sep='\s')
bb = aa.squeeze()
# %%
bb[1]
import re
# %%
# rex = '\bGLDAS[\w]*.A\d*.\d*.\d*.nc4'
# rex = r'\bGLDAS[\w]+.[\S]*.nc4.SUB.nc4'
rex = r'\bGLDAS[\w]*.[\S]*020.nc4'
pattern = re.compile(rex)
print(bb[1])
# pattern
# m = pattern.match(bb[1])
# m
# print(m)
# %%
# cc = re.match(rex, bb[1]).group()
m = re.search(rex, bb[1])
# print(m)
# type(m)
m.group()
# cc







# %%
# import iris
import netCDF4 as  nc
from pip import main
import wrf
from cv2 import transform
import xarray as xr
from cfgrib.xarray_to_grib import to_grib  # cfgrib库写grib文件还不太完善
import numpy as np
import pandas as pd
import xarray as xr
import eccodes as ec
# from iris_grib import save_grib2
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cnmaps import get_map, draw_map
from baobao.map import Map
import netCDF4 as nc
import os
import iris
import iris_grib
import gribapi
from iris_grib import save_grib2
import cfgrib
# import cf2cfm
from cfgrib.xarray_to_grib import to_grib



# %%

# flnm = '/home/fengxiang/TEST/transform_grib/old/era5.pl.20210719.grib'
# ds = xr.open_dataset(flnm, engine='cfgrib')
# # ds
# to_grib(ds, 'ttt.grib', grib_keys={'edition': 1})
# %%
def get_mask(ds):
    st = {
            'lat': 24.0,
            'lon': 131.1,
            'abbreviation':'烟花'
        }
    # st = station['YanHua']
    area = {
        'lon1':st['lon']-2,
        'lon2':st['lon']+2,
        'lat1':st['lat']-2,
        'lat2':st['lat']+2,
        'interval':0.05,
    }
    ## 制作二维的掩膜
    mask = (
        (ds.coords['latitude']>area['lat1'])
        &(ds.coords['latitude']<area['lat2'])
        &(ds.coords['longitude']<area['lon2'])
        &(ds.coords['longitude']>area['lon1'])
    )
    return mask




# %%




def pl(flnm):
    # flnm = '/home/fengxiang/TEST/transform_grib/old/era5.pl.20210719.grib'
    ds = xr.open_dataset(flnm, engine='cfgrib')
    mask = get_mask(ds)
    radd = (xr.where(mask, 1.5, 1)).astype('float32')
    rminus = (xr.where(mask, 0.9, 1)).astype('float32')
    # bb = aa.astype('float32')

    var_list_pl_upar = ['u','v']
    var_list_pl_minus = ['z']
    for var in var_list_pl_minus:
        attr = ds[var].attrs
        ds[var] = ds[var]*rminus
        ds[var].attrs = attr
        print(ds[var].dtype)

    for var in var_list_pl_upar:
        attr = ds[var].attrs
        ds[var] = ds[var]*radd
        ds[var].attrs = attr
        print(ds[var].dtype)

    # file_grib_to = 'test.grib'
    file_grib_to = flnm.replace('old','newnew')
    to_grib(ds, file_grib_to, grib_keys={'edition': 1})


def sl(flnm):
    # flnm = '/home/fengxiang/TEST/transform_grib/old/era5.sl.20210719.grib'
    ds = xr.open_dataset(flnm, engine='cfgrib')
    mask = get_mask(ds)
    radd = (xr.where(mask, 1.5, 1)).astype('float32')
    rminus = (xr.where(mask, 0.9, 1)).astype('float32')
    # bb = aa.astype('float32')

    # var_list_pl_upar = ['u','v']
    # var_list_pl_minus = ['z']
    var_list_sl_upar = ['u10','v10']
    var_list_sl_minus = ['sp']
    for var in var_list_sl_minus:
        attr = ds[var].attrs
        ds[var] = ds[var]*rminus
        ds[var].attrs = attr
        print(ds[var].dtype)

    for var in var_list_sl_upar:
        attr = ds[var].attrs
        ds[var] = ds[var]*radd
        ds[var].attrs = attr
        print(ds[var].dtype)

    # file_grib_to = 'test.grib'

    file_grib_to = flnm.replace('old','newnew')
    to_grib(ds, file_grib_to, grib_keys={'edition': 1})
# sl()


# %%
path='/home/fengxiang/TEST/transform_grib/old'
fl_list1 = os.popen('ls {}/*sl*.grib'.format(path))  # 打开一个管道
fl_list1 = fl_list1.read().split()
fl_list2 = os.popen('ls {}/*pl*.grib'.format(path))  # 打开一个管道
fl_list2 = fl_list2.read().split()

# for fl in fl_list1:
#     sl(fl)

for fl in fl_list2:
    pl(fl)


# var_list_pl_upar = ['u','v']
# var_list_pl_minus = ['z']

# var_list_sl_upar = ['u10','v10']
# var_list_sl_minus = ['sp']
# pl()

# attr = ds['u10'].attrs
# ds['u10'] = ds['u10']*aa  ## 改变其数据大小
# ds['u10'].attrs = attr

# var_list_pl = ['u','v']




# file_grib_to = 'test.grib'
# to_grib(ds, file_grib_to, grib_keys={'edition': 1})
# index
# dac.loc[:,index]
# aa.values
# %%
# xr.open_dataset(file_grib_to,engine='cfgrib')
# ds.step
# np.where(bb<70&bb>10,1,0)
# ddd*ds
# ds
# 改变区域的值
# ds['u10'].shape
# path='/home/fengxiang/TEST/transform_grib/met/old/met_em.d01.2021-07-19_00:00:00.nc'
# ds = xr.open_dataset(path)
# ds




# %%

# ## 构建插值系数
# # np.zeros(u.shape)
def create_coff2d(da, coff):
    """构建二维的插值系数矩阵
    da, DataArry, 水平2维
    coff:关键区域变化的系数
    """
    shape = da.shape
    dac = np.ones(shape)  # 构建二维系数矩阵
    ## 系数矩阵变成和数据一样的维度
    dac = xr.DataArray(
        dac,
        coords=da.coords,
        dims=da.dims,
    )
    ## 设置作用范围
    station = {
            'YanHua': {
                'lat': 24.0,
                'lon': 131.1,
                'abbreviation':'烟花'
            },
            'ChaPaka': {
                'lat': 21.1,
                'lon': 112.8,
                'abbreviation':'查帕卡'
            },
        }
    st = station['YanHua']
    area = {
        'lon1':st['lon']-2,
        'lon2':st['lon']+2,
        'lat1':st['lat']-2,
        'lat2':st['lat']+2,
        'interval':0.05,
    }
    ## 修改范围内的值
    dac = dac.rename({'latitude':'lat', 'longitude':'lon'})
    index = ((dac.lat<=area['lat2']) & (dac.lat>=area['lat1']) & (dac.lon>=area['lon1']) & (dac.lon<=area['lon2']))
    cc = dac.values
    cc[index.values] = coff   # 得到二维的系数矩阵

    cc1 = xr.DataArray(
        cc,
        coords=da.coords,
        dims = da.dims,
    )
    return cc1

def upgrade_dimesion_2(cc1, ds2):
    """升级两次维度

    Args:
        cc (_type_): _description_
        ds2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    ## 将二维升级到3维
    level = ds2['u'].isobaricInhPa.values
    dd_list1 = []
    for i in range(len(level)):
        # print(i, ll[i])
        dd_list1.append(cc1)
    cc2 = xr.concat(dd_list1, pd.Index(level, name='isobaricInhPa'))

    ## 将三维升级到4维度
    time_range = ds2['u'].time.values
    dd_list2 = []
    for i in range(len(time_range)):
        # print(i, ll[i])
        dd_list2.append(cc2)
    cc3 = xr.concat(dd_list2, pd.Index(time_range, name='time'))
    cc3
    return cc3

def upgrade_dimesion_1(cc, ds2):
    """升级一次维度, 没有level

    Args:
        cc (_type_): _description_
        ds2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    ## 将二维升级到3维

    # cc1 = xr.DataArray(
    #     cc,
    #     coords=da.coords,
    #     dims = da.dims,
    # )
    
    time_range = ds2['sp'].time.values
    dd_list2 = []
    for i in range(len(time_range)):
        # print(i, ll[i])
        dd_list2.append(cc)
    cc2 = xr.concat(dd_list2, pd.Index(time_range, name='time'))
    cc2
    return cc2

# # %%
# ## 改变高空数据
def change_pl():
    flnm = '/home/fengxiang/TEST/transform_grib/old/era5.pl.20210719.grib'
    ds = xr.open_dataset(flnm, engine='cfgrib')
    da = ds['u'].sel(isobaricInhPa=500).isel(time=0)

    u_attr = ds['u'].attrs
    v_attr = ds['v'].attrs
    z_attr = ds['z'].attrs
    

    cc1 = create_coff2d(da, 1.5)
    cc1 = upgrade_dimesion_2(cc1,ds)
    cc2 = create_coff2d(da, 0.98)
    cc2 = upgrade_dimesion_2(cc2,ds)
    ds['u'] = ds['u']*cc1
    ds['v'] = ds['v']*cc1
    ds['z'] = ds['z']*cc2

    ds['u'].attrs = u_attr
    ds['v'].attrs = v_attr
    ds['z'].attrs = z_attr

    file_grib_to = '/home/fengxiang/TEST/transform_grib/new/era5.pl.20210719.grib'
    # ds.to_netcdf(file_nc)
    to_grib(ds, file_grib_to, grib_keys={'edition': 1})


# # %%
# ## 改变地面数据
def change_sl():
    flnm = '/home/fengxiang/TEST/transform_grib/old/era5.sl.20210719.grib'
    ds = xr.open_dataset(flnm, engine='cfgrib')
    da = ds['sp'].isel(time=0)
    u10_attr = ds['u10'].attrs
    v10_attr = ds['v10'].attrs
    sp_attr = ds['sp'].attrs
    cc1 = create_coff2d(da, 1.5)
    cc1 = upgrade_dimesion_1(cc1,ds)
    cc2 = create_coff2d(da, 0.98)
    cc2 = upgrade_dimesion_1(cc2,ds)
    ds['sp'] = ds['sp']*cc2
    ds['u10'] = ds['u10']*cc1
    ds['v10'] = ds['v10']*cc1
    ds['sp'].attrs = sp_attr
    ds['u10'].attrs = u10_attr
    ds['v10'].attrs = v10_attr
    file_grib_to = '/home/fengxiang/TEST/transform_grib/new/era5.sl.20210719.grib'
    to_grib(ds, file_grib_to, grib_keys={'edition': 1})
    # ds.to_netcdf(file_nc)
# change_sl()
# change_pl()
# %%