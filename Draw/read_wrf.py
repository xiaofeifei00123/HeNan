#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取河南暴雨的数据
-----------------------------------------
Time             :2021/09/09 10:53:08
Author           :Forxd
Version          :1.0
'''


# %%
import xarray as xr
import os
import xesmf as xe
import numpy as np
import netCDF4 as nc
import wrf


# from metpy.units import units
# from metpy.calc import specific_humidity_from_dewpoint
# from metpy.calc import mixing_ratio_from_specific_humidity
# from metpy.calc import virtual_potential_temperature
# from metpy.calc import potential_temperature
# from metpy.calc import relative_humidity_from_dewpoint
# from metpy.calc import dewpoint_from_relative_humidity








path = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912/'
fl_list = os.popen('ls {}/wrfout*'.format(path))  # 打开一个管道
fl_list = fl_list.read().split()

fl_list = fl_list[0:1]
# %%
# fl = fl_list[0]
# ds = xr.open_dataset(fl)
# data_nc = nc.Dataset(fl)

# rain = ds['RAINNC']
# rh = wrf.getvar(data_nc, 'rh')

# %%
# u
# rh
# rain
# rh.Time.Time.values
# rh.set_dims(dim=rh.Time.Time.values)
# rh.assign_coords('Time':)
# rh.swap_dims({'time':'Time'})
# a = []
# a.append(rh.Time.values)
# a
# rh.Time
# rh2 = rh.expand_dims(dim='Time')
# rain
# ds = xr.Dataset()
# ds['rain'] = rain
# ds['rh'] = rh2
# ds
# ds['U']
# u = wrf.getvar(data_nc, 'ua')
# u


# %%
dds_list = []
r = 0   # 第0个时次就算降水为0
for fl in fl_list:
    
    ds = xr.open_dataset(fl)
    data_nc = nc.Dataset(fl_list[0])
    print(ds.Time.XTIME.values)
    # print(ds.Time.values)
    dds = xr.Dataset()
    dds['RAINNC'] = ds['RAINNC']-r
    r = ds['RAINNC']

    for var in ['ua', 'va', 'td', 'temp', 'theta_e', 'pressure', 'height_agl']:
        dds[var] = wrf.getvar(data_nc, var, squeeze=False)
        # dds[var] = da.expand_dims(dim='Time')

    dds_list.append(dds)
dds_concate = xr.concat(dds_list, dim='time')
dds_concate

# %%
# dds_concate['ua']
u = wrf.getvar(data_nc, 'ua', squeeze=False)
v = wrf.getvar(data_nc, 'va', squeeze=False)
p = wrf.getvar(data_nc, 'pressure', squeeze=False)
# dda
wrf.interplevel(u, p, [500,200], squeeze=False)
# ds.var



# %%

# def caculate_diagnostic(self, ds):
#     """计算比湿，位温等诊断变量
#     根据td或者rh计算q,theta_v
#     返回比湿q, 虚位温theta_v, 相对湿度rh

#     Args:
#         ds (Dataset): 包含有temp ,td的多维数据
#         这里传入Dataset合理一点
#     """
#     pass        
#     ## 获得温度和露点温度
#     dims_origin = ds['temp'].dims  # 这是一个tuple, 初始维度顺序
#     ds = ds.transpose(*(...,'pressure'))

#     var_list = ds.data_vars
#     t = ds['temp']

#     ## 转换单位
#     pressure = units.Quantity(t.pressure.values, "hPa")

#     ## 针对给的rh 或是td做不同的计算
#     ## 需要确定t的单位
#     if 'td' in var_list:
#         """探空资料的温度单位大多是degC"""
#         td = ds['td']
#         dew_point = units.Quantity(td.values, "degC")
#         temperature = units.Quantity(t.values, "degC")
#     elif 'rh' in var_list:
#         """FNL资料的单位是K"""
#         # rh = da.sel(variable='rh')
#         rh = ds['rh']
#         rh = units.Quantity(rh.values, "%")
#         temperature = units.Quantity(t.values, "K")
#         dew_point = dewpoint_from_relative_humidity(temperature, rh)
#     else:
#         print("输入的DataArray中必须要有rh或者td中的一个")
    
#     ## 记录维度坐标
#     # time_coord = t.time.values
#     # pressure_coord = t.pressure.values

#     ## 计算诊断变量
#     q = specific_humidity_from_dewpoint(pressure, dew_point)
#     w = mixing_ratio_from_specific_humidity(q)
#     theta_v = virtual_potential_temperature(pressure, temperature, w)

#     if 'td' in var_list:
#         rh = relative_humidity_from_dewpoint(temperature, dew_point)
#         var_name_list = ['q', 'rh', 'theta_v']
#         var_data_list = [q, rh, theta_v]
#     elif 'rh' in var_list:
#         pass
#         var_name_list = ['q', 'td', 'theta_v']
#         var_data_list = [q, dew_point, theta_v]

#     ## 融合各物理量为一个DataArray

#     ds_return = xr.Dataset()

#     for var_name, var_data in zip(var_name_list, var_data_list):
#         pass
#         ## 为了去除单位
#         dda = xr.DataArray(
#             var_data, 
#             # coords=[time_coord,pressure_coord],
#             coords = t.coords,
#             dims=t.dims)
#             # dims=['time', 'pressure'])

#         ds_return[var_name] = xr.DataArray(
#             dda.values, 
#             # coords=[time_coord,pressure_coord],
#             coords=t.coords,
#             dims=t.dims)
#     # da_return  = ds_return.to_array()

#     ## 转换维度顺序        
#     # da_return = da_return.transpose(*dims_origin)
#     ds_return = ds_return.transpose(*dims_origin)
#     return ds_return






# %%
# dds_concate['RAINNC'].max()

# %%
# dds_concate.to_netcdf('/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1912.nc')

# %%
aa = wrf.getvar('pressure')
aa
# ds['RAINNC']
# dds['ua']
# ds['pr']

# %%
## 插值
# ds_rain = dds_concate['RAINNC'].to_dataset()
# ds_rain = ds_rain.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'})
# ds_in = ds_rain.isel(Time=0)
# ## 插值目标区域
# ds_regrid = xe.util.grid_2d(108.975, 116, 0.05, 31.975, 37, 0.05)
# ## 插值器构建
# regridder = xe.Regridder(ds_in, ds_regrid, 'bilinear')  # 好像是创建了一个掩膜一样
# ## 插值 
# ds_out = regridder(ds_rain)  # 插值完成的数据依然是二维的，需要重构成一维的
# # dp_out
# # %%
# # ds_out['RAINNC'].max()

# # %%
# ### 重新构建经纬度坐标
# lat = ds_out.lat.sel(x=0).values
# lon = ds_out.lon.sel(y=0).values
# da = ds_out['RAINNC']

# # da.dims
# ds_rr = xr.Dataset(
#     {'RAINNC':(['time','lat', 'lon'], da.values)},
#     coords={
#         'time':da.time.values,
#         'lat':lat,
#         'lon':lon,
#     },
# )
# # %%
# # rain = ds_rr['RAINNC']

# # # rain.mean(dim=['lat', 'lon']).max()
# # # rain.sel(lat=)
# # # rain.max()




# lat = np.arange(24.875, 45.125+0.25,0.25)
# lon = np.arange(69.875, 105.125+0.25,0.25)



# ds_input = xr.Dataset(
#     {'RAINNC':(['time', 'lat', 'lon'], ds_rain.values)}
# )

# ds_input


# ds_out = xe.util.grid_2d(77.875, 105, 0.25, 25.875, 38, 0.25)
# ds_out = xe.util.grid_2d(109.875, 116, 0.25, 31.875, 37, 0.25)
# ds_out = xe.util.grid_2d(108.975, 116, 0.05, 31.975, 37, 0.05)
# regridder = xe.Regridder(ds_rain, ds_out, 'bilinear')  # 好像是创建了一个掩膜一样
# ds_out





