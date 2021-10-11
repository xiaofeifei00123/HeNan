#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
存放公共的计算函数
插值函数(水平插值，垂直插值)
计算诊断变量
-----------------------------------------
Time             :2021/09/22 11:01:19
Author           :Forxd
Version          :1.0
'''

# %%
# from SWV.draw_map_model import main
from metpy.units import units
from metpy.calc import specific_humidity_from_dewpoint
from metpy.calc import mixing_ratio_from_specific_humidity
from metpy.calc import virtual_potential_temperature
from metpy.calc import potential_temperature
from metpy.calc import relative_humidity_from_dewpoint
from metpy.calc import dewpoint_from_relative_humidity
import metpy.interpolate as interp
from wrf import projection
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np


def caculate_diagnostic(ds):
    """计算比湿，位温等诊断变量
    传入的温度和露点必须是摄氏温度
    根据td或者rh计算q,theta_v
    返回比湿q, 虚位温theta_v, 相对湿度rh

    Args:
        ds (Dataset): 包含有temp ,td的多维数据
        这里传入Dataset合理一点
    """
    pass        
    ## 获得温度和露点温度
    dims_origin = ds['temp'].dims  # 这是一个tuple, 初始维度顺序
    ds = ds.transpose(*(...,'pressure'))

    var_list = ds.data_vars
    t = ds['temp']

    ## 转换单位
    pressure = units.Quantity(t.pressure.values, "hPa")

    ## 针对给的rh 或是td做不同的计算
    ## 需要确定t的单位
    if 'td' in var_list:
        """探空资料的温度单位大多是degC"""
        td = ds['td']
        dew_point = units.Quantity(td.values, "degC")
        temperature = units.Quantity(t.values, "degC")
    elif 'rh' in var_list:
        """FNL资料的单位是K"""
        # rh = da.sel(variable='rh')
        rh = ds['rh']
        rh = units.Quantity(rh.values, "%")
        temperature = units.Quantity(t.values, "degC")
        dew_point = dewpoint_from_relative_humidity(temperature, rh)
    else:
        print("输入的DataArray中必须要有rh或者td中的一个")
    
    ## 记录维度坐标
    # time_coord = t.time.values
    # pressure_coord = t.pressure.values

    ## 计算诊断变量
    q = specific_humidity_from_dewpoint(pressure, dew_point)
    w = mixing_ratio_from_specific_humidity(q)
    theta_v = virtual_potential_temperature(pressure, temperature, w)

    if 'td' in var_list:
        rh = relative_humidity_from_dewpoint(temperature, dew_point)
        var_name_list = ['q', 'rh', 'theta_v']
        var_data_list = [q, rh, theta_v]
    elif 'rh' in var_list:
        pass
        var_name_list = ['q', 'td', 'theta_v']
        var_data_list = [q, dew_point, theta_v]

    ## 融合各物理量为一个DataArray

    ds_return = xr.Dataset()

    for var_name, var_data in zip(var_name_list, var_data_list):
        pass
        ## 为了去除单位
        dda = xr.DataArray(
            var_data, 
            # coords=[time_coord,pressure_coord],
            coords = t.coords,
            dims=t.dims)
            # dims=['time', 'pressure'])

        ds_return[var_name] = xr.DataArray(
            dda.values, 
            # coords=[time_coord,pressure_coord],
            coords=t.coords,
            dims=t.dims)
    ## 转换维度顺序        
    ds_return = ds_return.transpose(*dims_origin)
    return ds_return

def regrid_xesmf(dataset, area):
    """利用xESMF库，将非标准格点的数据，插值到标准格点上去
    Args:
        dataset ([type]): Dataset格式的数据, 多变量，多时次，多层次
    读的是80-102度的数据
        area, 需要插值的网格点范围, 即latlon坐标的经纬度范围
    """
    ## 创建ds_out, 利用函数创建,这个东西相当于掩膜一样
    ds_regrid = xe.util.grid_2d(area['lon1'], area['lon2'], area['interval'], area['lat1'], area['lat2'], area['interval'])
    # ds_regrid = xe.util.grid_2d(110, 116, 0.05, 32, 37, 0.05)
    regridder = xe.Regridder(dataset, ds_regrid, 'bilinear')  # 好像是创建了一个掩膜一样
    ds_out = regridder(dataset)  # 返回插值后的变量

    ### 重新构建经纬度坐标
    lat = ds_out.lat.sel(x=0).values.round(2)
    lon = ds_out.lon.sel(y=0).values.round(2)
    ds_1 = ds_out.drop_vars(['lat', 'lon'])  # 可以删除variable和coords

    ## 设置和dims, x, y相互依存的coords, 即lat要和y的维度一样
    ds2 = ds_1.assign_coords({'lat':('y',lat), 'lon':('x',lon)})
    # ## 将新的lat,lon, coords设为dims
    ds3 = ds2.swap_dims({'y':'lat', 'x':'lon'})
    ## 删除不需要的coords
    # ds_return = ds3.drop_vars(['XTIME'])
    ds_return = ds3
    return ds_return
    

def interp_metpy(sta):
    """
    站点插值到格点
    反距离权重插值

    Args:
        sta (DataFrame): [lon,lat,height]

    Returns:
        [type]: [description]
    """
    # h = sta['temperature']
    # import metpy.interpolate as interp
    h = sta['height']

    lon = sta['lon']
    lat = sta['lat']
    # x,y,z = interp.interpolate_to_grid(lon, lat, h, 'barnes', hres=0.5, minimum_neighbors=2)
    # grid1 = meb.grid([60,150,0.25],[10,60,0.25])
    # x0, x1 = 69.05, 150.1
    # y0, y1 = 0, 55.1 

    x0, x1 = 0, 160
    y0, y1 = 0, 80
    # res = 1 / 32.0
    # res = 1 / 10  # 0.1°
    res = 1
    mx, my = np.meshgrid(np.arange(x0, x1, res),
                            np.arange(y0, y1, res),
                            indexing="ij")
    # mx
    # grd2 = meb.interp_sg_idw(stb, grid1)
    # z = interp.inverse_distance_to_grid(lon, lat, h, mx, my, r=10)
    # x,y,z = interp.remove_nan_observations(lon,lat, h)

    z = interp.inverse_distance_to_grid(lon, lat, h, mx, my, r=5, min_neighbors=1)
    return mx,my,z
###  测试非均匀网格点插值程序
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d03.nc'
# ds_input = xr.open_dataset(flnm)
# ds_input.XTIME
###  测试结束

if __name__ == '__main__':
    # main()
    pass
    