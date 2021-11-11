import numpy as np
import xarray as xr
import pandas as pd
from metpy import calc as ca  # calc是一个文件夹
from metpy.units import units  # units是units模块中的一个变量名(函数名？类名？)
from metpy import constants  # constatns是一个模块名

# %%
class qv_div():
    """水汽通量散度的计算,
    这里分装这个类的意图是使得计算水汽通量散度更好区分，
    每个函数需要传递的ds没有相关关系
    """
    def __init__(self) -> None:
        pass
    def caculate_qv_div(self, ds):
        """计算单层, 单个时次水汽通量散度

        Args:
            ds (Dataset): 
                lat, lon坐标的多维数组, 包含有q,u,v等变量
                q: kg/kg
                u: m/s
                v: m/s

        Returns:
            qv_div : 计算好的水汽通量散度
        """
        lon = ds.lon
        lat = ds.lat
        u = ds.u*units('m/s')
        v = ds.v*units('m/s')
        q = ds.q*10**3*units('g/kg')
        qv_u = q*u/constants.g
        qv_v = q*v/constants.g

        dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
        qv_div = ca.divergence(u=qv_u, v=qv_v, dx=dx, dy=dy)
        return qv_div

    # dds.sel(pressure=500)
    def caculate_qv_div_integral(self, dds):
        """计算单个时次整层水汽通量散度, 多个层次的积分值
        其实就是个算梯度积分的函数

        Args:
            dds (Dataset): 包含u,v,q(kg/kg)变量, 维度是lat, lon, pressure三维, 单个时次的

        Example:
            flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
            ds = xr.open_dataset(flnm_wrf)
            dds = ds.sel(pressure=[1000, 925, 850, 700, 500]).isel(time=0)
        """
        lev = dds.pressure
        qv_div_list = []
        for i in lev.values:
            ds_one = dds.sel(pressure=i)
            qv_div = self.caculate_qv_div(ds_one)
            qv_div_list.append(qv_div)
        qv_div_integral = np.trapz(qv_div_list[::-1], lev.values[::-1], axis=0)
        da_qv_div_integral = xr.DataArray(qv_div_integral, 
                        coords={
                            'lat':qv_div.lat.values,
                            'lon':qv_div.lon.values,
                        }, 
                        dims=['lat', 'lon'])
        return da_qv_div_integral


    def caculate_qv_div_time(self, ds):
        """计算某一层所有时次的水汽通量散度
        例如500hPa多个时次的水汽通量散度
        把不同时次的值聚合在一起而已

        Args:
            ds ([type]): [description]

        Returns:
            [type]: [description]
        """
        pass
        tt = ds.time
        qv_div_list = []
        for t in tt:
            qvd = self.caculate_qv_div(ds.sel(time=t))
            qv_div_list.append(qvd)
        qv_div_dual_time = xr.concat(qv_div_list, dim='time')
        return qv_div_dual_time
            
    def caculate_qv_div_integral_time(self, ds):
        """计算多个时次的整层水汽通量散度
        把不同时次的值聚合在一起而已

        Args:
            ds (Dataset): 多时次，多层次的, 含有变量q,v,u

        Returns:
            [type]: [description]
        """
        pass
        tt = ds.time
        qv_div_integral_list = []
        for t in tt:
            qvdi = self.caculate_qv_div_integral(ds.sel(time=t))
            qv_div_integral_list.append(qvdi)
        qv_div_dual_time = xr.concat(qv_div_integral_list, dim='time')
        return qv_div_dual_time

    def test(self,):
        """用来测试上述函数的准确性的
        """
        flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
        ds = xr.open_dataset(flnm_wrf)
        # dds = ds.sel(pressure=500).isel(time=0)
        # dds = ds.sel(pressure=500)
        dds = ds.sel(pressure=[1000, 925, 850, 700, 500])
        cc = self.caculate_qv_div_integral_time(ds)
        # da = self.caculate_qv_div_time(dds)
        print(cc.min())


class qv():
    """水汽通量的计算,
    每个函数需要传递的ds没有相关关系, 只是都是Dataset这种格式,
    所以不能做到统一输入
    """
    def __init__(self) -> None:
        pass

    def caculate_qv(self, ds):
        """计算单层, 单个时次水汽通量

        Args:
            ds (Dataset): 
                lat, lon坐标的多维数组, 包含有q,u,v等变量
                q: kg/kg
                u: m/s
                v: m/s

        Returns:
            qv_div : 计算好的水汽通量散度
        """
        # lon = ds.lon
        # lat = ds.lat
        u = ds.u*units('m/s')
        v = ds.v*units('m/s')
        q = ds.q*10**3*units('g/kg')
        qv_u = q*u/constants.g
        qv_v = q*v/constants.g
        qf = xr.ufuncs.sqrt(qv_u**2+qv_v**2) # 水汽通量的大小,模

        dda = xr.concat([qv_u, qv_v, qf], pd.Index(['qv_u', 'qv_v', 'qv_f'], name='model'))
        dds = dda.to_dataset(dim='model')

        # dic_qv = {
        #     'qv_u':qv_u,
        #     'qv_v':qv_v,
        #     'qv_f':qf
        # }
        # return dic_qv
        return dds

        # dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
        # qv_div = ca.divergence(u=qv_u, v=qv_v, dx=dx, dy=dy)
        # return qv_div

    # dds.sel(pressure=500)
    def caculate_qv_integral(self, dds):
        """计算水汽通量整层积分

        Args:
            dds (Dataset): 包含u,v,q(kg/kg)变量, 维度是lat, lon, pressure三维, 单个时次的

        Example:
            flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
            ds = xr.open_dataset(flnm_wrf)
            dds = ds.sel(pressure=[1000, 925, 850, 700, 500]).isel(time=0)
        """
        lev = dds.pressure
        qv_list = []
        for i in lev.values:
            ds_one = dds.sel(pressure=i)
            qv_ds = self.caculate_qv(ds_one)  # qv, qu, qf三个变量同时处理了, 这样做可以吗
            qv = qv_ds.to_array(dim='model')  # 将Dataset 处理为DataArray
            qv_list.append(qv)  # 水汽通量大小, 以及两个矢量

        qv_integral = np.trapz(qv_list[::-1], lev.values[::-1], axis=0)
        da_qv_integral = xr.DataArray(qv_integral, 
                        coords={
                            'model':qv.model.values,
                            'lat':qv.lat.values,
                            'lon':qv.lon.values,
                        }, 
                        dims=['model','lat', 'lon'])
        ds_qv_integral = da_qv_integral.to_dataset(dim='model')
        
        return ds_qv_integral


    def caculate_qv_time(self, ds):
        """计算某一层所有时次的水汽通量
        例如500hPa多个时次的水汽通量
        把不同时次的值聚合在一起而已

        Args:
            ds ([type]): [description]

        Returns:
            [type]: [description]
        """
        pass
        tt = ds.time
        qv_list = []
        for t in tt:
            qvd = self.caculate_qv(ds.sel(time=t))
            qv_list.append(qvd)
        qv_div_dual_time = xr.concat(qv_list, dim='time')
        return qv_div_dual_time
            
    def caculate_qv_integral_time(self, ds):
        """计算多个时次的整层水汽通量散度
        把不同时次的值聚合在一起而已

        Args:
            ds (Dataset): 多时次，多层次的, 含有变量q,v,u

        Returns:
            [type]: [description]
        """
        pass
        tt = ds.time
        qv_div_integral_list = []
        for t in tt:
            qvdi = self.caculate_qv_div_integral(ds.sel(time=t))
            qv_div_integral_list.append(qvdi)
        qv_div_dual_time = xr.concat(qv_div_integral_list, dim='time')
        return qv_div_dual_time

    def test(self,):
        """用来测试上述函数的准确性的
        """
        pass
        flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
        ds = xr.open_dataset(flnm_wrf)
        # dds = ds.sel(pressure=500).isel(time=0)
        # dds = ds.sel(pressure=500)
        dds = ds.sel(pressure=[925, 850, 700, 500])  # 1000hPa的数据存在很多缺测，故舍弃
        cc = self.caculate_qv_div_integral_time(ds)
        # da = self.caculate_qv_div_time(dds)
        print(cc.min())

flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
ds = xr.open_dataset(flnm_wrf)
# dds = ds.sel(pressure=[1000, 925, 850, 700, 500]).isel(time=0)
dds = ds.sel(pressure=[925, 850, 700, 500]).isel(time=0)
# dds = ds.sel(pressure=[925]).isel(time=0)

cc = qv()
# da = cc.caculate_qv(dds)
ds = cc.caculate_integral(dds)
# da
# %%
# lev = dds.pressure
# lev
# len(li)
# li[::-1]
# li[4].mean()
# qv_integral = np.trapz(li[::-1], lev.values[::-1], axis=0)
# da = xr.DataArray(qv_integral)
# %%
# (li[0]+li[1]+li[2]+li[3]).plot()
# li[4].plot()

# ds['qv_v'].plot()
# ds.plot()
# ds.mean()
# ds.plot()
# ds.min()


# cc = self.caculate_qv_div_integral_time(ds)
# da = self.caculate_qv_div_time(dds)
# print(cc.min())