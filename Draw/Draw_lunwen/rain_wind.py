#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图, 小时降水
实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
改造，数据处理和画图完全分开
-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''

# %%
import sys,os
import xarray as xr
import numpy as np
import pandas as pd
import wrf
import netCDF4 as nc

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cmaps
from baobao.map import Map



# %%
class Draw(object):

    def __init__(self,):
        # self.ax
        ## 作为类的属性, 所有的值都可以在对象里面重新赋值, 这样的话，相当于把这个变量赋值了一个初始值，但是可以改变
        self.colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
        self.colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表

        # self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        # self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'

    def draw_quiver(self,u,v,ax, scale=80, ulength=10):
        """风矢

        Args:
            u (_type_): _description_
            v (_type_): _description_
            x (_type_): _description_
            y (_type_): _description_
            ax (_type_): _description_
        """
        # u = u[::3,::3]
        # v = v[::3,::3]
        u = u[::10,::10]
        v = v[::10,::10]
        # u = u[::20,::20]
        # v = v[::20,::20]
        y = u.lat.values
        x = u.lon.values

        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=80,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        qk = ax.quiverkey(Q,
                        X=0.05, Y=0.15, 
                        U=10 ,
                        label=r'$10 m/s$', 
                        labelpos='S',  # label在参考箭头的哪个方向
                        labelsep=0.05, # 箭头和标签之间的距离
                        coordinates='figure',  
                        fontproperties={'size':8}
                        )   # 设置参考风矢

    def draw_contourf(self, data, ax):
        """画填色图
        """

        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140,  700]#雨量等级
        # colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
        # colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        x = data.lon
        y = data.lat
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                          levels=self.colorlevel,
                          colors = self.colordict,
                          transform=ccrs.PlateCarree())
        return crx


def draw_single(da, u, v):
    """画单个的降水及其风场图

    Args:
        da (DataArray): 单个时次的降水
        ds: 单个时次的wrfout数据
    """
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    proj = ccrs.PlateCarree()  # 创建坐标系
    ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        'extent':[110.5, 116, 32, 36.5],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    station = {
        'ZhengZhou': {
            'abbreviation':'郑州',
            'lat': 34.76,
            'lon': 113.65
        },
    }
    mp.add_station(ax, station, justice=True)


    dr = Draw()
    cf = dr.draw_contourf_single(da, ax)
    dr.draw_quiver(u,v, ax)
    colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
    colorticks = colorlevel[1:-1]
    
    cb = fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.08,  #  色标和子图间距离
    )

    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    return fig
    # fig_name = 'rain1h'+picture_dic['type']
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    # fig.savefig(fig_path+fig_name)


def draw_one(model='gwd0'):
    """_summary_

    Args:
        model (str, optional): _description_. Defaults to 'gwd0'.
    """

    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
    da = xr.open_dataarray(flnm)
    t = '2021-07-20 12'
    rain1h = da.sel(time=t)
    flnm_upar = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
    ds = xr.open_dataset(flnm_upar)
    ds = ds.rename({'ua':'u', 'va':'v'})
    level =  900
    ds2 = ds.sel(time=t, pressure=level)
    u = ds2['u']
    v = ds2['v']
    fig = draw_single(rain1h, u, v)
    fig_name = 'rain_wind_hourly'+model
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig.savefig(fig_path+fig_name)


def draw_minus():
    pass

    t = '2021-07-20 12'
    flnm1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'
    flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain.nc'
    da1 = xr.open_dataarray(flnm1).sel(time=t)
    da2 = xr.open_dataarray(flnm2).sel(time=t)
    rain1h = da2-da1

    flnm_upar1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/upar.nc'
    flnm_upar2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
    level =  900
    ds1 = xr.open_dataset(flnm_upar1).sel(time=t, pressure=level)
    ds2 = xr.open_dataset(flnm_upar2).sel(time=t, pressure=level)
    ds = ds2-ds1
    ds = ds.rename({'ua':'u', 'va':'v'})
    # ds2 = ds.sel(time=t, pressure=level)
    u = ds['u']
    v = ds['v']
    fig = draw_single(rain1h, u, v)

    fig_name = 'rain_wind_hourly'+'minus'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig.savefig(fig_path+fig_name)




def draw_dual():
    model_list = ['gwd0', 'gwd3']
    # model_list = ['gwd0','gwd1', 'gwd3', 'gwd3-LS', 'gwd3-BL', 'gwd3-SS', 'gwd3-FD']
    # model_list = ['1912_90m_gwd3']
    for model in model_list:
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/'
        draw_one(model)


if __name__ == '__main__':

    # draw_dual()
    draw_minus()