#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图, 小时降水
实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)
改造，数据处理和画图完全分开

画图的类
获得数据的类
综合起来画图的函数

多子图
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
import matplotlib as mpl
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

    def draw_quiver(self,u,v,ax, scale=200, ulength=10):
        """风矢
        sacle 越小，箭头长度越大

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

        # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=scale,pivot='middle', transform=ccrs.PlateCarree(), zorder=3)  # 绘制风矢
        # Q = ax.quiver(x, y, u.values,v.values,units='width',scale=scale,pivot='middle', transform=ccrs.PlateCarree(), zorder=3)  # 绘制风矢
        Q = ax.quiver(x, y, u.values,v.values,units='width',scale=scale,pivot='middle', width=0.005, transform=ccrs.PlateCarree(), zorder=3)  # 绘制风矢
        qk = ax.quiverkey(Q,
                        X=0.8, Y=1.04, 
                        U=ulength ,
                        label=r'${} m/s$'.format(ulength), 
                        labelpos='E',  # label在参考箭头的哪个方向
                        labelsep=0.05, # 箭头和标签之间的距离
                        fontproperties={'size':8}, 
                        coordinates='axes', # 是相对于ax的位置，还是相对于figure的位置
                        edgecolor='white',
                        facecolor='black'
                        )   # 设置参考风矢
        return qk

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
    cf = dr.draw_contourf(da, ax)
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


class GetData():

    def get_data_hourly(self,model='gwd0'):
        """_summary_

        Args:
            model (str, optional): _description_. Defaults to 'gwd0'.
        """

        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
        da = xr.open_dataarray(flnm)
        t = '2021-07-20 12'
        # t = '2021-07-20 00'
        rain1h = da.sel(time=t)

        flnm_upar = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar.nc'
        ds = xr.open_dataset(flnm_upar)
        ds = ds.rename({'ua':'u', 'va':'v'})
        level =  900
        ds2 = ds.sel(time=t, pressure=level)
        u = ds2['u']
        v = ds2['v']
        return rain1h, u, v
        # fig = draw_single(rain1h, u, v)
        # fig_name = 'rain_wind_hourly'+model
        # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
        # fig.savefig(fig_path+fig_name)


    def get_data_minus(self,):
        pass

        # t = '2021-07-20 12'
        t = '2021-07-20 00'
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
        return rain1h, u, v
        # fig = draw_single(rain1h, u, v)

        # fig_name = 'rain_wind_hourly'+'minus'
        # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
        # fig.savefig(fig_path+fig_name)

def add_map(ax,):
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
    return ax

def main():
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    proj = ccrs.PlateCarree()  # 创建坐标系
    ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
## 画降水
    ## 获得降水数据
    gd = GetData()
    # rain, u, v = gd.get_data_hourly('gwd0')
    rain, u, v = gd.get_data_hourly('gwd3')

    ## 将降水画到图上
    ## 添加底图和站点标注
    ax= add_map(ax)
    dr = Draw()

    cf = dr.draw_contourf(rain, ax)
    dr.draw_quiver(u,v, ax)

    cb = fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorlevel[1:-1],
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.08,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小


def multi():
    cm = 1/2.54
    # fig = plt.figure(figsize=(19*cm, 8*cm), dpi=300)
    fig = plt.figure(figsize=(8*cm, 24*cm), dpi=600)
    proj = ccrs.PlateCarree()  # 创建坐标系
    grid = plt.GridSpec(3,
                        1,
                        figure=fig,
                        left=0.1,
                        right=0.95,
                        bottom=0.05,
                        top=0.98,
                        wspace=0.2, # 两子图之间的宽
                        # hspace=0.25)
                        hspace=0.1)  # 两子图之间的高

    num = 3
    axes = [None]*num
    for i in range(num):
        axes[i] = fig.add_subplot(grid[i], projection=proj)


    gd = GetData()
    rain1, u1, v1 = gd.get_data_hourly('gwd0')
    rain2, u2, v2 = gd.get_data_hourly('gwd3')
    rain3, u3, v3 = gd.get_data_minus()

    ax1= add_map(axes[0])
    ax2= add_map(axes[1])
    ax3= add_map(axes[2])

    ax1.set_title('(a)', loc='left', y=0.98, fontsize=9)
    ax2.set_title('(b)', loc='left', y=0.98, fontsize=9)
    ax3.set_title('(c)', loc='left', y=0.98, fontsize=9)

    dr1 = Draw()
    cf1 = dr1.draw_contourf(rain1, ax1)
    dr1.draw_quiver(u1,v1, ax1)

    dr2 = Draw()
    cf2 = dr2.draw_contourf(rain2, ax2)
    dr2.draw_quiver(u2,v2, ax2)

    # dr.colorlevel = []
    dr3 = Draw()
    dr3.colordict=['#0000fb','#3232fd','#6464fd','white','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    dr3.colorlevel= [-90,-20,-10,-3, 0, 3,10,20,90]  # 垂直速度的色标
    cf3 = dr3.draw_contourf(rain3, ax3)
    dr3.draw_quiver(u3,v3, ax3, scale=80, ulength=5)


    cb = fig.colorbar(
        cf1,
        ax=axes[0],
        orientation='horizontal',
        ticks=dr1.colorlevel[1:-1],
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.08,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小


    cb = fig.colorbar(
        cf2,
        ax=axes[1],
        orientation='horizontal',
        ticks=dr2.colorlevel[1:-1],
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.08,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小


    cb = fig.colorbar(
        cf3,
        ax=axes[2],
        orientation='horizontal',
        ticks=dr3.colorlevel[1:-1],
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.08,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小



    fig_name = 'rain1hhhhhhhhhhhh'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig.savefig(fig_path+fig_name)

# %%


if __name__ == '__main__':

    # draw_dual()
    # draw_minus()
    pass