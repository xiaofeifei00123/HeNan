#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
降水分布图
实况降水，站点插值出
模式降水，原始的wrfout网格点(未插值)


更改这个东西，留下绘制多图的接口
色标， colorlevel, colordict这些可以放在外面定义
这个东西比较麻烦的一个是传参， 传着传着就复杂了
内核需要什么东西要搞清楚，然后要传递啥东西, 不同的图哪些东西是需要变的，哪些是不变的
有些东西需要变，但是在这里不是一个常变的变量，可以在类的属性里面设置， 在类的属性里面设置了，也可以变的好像
类的属性有个默认值之后，还可以改变
对象的属性是可以重新定义的
dr.colorlevel = [0, 1, 3,]
一般而言，画降水嘛，就是数据不一样
不同的试验会是地图啥的不一样，重写类就可以了

东西一定要简洁，思路一定要清晰，让别人和未来的自己一看就明白
把数据和图耦合起来
函数不要过长，也不要过短
-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cmaps
from baobao.map import Map
from baobao.get_cmap import select_cmap
import pandas as pd
import os
import numpy as np
from baobao.map import get_rgb
import matplotlib.patches as patches
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
from common import Common

# %%

class Rain(Common):
    def __init__(self):
        pass
        super().__init__()

class GetData(Rain):

    def obs(self,):
        """这个主要是进行数据处理
        包括创建一个图片对象

        Args:
            dr (_type_): 子图的类
        """
        pass
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc'
        da = xr.open_dataarray(flnm)
        # da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        da = da.sel(time=slice('2021-07-17 01', '2021-07-23 00'))
        da = da.sum(dim='time') 
        return da

    def onemodel(self, model='gwd3'):
        pass

        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain_d03.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_d02_grd.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify/'+model+'/'+'rain_wrfout_d03.nc'
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'+model+'/wrfout/'+'rain_wrfout_d03.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify/CTRL/rain_d02.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Typhoon/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
        da = xr.open_dataarray(flnm)
        # da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        da = da.sum(dim='time') 
        return da

    def onemodel_newall(self, model='gwd3'):
        pass

        def acsum(flnm1, flnm2):
            # flnm1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/wrfout_d03_2021-07-17_00:00:00'
            # flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/wrfout_d03_2021-07-17_23:00:00'
            ds1 = xr.open_dataset(flnm1)
            ds2 = xr.open_dataset(flnm2)
            da1 = ds1['RAINNC']+ds1['RAINC']+ds1['RAINSH']
            da2 = ds2['RAINNC']+ds2['RAINC']+ds2['RAINSH']

            ## 对流降水
            # da1 = ds1['RAINC']+ds1['RAINSH']   #  深对流+浅对流
            # da2 = ds2['RAINC']+ds2['RAINSH']   #  深对流+浅对流
            # 格点降水
            # da1 = ds1['RAINNC']   #  
            # da2 = ds2['RAINNC']   #  
            
            r = da2-da1
            r = r.squeeze()  # 该是几维的就是几维的
            r = r.rename({'XLAT':'lat', 'XLONG':'lon'})
            return r

        tt = pd.date_range('2021-07-17 00', '2021-07-22 00', freq='24H')
        # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/'
        path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'+model+'/wrfout/'
        rain = 0.
        for t in tt:
            t1 = t
            t2 = t+pd.Timedelta('23H')
            # print(t2)
            strt1 = 'wrfout_d03_'+t1.strftime('%Y-%m-%d_%H:%M:%S')
            strt2 = 'wrfout_d03_'+t2.strftime('%Y-%m-%d_%H:%M:%S')
            flnm1 = os.path.join(path, strt1)
            flnm2 = os.path.join(path, strt2)
            rain = rain+acsum(flnm1, flnm2) 
        # print(rain)
        return rain

    def EC(self,):
        """绘制EC细网格降水

        Args:
            model (str, optional): _description_. Defaults to 'gwd3'.
        """
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_ec.nc'
        da = xr.open_dataarray(flnm)
        da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
        da = da.sum(dim='time').squeeze()
        return da


class Draw(Rain):
    """画单张降水图的
    只管画图，不管标注
    Args:
        object (_type_): _description_
    """

    def __init__(self,):
        """存储色标和地图文件的路径

        Args:
            fig (_type_): _description_
            ax (_type_): _description_
        """
        super().__init__()
        self.map_dic = {
                'proj':ccrs.PlateCarree(),
                'extent':[110.5, 116, 32, 36.5],
                'extent_interval_lat':1,
                'extent_interval_lon':1,
            } 

        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'

    def set_colorbar(self,cf, ax, fig):    
        ## 设置colorbar
        levels = cf.levels
        colorticks = levels[1:-1]
        cb = fig.colorbar(
            cf,
            # cax=ax,
            # ax,
            orientation='horizontal',
            ticks=colorticks,
            fraction = 0.06,  # 色标大小,相对于原图的大小
            pad=0.1,  #  色标和子图间距离
            )
        ax.tick_params(labelsize=10)  # 设置色标标注的大小
        tic = cb.get_ticks()
        labels = list(map(lambda x: str(x) if x<1 else str(int(x)), tic))  # 将colorbar的标签变为字符串
        cb.set_ticklabels(labels)

    def add_patch(self, area, ax, **kw):
            xy = (area['lon1'], area['lat1'])
            width = area['lon2']-area['lon1']
            height = area['lat2']-area['lat1']
            rect = patches.Rectangle(xy=xy, width=width, height=height, fill=False, lw=1.5, **kw) # 左下角的点的位置
            ax.add_patch(rect)

    def draw_rain_24h(self, da, ax):
        """单张图的功能基本都在这个函数里面进行
        除了保存图片, fig和ax之外

        Args:
            da (DataArray): 单个时次的降水
        """
        ## 画图
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])

        x = da.lon
        y = da.lat

        colorlevel=[0, 1, 10, 25, 50, 100, 250, 400,600,800,1000, 2000]#雨量等级
        rgbtxt = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_6d/11colors.txt'
        rgb = get_rgb(rgbtxt)
        colordict = rgb        
        colorticks = colorlevel[1:-1]
        crx = ax.contourf(x,
                          y,
                          da,
                          corner_mask=False,
                          levels=colorlevel,
                          colors = colordict,
                          transform=ccrs.PlateCarree()
                          )

        ## 设置标注
        if 'south_north' in da.dims:
            rain_max = da.max(dim=['south_north', 'west_east'])        
        elif 'lat' in da.dims:
            rain_max = da.max(dim=['lat', 'lon'])        
        else:
            print("出错啦")
        mp.add_station(ax, self.station_zz, justice=True, delx=-0.1)
        ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='right')

        # ax.text(self.areaA['lon1']-0.3, self.areaA['lat1']+0.2, 'A')
        # ax.text(self.areaB['lon1']-0.3, self.areaB['lat1']+0.2, 'B')
        ax.text(self.areaD['lon1']-0.3, self.areaD['lat1']+0.2, 'B')
        # ax.text(self.areaC['lon1']-0.3, self.areaC['lat1']+0.2, 'C')
        # self.add_patch(self.areaA, ax, color='blue')
        self.add_patch(self.areaD, ax, color='blue')
        # self.add_patch(self.areaC, ax, color='blue')
        
        

        return crx   # 可以从crx里面获得colorticks, 色标相关的变量
        
    def draw_rain_1h(self, data, ax):
        """1小时降水量
        """

        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])

        # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140,  700]#雨量等级
        colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
        colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        x = data.lon
        y = data.lat
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                          levels=colorlevel,
                          colors = colordict,
                          transform=ccrs.PlateCarree())

        # mp.add_station(ax, self.station, justice=True, ssize=10, marker='o')

        return crx

    def draw_rain_minus_24h(self, da,ax):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        # ax = self.ax
        ## 给图像对象叠加地图
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])

        
        if 'south_north' in da.dims:
            rain_max = da.max(dim=['south_north', 'west_east'])        
        elif 'lat' in da.dims:
            rain_max = da.max(dim=['lat', 'lon'])        
        else:
            print("出错啦")
            
        ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='right')

        x = da.lon
        y = da.lat
        crx = ax.contourf(x,
                          y,
                          da,
                          corner_mask=False,
                          levels=self.colorlevel,
                          colors = self.colordict,
                          transform=ccrs.PlateCarree()
                          )
        # mp.add_station(ax, self.station, justice=True, delx=-0.1)
        # areaA = self.areaA
        # areaB = self.areaB
        ax.text(self.areaA['lon1']+0.3, self.areaA['lat1']+0.3, 'A')
        ax.text(self.areaB['lon1']+0.3, self.areaB['lat1']+0.3, 'B')
        ax.text(self.areaC['lon1']+0.3, self.areaB['lat1']+0.3, 'C')
        self.add_patch(self.areaA, ax)
        self.add_patch(self.areaB, ax)
        self.add_patch(self.areaC, ax)
        # ax.text(114.4, 33.7, 'D', transform=ccrs.PlateCarree())
        com = Common()
        ax.plot(np.linspace(com.cross_start[0], com.cross_end[0], 10), np.linspace(com.cross_start[1], com.cross_end[1], 10), color='black')
        return crx

    def draw_tricontourf(self, rain, ax):
        """rain[lon, lat, data],离散格点的DataArray数据
        由离散格点的数据绘制降水
        Args:
            rain ([type]): [description]
        Example:
        da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
        da.max()
        rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
        """
        # ax = self.ax
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        # mp.add_station(ax, self.station, justice=True)

        ax.set_extent(self.map_dic['extent'])
        cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=self.colorlevel,colors=self.colordict, transform=ccrs.PlateCarree())
        return cs


# %%

if __name__ == '__main__':
    pass
