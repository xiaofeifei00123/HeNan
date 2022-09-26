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

# %%

class Rain():
    pass
    def __init__(self):

        self.cross_start = [111.5, 33]
        self.cross_end = [114.3, 36]

        self.station = {
                'A': {
                    'abbreviation':'A',
                    'lat': 33.6,
                    'lon': 112.1,
                },
                'B': {
                    'abbreviation':'B',
                    'lat': 34.5,
                    'lon': 112.9,
                },
                'C': {
                    'abbreviation':'C',
                    'lat': 35.15,
                    'lon': 113.5,
                },
                'D': {
                    'abbreviation':'D',
                    'lat': 35.75,
                    'lon': 114.0,
                },
            }

class Draw(Rain):
    """画单张降水图的
    只管画图，不管标注
    Args:
        object (_type_): _description_
    """

    def __init__(self, fig, ax) -> None:
        super().__init__()
        self.fig = fig
        self.ax = ax
        # self.colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
        # self.colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表


        # self.colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250, 400,600, 1000]#雨量等级
        # self.colordict = select_cmap('rain9')

        # self.colorlevel=[0, 1, 10, 25, 50, 100, 250, 400,600,800,1000, 2000]#雨量等级
        # rgbtxt = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_6d/11colors.txt'
        # rgb = get_rgb(rgbtxt)
        # self.colordict = rgb        
        
        self.colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
        self.colorticks = self.colorlevel[1:-1]
        self.colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
        

        self.colorticks = self.colorlevel[1:-1]
        self.map_dic = {
                'proj':ccrs.PlateCarree(),
                'extent':[110.5, 116, 32, 36.5],
                'extent_interval_lat':1,
                'extent_interval_lon':1,
            } 

        # self.cross_start = [111.2, 32]
        # self.cross_end = [114.8, 36.5]
        # self.cross_start = [111.5, 33]
        # self.cross_end = [114.3, 36]

        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'

        # self.station = {
        #         'ZhengZhou': {
        #             'abbreviation':'郑州',
        #             'lat': 34.76,
        #             'lon': 113.65
        #         },
        #         'NanYang': {
        #             'abbreviation':'南阳',
        #             'lat': 33.1,
        #             'lon': 112.49,
        #         },
        #         'LuShi': {
        #             'abbreviation':'卢氏',
        #             'lat': 34.08,
        #             'lon': 111.07,
        #         },
        #     }
        # loc1 = (33.7, 112.55)
        # loc2 = (34.7, 113.35)
        # loc3 = (35.5, 114.0)

        # self.station = {
        #         'A': {
        #             'abbreviation':'A',
        #             'lat': 33.5,
        #             'lon': 112.55,
        #         },
        #         'B': {
        #             'abbreviation':'B',
        #             'lat': 34.7,
        #             'lon': 113.35,
        #         },
        #         'C': {
        #             'abbreviation':'C',
        #             'lat': 35.5,
        #             'lon': 114.0,
        #         },
        #     }



    def add_patch(self, area, ax):
            xy = (area['lon1'], area['lat1'])
            width = area['lon2']-area['lon1']
            height = area['lat2']-area['lat1']
            rect = patches.Rectangle(xy=xy, width=width, height=height, edgecolor='black', fill=False, lw=1.5, ) # 左下角的点的位置
            ax.add_patch(rect)

    def draw_single(self, da,):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        ax = self.ax
        ## 给图像对象叠加地图
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])
        # mp.add_station(ax, self.station, justice=True, delx=-0.1)

        
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
        # mp.add_station(ax, self.station, justice=True)
        # plt.plot(self.cross_start, self.cross_end, transform=ccrs.PlateCarree(),zorder=5)
        # plt.plot(112, 35, 113,36,transform=ccrs.PlateCarree(),zorder=5)
        # ax.scatter(np.linspace(self.cross_start[0], self.cross_end[0], 10), np.linspace(self.cross_start[1], self.cross_end[1], 10), transform=ccrs.PlateCarree())
        # ax.plot(np.linspace(self.cross_start[0], self.cross_end[0], 10), np.linspace(self.cross_start[1], self.cross_end[1], 10), color='black')

        areaA = {
            'lat1':33.6,
            'lat2':34.0,
            'lon1':111.8,
            'lon2':112.2,
        }        
        areaB = {
            'lat1':33.4,
            'lat2':33.8,
            'lon1':112.4,
            'lon2':112.8,
        }        
        self.add_patch(areaA, ax)
        self.add_patch(areaB, ax)
        # ax.text(114.4, 33.7, 'D', transform=ccrs.PlateCarree())
        return crx
        
    def draw_tricontourf(self, rain):
        """rain[lon, lat, data],离散格点的DataArray数据
        由离散格点的数据绘制降水
        Args:
            rain ([type]): [description]
        Example:
        da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
        da.max()
        rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
        """
        ax = self.ax
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        # mp.add_station(ax, self.station, justice=True)

        ax.set_extent(self.map_dic['extent'])
        cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=self.colorlevel,colors=self.colordict, transform=ccrs.PlateCarree())

        rain_max = rain.max()        
        # ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='right')
        # ax.set_title('2021-07 20/00--21/00', fontsize=35,)
        # ax.set_title('OBS', fontsize=10,loc='right')
        # ax.plot(np.linspace(self.cross_start[0], self.cross_end[0], 10), np.linspace(self.cross_start[1], self.cross_end[1], 10), color='black')
        return cs


class GetData():
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

        
        

        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify/CTRL/rain_d02.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Typhoon/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
        # da = xr.open_dataarray(flnm)
        # # da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        # da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        # da = da.sum(dim='time') 
        # return da


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

def get_dr():
    cm = round(1/2.54, 2)
    proj = ccrs.PlateCarree()  # 创建坐标系
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    ax = fig.add_axes([0.13,0.1,0.82,0.8], projection=proj)
    dr = Draw(fig, ax)
    return dr

def draw_obs(tl):
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_obs = xr.open_dataset(flnm_obs)
    dr = get_dr()  # 画图的对象
    # ds_obs = ds_obs.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    ds_obs = ds_obs.sel(time=tl)
    da_obs = ds_obs.sum(dim='time')['PRCP']
    cf = dr.draw_single(da_obs)    

    cb = dr.fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorticks,
        fraction = 0.06,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    labels = list(map(lambda x: str(x) if x<1 else str(int(x)), dr.colorticks))  # 将colorbar的标签变为字符串
    cb.set_ticklabels(labels)
    dr.ax.set_title('OBS', fontsize=10,loc='left')
    t1 = str(ds_obs.time.dt.strftime('%d/%H')[0].values)
    t2 = str(ds_obs.time.dt.strftime('%d/%H')[-1].values)
    tfig = str(ds_obs.time.dt.strftime('%d%H')[-1].values)
    tt = t1+'-'+t2
    # dr.ax.set_title(tt, fontsize=10,loc='center')
    
    fig_name = 'obs'+tfig+'aa'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_new/'
    dr.fig.savefig(fig_path+fig_name)

def draw_model(tl):
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    ds_model = xr.open_dataset(flnm_model)
    # for model in ['gwd3','gwd1', 'gwd0']:
    # model_list = ['SS2']
    # model_list = ['CTRL','FD','GWD3','SS']
    # model_list = ['GWD3']
    # for model in model_list:
    dr = get_dr()  # 画图的对象
    gd = GetData()  # 数据的对象
    # da = gd.onemodel(model)
    # da = gd.onemodel_newall(model)
    # da
    # ds = ds_model.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    ds = ds_model.sel(time=tl)
    da1 = ds['CTRL'].sum(dim='time')
    da2 = ds['GWD3'].sum(dim='time')
    da = da2 - da1
    cf = dr.draw_single(da)    


    cb = dr.fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorticks,
        fraction = 0.06,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    labels = list(map(lambda x: str(x) if x<1 else str(int(x)), dr.colorticks))  # 将colorbar的标签变为字符串
    cb.set_ticklabels(labels)
    model = 'minus'
    dr.ax.set_title(model, fontsize=10,loc='left')
    t1 = str(ds.time.dt.strftime('%d/%H')[0].values)
    t2 = str(ds.time.dt.strftime('%d/%H')[-1].values)
    tt = t1+'-'+t2
    tfig = str(ds.time.dt.strftime('%d%H')[-1].values)
    dr.ax.set_title(tt, fontsize=10,loc='center')
    
    fig_name = model+tfig+'minus'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_new/'
    dr.fig.savefig(fig_path+fig_name)

if __name__ == '__main__':
    pass
    # draw_model()
    # draw_obs()



    ## 画模式降水
    
    
    time1=slice('2021-07-17 01', '2021-07-18 00')
    time2=slice('2021-07-18 01', '2021-07-19 00')
    time3=slice('2021-07-19 01', '2021-07-20 00')
    time4=slice('2021-07-20 01', '2021-07-21 00')
    time5=slice('2021-07-21 01', '2021-07-22 00')
    # time6=slice('2021-07-22 01', '2021-07-23 00')
    time7=slice('2021-07-17 01', '2021-07-23 00')
    # time8=slice('2021-07-19 13', '2021-07-20 00')
    # time_list = [time1, time2, time3, time4, time5, time6]
    # time_list = [time4,time7]
    time_list = [time7]
    # for tl in time_list[0:1]:
    for tl in time_list:
        draw_model(tl)
        # draw_obs(tl)

# %%