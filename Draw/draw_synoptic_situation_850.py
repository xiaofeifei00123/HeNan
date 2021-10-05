#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
200hPa的急流
形势场图， 风场、等高线、等温线



处理了一个问题是， 图片是不是多张合成一张

画单个图，还是多个图
-----------------------------------------
Time             :2021/09/13 11:39:04
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import cmaps
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from multiprocessing import Pool
from netCDF4 import Dataset
from wrf import get_cartopy, smooth2d, getvar
import metpy.calc as ca
import metpy
from metpy.units import units


plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# %%
class Data():
    """数据处理相关的类
    """
    def get_data(self,):
        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/ERA5_202107.nc'
        ds = xr.open_dataset(flnm)
        ds = ds.rename({'level':'pressure', 't':'temperature', 'longitude':'lon', 'latitude':'lat'})
        # ds.level
        return ds

# %%
class Map():
    """控制地图的类
    """
    def __init__(self, ax):
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.ax = ax

    def draw_station(self, ):
        """在地图上标记站点
        """
        pass
        # station = station_dic
        station = {
            'ZhengZhou': {
                'abbreviation':'郑州',
                'lat': 34.75,
                'lon': 113.62
            },
        }
        # station = station_dic
        values = station.values()
        station_name = list(station.keys())
        station_name = []
        x = []
        y = []
        for i in values:
            y.append(float(i['lat']))
            x.append(float(i['lon']))
            station_name.append(i['abbreviation'])

        # 标记出站点
        self.ax.scatter(x,
                   y,
                #    color='black',
                   color='black',
                #    color='red',
                   transform=ccrs.PlateCarree(),
                   alpha=1.,
                   linewidth=5,
                   s=20,
                   )
        # 给站点加注释
        for i in range(len(x)):
            # print(x[i])
            self.ax.text(x[i]-0.4,
                    y[i] + 0.4,
                    station_name[i],
                    transform=ccrs.PlateCarree(),
                    alpha=1.,
                    fontdict={
                    'size': 28,
            })

    def draw_patch(self, ):
        """在地图上绘制画矩形框
        """
        area = [None] * 3  # 设置一个维度为8的空列表
        area[0] = {"lat1": 33.5, "lat2": 40, "lon1": 80, "lon2": 90}  # north
        area[1] = {"lat1": 28, "lat2": 33,
                   "lon1": 83, "lon2": 94}  # south left
        area[2] = {"lat1": 26, "lat2": 33,
                   "lon1": 95, "lon2": 103}  # south right
        for i in range(3):
            lon = np.empty(4)
            lat = np.empty(4)
            lon[0], lat[0] = area[i]['lon1'], area[i]['lat1']
            lon[1], lat[1] = area[i]['lon2'], area[i]['lat1']
            lon[2], lat[2] = area[i]['lon2'], area[i]['lat2']
            lon[3], lat[3] = area[i]['lon1'], area[i]['lat2']
            x, y = lon, lat
            xy = list(zip(x, y))
            poly = plt.Polygon(xy, edgecolor="red", fc="none", lw=.9, alpha=1)
            self.ax.add_patch(poly)

    def create_map(self,):
        """在底图上添加底图的要素, 省界等
        """
        proj = ccrs.PlateCarree()
        # proj = ccrs.LambertConformal(central_latitude=33, central_longitude=120)  # 创建坐标系
        # --设置地图属性
        # 画省界
        provinces = cfeat.ShapelyFeature(
            Reader(self.path_province).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')
        
        city = cfeat.ShapelyFeature(
            Reader(self.path_city).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')

        Tibet = cfeat.ShapelyFeature(
            Reader(self.path_tibet).geometries(),
            proj,
            edgecolor='k',
            facecolor='none')

        # --设置图像刻度
        ax = self.ax
        
        # # -- 设置图像范围
        # ax.set_extent([108, 117, 31, 38], crs=ccrs.PlateCarree())
        
        ## --- 添加边界线
        # ax.add_feature(provinces, linewidth=1, zorder=2)
        # # ax.add_feature(Tibet, linewidth=2, zorder=2)  # 添加青藏高原区域
        ax.add_feature(provinces, linewidth=0.7, zorder=10)
        # # ax.add_feature(city, linewidth=1, zorder=2)  # 添加青藏高原区域

        ## -- 画海岸线
        # ax.coastlines(resolution='110m')

        
        ## --设置网格属性, 不画默认的标签
        gl = ax.gridlines(draw_labels=True,
                        dms=True,
                        linestyle=":",
                        linewidth=0.3,
                        x_inline=False,
                        y_inline=False,
                        color='k')
        # # gl=ax.gridlines(draw_labels=True,linestyle=":",linewidth=0.3 , auto_inline=True,x_inline=False, y_inline=False,color='k')

        ## 关闭上面和右边的经纬度显示
        gl.top_labels = False  #关闭上部经纬标签
        # gl.bottom_labels = False
        # # gl.left_labels = False
        gl.right_labels = False
        ## 这个东西还挺重要的，对齐坐标用的
        gl.rotate_labels = None

        gl.xformatter = LONGITUDE_FORMATTER  #使横坐标转化为经纬度格式
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator(np.arange(80, 130 + 1, 5))
        gl.ylocator = mticker.FixedLocator(np.arange(10, 50+1, 5))
        gl.xlabel_style = {'size': 16}  #修改经纬度字体大小
        gl.ylabel_style = {'size': 16}
        ax.spines['geo'].set_linewidth(0.7)  #调节边框粗细
        return ax
        

class DataShow():
    """绘制等值线图等各式图形
    """

    def draw_contourf(self, da, ax):
        """在地图上绘制填色图
        """
        x = da.lon
        y = da.lat
        # levels = [205, 210, 215, 220, 225, 235, 245]  # 需要画出的等值线
        # contour_levels = np.arange(30,61,5)
        # contour_levels = np.arange(-0.004,0.004+0.001,0.001)
        # contour_levels = [-0.004, -0.003, -0.002,-0.001, 0.001, 0.002, 0.003, 0.004]
        # contour_levels = [-0.003, -0.002,-0.001, 0.001, 0.002, 0.003]
        contour_levels = [-0.003, -0.002,-0.001,-0.0005, 0.0005, 0.001, 0.002, 0.003]
        
        # colormap = cmaps.precip3_16lev
        colormap = cmaps.ViBlGrWhYeOrRe
        # colormap = cmaps.precip3_16lev_r  # 反转色标
        crx = ax.contourf(x,
                          y,
                          da,
                          cmap=colormap,
                        #   norm=norm,
                          extend='both',
                        #   extend='max',
                          levels=contour_levels,
                          transform=ccrs.PlateCarree())
        return crx

    def draw_contour_temp(self, da, ax, **kw):
        """在地图上绘制等温线
        """
        x = da.lon
        y = da.lat
        # level_contour = np.arange(-32, 36, 4)
        level_contour = np.arange(16, 25, 4)
        da = smooth2d(field=da, passes=12)
        
        crx = ax.contour(x,
                          y,
                          da,
                          colors = 'red',
                          levels=level_contour,
                          linestyles = 'dashed',
                          transform=ccrs.PlateCarree())
        ax.clabel(crx,inline=1, fontsize=15, colors='red') # 等值线的标注
        return crx

    def draw_contour_height(self, da, ax, **kw):
        """在地图上绘制等值线
        """
        x = da.lon
        y = da.lat
        if da.pressure.values == 500:
            level_contour = np.arange(556, 596,4) # 500hPa
        elif da.pressure.values == 700:
            level_contour = np.arange(3040, 3100,10) # 700hPa
        elif da.pressure.values == 850:
            level_contour = np.arange(140, 156,4) # 850hPa
        else:
            print(da.pressure.values)

        da = smooth2d(field=da, passes=8)
        
        crx = ax.contour(x,
                          y,
                          da,
                          colors = 'blue',
                          levels=level_contour,
                          linestyles = 'solid',
                        #   levels=self.levels,
                          transform=ccrs.PlateCarree())
        ax.clabel(crx,inline=1, fontsize=15, colors='blue') # 等值线的标注
        return crx
        


    def draw_barbas(self,u,v):
        '''
        绘制风羽图图
        '''
        u = u[::12,::12]
        v = v[::12,::12]
        # y = u.coords['lat']
        # x = u.coords['lon']
        y = u.lat.values
        x = u.lon.values
        
        # emptyarb设置风旋转的那一点的大小;spacing是F的两条线之间的距离;height:F的横线长短
        self.ax.barbs(x, y, u.values, v.values,length=6,pivot='middle',
            sizes= {'emptybarb':0, 'spacing':0.2, 'height':0.5}, transform=ccrs.PlateCarree())
        # self.ax.barbs(x, y, u.values, v.values)


    def draw_quiver(self,u,v):
        '''
        绘制风矢图
        '''
        u = u[::12,::12]
        v = v[::12,::12]
        # y = u.coords['lat']
        y = u.lat.values
        x = u.lon.values
        # print(y)
        # print(type(u))
        # print(type(y))
        # x = u.coords['lon']


        ax = self.ax
        # Q = ax.quiver(x,y,u,v,units='inches',scale=18,pivot='middle')  # 绘制风矢
        # Q = ax.quiver(x, y, u,v,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=30*10**1,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        qk = ax.quiverkey(Q, X=0.8, Y=0.15, U=300, label=r'$300\ g/kg \cdot m/s$', labelpos='E',coordinates='figure')   # 设置参考风矢

    def draw_streamplot(self,u,v, **kw):
        '''
        绘制风矢图
        '''
        u = u[::12,::12]
        v = v[::12,::12]
        # y = u.coords['lat']
        y = u.lat.values
        x = u.lon.values
        if kw['color']:
            color = kw['color']
        else:
            color = 'blue'

        ax = self.ax
        sr = ax.streamplot(x, y, u.values, v.values, density=[1, 1], 
                           arrowsize = 2.5,
                           arrowstyle = '->',
                           color=color, 
                           transform=ccrs.PlateCarree())  # 绘制风矢
        # Q = ax.quiver(x,y,u,v,units='inches',scale=18,pivot='middle')  # 绘制风矢
        # Q = ax.quiver(x, y, u,v,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        # qk = ax.quiverkey(Q, X=0.7, Y=0.15, U=10, label=r'$10\ m/s$', labelpos='E',coordinates='figure')   # 设置参考风矢

class Axes(Map, DataShow):
    """图片的子类
    所有的操作都是针对一个ax的
    针对一个ax的
    所有的都和ax有关
    所有图片的基类
    控制地图还有DataShow的ax， 在这里集中
    """

    def __init__(self, ax):
        super().__init__(ax)

    def draw_picture(self, fig, dict, picture_dic):
        """绘制单幅图
        设置它的属性之类的

        Args:
            ax ([type]): [画纸对象]
            ds ([Dataset]): [数据]
            picture_dic ([dict]): [图片属性，标注的字典]
        """

        ## 设置地图的范围
        ax = self.ax
        ## 设置图片标题
        self.create_map()
        ax.set_title(picture_dic['title'], fontsize=30)
        ax.set_title(picture_dic['pressure'], fontsize=20,loc='right')
        ## 画填色图 
        cs = self.draw_contourf(dict['div_q'], ax)
        # cs = self.draw_contourf(dict['speed'], ax,)
        ## 画等值线图
        ccc = self.draw_contour_height(dict['geopt']/10, ax)
        # ccc = self.draw_contour_temp(dict['temperature'], ax)
        ## 画风场
        # cccc = self.draw_barbas(dict['u'], dict['v'])
        # cccc = self.draw_quiver(dict['u'], dict['v'])
        cccc = self.draw_quiver(dict['qu'], dict['qv'])
        ## 画流量图
        # cccc = self.draw_streamplot(dict['u'], dict['v'], color='black')
        ## 标注站点
        self.draw_station()
        cb = fig.colorbar(
            cs,
            orientation='horizontal',
            fraction=0.04,  # 色标大小
            pad=0.08,  # colorbar和图之间的距离
        )
        font = {'size':15}
        # qk = ax.quiverkey(Q, X=0.8, Y=0.15, U=300, label=r'$300\ g/kg \cdot m/s$', labelpos='E',coordinates='figure')   # 设置参考风矢
        # cb.set_label('water vapor flux divergence ($\nabla \cdot (q	\overrightarrow{v})$,g/(s*hPa*cm))', fontdict=font)
        cb.set_label(r'water vapor flux divergence ,$\nabla \cdot (q \overrightarrow{v})$  (g/(s*hPa*cm))', fontdict=font)
        # return cs


class Figure():
    """
    所有的操作都是针对的一个figure的
    画图这个动作的类
        创建画纸
        保存图片
        主要是提供ax
        对图片做一些标注
    """
    def one_subplot(self,dict):
        """创建ax, 设置投影
        画单张图

        Args:
            dict ([type]): [description]
        """
        pass 
        fig = plt.figure(figsize=(12, 10), dpi=300)
        ## 手动设置投影格式
        proj = ccrs.LambertConformal(
            central_longitude=105.0,
            central_latitude=32,
            standard_parallels=(30, 60),
            cutoff=-30,)
        # proj = ccrs.LambertConformal(central_longitude=120, central_latitude=33)
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        extent= [75, 135, 10, 60]
        # extent= [88, 110, 20, 40]
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        # ax.set_extent(extent, crs=proj)
        
        tt = dict['q'].time
        titile = str(tt.dt.strftime('%Y-%m-%d %H').values)
        pressure = str(int(dict['u'].pressure.values))+' hPa'
        picture_dic = {'pressure':pressure, 'title':titile}
        
        pr = Axes(ax, )
        pr.draw_picture(fig, dict, picture_dic)
        figpath = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar_OBS/'
        fig_name = str(tt.dt.strftime('%Y-%m-%d_%H').values)+"_"+pressure+'_ERA5_'
        fig.savefig(figpath+fig_name)

    def multi_subplot():
        pass


##### 后面这些是调用前面的类画图, 不同时次，不同层次，使用多线程之类的
##### 然后画图又包括获得数据，切分处理成需要的数据，画图, 保存图片
##### 很多时候画图是根据数据来的，多个数据做循环肯定是画多个图了



class Draw():
    """画图，对单个或多个数据进行画图
    """
    def draw_single(self, ds, kw):
        """
        画单个时次，单层的图
        """
        ## 处理数据
        pre = kw['pre']
        selt = kw['select_time']
        # print(kw)
        # title = str(selt.dt.strftime('%Y-%m-%d %H').values)
        title = str(selt.strftime('%Y-%m-%d %H'))
        print('画%s时刻的图'%title)
        # print(title)
        dict = {}
        ds.sel(pressure=pre, time=selt)
        # dict['q'] = ds_fnl['q'].sel(time=t, pressure=pre).T*10**3
        dict['q'] = ds['q'].sel(time=selt, pressure=pre)*10**3
        # dict['height_agl'] = ds['height_agl'].sel(pressure=pre, time=t)
        dict['temperature'] = ds['temperature'].sel(pressure=pre, time=selt)-273.15
        dict['geopt'] = ds['z'].sel(pressure=pre, time=selt)/10
        dict['u'] = ds['u'].sel(pressure=pre, time=selt)
        dict['v'] = ds['v'].sel(pressure=pre, time=selt)
        dds = ds.sel(pressure=pre, time=selt)
        u = dds.metpy.parse_cf('u')
        v = dds.metpy.parse_cf('v')
        qu = dict['q']*u
        qv = dict['q']*v
        dict['qu'] = qu
        dict['qv'] = qv
        div = ca.divergence(u,v)
        div_q = ca.divergence(qu, qv)
        dict['div'] = div
        dict['div_q'] = div_q
        # dict['speed'] = xr.ufuncs.sqrt(dict['u']**2+dict['v']**2)
        ## 画图
        fi = Figure()
        fi.one_subplot(dict)
    

    def draw_multi(self,):
        """单进程绘图
        """
        # fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar.nc'
        # ds_fnl = xr.open_dataset(fl_wrf)
        gd = Data()
        ds = gd.get_data()
        pre = 200
        for t in ds.time:
            pass
            kw = {'pre':200, 'select_time':t}
            self.draw_single(ds, kw)

    def draw_multi_pool(self,):
        """多进程绘图
        """
        # fl = '/mnt/zfm_18T/fengxiang/DATA/FY_TBB/TBB_FY2G_201607.nc'
        # ds = xr.open_dataset(fl)
        # fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar.nc'
        # fl = '/mnt/zfm_18T/fengxiang/HeNan/Data/GDAS/'
        # ti_list = ['1800', '1812', '1900', '1912']
        
        # for i in ti_list:
            # fl_wrf = fl+'YSU_'+i+'_upar_d02.nc'
            # ds_fnl = xr.open_dataset(fl_wrf)
            # pre = 500
            # pre_list = [500, 700, 850]
            # pre_list = [500]
        gd = Data()
        ds = gd.get_data()
        # for pre in pre_list:
        pool = Pool(10)
        time_index = pd.date_range(start='2021-07-18 00', end='2021-07-22 00', freq='1H')
        # for t in ds.time:
        for t in time_index:
            # draw_single(t, ds_fnl, pre)
            kw = {'pre':850, 'select_time':t}
            tr = pool.apply_async(self.draw_single, args=(ds,kw,))
        pool.close()
        pool.join()


### 测试开始  ###
def test():
    pass
    gd = Data()
    ds = gd.get_data()
    fltime = '2021-07-20 1500'
    # tt = ds.time.sel(time=fltime)
    tt = pd.Timestamp('2021-07-20 1700')

    import metpy.calc as ca
    import metpy
    from metpy.units import units
    u = ds['u'].sel(pressure=500, time=tt)
    v = ds['v'].sel(pressure=500, time=tt)
    # v
    dx = units.Quantity(3, "km")
    dy = units.Quantity(3, "km")
    # dy = units.Quantity(3000, "m")
    # temperature = units.Quantity(t.values, "degC")
    # ca.divergence(u,v)
    # u = u.rename({'lat':'latitude','lon':'longitude'})
    # v = v.rename({'lat':'latitude','lon':'longitude'})
    # proj = ccrs.LambertConformal(
    #     central_longitude=105.0,
    #     central_latitude=32,
    #     standard_parallels=(30, 60),
    #     cutoff=-30,)
    # proj = ccrs.PlateCarree()
    # ds.metpy.parse_cf('pv')
    dds = ds.sel(pressure=500, time='2021-07-20 1200')
    # ddds = dds.rename({'lat':'latitude', 'lon':'longitude'})
    # dds
    bds = dds['u']
    # bbds = bds.assign_attrs({'metpy_crs': ccrs.PlateCarree()})
    # bbds
    u = dds.metpy.parse_cf('u')
    v = dds.metpy.parse_cf('v')
    ca.divergence(u,v)
    # ca.first_derivative(u)
    # cds.metpy_crs
    # bbs = bds.assign_coords({'metpy_crs':ccrs.PlateCarree()})
    # bbs

# ds
gd = Data()
ds = gd.get_data()
tt = pd.Timestamp('2021-07-20 1200')
kw = {'pre':850, 'select_time':tt}
dr = Draw()
dict = dr.draw_single(ds, kw)
# dr.draw_single(dict)
### 测试结束  ###
# %%
if __name__ == "__main__":
    pass
    # dr = Draw()
    # dr.draw_multi()
    # dr.draw_multi_pool()
    # single_process_draw()
    # multi_process_draw()