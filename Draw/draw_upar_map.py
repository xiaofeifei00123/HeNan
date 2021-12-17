#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
在地图上画图的模板
cartopy库
只用改ax的投影格式为lambert即可
所有的数据用platecaree
地图上的等值线、填色、风矢图(变为np.dataarray)
-----------------------------------------
Time             :2021/09/13 11:39:04
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import cmaps
import numpy as np
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



plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['SimHei']

# %%
# fl_fnl = '/mnt/zfm_18T/fengxiang/DATA/FNL/fnl_2016.nc'
# fl_gdas = '/mnt/zfm_18T/fengxiang/DATA/FNL/gdas_2016_0710_20.nc'
# fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar.nc'
# ds_fnl = xr.open_dataset(fl_wrf)

# %%
# ds_fnl


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
                   transform=ccrs.PlateCarree(),
                   alpha=1.,
                   linewidth=5,
                   s=30,
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
        # # ax.add_feature(provinces, linewidth=1, zorder=2)
        # # ax.add_feature(Tibet, linewidth=2, zorder=2)  # 添加青藏高原区域
        ax.add_feature(provinces, linewidth=1, zorder=10)
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
        gl.xlocator = mticker.FixedLocator(np.arange(106, 138 + 1, 1))
        gl.ylocator = mticker.FixedLocator(np.arange(20, 41, 1))
        gl.xlabel_style = {'size': 16}  #修改经纬度字体大小
        gl.ylabel_style = {'size': 16}
        ax.spines['geo'].set_linewidth(0.7)  #调节边框粗细
        return ax
        

class DataShow():

    def draw_contourf(self, da, ax):
        """在地图上绘制填色图
        """
        x = da.lon
        y = da.lat
        # levels = [205, 210, 215, 220, 225, 235, 245]  # 需要画出的等值线
        contour_levels = np.arange(2,19,2)
        
        colormap = cmaps.precip3_16lev
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

    def draw_contour(self, da, ax):
        """在地图上绘制等值线
        """
        x = da.lon
        y = da.lat
        # levels = [220,230, 240, 250, 260, 270, 280, 290, 300]  # 需要画出的等值线
        # levels = [205,210, 215,220, 225,235,245]  # 需要画出的等值线
        # levels = [205, 210, 215, 220, 225, 235, 245]  # 需要画出的等值线
        # levels = np.arange(5000, 5800, 10)
        if da.pressure.values == 500:
            level_contour = np.arange(5720, 5770,10) # 500hPa
        elif da.pressure.values == 700:
            level_contour = np.arange(3040, 3100,10) # 700hPa
        elif da.pressure.values == 850:
            level_contour = np.arange(1410, 1470,10) # 850hPa
        else:
            print(da.pressure.values)
        # level_contour = np.arange(5720, 5770,10) # 500hPa
        # level_contour = np.arange(3040, 3100,10) # 700hPa
        # level_contour = np.arange(1410, 1470,10) # 850hPa
        # cmap = cmaps.precip_16lev 
        
        crx = ax.contour(x,
                          y,
                          da/10,
                          colors = 'red',
                        #   cmap=cmap,
                        #   norm=norm,
                        #   extend='both',
                        #   extend='max',
                          levels=level_contour,
                        #   levels=self.levels,
                          transform=ccrs.PlateCarree())
        ax.clabel(crx,inline=1, fontsize=20, colors='black') # 等值线的标注
        return crx

    def draw_barbas(self,u,v):
        '''
        绘制风羽图图
        '''
        u = u[::8,::8]
        v = v[::8,::8]
        y = u.coords['lat']
        x = u.coords['lon']
        
        # emptyarb设置风旋转的那一点的大小;spacing是F的两条线之间的距离;height:F的横线长短
        self.ax.barbs(x, y, u, v,length=5,pivot='middle',
            sizes=dict(emptybarb=0, spacing=0.3, height=0.5))


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
        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        qk = ax.quiverkey(Q, X=0.7, Y=0.15, U=10, label=r'$10\ m/s$', labelpos='E',coordinates='figure')   # 设置参考风矢


class Picture(Map, DataShow):
    """绘图的子类
    所有图片的基类
    """

    def __init__(self, ax):
        super().__init__(ax)

    def draw_single(self, fig, dict, picture_dic):
        """绘制单幅图
        设置它的属性之类的

        Args:
            ax ([type]): [画纸对象]
            ds ([Dataset]): [数据]
            picture_dic ([dict]): [图片属性，标注的字典]
        """

        ## 设置地图的范围
        # lat = ds_fnl.lat
        # lon = ds_fnl.lon
        ax = self.ax
        # ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

        # box = [110, 116, 32, 37]
        # ax.set_extent(box, crs=ccrs.PlateCarree())
        # ax.set_extent(box)

        ## 设置图片标题
        self.create_map()
        ax.set_title(picture_dic['title'], fontsize=30)
        ax.set_title(picture_dic['pressure'], fontsize=20,loc='right')
        ax.set_title(picture_dic['initial_time'], fontsize=20,loc='left')
        # cf = self.draw_contourf_single(da, ax, dic)
        
        ## 绘制云顶亮温
        cs = self.draw_contourf(dict['q'], ax)

        # ccc = self.draw_contour(dict['geopt'], ax)
        cccc = self.draw_quiver(dict['u'], dict['v'])
        self.draw_station()
        return cs

class Draw():
    """绘图的类
    负责创建画纸, 一张图或多张图的区别
    保存图片
    """

    def draw_one_subplot(self,dict, initialtime):
        pass 
        fig = plt.figure(figsize=(12, 10), dpi=600)
        # proj = ccrs.PlateCarree()  # 创建坐标系
        # flnm_nc = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800/wrfout_d02_2021-07-18_00:00:00'
        # ncfile = Dataset(flnm_nc)
        # slp = getvar(ncfile, 'slp')
        # proj = get_cartopy(slp)
        ## 手动设置投影格式
        proj = ccrs.LambertConformal(
            central_longitude=120.0,
            central_latitude=32,
            standard_parallels=(30, 60),
            cutoff=-30,)
        # proj = ccrs.LambertConformal(central_longitude=120, central_latitude=33)
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        pr = Picture(ax, )
        
        tt = dict['q'].time
        titile = str(tt.dt.strftime('%Y-%m-%d %H').values)
        pressure = str(int(dict['u'].pressure.values))+' hPa'
        picture_dic = {'pressure':pressure, 'initial_time':initialtime, 'title':titile}
        
        cs = pr.draw_single(fig, dict, picture_dic)
        cb = fig.colorbar(
            cs,
            # cax=ax6,
            orientation='horizontal',
            # ticks=bounds,
            fraction=0.04,  # 色标大小
            pad=0.08,  # colorbar和图之间的距离
            # label='q (g/kg)',
        )
        font = {'size':15}
        cb.set_label('q (g/kg)', fontdict=font)
        figpath = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar_GDAS/'
        fig_name = str(tt.dt.strftime('%Y-%m-%d_%H').values)+"_"+pressure+'_GDAS_'+picture_dic['initial_time']+'do2'
        fig.savefig(figpath+fig_name)




def draw_single(t, ds_fnl, pre, initial_time):
    """画单个时次，单层的图

    Args:
        t ([type]): 时间
        ds_fnl ([type]): 多维数据
        pre ([type]): 气压层
    """
    title = str(t.dt.strftime('%Y-%m-%d %H').values)
    print(title)
    dict = {}
    ds_fnl.sel(pressure=pre, time=t)
    # dict['q'] = ds_fnl['q'].sel(time=t, pressure=pre).T*10**3
    dict['q'] = ds_fnl['q'].sel(time=t, pressure=pre)*10**3
    dict['height_agl'] = ds_fnl['height_agl'].sel(pressure=pre, time=t)
    dict['geopt'] = ds_fnl['geopt'].sel(pressure=pre, time=t)
    dict['u'] = ds_fnl['ua'].sel(pressure=pre, time=t)
    dict['v'] = ds_fnl['va'].sel(pressure=pre, time=t)
    dr = Draw()
    dr.draw_one_subplot(dict, initial_time)
    



def single_process_draw():
    """单进程绘图
    """
    fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar.nc'
    ds_fnl = xr.open_dataset(fl_wrf)
    pre = 500
    for t in ds_fnl.time:
        pass
        draw_single(t, ds_fnl, pre)

def multi_process_draw():
    """多进程绘图
    """
    # fl = '/mnt/zfm_18T/fengxiang/DATA/FY_TBB/TBB_FY2G_201607.nc'
    # ds = xr.open_dataset(fl)
    # fl_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar.nc'
    fl = '/mnt/zfm_18T/fengxiang/HeNan/Data/GDAS/'
    ti_list = ['1800', '1812', '1900', '1912']
    
    for i in ti_list:
        fl_wrf = fl+'YSU_'+i+'_upar_d02.nc'
        ds_fnl = xr.open_dataset(fl_wrf)
        # pre = 500
        # pre_list = [500, 700, 850]
        pre_list = [500]
        for pre in pre_list:
            pool = Pool(10)
            for t in ds_fnl.time:
                # draw_single(t, ds_fnl, pre)
                tr = pool.apply_async(draw_single, args=(t,ds_fnl,pre,i,))
            pool.close()
            pool.join()


if __name__ == "__main__":
    pass
    # main()
    # single_process_draw()
    multi_process_draw()



