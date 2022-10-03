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
# from baobao.get_cmap

# %%
class Draw(object):
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

        from baobao.get_cmap import select_cmap
        self.colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250, 400,600, 1000]#雨量等级
        self.colordict = select_cmap('rain9')
        

        self.colorticks = self.colorlevel[1:-1]
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
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'
        self.station = {
                'ZhengZhou': {
                    'abbreviation':'郑州',
                    'lat': 34.76,
                    'lon': 113.65
                },
                'NanYang': {
                    'abbreviation':'南阳',
                    'lat': 33.1,
                    'lon': 112.49,
                },
                'LuShi': {
                    'abbreviation':'卢氏',
                    'lat': 34.08,
                    'lon': 111.07,
                },
            }


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
        mp.add_station(ax, {'ZhengZhou':self.station['ZhengZhou']}, justice=False, delx=-0.2, marker='o', ssize=15)
        mp.add_station(ax, {'NanYang':self.station['NanYang']}, justice=False, delx=-0.2, marker='s', ssize=15)
        mp.add_station(ax, {'LuShi':self.station['LuShi']}, justice=False, delx=-0.2, marker='x', ssize=15)

        
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
        mp.add_station(ax, {'ZhengZhou':self.station['ZhengZhou']}, justice=False, delx=-0.2, marker='o', ssize=15)
        mp.add_station(ax, {'NanYang':self.station['NanYang']}, justice=False, delx=-0.2, marker='s', ssize=15)
        mp.add_station(ax, {'LuShi':self.station['LuShi']}, justice=False, delx=-0.2, marker='x', ssize=15)

        ax.set_extent(self.map_dic['extent'])
        cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=self.colorlevel,colors=self.colordict, transform=ccrs.PlateCarree())

        rain_max = rain.max()        
        ax.set_title('Max = %s'%(rain_max.values.round(1)), fontsize=10,loc='right')
        # ax.set_title('2021-07 20/00--21/00', fontsize=35,)
        # ax.set_title('OBS', fontsize=10,loc='right')
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
        da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        # da = da.sel(time=slice('2021-07-19 13', '2021-07-21 00'))
        da = da.sum(dim='time') 
        return da

    def onemodel(self, model='gwd3'):
        pass

        flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/Typhoon/'+model+'/'+'rain.nc'
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
        da = xr.open_dataarray(flnm)
        da = da.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
        da = da.sum(dim='time') 
        return da


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

if __name__ == '__main__':

    def get_dr():
        cm = round(1/2.54, 2)
        proj = ccrs.PlateCarree()  # 创建坐标系
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
        ax = fig.add_axes([0.1,0.08,0.85,0.85], projection=proj)
        dr = Draw(fig, ax)
        return dr

    ## 画观测降水
    dr = get_dr()
    gd = GetData()
    da = gd.obs()
    cf = dr.draw_tricontourf(da)    
    cb = dr.fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    # cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    labels = list(map(lambda x: str(x) if x<1 else str(int(x)), dr.colorticks))  # 将colorbar的标签变为字符串
    cb.set_ticklabels(labels)
    fig_name = 'OBS'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    dr.fig.savefig(fig_path+fig_name)

    ## 画EC降水
    dr = get_dr()
    gd = GetData()
    da = gd.EC()
    cf = dr.draw_single(da)    
    cb = dr.fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=dr.colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    fig_name = 'EC'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    dr.fig.savefig(fig_path+fig_name)

    ## 画模式降水
    for model in ['gwd3', 'gwd0']:
        dr = get_dr()  # 画图的对象
        gd = GetData()  # 数据的对象
        da = gd.onemodel(model)
        cf = dr.draw_single(da)    
        cb = dr.fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=dr.colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.1,  #  色标和子图间距离
            )
        cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
        dr.ax.set_title(model, fontsize=10,loc='right')
        fig_name = model
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
        dr.fig.savefig(fig_path+fig_name)
