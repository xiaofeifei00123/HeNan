#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制地形高度
-----------------------------------------
Time             :2021/11/11 19:34:33
Author           :Forxd
Version          :1.0
'''
# %%
import wrf
import xarray as xr
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import cmaps
from baobao.map import Map

from draw_rain_distribution_24h import Draw
# %%

def get_hgt_met(flnm):
    # 从met文件中获取海拔高度数据
    ds = xr.open_dataset(flnm)
    hgt_m = ds['HGT_M'].squeeze()
    lat = ds['XLAT_M'].squeeze()
    lon = ds['XLONG_M'].squeeze()

    hgt = hgt_m.assign_coords({'lat':(['south_north', 'west_east'],lat.values),
                        'lon':(['south_north', 'west_east'],lon.values)})
    return hgt


def get_hgt(flnm):
    """从wrfout文件中获取海拔高度数据
    """
    wrfnc = nc.Dataset(flnm)
    h = wrf.getvar(wrfnc, 'ter')
    h = h.rename({'XLAT':'lat', 'XLONG':'lon'})
    return h
# %%
def draw_contourf_lambert(rain, pic_dic):

    """rain[lon, lat, data],离散格点的DataArray数据
    使用lambert投影画这个地形图

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    # from nmc_met_graphics.plot import mapview
    # mb = mapview.BaseMap()
    cm = round(1/2.54, 2)
    fig = plt.figure(figsize=[8*cm, 8*cm], dpi=300)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.LambertConformal(central_latitude=34, central_longitude=113))
    mp = Map()

    # colorlevel = np.arange(0, 2300, 100)
    colorlevel = np.arange(0, 4000, 200)
    cmap = cmaps.MPL_terrain
    cs = ax.contourf(rain.lon, 
                    rain.lat,
                    rain,
                    levels=colorlevel,
                    # colors=colordict,
                    cmap=cmap,
                    transform=ccrs.PlateCarree())
    colorticks = colorlevel[1:-1][::4]
    
    cb = fig.colorbar(
        cs,
        # cax=ax6,
        # orientation='horizontal',
        orientation='vertical',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.01,  #  色标和子图间距离
    )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小

    station = {
        'ZhengZhou': {
            'abbreviation':'郑州',
            'lat': 34.76,
            'lon': 113.65
        },
    }
    proj = ccrs.PlateCarree()  # 创建坐标系
    mp.add_station(ax, station, justice=True, fontsize=10, ssize=8, dely=0.2)

    Henan = cfeat.ShapelyFeature(
        Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp').geometries(),
        # Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp').geometries(),
        proj,
        edgecolor='black',
        lw=1.,
        linewidth=1.,
        facecolor='none',
        alpha=1.)
    ax.add_feature(cfeature.RIVERS, lw=1)
    ax.add_feature(Henan, linewidth=1, zorder=2)
    ax.add_feature(cfeature.LAKES, lw=1)
    
    gl = ax.gridlines(draw_labels=True,
                      dms=True,
                      linestyle=":",
                      linewidth=0.2,
                      x_inline=False,
                      y_inline=False,
                      color='k',)
    
    gl.top_labels = False  #关闭上部经纬标签
    gl.right_labels = False
    ## 这个东西还挺重要的，对齐坐标用的
    gl.rotate_labels = None
    ## 坐标的范围
    gl.xlocator = mticker.FixedLocator(np.arange(100, 120, 2))
    gl.ylocator = mticker.FixedLocator(np.arange(20, 40, 2))
    ## 坐标标签的大小
    gl.xlabel_style = {'size': 10}  #修改经纬度字体大小
    gl.ylabel_style = {'size': 10}
    ## 坐标标签样式
    gl.xformatter = LongitudeFormatter(degree_symbol="${^\circ}$")
    gl.yformatter = LatitudeFormatter(degree_symbol="${^\circ}$")
    ax.spines['geo'].set_linewidth(1.0)  #调节图片边框粗细
    
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_terrain/'
    fig_name = fig_path+pic_dic['title']+'_lambert'
    fig.savefig(fig_name, bbox_inches = 'tight')


def draw_contourf_latlon(rain, pic_dic):

    """rain[lon, lat, data],离散格点的DataArray数据
    使用等经纬线投影画这地形图

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    # from nmc_met_graphics.plot import mapview
    # mb = mapview.BaseMap()
    # cm = round(1/2.54, 2)
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    proj = ccrs.PlateCarree()  # 创建坐标系
    ax = fig.add_axes([0.105,0.1,0.85,0.85], projection=ccrs.PlateCarree())
    

    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        'extent':[110.5, 116, 32, 36.5],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])
    


    colorlevel = np.arange(0, 2300, 100)
    # colorlevel = np.arange(0, 4000, 200)
    cmap = cmaps.MPL_terrain
    cs = ax.contourf(rain.lon, 
                    rain.lat,
                    rain,
                    levels=colorlevel,
                    # colors=colordict,
                    cmap=cmap,
                    transform=ccrs.PlateCarree())
    colorticks = colorlevel[1:-1][::4]

    cb = fig.colorbar(
        cs,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
    )

    station = {
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
    mp.add_station(ax, station, justice=True)

    dr = Draw(fig, ax)
    # areaC = {
    #     'lat1':33.5,
    #     'lat2':36,
    #     'lon1':112.2,
    #     'lon2':114.8,
    # }        
    # dr.add_patch(areaC, ax)
    ## 画斜线
    # ax.plot(np.linspace(dr.cross_start[0], dr.cross_end[0], 10), np.linspace(dr.cross_start[1], dr.cross_end[1], 10), color='black')
    # mp.add_station(ax, dr.station, justice=True, ssize=30)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_terrain/'
    fig_name = fig_path+pic_dic['title']+'_latlon'
    # fig.savefig(fig_name, bbox_inches = 'tight')
    fig.savefig(fig_name)


if __name__ == '__main__':
    flnm_90m_met = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/geo_em.d04.nc'
    # flnm_90m_met = '/home/fengxiang/HeNan/Data/GWD/d03/new_modify/MET/met_em.d03.2021-07-20_06:00:00.nc'
    met_h90 = get_hgt_met(flnm_90m_met)
    # draw_contourf_lambert(met_h90, {'title':'geo_90m'})
    draw_contourf_latlon(met_h90, {'title':'met_90m'})