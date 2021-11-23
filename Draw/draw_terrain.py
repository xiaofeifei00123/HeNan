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
import matplotlib.pyplot as plt
import numpy as np
import cmaps
from baobao.map import Map
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
def draw_contourf(rain, pic_dic):

    """rain[lon, lat, data],离散格点的DataArray数据

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    from nmc_met_graphics.plot import mapview
    mb = mapview.BaseMap()
    fig = plt.figure(figsize=[12, 10], dpi=600)
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.LambertConformal(central_latitude=34, central_longitude=113))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
    ax = fig.add_axes([0.1, 0.1, 0.85, 0.8], projection=ccrs.PlateCarree())
    # ax.plot(rain.lon.values, rain.lat.values, 'ko',ms=3,zorder=2, transform=ccrs.PlateCarree())
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        # 'extent':[110.5, 116, 32, 36.5],
        # 'extent':[111, 112, 33.5, 34.5],
        'extent':[109, 117, 31, 37],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
    colorlevel = np.arange(0, 2300, 100)
    # colorlevel = np.arange(-50, 50, 1)
    colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000',]#颜色列表
    # cmap = cmaps.precip3_16lev
    cmap = cmaps.MPL_terrain
    # cmap = cmaps.ViBlGrWhYeOrRe
    # cmap = cmaps.MPL_gist_yarg
    # cs = ax.tricontour(rain.lon, rain.lat, rain, levels=colorlevel, transform=ccrs.PlateCarree())
    # cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=colorlevel,colors=colordict, transform=ccrs.PlateCarree())
    cs = ax.contourf(rain.lon, 
                    rain.lat,
                    rain,
                    levels=colorlevel,
                    # colors=colordict,
                    cmap=cmap,
                    transform=ccrs.PlateCarree())
    # colorticks=[0.1,5,15,30.0,70,140]#雨量等级
    colorticks = colorlevel[1:-1][::4]
    
    cb = fig.colorbar(
        cs,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.07,  #  色标和子图间距离
    )

    station = {
        'ZhengZhou': {
            'abbreviation':'ZZ',
            'lat': 34.76,
            'lon': 113.65
        },
    }
    mp.add_station(ax, station, justice=True)
    ax.add_feature(cfeature.RIVERS, lw=2.5)
    # ax.add_feature(cfeature.LAKES, lw=1)
    
    
    

    # mb.cities(ax, city_type='base_station', color_style='black', 
    #             marker_size=5, font_size=16)
    # mb.gridlines()
    cb.ax.tick_params(labelsize=30)  # 设置色标标注的大小
    ax.set_title(pic_dic['title'], loc='left', fontsize=30)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_terrain/'
    fig_name = fig_path+pic_dic['title']
    fig.savefig(fig_name)

def draw_contourf_minus(rain, pic_dic):
    """地形高度的差

    Args:
        rain ([type]): [description]
    Example:
    da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_station.nc')
    da.max()
    rain = da.sel(time=slice('2021-07-20 00', '2021-07-20 12')).sum(dim='time')
    """
    from nmc_met_graphics.plot import mapview
    mb = mapview.BaseMap()
    fig = plt.figure(figsize=[12, 12], dpi=600)
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.LambertConformal(central_latitude=34, central_longitude=113))
    # ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
    ax = fig.add_axes([0.1, 0.1, 0.85, 0.8], projection=ccrs.PlateCarree())
    # ax.plot(rain.lon.values, rain.lat.values, 'ko',ms=3,zorder=2, transform=ccrs.PlateCarree())
    mp = Map()
    map_dic = {
        'proj':ccrs.PlateCarree(),
        # 'extent':[110.5, 116, 32, 36.5],
        # 'extent':[111, 112, 33.5, 34.5],
        'extent':[109, 117, 31, 37],
        'extent_interval_lat':1,
        'extent_interval_lon':1,
    }

    ax = mp.create_map(ax, map_dic)
    ax.set_extent(map_dic['extent'])

    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 700]#雨量等级
    # colorlevel = np.arange(0, 2300, 100)
    # colorlevel = np.arange(-60, 60+5, 5)
    colorlevel = np.arange(-160, 160+20, 20)
    colordict=['#F0F0F0','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000',]#颜色列表
    # cmap = cmaps.precip3_16lev
    # cmap = cmaps.MPL_terrain
    cmap = cmaps.ViBlGrWhYeOrRe
    # cmap = cmaps.MPL_gist_yarg
    # cs = ax.tricontour(rain.lon, rain.lat, rain, levels=colorlevel, transform=ccrs.PlateCarree())
    # cs = ax.tricontourf(rain.lon, rain.lat, rain, levels=colorlevel,colors=colordict, transform=ccrs.PlateCarree())
    cs = ax.contourf(rain.lon, 
                    rain.lat,
                    rain,
                    levels=colorlevel,
                    # colors=colordict,
                    cmap=cmap,
                    transform=ccrs.PlateCarree())
    # colorticks=[0.1,5,15,30.0,70,140]#雨量等级
    colorticks = colorlevel[1:-1][::2]
    
    cb = fig.colorbar(
        cs,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.05,  # 色标大小,相对于原图的大小
        pad=0.05,  #  色标和子图间距离
    )
    # mb.cities(ax, city_type='base_station', color_style='black', 
    #             marker_size=5, font_size=16)
    # mb.gridlines()
    cb.ax.tick_params(labelsize=30)  # 设置色标标注的大小
    ax.set_title(pic_dic['title'], loc='left', fontsize=30)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_terrain/'
    fig_name = fig_path+pic_dic['title']
    fig.savefig(fig_name)






if __name__ == '__main__':
    flnm_900m = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_900m/wrfout_d04_2021-07-19_01:00:00'
    flnm_90m = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d04_2021-07-19_01:00:00'
    flnm_900m_met = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_900m/geo_em.d03.nc'
    flnm_90m_met = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/geo_em.d03.nc'

    # h90 = get_hgt(flnm_90m)
    # h900 = get_hgt(flnm_900m)
    # h99 = h90-h900

    met_h90 = get_hgt_met(flnm_90m_met)
    met_h900 = get_hgt_met(flnm_900m_met)
    met_h99 = met_h90-met_h900


    # draw_tricontourf(h90, {'title':'90m'})
    # draw_tricontourf(h900, {'title':'900m'})
    # draw_tricontourf_minus(h99, {'title':'90m-900m'})

    draw_contourf(met_h90, {'title':'geo_90m'})
    draw_contourf(met_h900, {'title':'geo_900m'})
    draw_contourf_minus(met_h99, {'title':'geo_90m-900m'})