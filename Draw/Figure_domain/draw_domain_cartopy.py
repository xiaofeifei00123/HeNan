#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
read in namelist.wps , draw wrf domain and plot some station
-----------------------------------------
Time             :2021/03/28 17:28:59
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import cartopy.feature as cf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import geopandas
import cmaps
import re

from matplotlib.path import Path
import matplotlib.patches as patches

# %
# from matplotlib import rcParams

# config = {
#     "font.family": 'serif', # 衬线字体
#     "font.size": 12, # 相当于小四大小
#     "font.serif": ['SimSun'], # 宋体
#     "mathtext.fontset": 'stix', # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
#     'axes.unicode_minus': False # 处理负号，即-号
# }
# rcParams.update(config)

# %%
def draw_screen_poly(lats, lons):
    '''
    lats: 纬度列表
    lons: 经度列表
    purpose:  画区域直线
    '''
    x, y = lons, lats
    xy = list(zip(x, y))
    # print(xy)
    poly = plt.Polygon(xy, edgecolor="black", fc="none", lw=1, alpha=1)
    plt.gca().add_patch(poly)


def create_map(info):
    """创建一个包含青藏高原区域的Lambert投影的底图

    Returns:
        ax: 坐标图对象
    """
    ## --创建画图空间

    ref_lat = info['ref_lat']
    ref_lon = info['ref_lon']
    true_lat1 = info['true_lat1']
    true_lat2 = info['true_lat2']
    false_easting = (info['e_we'][0] - 1) / 2 * info['dx']
    false_northing = (info['e_sn'][0] - 1) / 2 * info['dy']

    proj_lambert = ccrs.LambertConformal(
        central_longitude=ref_lon,
        central_latitude=ref_lat,
        standard_parallels=(true_lat1, true_lat2),
        cutoff=-30,
        false_easting=false_easting,
        false_northing=false_northing,
    )
    # proj = ccrs.PlateCarree(central_longitude=ref_lon)  # 创建坐标系
    proj = ccrs.PlateCarree()  # 创建坐标系
    ## 创建坐标系
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 6*cm), dpi=600)  # 创建页面
    ax = fig.add_axes([0.12, 0.08, 0.85, 0.9], projection=proj_lambert)

    ## 读取青藏高原地形文件
    Province = cfeat.ShapelyFeature(
        Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp').geometries(),
        proj,
        edgecolor='black',
        lw=1.,
        linewidth=1.,
        facecolor='none',
        alpha=1.)

    Henan = cfeat.ShapelyFeature(
        Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp').geometries(),
        # Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp').geometries(),
        proj,
        edgecolor='black',
        lw=1.,
        linewidth=1.,
        facecolor='none',
        alpha=1.)
    city = cfeat.ShapelyFeature(
        Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp').geometries(),
        proj,
        edgecolor='black',
        lw=2.,
        linewidth=2.,
        facecolor='none',
        alpha=1.)
    ## 将青藏高原地形文件加到地图中区
    # ax.add_feature(Province, linewidth=0.5, zorder=2, alpha=0.3)

    ax.add_feature(Henan, linewidth=0.6, zorder=2)
    # ax.add_feature(city, linewidth=0.5, zorder=2)
    ax.coastlines(linestyle='--', linewidth=0.5, alpha=0.7)
    # ax.coastlines()
    # import cartopy.feature as cfeature
    # ax.add_feature(cfeature.BORDERS, linestyle=':')

    ## --设置网格属性, 不画默认的标签
    gl = ax.gridlines(draw_labels=True,
                      dms=True,
                      linestyle=":",
                      linewidth=0.2,
                      x_inline=False,
                      y_inline=False,
                      color='k',)
    # # gl=ax.gridlines(draw_labels=True,linestyle=":",linewidth=0.3 , auto_inline=True,x_inline=False, y_inline=False,color='k')

    ## 关闭上面和右边的经纬度显示
    gl.top_labels = False  #关闭上部经纬标签
    gl.right_labels = False
    ## 这个东西还挺重要的，对齐坐标用的
    gl.rotate_labels = None
    ## 坐标的范围
    # gl.xlocator = mticker.FixedLocator(np.arange(90, 140, 5))
    gl.xlocator = mticker.FixedLocator(np.arange(70, 150, 10))
    
    gl.ylocator = mticker.FixedLocator(np.arange(10, 50, 5))
    ## 坐标标签的大小
    gl.xlabel_style = {'size': 8}  #修改经纬度字体大小
    gl.ylabel_style = {'size': 8}
    ## 坐标标签样式
    gl.xformatter = LongitudeFormatter(degree_symbol="${^\circ}$")
    gl.yformatter = LatitudeFormatter(degree_symbol="${^\circ}$")

    ax.spines['geo'].set_linewidth(1.0)  #调节图片边框粗细

    ## 画d02
    ax.set_extent([0, false_easting * 2, 0, false_northing * 2],
                  crs=proj_lambert)

    # 标注d01, 这个位置需要根据经纬度手动调整
    ax.text(78, # 经度
            45,  # 纬度
            'd01',
            transform=ccrs.PlateCarree(),
            fontdict={
                'size': 10,
            })
    return ax


def get_information(flnm):
    """根据namelist.wps文件，获取地图的基本信息

    Args:
        flnm ([type]): [description]

    Returns:
        [type]: [description]
    """

    ## 设置正则表达式信息
    pattern = {}
    # pattern['dx'] = 'dx\s*=\s*\d*,'
    # pattern['dy'] = 'dy\s*=\s*\d*,'
    pattern['dx'] = 'dx\s*=\s*[1-9].?[1-9]*'
    pattern['dy'] = 'dy\s*=\s*[1-9].?[1-9]*'
    pattern['max_dom'] = 'max_dom\s*=\s*\d\s*,'
    pattern[
        'parent_grid_ratio'] = 'parent_grid_ratio\s*=\s*\d,\s*\d,\s*\d+,\s*\d,'
    pattern['j_parent_start'] = 'j_parent_start\s*=\s*\d,\s*\d*,\s*\d*,\s*\d*,'
    pattern['i_parent_start'] = 'i_parent_start\s*=\s*\d,\s*\d*,\s*\d*,\s*\d*,'
    pattern['e_sn'] = 'e_sn\s*=\s*\d*,\s*\d*,\s*\d*,\s*\d*'
    pattern['e_we'] = 'e_we\s*=\s*\d*,\s*\d*,\s*\d*,\s*\d*'
    pattern['ref_lat'] = 'ref_lat\s*=\s*\d*.?\d*,'
    pattern['ref_lon'] = 'ref_lon\s*=\s*\d*.?\d*,'
    pattern['true_lat1'] = 'truelat1\s*=\s*\d*.?\d*,'
    pattern['true_lat2'] = 'truelat2\s*=\s*\d*.?\d*,'

    f = open(flnm)
    fr = f.read()

    def get_var(var, pattern=pattern, fr=fr):
        """处理正则表达式得到的数据"""
        ff1 = re.search(pattern[var], fr, flags=0)
        # print(ff1)
        str_f1 = ff1.group(0)

        str1 = str_f1.replace('=', ',')
        aa = str1.split(',')
        bb = []
        for i in aa[1:]:
            if i != '':
                bb.append(i.strip())
        return bb

    dic_return = {}
    aa = get_var('parent_grid_ratio')

    var_list = [
        'dx',
        'dy',
        'max_dom',
        'parent_grid_ratio',
        'j_parent_start',
        'i_parent_start',
        'e_sn',
        'e_we',
        'ref_lat',
        'ref_lon',
        'true_lat1',
        'true_lat2',
    ]

    for i in var_list:
        aa = get_var(i)
        if i in [
                'parent_grid_ratio',
                'j_parent_start',
                'i_parent_start',
                'e_we',
                'e_sn',
        ]:
            bb = aa
            bb = [float(i) for i in bb]
        else:
            bb = float(aa[0])
        dic_return[i] = bb

    return dic_return


def draw_d02(info):
    """绘制domain2

    Args:
        info ([type]): [description]
    """
    max_dom = info['max_dom']
    dx = info['dx']
    dy = info['dy']
    i_parent_start = info['i_parent_start']
    j_parent_start = info['j_parent_start']
    parent_grid_ratio = info['parent_grid_ratio']
    e_we = info['e_we']
    e_sn = info['e_sn']

    if max_dom >= 2:
        ### domain 2
        # 4 corners 找到四个顶点和距离相关的坐标
        ll_lon = dx * (i_parent_start[1] - 1)
        ll_lat = dy * (j_parent_start[1] - 1)
        ur_lon = ll_lon + dx / parent_grid_ratio[1] * (e_we[1] - 1)
        ur_lat = ll_lat + dy / parent_grid_ratio[1] * (e_sn[1] - 1)

        lon = np.empty(4)
        lat = np.empty(4)

        lon[0], lat[0] = ll_lon, ll_lat  # lower left (ll)
        lon[1], lat[1] = ur_lon, ll_lat  # lower right (lr)
        lon[2], lat[2] = ur_lon, ur_lat  # upper right (ur)
        lon[3], lat[3] = ll_lon, ur_lat  # upper left (ul)

        draw_screen_poly(lat, lon)  # 画多边型

        ## 标注d02
        # plt.text(lon[0] * 1+100000, lat[0] * 1. - 225000, "d02", fontdict={'size':14})
        # plt.text(lon[2] * 1 - 440000,
        plt.text(lon[2] * 1 - 540000,
                #  lat[2] * 1. - 200000,
                 lat[2] * 1. - 300000,
                 "d02",
                 fontdict={'size': 10})

    if max_dom >= 3:
        ### domain 3
        ## 4 corners
        ll_lon += dx / parent_grid_ratio[1] * (i_parent_start[2] - 1)
        ll_lat += dy / parent_grid_ratio[1] * (j_parent_start[2] - 1)
        ur_lon = ll_lon + dx / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            e_we[2] - 1)
        ur_lat = ll_lat + dy / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            e_sn[2] - 1)

        ## ll
        lon[0], lat[0] = ll_lon, ll_lat
        ## lr
        lon[1], lat[1] = ur_lon, ll_lat
        ## ur
        lon[2], lat[2] = ur_lon, ur_lat
        ## ul
        lon[3], lat[3] = ll_lon, ur_lat

        draw_screen_poly(lat, lon)

        ## 标注d03
        # plt.text(lon[0] * 1+100000, lat[0] * 1. - 225000, "d02", fontdict={'size':14})
        plt.text(lon[3] * 1+50000 ,
                #  lat[3] * 1. - 200000,
                 lat[3] * 1. - 300000,
                 "d03",
                 fontdict={'size': 10})
        
    if max_dom >= 4:

        ### domain 4
        ## 4 corners
        ll_lon += dx / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            i_parent_start[3] - 1)
        ll_lat += dy / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            j_parent_start[3] - 1)
        ur_lon = ll_lon + dx / parent_grid_ratio[1] / parent_grid_ratio[
            2] / parent_grid_ratio[3] * (e_we[3] - 1)
        ur_lat = ll_lat + dy / parent_grid_ratio[1] / parent_grid_ratio[
            2] / parent_grid_ratio[3] * (e_sn[3] - 1)

        ## ll
        lon[0], lat[0] = ll_lon, ll_lat
        ## lr
        lon[1], lat[1] = ur_lon, ll_lat
        ## ur
        lon[2], lat[2] = ur_lon, ur_lat
        ## ul
        lon[3], lat[3] = ll_lon, ur_lat
        draw_screen_poly(lat, lon)

    if max_dom >= 5:

        ### domain 5
        ## 4 corners
        ll_lon += dx / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            i_parent_start[4] - 1)
        ll_lat += dy / parent_grid_ratio[1] / parent_grid_ratio[2] * (
            j_parent_start[4] - 1)
        ur_lon = ll_lon + dx / parent_grid_ratio[1] / parent_grid_ratio[
            2] / parent_grid_ratio[3] * (e_we[4] - 1)
        ur_lat = ll_lat + dy / parent_grid_ratio[1] / parent_grid_ratio[
            2] / parent_grid_ratio[3] * (e_sn[4] - 1)

        ## ll
        lon[0], lat[0] = ll_lon, ll_lat
        ## lr
        lon[1], lat[1] = ur_lon, ll_lat
        ## ur
        lon[2], lat[2] = ur_lon, ur_lat
        ## ul
        lon[3], lat[3] = ll_lon, ur_lat
        draw_screen_poly(lat, lon)

def draw_station(ax):
    """画站点
    """
    station = {
        'ZhengZhou': {
            'lat': 34.76,
            'lon': 113.65,
            'abbreviation':'郑州'
        },
        'YanHua': {
            'lat': 22.4,
            'lon': 132.5,
            'abbreviation':'烟花'
        },
        'ChaPaka': {
            'lat': 21.1,
            'lon': 112.9,
            'abbreviation':'查帕卡'
        },
    }
    values = station.values()
    # station_name = list(station.keys())
    # station_name = list(station[])
    # print(type(station_name[0]))
    # print(station_name[0])
    x = []
    y = []
    station_name = []
    for i in values:
        y.append(float(i['lat']))
        x.append(float(i['lon']))
        station_name.append(i['abbreviation'])

    ## 标记出站点
    ax.scatter(x,
               y,
               color='black',
               transform=ccrs.PlateCarree(),
               linewidth=0.2,
               s=12)
    ## 给站点加注释
    for i in range(len(x)):
    # for i,j in zip(len(x), values):
        ax.text(x[i] - 1.,
                 y[i] + 0.8,
                 station_name[i],
                 transform=ccrs.PlateCarree(),
                 fontdict={
                     'size': 8,
                 })
def draw():
    pass
    file_folder = "/mnt/zfm_18T/fengxiang/HeNan/Draw/Figure_domain/"
    file_name = "namelist.wps_LES"
    # file_name = "namelist.wps"
    flnm = file_folder + file_name
    print(flnm)

    info = get_information(flnm)  # 获取namelist.wps文件信息
    # print(info['ref_lat'])
    ax = create_map(info)  # 在domain1区域内，添加地理信息，创建底图
    print("创建地图完毕")
    draw_d02(info)  # 绘制domain2区域
    # print("绘制完毕")

    draw_station(ax)
    # print("标注站点完毕")
    fig_name = file_folder+'domain_LES.png'
    plt.savefig(fig_name)

if __name__ == '__main__':
    draw()

# %%
