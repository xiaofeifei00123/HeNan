#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取micaps高空形势场数据
使用nwc(国家气象中心的python包)
绘制高空形势场
-----------------------------------------
Time             :2021/10/03 14:33:01
Author           :Forxd
Version          :1.0
'''


# %%
from time import strftime
from cartopy.crs import Projection
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import xarray as xr
import meteva.base as meb
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# from metdig.io.cassandra import get_obs_station

# %%
# def get_analaysis(dic):
def get_analysis(dic={'var':'height', 'level':'500', 'time':pd.Timestamp('2021-07-20 0800') }):
    """读micaps 14类数据

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/HGT/500/20210720080000.000'
    print(dic['time'].strftime('%Y-%m-%d %H%M'), dic['level'])
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/MANUAL_ANALYSIS/'
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/'
    if dic['var'] == 'height':
        flnm = path+'HGT/'+dic['level']+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')
        # print(flnm)
    elif dic['var'] == 'temp':
        flnm = path+'TMP/'+dic['level']+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')

    # print(flnm)
    df = read_micaps_14(flnm)
    dc = df['lines']
    loc = dc['line_xyz']
    label = dc['line_label']
    label_loc = dc['line_label_xyz']

    num = len(loc)
    ## 处理line
    line_list = []
    for i in range(num):
        aa = loc[i]
        bb = label[i]
        cc = label_loc[i]
        ## 处理line
        df1 = pd.DataFrame(aa[:,0:2], columns=['lon', 'lat'])
        ## 处理label
        df1[dic['var']] = bb
        df1['label_lon'] = cc[0][0]
        df1['label_lat'] = cc[0][1]
        line_list.append(df1)
    return line_list

# dic = {
#     'var':'temp',
#     'level':'500',
#     'time':pd.Timestamp('2021-07-20 0800')
# }
        
# aa = get_analysis(dic)


# %%
def get_plot(dic):
    """读取micaps 2类数据，高空填图数据"""
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/'
    flnm = path+dic['level']+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/500/20210720080000.000'
    df = read_micaps_2(flnm)
    sta = meb.sta_data(df, columns = [
                'id', 'lon', 'lat', 'alt', 'grade', 'height', 'temperature', 'dewpoint',
                'wind_angle', 'wind_speed', 'time', 'level'])

    return sta




# %%
def interp_metpy(sta):
    """
    站点插值到格点
    反距离权重插值

    Args:
        sta (DataFrame): [lon,lat,height]

    Returns:
        [type]: [description]
    """
    # h = sta['temperature']
    import metpy.interpolate as interp
    h = sta['height']

    lon = sta['lon']
    lat = sta['lat']
    # x,y,z = interp.interpolate_to_grid(lon, lat, h, 'barnes', hres=0.5, minimum_neighbors=2)
    # grid1 = meb.grid([60,150,0.25],[10,60,0.25])
    # x0, x1 = 69.05, 150.1
    # y0, y1 = 0, 55.1 

    x0, x1 = 0, 160
    y0, y1 = 0, 80
    # res = 1 / 32.0
    # res = 1 / 10  # 0.1°
    res = 1
    mx, my = np.meshgrid(np.arange(x0, x1, res),
                            np.arange(y0, y1, res),
                            indexing="ij")
    # mx
    # grd2 = meb.interp_sg_idw(stb, grid1)
    # z = interp.inverse_distance_to_grid(lon, lat, h, mx, my, r=10)
    # x,y,z = interp.remove_nan_observations(lon,lat, h)

    z = interp.inverse_distance_to_grid(lon, lat, h, mx, my, r=5, min_neighbors=1)
    return mx,my,z



def interp_pyinterp():
    import pyinterp
    import pyinterp.tests
    import pyinterp.backends.xarray

    lons = sta['lon']
    lats = sta['lat']
    h = sta['height']
    # lons = lon
    # lats = lat
    # da = ds['tbb']
    # tt = ds.time
    # da = ds['tbb'].isel(time=0)
    mesh = pyinterp.RTree()
    mesh.packing(
        np.vstack((lons.values, lats.values)).T,
        h)
    x0, x1 = 69.05, 150.1
    # y0, y1 = -40, 40
    y0, y1 = 0, 55.1 
    # res = 1 / 32.0
    # res = 1 / 10  # 0.1°
    res = 1
    mx, my = np.meshgrid(np.arange(x0, x1, res),
                            np.arange(y0, y1, res),
                            indexing="ij")
    idw_eta, neighbors = mesh.inverse_distance_weighting(
        np.vstack((mx.flatten(), my.flatten())).T,
        within=True,  # Extrapolation is forbidden
        # radius=55000,  # In a radius of 5.5 Km
        radius=5000000,  # In a radius of 5.5 Km
        k=8,  # We are looking for at most 8 neighbours
        num_threads=0)
    idw_eta = idw_eta.reshape(mx.shape)

    lon = np.arange(x0,x1,res)
    lat = np.arange(y0,y1,res)

def draw_station(ax):
    pass
    # station = station_dic
    station = {
        'ZhengZhou': {
            'abbreviation':'郑州',
            'lat': 34.75,
            'lon': 113.62
        },
    }
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
    ax.scatter(x,
                y,
            #    color='black',
                color='green',
                transform=ccrs.PlateCarree(),
                alpha=1.,
                linewidth=5,
                s=20,
                marker='o',
                )
    # 给站点加注释
    # for i in range(len(x)):
    #     # print(x[i])
    #     ax.text(x[i]-0.2,
    #             y[i] + 0.2,
    #             station_name[i],
    #             transform=ccrs.PlateCarree(),
    #             alpha=1.,
    #             fontdict={
    #             'size': 28,
    #     })

def add_ticks(ax,):
    """添加坐标标签"""
    # ax.set_yticks(np.arange(10, 60 + 1, 10))
    # ax.set_xticks(np.arange(70, 140 + 1, 10))
    ax.set_yticks(np.arange(10, 70 + 1, 10))
    # ax.set_xticks(np.arange(60, 140 + 1, 10))
    ax.set_xticks(np.arange(60, 150 + 1, 20))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.tick_params(which='major',length=10,width=2.0) # 控制标签大小 
    ax.tick_params(which='minor',length=5,width=1.0)  #,colors='b')
    ax.tick_params(axis='both', labelsize=25, direction='out')

def draw_south_sea(fig,):
    pass
    ax2 = fig.add_axes([0.102, 0.138, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2.set_extent([105.8, 122,0,25])
    # ax2.set_extent([60, 130,0,75])
    # ax2.add_feature(cfeature.LAKES.with_scale('50m'))
    ax2.add_geometries(Reader('/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='black',linewidth=0.8)
    # ax2.add_geometries(Reader(r'F:/Rpython/lp27/data/china1.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='k',linewidth=0.2)
    # ax2.add_geometries(Reader(r'F:/Rpython/lp27/data/china2.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='k',linewidth=0.1)
    # ax2.add_geometries(Reader(r'F:/Rpython/lp27/data/ne_10m_land.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='k',linewidth=0.2)
    # ax2.add_geometries(Reader(r'F:/Rpython/lp27/data/ne_50m_lakes.shp').geometries(),ccrs.PlateCarree(),facecolor='none',edgecolor='k',linewidth=0.2)
    # c21=ax2.contourf(lon0,lat0,z,np.arange(-30,31,5),extend='both',transform=ccrs.PlateCarree(),cmap='gist_rainbow')
    # clip=maskout31.shp2clip(c21,ax2,'F:/Rpython/lp27/data/china0')
    # ax0 = plt.gca()   #获取边框
    # ax0.outline_patch.set_linewidth(0.5)    #修改边框粗细
    # plt.savefig('F:/Rpython/lp28/plot26.png',dpi=600)

def draw(hgt_list, tmp_list, ddf, dic):
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地及周边')
    ax.set_extent([60, 150,10,70])
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    tt = (dic['time']+pd.Timedelta('-8H')).strftime('%Y-%m-%d %H%M')
    ax.set_title(tt, loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=25)

    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    # mb.set_extent('中国陆地')
    # mb.southsea(zoom=0.3, loc='left_bottom')

    for df in hgt_list:
        # print(df)
        y = df['label_lat'][0]
        x = df['label_lon'][0]
        if x>73 and x<136 and y>15 and y<56:
            if df['height'][0] in ['140', '144', '148',  '584', '588', '1240', '1248', '1256', '1264']:
                ax.text(df['label_lon'][0], df['label_lat'][0], df['height'][0], fontsize=15,color='blue', transform=ccrs.PlateCarree())
        c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='blue')
        
    for df in tmp_list:
        # print(df)
        y = df['label_lat'][0]
        x = df['label_lon'][0]
        if x>73 and x<136 and y>15 and y<56:
            ax.text(df['label_lon'][0], df['label_lat'][0], df['temp'][0], fontsize=15,color='red', transform=ccrs.PlateCarree())
        c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='red', linewidth=0.8, alpha=0.8)
    # gl = mb.gridlines(font_size=25, alpha=0)
    add_ticks(ax)
    draw_station(ax)

    x = ddf['lon']
    y = ddf['lat']
    u = -1*np.sin(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    v = -1*np.cos(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    # ax.barbs(x,y,u,v)
    # Q = ax.quiver(x, y, u,v,units='inches',scale=30,pivot='tip',width=0.04, transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u,v,headwidth=3, headlength=4,width=0.03, units='inches',scale=,pivot='tip', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.8, Y=0.05, U=10, label=r'$10\ m/s$', labelpos='E',coordinates='figure', fontproperties={'size':25})   # 设置参考风矢
    # mb.southsea(zoom=0.3, loc='left_bottom')
    draw_south_sea(fig)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar_micaps/'
    fig_name = fig_path+dic['level']+'_'+(dic['time']+pd.Timedelta('-8H')).strftime('%Y%m%d%H')
    # fig.savefig('test.png')
    fig.savefig(fig_name)


def draw_all(t, level='500'):
    # t = pd.Timestamp('2021-07-19 2000')
    dic_t = {
        'var':'temp',
        'level':level,
        'time':t
    }
    dic_h = {
        'var':'height',
        'level':level,
        'time':t
    }
    hgt_list = get_analysis(dic_h)
    tmp_list = get_analysis(dic_t)
    ddf = get_plot(dic_h)
    draw(hgt_list, tmp_list, ddf, dic_t)

### 单个时次，单个层次测试
# t = pd.Timestamp('2021-07-19 2000')
# draw_all(t, '500')
### 测试结束
if __name__ == '__main__':
    pass
    ttt = pd.date_range(start='2021-07-18 08', end='2021-07-20 20', freq='12H')
    # for level in ['200', '500', '850']:
    for level in ['500', '700', '850']:   # 人工分析的只有这几个高度的
        for t in ttt:
            draw_all(t, level)
    



# %%