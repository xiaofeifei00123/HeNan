#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
不同模式d04区域， 850hPa水汽通量的比较
各模式的图
各模式的差值的图
绘制的要素:
    850hPa
    水汽通量(填色)
    风场(矢量)
使用的数据：
    模式预报插值的格点数据
-----------------------------------------
Time             :2021/10/28 
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
import matplotlib as mpl
import cartopy.crs as ccrs
import cmaps
from wrf import get_cartopy, smooth2d, getvar
import metpy.calc as ca
import metpy
from metpy.units import units
from pyproj import CRS

# from metdig.io.cassandra import get_obs_station


### 测试
# %%
flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/high_resolution_high_hgt_upar_d04_latlon.nc'
ds = xr.open_dataset(flnm_wrf)
ds = ds.isel(time=0).isel(pressure=2)
# da = ds['q'].sel(time='2021-07-20 00').sel(pressure=850)
# da
# %%
lon = ds.lon
lat = ds.lat
# lon
# lat
u = ds.u*units('m/s')
v = ds.v*units('m/s')
q = ds.q*

# ca.divergence(u, v)
# crs = CRS.from_proj4("+proj=latlon")
# lon
dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
dx
# dx
# %%
# ds.lon
# dx.shape
# u
ca.divergence(u=u, v=v, dx=dx, dy=dy)

# dx
# dy.max()
# help(ca.lat_lon_grid_deltas)
# lon
# help(ca.divergence)






# %%
def caculate_data(dds):
    """ 抽取并计算出需要的数据
    u, v, qu, qv

    Args:
        dds ([type]): [description]

    Returns:
        [type]: [description]
    """
    dict = {}
    # u = dds.metpy.parse_cf('u')
    # v = dds.metpy.parse_cf('v')
    # q = dds.metpy.parse_cf('q')*10**3/9.8
    u = dds['u']
    v = dds['v']
    q = dds['q']*10**3
    qu = q*u/9.8
    qv = q*v/9.8
    dict['q'] = q
    dict['qu'] = qu
    dict['qv'] = qv
    dict['qf'] = xr.ufuncs.sqrt(qu**2+qv**2)
    dict['u'] = u
    dict['v'] = v
    # div = ca.divergence(u,v)
    # div_q = ca.divergence(qu, qv)
    # dict['div'] = div
    # dict['div_q'] = div_q
    return dict

def get_data(dic):
    """获得模式和观测
    水汽和风场数据

    Args:
        dic ([dict]): 时间和层次信息

    Returns:
        [type]: [description]
    """
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_upar_latlon1.nc' 
    ds_obs = xr.open_dataset(flnm_obs)
    # flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GDAS/YSU_1900_upar_d02_latlon.nc'
    flnm_wrf = dic['flnm']
    ds_wrf = xr.open_dataset(flnm_wrf)
    # ds_wrf.time
    t = dic['time']
    level = dic['level']
    # ds1 = ds_obs.sel(time=t, pressure=level)
    ds2 = ds_wrf.sel(time=t, pressure=level)
    # dic1 = caculate_data(ds1)
    dic2 = caculate_data(ds2)
    # qf_dif = dic2['qf'] - dic1['qf']  # 水汽通量大小差
    # qu_dif = dic2['qu'] - dic1['qu']  # 水汽通量矢量差
    # qv_dif = dic2['qv'] - dic1['qv']
    # u_dif = dic2['u'] - dic1['u']  # 风场差
    # v_dif = dic2['v'] - dic1['v']

    dic_return = {
        # 'u_obs':dic1['u'],
        # 'v_obs':dic1['v'],
        'u_model':dic2['u'],
        'v_model':dic2['v'],
        # 'qf_obs':dic1['qf'],
        'qf_model':dic2['qf'],
        # 'qf_dif':qf_dif,
        # 'qu_dif':qu_dif,
        # 'qv_dif':qv_dif,
        # 'v_dif':v_dif,
        # 'u_dif':u_dif,
    }

    return dic_return 



def get_analysis(dic={'var':'height', 'level':'500', 'time':pd.Timestamp('2021-07-20 0800') }):
    """读micaps 14类数据

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/HGT/500/20210720080000.000'
    # print(dic['time'].strftime('%Y-%m-%d %H%M'), dic['level'])
    # path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/MANUAL_ANALYSIS/'
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/'
    if dic['var'] == 'height':
        flnm = path+'HGT/'+str(dic['level'])+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')
    elif dic['var'] == 'temp':
        flnm = path+'TMP/'+str(dic['level'])+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')

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

def get_plot(dic):
    """读取micaps 2类数据，高空填图数据"""
    path = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/'
    flnm = path+str(dic['level'])+'/'+dic['time'].strftime('%Y%m%d%H%M%S.000')
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/500/20210720080000.000'
    df = read_micaps_2(flnm)
    sta = meb.sta_data(df, columns = [
                'id', 'lon', 'lat', 'alt', 'grade', 'height', 'temperature', 'dewpoint',
                'wind_angle', 'wind_speed', 'time', 'level'])

    return sta

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
               color='black',
                # color='green',
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
    # mb.set_extent([110, 116, 32, 36])
    ax.set_yticks(np.arange(32, 36 + 1, 1))
    ax.set_xticks(np.arange(110, 116 + 1, 1))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.tick_params(which='major',length=10,width=2.0) # 控制标签大小 
    ax.tick_params(which='minor',length=5,width=1.0)  #,colors='b')
    ax.tick_params(axis='both', labelsize=25, direction='out')

def draw_south_sea(fig,):
    pass
    # ax2 = fig.add_axes([0.102, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2 = fig.add_axes([0.798, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    # ax2.set_extent([105.8, 122,0,25])
    # ax2.set_extent([110, 116, 32, 37], crs=ccrs.PlateCarree())
    # ax2.set_extent([110, 116, 32, 37])
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
# def draw_contourf(ax, u):
#     ax.contourf(u.lon, u.lat, u)
#     pass

def draw_contourf(ax, da):
    """在地图上绘制填色图
    """
    x = da.lon
    y = da.lat
    # levels = [205, 210, 215, 220, 225, 235, 245]  # 需要画出的等值线
    # contour_levels = np.arange(30,61,5)
    # contour_levels = np.arange(-24,24+1,6)
    # contour_levels = [-30, -24, -18, -12,-6,-3, 3, 6, 12, 18, 24, 30]
    # contour_levels = [-30,-22,-18, -14, -10,-6, -2, 2,  6, 10, 14, 18, 22, 30]
    # contour_levels = [-30,-22,-18, -14, -10,-6, -2, 2,  6, 10, 14, 18, 22, 30]
    # contour_levels = [-15,-11,-9, -7, -5,-3, -1, 1,  3, 5, 7, 9,11,15]
    # contour_levels = [-16, -12, -8, -6, -4,-2,-1, 1, 2, 4, 6, 8, 12, 16]
    # contour_levels = [-0.020, -0.016, -0.012, -0.008, -0.004, -0.002, 0.002, 0.004, 0.008, 0.012, 0.016, 0.020]
    # contour_levels = [-50, -25, -16, -8, -4, 4, 16, 20, 25, 50]
    # contour_levels = [0, 10, 12, 14, 16, 18, 20, 22, 24, 30, 100]
    contour_levels = [0, 10, 12, 14, 16, 18, 20, 22,24,30, 100]
    # contour_levels = np.arange(0, 17, 1)
    color_li = ['white', '#6CA6CD', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','#EE7600','red']
    # color_li = ['white', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','red']
    
    # contour_levels = [2, 4, 6, 8, 10, 12, 14, 16]
    # contour_levels = np.arange(0, 17, 1)
    
    # colormap = cmaps.precip3_16lev
    # colormap = cmaps.ViBlGrWhYeOrRe
    # colormap = cmaps.precip3_16lev_r  # 反转色标
    crx = ax.contourf(x,
                        y,
                        da,
                        # cmap=colormap,
                        colors=color_li,
                    #   norm=norm,
                    #   extend='both',
                    #   extend='max',
                        levels=contour_levels,
                        transform=ccrs.PlateCarree())
    # crx.cmap.set_over('red')
    # crx.cmap.set_under('white')
    return crx

def draw_contour(ax, da, **kw):
    """在地图上绘制等温线
    """
    x = da.lon
    y = da.lat
    # level_contour = np.arange(0, 12.1, 0.5)
    # level_contour = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    # level_contour = [5, 7, 9, 11, 13]
    level_contour = [6, 8, 10, 12, 14, 16] 
    # da = smooth2d(field=da, passes=8)
    
    da = smooth2d(field=da, passes=2)
    crx = ax.contour(x,
                        y,
                        da,
                        colors = 'red',
                        levels=level_contour,
                        linestyles = 'dashed',
                        transform=ccrs.PlateCarree())

    ax.clabel(crx,inline=1, fontsize=20, colors='red') # 等值线的标注
    return crx

# def draw_contour(ax, da):
#     ax.contour(da.lon, da.lat, da)
def draw_quiver(u,v, ax):
    '''
    绘制风矢图
    '''
    u = u[::3,::3]
    v = v[::3,::3]
    # y = u.coords['lat']
    y = u.lat.values
    x = u.lon.values
    # print(y)
    # print(type(u))
    # print(type(y))
    # x = u.coords['lon']


    # ax = self.ax
    # Q = ax.quiver(x,y,u,v,units='inches',scale=18,pivot='middle')  # 绘制风矢
    # Q = ax.quiver(x, y, u,v,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=20,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=25,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    # qk = ax.quiverkey(Q, X=0.55, Y=0.15, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=0.55, Y=0.12, U=10, label=r'($10 ,g \cdot cm^{-1} \cdot hPa^{-1} \cdot s^{-1}$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=1.55, Y=0.05, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢

def draw(qdif, qu,qv, dic):
    """画

    Args:
        hgt_list ([type]): [description]
        tmp_list ([type]): [description]
        ddf ([type]): [description]
        dic ([type]): [description]
    """
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地')
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    tt = (dic['time']).strftime('%Y-%m-%d %H%M')
    print(tt)
    ax.set_title(tt, loc='center', fontsize=25)
    ax.set_title(dic['model'], loc='left', fontsize=25)
    # ax.set_title('OBS', loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=25)

    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    # mb.set_extent('中国陆地')
    # mb.set_extent([107, 135, 20,40])
    mb.set_extent([110, 116, 32, 36])
    # mb.southsea(zoom=0.3, loc='left_bottom')

    ## 等值线(等高线)
    # for df in hgt_list:
    #     # print(df)
    #     y = df['label_lat'][0]
    #     x = df['label_lon'][0]
    #     if x>107 and x<135 and y>20 and y<40:
    #         if df['height'][0] in ['140', '144', '148',  '584', '588', '1240', '1248', '1256', '1264']:
    #             ax.text(df['label_lon'][0], df['label_lat'][0], df['height'][0], fontsize=15,color='blue', transform=ccrs.PlateCarree())
    #     c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='blue')
        
    # for df in tmp_list:
    #     # print(df)
    #     y = df['label_lat'][0]
    #     x = df['label_lon'][0]
    #     if x>107 and x<135 and y>20 and y<40:
    #         ax.text(df['label_lon'][0], df['label_lat'][0], df['temp'][0], fontsize=15,color='red', transform=ccrs.PlateCarree())
    #     c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='red', linewidth=0.8, alpha=0.8)

    # wdif, wobs, wf = get_dif(dic)

    

    cs = draw_contourf(ax,qdif)
    # # contour_levels = [-0.020, -0.016, -0.012, -0.008, -0.004, -0.002, 0.002, 0.004, 0.008, 0.012, 0.016, 0.020]
    # contour_levels = [10, 12, 14, 16, 18, 20, 22, 24, 30]
    # contour_levels = [-50, -25, -16, -8, -4, 4, 16, 20, 25, 50]
    # # contour_levels = [-20, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 20]
    contour_levels = [10, 12, 14, 16, 18, 20, 22, 24, 30]
    # contour_levels = np.arange(0, 17, 1)
    # contour_levels = [-16, -12, -8, -6, -4,-2,-1, 1, 2, 4, 6, 8, 12, 16]
    # contour_levels = [-24, -18, -12,-6,-3, 3, 6, 12, 18, 24]
    # contour_levels = [-22,-18, -14, -10,-6, -2, 2,  6, 10, 14, 18, 22]
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.05,  # 色标大小
        pad=0.12,  # colorbar和图之间的距离
        ticks=contour_levels,

    )
    # cmp = cmaps.precip3_16lev
    # norm = mpl.colors.BoundaryNorm([0, 2, 4, 6, 8, 12, 14, 16, 20 , 30], cmp.N)
    # # norm = mpl.colors.Normalize(vmin=2, vmax=16)
    # cax, _ = mpl.colorbar.make_axes(ax, location='right')
    # cb3 = mpl.colorbar.ColorbarBase(
    #     norm = norm,
    #     cmap = cmp,
    #     ax = cax,
    # )
    # font = {'size':18}
    # cb.set_label(r'$q_{fo}-q_{obs}$  $(g \cdot cm^{-1} \cdot s^{-1})$', fontdict=font)
    cb.ax.tick_params(labelsize=24)  # 设置色标标注的大小
    # cb.ax.ticks(contour_levels)
    # cb.ax.set_xticks(contour_levels)
    # draw_contour(ax, wobs)
    # ax.contourf(u.lon, u.lat, u)
    # gl = mb.gridlines(font_size=25, alpha=0)
    add_ticks(ax)
    draw_station(ax)

    # x = ddf['lon']
    # y = ddf['lat']
    # u = -1*np.sin(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    # v = -1*np.cos(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    # ax.barbs(x,y,u,v)
    # Q = ax.quiver(x, y, u,v,units='inches',scale=30,pivot='tip',width=0.04, transform=ccrs.PlateCarree())  # 绘制风矢
    draw_quiver(qu,qv,ax)
    # x = qu.lon
    # y = qu.lat
    # Q = ax.quiver(x, y, qu,qv, headwidth=3, headlength=4,width=0.03, units='inches',scale=30,pivot='tip', transform=ccrs.PlateCarree())  # 绘制风矢
    # qk = ax.quiverkey(Q, X=0.8, Y=0.27, U=10, label=r'$10\ m/s$', labelpos='E',coordinates='figure', fontproperties={'size':25})   # 设置参考风矢
    # mb.southsea(zoom=0.3, loc='left_bottom')
    # draw_south_sea(fig)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/850/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')
    # fig.savefig('test.png')
    fig.savefig(fig_name)


def draw_all(dic_model):
    """
    在一张图上画所有需要的要素
    获取数据，调用画图程序
    高度场
    温度场
    风场差

    Args:
        t ([type]): [description]
        level (str, optional): [description]. Defaults to '500'.
    """
    # t = pd.Timestamp('2021-07-19 2000')
    level = dic_model['level']
    t = dic_model['time']
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
    dic_v = {
        'var':'height',
        'level':level,
        'time':t-pd.Timedelta('8H')
    }
    # hgt_list = get_analysis(dic_h)  # 人工分析高度场
    # tmp_list = get_analysis(dic_t)  # 人工分析温度场
    # ddf = get_plot(dic_h)  # 高空填图信息
    # wdif, wobs, wf = get_dif(dic_v)  # 模式和观测的风速差
    # draw(hgt_list, tmp_list, ddf, wdif, wobs, dic_t)
    dic_model['time'] = dic_model['time']-pd.Timedelta('8H')
    dic_data = get_data(dic_model)
    # qdif = dic_data['qf_obs']
    # qu = dic_data['u_obs']
    # qv = dic_data['v_obs']
    qdif = dic_data['qf_model']
    qu = dic_data['u_model']
    qv = dic_data['v_model']
    draw(qdif, qu, qv, dic_model)

def draw_one_model(dic_model):
    pass
    if dic_model['model'][-4:] == '1800':
        ## 这里写北京时的目的是为了考虑到观测的分析场(high/Analysis)是北京时
        # ttt = pd.date_range(start='2021-07-18 08', end='2021-07-20 20', freq='12H')
        ttt = pd.DatetimeIndex(['2021-07-18 08', '2021-07-20 08'])
    if dic_model['model'][-4:] == '1812':
        # ttt = pd.date_range(start='2021-07-18 20', end='2021-07-20 20', freq='12H')
        ttt = pd.DatetimeIndex(['2021-07-18 20', '2021-07-20 08'])
    if dic_model['model'][-4:] == '1900':
        # ttt = pd.date_range(start='2021-07-19 08', end='2021-07-20 20', freq='12H')
        ttt = pd.DatetimeIndex(['2021-07-19 08', '2021-07-20 08'])
    if dic_model['model'][-4:] == '1912':
        # ttt = pd.date_range(start='2021-07-19 20', end='2021-07-20 20', freq='12H')
        ttt = pd.DatetimeIndex(['2021-07-19 20', '2021-07-20 08'])
    for level in [850]:
        for t in ttt:
            print("画 [%s] 时 [%s] 的图"%(t.strftime('%Y-%m-%d %H'), dic_model['model']))
            dic_model['time'] = t
            dic_model['level'] = level
            draw_all(dic_model)

def draw_model():
    """画所有模式的数据
    """
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    for f in initial_file_list:
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
        for t in time_list:
            flnm = 'YSU_'+t
            path_out = path_main+flnm+'_upar_d02_latlon.nc'
            # print(path_out)
            dic_model = {'model': f+t, 'flnm':path_out}
            # print(dic_model['flnm'])
            draw_one_model(dic_model)


def draw_model_once():
    """画所有模式的数据
    """
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    # flnm = 'high_resolution_high_hgt_upar_d04_latlon.nc'
    flnm = 'high_resolution_upar_d04_latlon.nc'
    path_out = path_main+flnm
    dic_model = {'model': '1km', 'flnm':path_out}
    # dic_model = {'model': '1km_hr', 'flnm':path_out}

    # for t in ttt:
    print("画 [%s] 的图"%(dic_model['model']))
    # dic_model['time'] = pd.Timestamp('2021-07-20 00')
    ### 北京时
    # ttt = pd.DatetimeIndex(['2021-07-19 08', '2021-07-20 08'])
    ttt = pd.DatetimeIndex(['2021-07-19 08', '2021-07-20 16'])
    dic_model['time'] = ttt[1]
    dic_model['level'] = 925
    # draw_one_model(dic_model)
    # print(dic_model)
    
    
    draw_all(dic_model)
    


### 测试结束
# %%
if __name__ == '__main__':
    pass
    draw_model_once()
    # draw_model()
    # draw_obs()



# %%