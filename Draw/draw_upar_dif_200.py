#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画高空模式和观测的差值填色图
并添加高空观测数据的形势场
200hPa， 作风速差
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
import cmaps
from wrf import get_cartopy, smooth2d, getvar
# from metdig.io.cassandra import get_obs_station

# %%
def get_dif(dic):
    """获得模式和观测的差值, 200,500,850各不相同

    Args:
        dic ([dict]): 时间和层次信息

    Returns:
        [type]: [description]
    """
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/obs_upar_latlon.nc' 
    ds_obs = xr.open_dataset(flnm_obs)
    # flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d02_latlon.nc'
    flnm_wrf = dic['flnm']
    ds_wrf = xr.open_dataset(flnm_wrf)
    # ds_wrf.time
    t = dic['time']
    level = dic['level']
    # t = pd.Timestamp('2021-07-19 00')
    ds1 = ds_obs.sel(time=t, pressure=level)
    ds2 = ds_wrf.sel(time=t, pressure=level)
    u1 = ds1['u'] 
    v1 = ds1['v']
    u2 = ds2['ua']
    v2 = ds2['va']
    ws1 = xr.ufuncs.sqrt(u1**2+v1**2)
    ws2 = xr.ufuncs.sqrt(u2**2+v2**2)
    # u = u2-u1
    w_dif = ws2-ws1
    # u.plot()
    return w_dif, ws1, ws2
# TODO 数据都处理好了，明天就是画图了

# ds1['q'].sel(pressure=200)
# ds.time
# %%
# wd, wo, wf = get_dif()
# wf.max()

# wo.max()
# wd.max()

# %%





# def get_analaysis(dic):
def get_analysis(dic={'var':'height', 'level':'500', 'time':pd.Timestamp('2021-07-20 0800') }):
    """读micaps 14类数据

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/HGT/500/20210720080000.000'
    print(dic['time'].strftime('%Y-%m-%d %H%M'), dic['level'])
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
    ax.set_yticks(np.arange(20, 40 + 1, 5))
    ax.set_xticks(np.arange(107, 135 + 1, 5))
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    # ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax.tick_params(which='major',length=10,width=2.0) # 控制标签大小 
    ax.tick_params(which='minor',length=5,width=1.0)  #,colors='b')
    ax.tick_params(axis='both', labelsize=25, direction='out')

def draw_south_sea(fig,):
    pass
    # ax2 = fig.add_axes([0.102, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2 = fig.add_axes([0.798, 0.145, 0.2, 0.2],projection=ccrs.PlateCarree())
    ax2.set_extent([105.8, 122,0,25])
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
    contour_levels = [-24, -18, -12,-6,-3, 3, 6, 12, 18, 24]
    
    # colormap = cmaps.precip3_16lev
    colormap = cmaps.ViBlGrWhYeOrRe
    # colormap = cmaps.precip3_16lev_r  # 反转色标
    crx = ax.contourf(x,
                        y,
                        da,
                        cmap=colormap,
                    #   norm=norm,
                    #   extend='both',
                    #   extend='max',
                        levels=contour_levels,
                        transform=ccrs.PlateCarree())
    return crx

def draw_contour(ax, da, **kw):
    """在地图上绘制等温线
    """
    x = da.lon
    y = da.lat
    level_contour = np.arange(10, 36, 4)
    # da = smooth2d(field=da, passes=8)
    
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

def draw(hgt_list, tmp_list, ddf, wdif, wobs, dic):
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
    ax.set_title(tt, loc='center', fontsize=25)
    ax.set_title(dic['model'], loc='left', fontsize=25)
    ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=25)

    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    # mb.set_extent('中国陆地')
    mb.set_extent([107, 135, 20,40])
    # mb.southsea(zoom=0.3, loc='left_bottom')

    for df in hgt_list:
        # print(df)
        y = df['label_lat'][0]
        x = df['label_lon'][0]
        if x>107 and x<135 and y>20 and y<40:
            if df['height'][0] in ['140', '144', '148',  '584', '588', '1240', '1248', '1256', '1264']:
                ax.text(df['label_lon'][0], df['label_lat'][0], df['height'][0], fontsize=15,color='blue', transform=ccrs.PlateCarree())
        c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='blue')
        
    # for df in tmp_list:
    #     # print(df)
    #     y = df['label_lat'][0]
    #     x = df['label_lon'][0]
    #     if x>107 and x<135 and y>20 and y<40:
    #         ax.text(df['label_lon'][0], df['label_lat'][0], df['temp'][0], fontsize=15,color='red', transform=ccrs.PlateCarree())
    #     c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='red', linewidth=0.8, alpha=0.8)

    # wdif, wobs, wf = get_dif(dic)
    cs = draw_contourf(ax,wdif)
    cb = fig.colorbar(
        cs,
        orientation='horizontal',
        fraction=0.05,  # 色标大小
        pad=0.1,  # colorbar和图之间的距离
    )
    draw_contour(ax, wobs)
    # ax.contourf(u.lon, u.lat, u)
    # gl = mb.gridlines(font_size=25, alpha=0)
    add_ticks(ax)
    draw_station(ax)

    x = ddf['lon']
    y = ddf['lat']
    u = -1*np.sin(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    v = -1*np.cos(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    # ax.barbs(x,y,u,v)
    # Q = ax.quiver(x, y, u,v,units='inches',scale=30,pivot='tip',width=0.04, transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u,v,headwidth=3, headlength=4,width=0.03, units='inches',scale=30,pivot='tip', transform=ccrs.PlateCarree())  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.8, Y=0.21, U=10, label=r'$10\ m/s$', labelpos='E',coordinates='figure', fontproperties={'size':25})   # 设置参考风矢
    # mb.southsea(zoom=0.3, loc='left_bottom')
    # draw_south_sea(fig)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar_micaps/'
    fig_name = str(fig_path)+str(dic['model'])+'_'+str(dic['level'])+'_'+(dic['time']).strftime('%Y%m%d%H')+'dif'
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
    hgt_list = get_analysis(dic_h)  # 人工分析高度场
    tmp_list = get_analysis(dic_t)  # 人工分析温度场
    ddf = get_plot(dic_h)  # 高空填图信息
    dic_model['time'] = dic_model['time']-pd.Timedelta('8H')
    wdif, wobs, wf = get_dif(dic_model)  # 模式和观测的风速差
    draw(hgt_list, tmp_list, ddf, wdif, wobs, dic_model)

### 单个时次，单个层次测试
# t = pd.Timestamp('2021-07-19 0800')
# draw_all(t, 200)
def draw_one_model(dic_model):
    pass
    if dic_model['model'][-4:] == '1800':
        ttt = pd.date_range(start='2021-07-18 08', end='2021-07-20 20', freq='12H')
    if dic_model['model'][-4:] == '1812':
        ttt = pd.date_range(start='2021-07-18 20', end='2021-07-20 20', freq='12H')
    if dic_model['model'][-4:] == '1900':
        ttt = pd.date_range(start='2021-07-19 08', end='2021-07-20 20', freq='12H')
    if dic_model['model'][-4:] == '1912':
        ttt = pd.date_range(start='2021-07-19 20', end='2021-07-20 20', freq='12H')
    for level in [200]:
        for t in ttt:
            dic_model['time'] = t
            dic_model['level'] = level
            draw_all(dic_model)
### 测试结束
# %%
if __name__ == '__main__':
    pass
    # ttt = pd.date_range(start='2021-07-18 08', end='2021-07-20 20', freq='12H')
    # for level in [200]:
    #     for t in ttt:
    #         draw_all(t, level)
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
    



# %%