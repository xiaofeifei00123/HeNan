#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
读取micaps高空形势场数据
第14类资料
使用nwc(国家气象中心的python包)
-----------------------------------------
Time             :2021/10/03 14:33:01
Author           :Forxd
Version          :1.0
'''


# %%
from cartopy.crs import Projection
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

def get_analaysis_h():
    """读micaps14类数据

    Returns:
        [type]: [description]
    """
    
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/HGT/500/20210720080000.000'
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
        df1['height'] = bb
        df1['label_lon'] = cc[0][0]
        df1['label_lat'] = cc[0][1]
        line_list.append(df1)
        
    return line_list


def get_analaysis_t():
    """读micaps14类数据

    Returns:
        [type]: [description]
    """
    
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/ANALYSIS/TMP/500/20210720080000.000'
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
        df1['temp'] = bb
        df1['label_lon'] = cc[0][0]
        df1['label_lat'] = cc[0][1]
        line_list.append(df1)
        
    return line_list


li = get_analaysis_t()
# %%
# li[0]
# for line in li:
    # print(line['label_lat'].max())


    # num = len(dc['line_xyz'])
    # loc_list = []
    # for i in range(num):
    #     aa = dc['line_xyz'][i]
    #     cc = float(dc['line_label'][i])
    #     aa[:,2] = cc
    #     loc_list.append(aa)
    # dd = np.concatenate(loc_list, axis=0)
    # df = pd.DataFrame(aa, columns=['lon', 'lat', 'height'])
    # df
    # return df


# %%
def get_plot():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/Micaps/high/PLOT/850/20210720080000.000'
    df = read_micaps_2(flnm)
    sta = meb.sta_data(df, columns = [
                'ID', 'lon', 'lat', 'alt', 'grade', 'height', 'temperature', 'dewpoint_depression',
                'wind_angle', 'wind_speed', 'time', 'level'])

    cols = ['height', 'ID', 'lon', 'lat', 'time', 'level']
    aa = sta[cols]
    # stb = meb.sta_data(aa, columns=['temperature', 'id','lon', 'lat','time', 'level'])
    stb = meb.sta_data(aa, columns=['height', 'id','lon', 'lat','time', 'level'])
    grid1 = meb.grid([60,150,0.25],[10,60,0.25])
    grd2 = meb.interp_sg_idw(stb, grid1)
    return grd2.squeeze(), sta
    # meb.plot_tools.contourf_2d_grid(grd2)

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

    # lona
    # np.nanmin(idw_eta)
    # tbb_return = xr.Dataset(
    #                 {
    #                     'tbb':(['lon', 'lat'], idw_eta)
    #                 },
    #                 coords={
    #                     'lat':lat,
    #                     'lon':lon,
    #                     'time':tt,
    #                     },
    #                 attrs={'var_name':'tbb' },)
    # # return tbb_return

# %%
# da, sta = get_plot()
# x,y, z = interp_metpy(sta)
# 
# pp  = get_analaysis()
# pp
# ll = get_analaysis()
# ddd = ll[0]
# ddd




# %%
def draw(hgt_list, tmp_list):
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.05,0.05,0.9,0.9], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地')
    # mb.drawcoastlines()

    mb.drawstates(linewidths=1.2) # 省界

    # c = mb.contourf(da.lon, da.lat, da.squeeze(),transform=ccrs.PlateCarree(), levels=10)
    # ax.scatter(lon, lat, h,transform=ccrs.PlateCarree())
    # c = ax.contourf(lon,lat,idw_eta.T, transform=ccrs.PlateCarree(),)
    # c = ax.pcolormesh(lon,lat,idw_eta.T, transform=ccrs.PlateCarree(),)
    # c = ax.pcolormesh(x,y,z, transform=ccrs.PlateCarree(),)
    # c = ax.scatter(sta['lon'], sta['lat'], transform=ccrs.PlateCarree())
    for df in hgt_list:
        print(df)
        y = df['label_lat'][0]
        x = df['label_lon'][0]
        if x>73 and x<136 and y>15 and y<56:
            ax.text(df['label_lon'][0], df['label_lat'][0], df['height'][0], fontsize=15,color='blue', transform=ccrs.PlateCarree())
        # ax.text(df['label_lon'][0], df['label_lat'][0], df['height'][0], transform=ccrs.PlateCarree())
        c = ax.plot(df['lon'], df['lat'].values,transform=ccrs.PlateCarree(), color='blue')

    # c = ax.plot(df['lon'].values, df['lat'].values,transform=ccrs.PlateCarree())
    mb.set_extent('中国陆地')

    fig.savefig('test.png')
    # ax.set_extent([73, 136, 15, 56], crs=ccrs.PlateCarree())
    # ax.set_extent([73, 136, 15, 56])
    # plt.show()
        # ax.text(99, 40, 'aa', transform=ccrs.PlateCarree())
    # c = ax.plot(x,y, transform=ccrs.PlateCarree())
    # c = ax.contour(df['lon'],df['lat'],df['height'], transform=ccrs.PlateCarree(), label='aa')
    # ax.clabel(c,inline=True)
    # c = ax.contour(x,y,z, transform=ccrs.PlateCarree(), levels=5)
    # c = ax.pcolormesh(lon,lat,filled.T, transform=ccrs.PlateCarree(),)
    # c = mb.contourf(mx,my,z, transform=ccrs.PlateCarree())

    # mb.gridlines(font_size=25)

    # fig.colorbar(c, pad=0.05, fraction=0.05, orientation='horizontal',)

hgt_list = get_analaysis_h()
tmp_list = get_analaysis_t()
draw(hgt_list, tmp_list)



# %%