#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制剖面图
散度场
风场
相当位温
-----------------------------------------
Time             :2021/11/09 20:51:01
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
# from wrf import xy, interp2dxy, to_np, getvar, CoordPair, vertcross, get_cartopy
from wrf import smooth2d
import wrf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from geopy.distance import distance  # 根据经纬度计算两点距离
# import cmaps
plt.rcParams['axes.unicode_minus']=False 
# from ..Data.read_vertical_cross import CrossData

# %%

flnm ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross3.nc'
# ds = xr.open_dataset(flnm)
# ds.sel(vertical=3000, method='nearest')
# flnm = flpath+'cross2.nc'
ds = xr.open_dataset(flnm)
# ds = ds.sel(time=t)
ds = ds.sel(vertical=3000, method='ffill',)
# print(ds)
# ds = ds[:,0:70]
# ds = ds.isel(cross_line_idx=np.arange(0,100,1))

w = ds['wa_cross']
ter_line = ds['ter']
u = ds['ua_cross']
v = ds['va_cross']
theta_e = ds['theta_e_cross']
div = ds['div_cross']
# %%
# w.T.plot()
div.plot()

# w.max()
# w.interpolate_na(dim='time', method='linear', fill_value='extrapolate')
# bb = w.interpolate_na(dim='cross_line_idx', method='linear',  fill_value="extrapolate")
# bb.plot()

# %%
def drop_na(da):
    """处理数据, 这一步是必须要的，不然好像画不出来图
    """
    for i in range(da.shape[-1]):
        column_vals = da[:,i].values
        # Let's find the lowest index that isn't filled. The nonzero function
        # finds all unmasked values greater than 0. Since 0 is a valid value
        # for dBZ, let's change that threshold to be -200 dBZ instead.
        first_idx = int(np.transpose((column_vals > -1).nonzero())[0])
        da[0:first_idx, i] = da[first_idx, i]
    da = da.dropna(dim='time')
    return da
# w.T.plot()
aa = drop_na(div)
aa.T.plot()










# %%
def draw_quiver(ax, u,v):
    '''
    绘制风矢图
    '''
    x = u.cross_line_idx.values[::3]
    # x = u.cross_line_idx.values
    # u = u[::3,::10]
    # v = v[::3,::10]
    # u = u[::2,::10]
    # v = v[::2,::10]
    u = u[::3,::3]
    v = v[::3,::3]
    # u = u[::9,::30]
    # v = v[::9,::30]
    # u = u[::3,::10]
    # v = v[::3,::10]
    y = u.coords['vertical'].values
    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=10,pivot='middle', zorder=2)  # 绘制风矢
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=40,pivot='tip',minlength=0.001, width=0.02,zorder=2)  # 绘制风矢
    # Q = ax.quiver(x, y, u.values,v.values,units='width',scale=20,pivot='tip', width=0.03,zorder=2)  # 绘制风矢

def draw_contour(ax, da):
    pass

    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    # levels=np.arange(342, 372, 2)
    levels=np.arange(336, 372, 4)
    cs = ax.contour(xs, ys, smooth2d(da.values, passes=16), levels=levels, colors='black')
    ax.clabel(cs, inline=True, fontsize=18)

def draw_contour2(ax,da):
    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    cs = ax.contour(xs, ys, smooth2d(da.values*10**4, passes=16), colors='red')
    ax.clabel(cs, inline=True, fontsize=10)



def draw_contourf(fig, ax_cross, da, ter_line):


    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
    # colorlevel=[-280, -60, -30, -10, -5, -1, 1, 5, 10, 30, 60, 280]#雨量等级
    colorlevel=[-680, -100, -40, -15, -5, -1, 1, 5, 15, 40, 100, 680]#雨量等级
    # colorticks=[-30, -20, -10, -5, -1, 1, 5, 10, 20, 30]#雨量等级
    colorticks = colorlevel[1:-1]
    dbz_contours = ax_cross.contourf(xs,
                                    ys,
                                    da.values*10,
                                    colors=colordict,
                                    levels=colorlevel
    )
    ax_cross.set_ylim(0, 10000)
    ax_cross.set_yticks(np.arange(0, 10000+1000, 1000))
    ax_cross.tick_params(axis='both', labelsize=18, direction='out')
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks)
    cb_dbz.ax.tick_params(labelsize=12)

    ## Set the x-ticks to use latitude and longitude labels
    coord_pairs = da.coords["xy_loc"].values
    x_ticks = np.arange(coord_pairs.shape[0])

    ## set tick labels use distance
    da1 = latlon2distance(da)    
    x_labels = da1.distance.values.astype(int)
    ax_cross.set_xticks(x_ticks[::8])
    ax_cross.set_xticklabels(x_labels[::8], rotation=30, fontsize=18)

    # ax_cross.set_ylim(0,10000)
    

    ## Set the desired number of x ticks below
    # x_labels = da.coords['xy_loc'].values
    # num_ticks = 6
    # thin = int((len(x_ticks) / num_ticks) + .5)
    # ax_cross.set_xticks(x_ticks[::thin])
    # ax_cross.set_xticklabels(x_labels[::thin], rotation=30, fontsize=18)

    ## Set the x-axis and  y-axis labels
    ax_cross.set_xlabel("Latitude, Longitude", fontsize=22)
    ax_cross.set_ylabel("Height (m)", fontsize=22)


    # ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
    #                                 facecolor="saddlebrown", zorder=2)
    ## 画山
    ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
                                    facecolor="#434343", zorder=3)



def drop_na(da):
    """处理数据, 这一步是必须要的，不然好像画不出来图
    """
    for i in range(da.shape[-1]):
        column_vals = da[:,i].values
        # Let's find the lowest index that isn't filled. The nonzero function
        # finds all unmasked values greater than 0. Since 0 is a valid value
        # for dBZ, let's change that threshold to be -200 dBZ instead.
        first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
        da[0:first_idx, i] = da[first_idx, i]
    da = da.dropna(dim='vertical')
    return da


def latlon2distance(da2):
    """将剖面数据的经纬度横坐标变为距离坐标

    Args:
        da2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
    # ds = xr.open_dataset(flnm)
    # ds1 = ds.sel(time='2021-07-20 08')
    # da = ds1['wa_cross']
    # da1 = da.interpolate_na(dim='vertical', method='linear',  fill_value="extrapolate")
    # da2 = da1.sel(vertical=2000, method='nearest')
    dd = da2.xy_loc
    def str_latlon(string):
        # d1 = dd.values[0]
        lat = float(string.split(',')[0])
        lon = float(string.split(',')[1])
        return lat, lon

    d2 = dd.values
    lat_list = []
    lon_list = []
    for i in d2:
        # print(i)
        lat, lon = str_latlon(i)
        lat_list.append(lat)
        lon_list.append(lon)

    dis_list = [0]
    di = 0
    for i in range(len(lat_list)-1):
        # print(i)
        lat1 = lat_list[i]
        lon1 = lon_list[i]
        loc1 = (lat1, lon1)
        lat2 = lat_list[i+1]
        lon2 = lon_list[i+1]
        loc2 = (lat2, lon2)
        dist = distance(loc1,loc2).km
        di = di+dist
        dis_list.append(di)
    dis_list
    dis_array = (np.array(dis_list)).round(1)
    dis_array
    da2 = da2.assign_coords({'distance':('cross_line_idx',dis_array)})
    da3 = da2.swap_dims({'cross_line_idx':'distance'})
    da3
    return da3





def draw(t='2021-07-20 08', flpath='/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/cross.nc'
    flnm = flpath+'cross2.nc'
    ds = xr.open_dataset(flnm)
    # ds = ds.sel(time=t)
    ds = ds.sel(vertical=3000, method='nearest')
    # print(ds)
    # ds = ds[:,0:70]
    # ds = ds.isel(cross_line_idx=np.arange(0,100,1))

    w = ds['wa_cross']
    ter_line = ds['ter']
    u = ds['ua_cross']
    v = ds['va_cross']
    theta_e = ds['theta_e_cross']
    div = ds['div_cross']

    w = drop_na(w)
    u = drop_na(u)
    v = drop_na(v)
    theta_e = drop_na(theta_e)
    div = drop_na(div)

    # fig = plt.figure(figsize=(10,8), dpi=400)
    # ax_cross = fig.add_axes([0.15, 0.2, 0.8, 0.7])

    # title_t = ds.time.dt.strftime('%d-%H').values
    # title_model = flnm.split('/')[-2]
    # print('画[{}]模式[{}]时刻的图'.format(title_model,title_t))
    # # ax_cross.set_title(title_t, loc='left', fontsize=26)
    # ax_cross.set_title(title_model, loc='right', fontsize=26)

    """
    draw_contour(ax_cross, theta_e)
    draw_contourf(fig, ax_cross, w*10, ter_line)


    ## 计算剖面风
    # hor = np.sqrt(u**2+v**2)    
    # ver = w
    ## 1. 计算切角
    ldic = {
        'lat1':33,
        'lon1':111,
        'lat2':36,
        'lon2':115.5,
    }
    dy = (ldic['lat2']-ldic['lat1'])
    dx = (ldic['lon2']-ldic['lon1'])
    angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度
    hor = u*np.cos(angle)+v*np.sin(angle)
    ver = w
    # draw_quiver(ax_cross,u,v)
    draw_quiver(ax_cross,hor,ver*10)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_cross/'
    fig_name = title_model+'_'+title_t
    fig.savefig(fig_path+fig_name+'test5.png')
    plt.clf()
    """

def draw_1time(t='2021-07-20 09'):
    path_main ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # model_list = ['1900_90m', '1900_900m','1912_90m', '1912_900m']
    # model_list = ['1912_90m', '1912_90m_OGWD']
    model_list = ['gwd0', 'gwd3']
    for model in model_list:
        fl = path_main+model+'/'
        draw(t=t, flpath=fl)

def draw_mtime():
    time_list = pd.date_range('2021-07-20 00', '2021-07-20 03', freq='3H')
    # time_list = pd.date_range('2021-07-20 06', '2021-07-20 06', freq='1H')
    # time_list = pd.date_range('2021-07-20 08', '2021-07-20 08', freq='1H')
    for t in time_list:
        draw_1time(t)
    pass

if __name__ == '__main__':
    # draw()
    draw_mtime()






# %%

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross2.nc'
ds = xr.open_dataset(flnm)
ds