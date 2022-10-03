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

# %%
def draw_quiver(ax, u,v):
    '''
    绘制风矢图
    '''
    # da = latlon2distance(u)
    # x = u.cross_line_idx.values[::10]

    # u = u[::3,::10]
    # v = v[::3,::10]

    
    # x = u.cross_line_idx.values[::5]
    # da = latlon2distance(u)
    # x = da.distance.values[::5]
    # u = u[::1,::5]
    # v = v[::1,::5]
    
    
    da = latlon2distance(u)
    x = da.distance.values[::10]
    u = u[::3,::10]
    v = v[::3,::10]
    

    # u = u[::4,::4]
    # v = v[::4,::4]
    # u = u[::9,::30]
    # v = v[::9,::30]
    # u = u[::3,::10]
    # v = v[::3,::10]
    y = u.coords['vertical'].values
    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=10,pivot='middle', zorder=2)  # 绘制风矢
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=50,pivot='tip',minlength=0.001, width=0.01,zorder=2)  # 绘制风矢
    # Q = ax.quiver(x, y, u.values,v.values,units='width',scale=20,pivot='tip', width=0.03,zorder=2)  # 绘制风矢
    qk = ax.quiverkey(Q,
                      X=0.05, Y=0.15, 
                      U=10 ,
                      label=r'$10 m/s$', 
                      labelpos='S',  # label在参考箭头的哪个方向
                      labelsep=0.05, # 箭头和标签之间的距离
                      coordinates='figure',  
                      fontproperties={'size':8}
                      )   # 设置参考风矢


def draw_contour(ax, da):
    pass

    da = latlon2distance(da)
    xs = da.distance.values
    # xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    # levels=np.arange(342, 372, 2)
    # levels=np.arange(336, 372, 4)
    levels=np.arange(336, 357, 4)
    plt.contour
    # cs = ax.contour(xs, ys, smooth2d(da.values, passes=16), levels=levels, colors='black', linewidths=0.5)
    # cs = ax.contour(xs, ys, smooth2d(da.values, passes=4), levels=levels, colors='black', linewidths=0.5)
    cs = ax.contour(xs, ys, smooth2d(da.values, passes=1), levels=levels, colors='black', linewidths=0.5)
    ax.clabel(cs, inline=True, fontsize=10)

def draw_contour2(ax,da):
    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    cs = ax.contour(xs, ys, smooth2d(da.values*10**4, passes=16), colors='red')
    ax.clabel(cs, inline=True, fontsize=10)



def draw_contourf(fig, ax_cross, da, ter_line):


    # colorlevel=[-168, -10, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 10, 168]#垂直速度
    # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
    # colorlevel=[-280, -60, -30, -10, -5, -1, 1, 5, 10, 30, 60, 280]#雨量等级
    # colorlevel=[-680, -100, -40, -15, -5, -1, 1, 5, 15, 40, 100, 680]#雨量等级
    ### 散度
    colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    # colorlevel=[-168, -10, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 10, 168]#垂直速度
    colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
    # colorlevel=[0., 0.1, 0.5, 1, 2,3, 4,6, 8,10,20,  168]#垂直速度

    ### 垂直速度
    # colordict= ['white', '#6CA6CD', '#436EEE', '#66CD00', '#7FFF00','#cdfd02', 'yellow','#fdaa48','#EE7600','red']
    # colorlevel=[0, 0.1, 0.5, 1, 2, 4, 8, 10,15, 20, 168]#垂直速度

    # colorlevel=[-1680, -100, -40, -20, -10, -1, 1, 10, 20, 40, 100, 1680]#雨量等级
    # colorticks=[-30, -20, -10, -5, -1, 1, 5, 10, 20, 30]#雨量等级
    colorticks = colorlevel[1:-1]

    
    da1 = latlon2distance(da)    
    xs = da1.distance.values
    # xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    

    dbz_contours = ax_cross.contourf(xs,
                                    ys,
                                    da.values,
                                    colors=colordict,
                                    levels=colorlevel
    )
    """
    ax_cross.set_ylim(0, 20000)
    ax_cross.set_yticks(np.arange(0, 20000+1, 2000))
    # ax_cross.set_ylim(0, 2000)
    # ax_cross.set_yticks(np.arange(0, 2000+1, 200))
    # y_labels = np.arange(0, 2+0.1, 0.2).astype(int)
    # ax_cross.set_yticklabels(y_labels)
    ax_cross.tick_params(axis='both', labelsize=10, direction='out')
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks, orientation='vertical', fraction=0.06, pad=0.02)

    labels = list(map(lambda x: str(x) if x<1 else str(int(x)), colorticks))  # 将colorbar的标签变为字符串
    cb_dbz.set_ticklabels(labels)  # 改变标签的格式
    
    

    cb_dbz.ax.tick_params(labelsize=10)
    ax_cross.xaxis.set_minor_locator(plt.MultipleLocator(10))
    """
    # print(xs)

    ## Set the x-ticks to use latitude and longitude labels
    # coord_pairs = da.coords["xy_loc"].values
    # x_ticks = np.arange(coord_pairs.shape[0])
    # ax_cross.set_xticks(x_ticks[::12])
    # x_labels = coord_pairs
    # # ax_cross.set_xticklabels(x_labels[::12] ,fontsize=10, rotation=45)
    # num_ticks = 10
    # thin = int((len(x_ticks) / num_ticks) + .5)
    # # 设置xticks
    # xticks = x_ticks[::thin]
    # ax_cross.set_xticks(xticks)
    # # 设置xticklabels
    # latlonlist = x_labels[::thin]
    # latlonlist_new = []
    # for latlon in latlonlist:
    #     latlon = latlon.split(',')
    #     lat = latlon[0]
    #     lon = latlon[1]
    #     latlon_new = lat+'$^{\circ}$'+'N'+'\n'+lon+'$^{\circ}$'+'E'
    #     latlonlist_new.append(latlon_new)
    # ax_cross.set_xticklabels(latlonlist_new)

    ## set tick labels use distance
    coord_pairs = da.coords["xy_loc"].values
    x_ticks = np.arange(coord_pairs.shape[0])
    da1 = latlon2distance(da)    
    x_labels = da1.distance.values.astype(int)

    # ax_cross.set_xticks(x_ticks[::12])
    # ax_cross.set_xticklabels(x_labels[::12] ,fontsize=10)

    # ax_cross.set_ylim(0,10000)

    def add_anotate(lat2, lon2, text):
        lat1, lon1 = coord_pairs[0].split(',')
        lat3, lon3 = coord_pairs[-1].split(',')
        # print(lat1, lon1)
        # print(lat3, lon3)
        # lat2, lon2 = 33.6, 112.6
        loc1 = (lat1, lon1)
        loc2 = (lat2, lon2)
        a = distance(loc1, loc2).km
        # print(type(a))
        idx = np.argmin(np.abs(a-x_labels))
        # print(idx)
        # print(x_labels[idx])
        # print(coord_pairs[idx])
        ax_cross.annotate(text, xy=(idx, 0), xytext=(idx, -5000),
                arrowprops=dict(facecolor='black',  headwidth=8, headlength=8),
                )
    # loc1 = (33.7, 112.55)
    # loc2 = (34.7, 113.35)
    # loc3 = (35.5, 114.0)
    loc1 = (34, 112.45)
    loc2 = (35.15, 113.5)
    loc3 = (35.7, 114.0)
    text_list = ['A', 'B', 'C']
    i = 0
    for loc in [loc1, loc2, loc3]:
        pass
        # lat = loc[0]
        # lon = loc[1]
        # add_anotate(lat, lon,text_list[i])
        # i+=1



    
    

    

    ## Set the desired number of x ticks below
    # x_labels = da.coords['xy_loc'].values
    # num_ticks = 6
    # thin = int((len(x_ticks) / num_ticks) + .5)
    # ax_cross.set_xticks(x_ticks[::thin])
    # ax_cross.set_xticklabels(x_labels[::thin], rotation=30, fontsize=18)

    ## Set the x-axis and  y-axis labels
    # ax_cross.set_xlabel("Latitude, Longitude", fontsize=10)
    ax_cross.set_xlabel("Distance (km)", fontsize=10)
    ax_cross.set_ylabel("Height (km)", fontsize=10)


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


def draw(t='2021-07-20 08', flpath='/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/', model='GWD3'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/cross.nc'
    # flnm = flpath+'cross2.nc'
    # flnm = flpath+'cross3.nc'
    # flnm = flpath+'cross3.nc'
    # flnm = flpath + 'cross4_1time.nc'
    flnm = flpath + 'cross8_1time.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=t)
    # print(ds)
    # ds = ds[:,0:70]
    # ds = ds.isel(cross_line_idx=np.arange(0,100,1))

    w = ds['wa_cross']
    ter_line = ds['ter']
    u = ds['ua_cross']
    v = ds['va_cross']
    theta_e = ds['theta_e_cross']
    theta_e = ds['theta_e_cross']
    div = ds['div_cross']
    drag = ds['drag_cross']
    ws = ds['ws_cross']

    w = drop_na(w)
    u = drop_na(u)
    v = drop_na(v)
    theta_e = drop_na(theta_e)
    div = drop_na(div)
    drag = drop_na(drag)
    ws = drop_na(ws)

    cm = round(1/2.54, 2)
    fig = plt.figure(figsize=(16*cm,6*cm), dpi=600)
    # fig = plt.figure(figsize=(8*cm,2*cm), dpi=600)
    # ax_cross = fig.add_axes([0.08, 0.2, 0.9, 0.7])
    ax_cross = fig.add_axes([0.1, 0.2, 0.88, 0.7])

    title_t = ds.time.dt.strftime('%d-%H').values
    # title_model = flnm.split('/')[-2]
    title_model = model
    print('画[{}]模式[{}]时刻的图'.format(title_model,title_t))
    ax_cross.set_title(title_t, loc='left', fontsize=10)
    ax_cross.set_title(title_model, loc='right', fontsize=10)

    draw_contour(ax_cross, theta_e)
    # draw_contour2(ax_cross, div)
    draw_contourf(fig, ax_cross, div*10**4, ter_line)
    # draw_contourf(fig, ax_cross, w*10, ter_line)
    # draw_contourf(fig, ax_cross, drag*10**5, ter_line)
    # draw_contourf(fig, ax_cross, ws, ter_line)
    # return u,v,w


    ## 计算剖面风
    # hor = np.sqrt(u**2+v**2)    
    # ver = w
    ## 1. 计算切角
    # self.cross_start= CoordPair(lat=33.5, lon=114.5)
    # self.cross_end= CoordPair(lat=35.5, lon=112.5)
    ldic = {
        'lat1':float(ds.attrs['cross_start'][0]),
        'lon1':float(ds.attrs['cross_start'][1]),
        'lat2':float(ds.attrs['cross_end'][0]),
        'lon2':float(ds.attrs['cross_end'][1]),
    }
        # self.cross_start = CoordPair(lat=32, lon=111.2)
        # self.cross_end = CoordPair(lat=36.5, lon=114.7)
    dy = (ldic['lat2']-ldic['lat1'])
    dx = (ldic['lon2']-ldic['lon1'])
    angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度
    hor = u*np.cos(angle)+v*np.sin(angle)
    ver = w
    # draw_quiver(ax_cross,u,v)
    draw_quiver(ax_cross,hor,ver*10)
    # print(hor)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_cross/'
    fig_name = title_model+'_'+title_t+'div'
    # fig_name = title_model+'_'+title_t+'ws'
    # fig.savefig(fig_path+fig_name+'test5.png', bbox_inches = 'tight')
    fig.savefig(fig_path+fig_name+'.png')
    # plt.clf()

def draw_1time(t='2021-07-20 00'):
    # path_main ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    path_main ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'
    # model_list = ['1900_90m', '1900_900m','1912_90m', '1912_900m']
    # model_list = ['1912_90m', '1912_90m_OGWD']
    # model_list = ['gwd0', 'gwd3']
    # model_list = ['CTRL', 'GWD3', 'FD', 'SS']
    model_list = ['CTRL', 'FD', 'SS', 'GWD3']
    # model_list = ['CTRL']
    # model_list = ['GWD3']
    for model in model_list:
        fl = path_main+model+'/wrfout/'
        draw(t=t, flpath=fl, model=model)

def draw_mtime():
    # time_list = pd.date_range('2021-07-17 00', '2021-07-23 00', freq='3H')
    # time_list = pd.date_range('2021-07-20 00', '2021-07-20 00', freq='3H')
    time_list = pd.date_range('2021-07-20 18', '2021-07-20 18', freq='3H')
    # time_list = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='3H')
    # time_list = pd.date_range('2021-07-20 06', '2021-07-20 06', freq='1H')
    # time_list = pd.date_range('2021-07-20 08', '2021-07-20 08', freq='1H')
    for t in time_list:
        draw_1time(t)
    pass

if __name__ == '__main__':
    # draw()
    draw_mtime()
    # draw_1time()

# %%





# %%

# flnm ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/SS/wrfout/cross2.nc'
# ds = xr.open_dataset(flnm)
# ds
# # %%
# da = ds['ua_cross']

# # %%
# #    coord_pairs = da.coords["xy_loc"].values
# coord_pairs = da.coords["xy_loc"].values
# x_ticks = np.arange(coord_pairs.shape[0])

# ## set tick labels use distance
# da1 = latlon2distance(da)    
# x_labels = da1.distance.values.astype(int)
# # %%
# x_labels
# # ax_cross.set_xticks(x_ticks[::12]) ds.xy_loc
