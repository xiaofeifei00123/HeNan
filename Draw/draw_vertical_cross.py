#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制剖面图
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
import cmaps
plt.rcParams['axes.unicode_minus']=False 
# from ..Data.read_vertical_cross import CrossData

# %%

# %%
def draw_quiver(ax, u,v):
    '''
    绘制风矢图
    '''
    x = u.cross_line_idx.values[::10]
    # u = u[::3,::10]
    # v = v[::3,::10]
    # u = u[::9,::30]
    # v = v[::9,::30]
    u = u[::3,::10]
    v = v[::3,::10]
    y = u.coords['vertical'].values
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=30,pivot='middle', zorder=2)  # 绘制风矢

def draw_contour(ax, da):
    pass

    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    levels=np.arange(342, 372, 4)
    cs = ax.contour(xs, ys, smooth2d(da.values, passes=16), levels=levels, colors='black')
    ax.clabel(cs, inline=True, fontsize=10)

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
    colorlevel=[-120, -50, -30, -10, -5, -1, 1, 5, 10, 30, 50, 120]#雨量等级
    # colorticks=[-30, -20, -10, -5, -1, 1, 5, 10, 20, 30]#雨量等级
    colorticks = colorlevel[1:-1]
    dbz_contours = ax_cross.contourf(xs,
                                    ys,
                                    da.values*10,
                                    colors=colordict,
                                    levels=colorlevel
    )
    ax_cross.set_ylim(0, 12000)
    ax_cross.set_yticks(np.arange(0, 12000+1000, 1000))
    ax_cross.tick_params(axis='both', labelsize=18, direction='out')
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks)
    cb_dbz.ax.tick_params(labelsize=12)

    ## Set the x-ticks to use latitude and longitude labels
    coord_pairs = da.coords["xy_loc"].values
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = da.coords['xy_loc'].values
    ## Set the desired number of x ticks below
    num_ticks = 6
    thin = int((len(x_ticks) / num_ticks) + .5)
    ax_cross.set_xticks(x_ticks[::thin])
    ax_cross.set_xticklabels(x_labels[::thin], rotation=30, fontsize=18)

    ## Set the x-axis and  y-axis labels
    ax_cross.set_xlabel("Latitude, Longitude", fontsize=22)
    ax_cross.set_ylabel("Height (m)", fontsize=22)


    # ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
    #                                 facecolor="saddlebrown", zorder=2)
    ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
                                    facecolor="#434343", zorder=2)



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



def draw(t='2021-07-20 09', flpath='/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/'):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/cross.nc'
    flnm = flpath+'cross.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=t)

    wa = ds['wa_cross']
    ter_line = ds['ter']
    u = ds['ua_cross']
    v = ds['va_cross']
    theta_e = ds['theta_e_cross']
    div = ds['div_cross']

    wa = drop_na(wa)
    theta_e = drop_na(theta_e)
    div = drop_na(div)

    fig = plt.figure(figsize=(10,8), dpi=400)
    ax_cross = fig.add_axes([0.2, 0.2, 0.75, 0.7])

    title_t = ds.time.dt.strftime('%d-%H').values
    title_model = flnm.split('/')[-2]
    print('画[{}]模式[{}]时刻的图'.format(title_model,title_t))
    ax_cross.set_title(title_t, loc='left', fontsize=18)
    ax_cross.set_title(title_model, loc='right', fontsize=18)

    draw_contour(ax_cross, theta_e)
    # draw_contour2(ax_cross, div)
    draw_contourf(fig, ax_cross, div*10**4, ter_line)
    draw_quiver(ax_cross,u,v)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_cross/'
    fig_name = title_model+'_'+title_t
    fig.savefig(fig_path+fig_name+'.png')
    plt.clf()

def draw_1time(t='2021-07-20 09'):
    path_main ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # model_list = ['1900_90m', '1900_900m','1912_90m', '1912_900m']
    # model_list = ['1912_90m', '1912_90m_OGWD']
    model_list = ['gwd0', 'gwd1', 'gwd3']
    for model in model_list:
        fl = path_main+model+'/'
        draw(t=t, flpath=fl)

def draw_mtime():
    time_list = pd.date_range('2021-07-20 00', '2021-07-20 00', freq='1H')
    # time_list = pd.date_range('2021-07-20 06', '2021-07-20 06', freq='1H')
    # time_list = pd.date_range('2021-07-20 08', '2021-07-20 08', freq='1H')
    for t in time_list:
        draw_1time(t)
    pass

if __name__ == '__main__':
    # draw()
    draw_mtime()






# %%
