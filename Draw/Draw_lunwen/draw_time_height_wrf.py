#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画时间高度廓线, 值是区域平均值(33-34N,111.5-113E)
-----------------------------------------
Time             :2022/01/05 17:40:10
Author           :Forxd
Version          :1.0
'''

# %%
import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from wrf import smooth2d
# matplotlib.rcParams['axes.unicode_minus']=False
# %%
# flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0.nc'
# flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/upar.nc'
# ds = xr.open_dataset(flpath)
# ds
# ds['q'].mean()*10**3
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3north.nc'
# ds = xr.open_dataset(flnm)
# ds['div'].dropna(dim='pressure')
# ds
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3right.nc'
# ds = xr.open_dataset(flnm)
# ds['w'].max()
# %%
# # flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0south.nc'
# # def get_data(flnm):
# ds = xr.open_dataset(flnm)
# ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
# ds
# he = np.arange(0,20000, 100)
# dds = ds.interpolate_na(dim='height')
# dds.interp(height=he, method='nearest', kwargs={'fill_value':'extrapolate'})
# dds
# dda = (dds['wind_speed']).T
# %%
# ds['u'].dropna(dim='pressure')
# ds['div'].magiitude()
# ds['u'].max()


# %%
def draw_contour(ax, x, y,da, **kw):
    """在地图上绘制等温线
    """
                        # levels=contour_levels,
    # da = smooth2d(field=da.T, passes=3)                
    # da = da.T
    da = smooth2d(field=da, passes=8)                
    # da = da.T
    crx = ax.contour(x,
                        y,
                        da,
                        colors = 'black',

                        # levels=[-16,-8,0, 8, 16],
                        levels=[-16,-8,-4, 4, 8, 16],
                        # levels=[-20,-10,0, 10, 20],
                        # levels = 4,
                        linestyles = 'solid',
                        # linestyles = '--',
                        # transform=ccrs.PlateCarree(),
                    
                        linewidths = 0.8,
                        alpha=1)

    ax.clabel(crx, inline=1, fontsize=12, colors='blue') # 等值线的标注
    return crx

def draw_quiver(ax,x,y,u,v, scale=60):
    '''
    绘制风矢图
    '''

    # u = u[::80,::2]
    # v = v[::80,::2]
    # x = x[::2]
    # y = y[::80]
    u = u[::8,::2]
    v = v[::8,::2]
    x = x[::2]
    y = y[::8]
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=scale,pivot='middle')  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.69, Y=0.05, U=10, label=r'$(v, 10 m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢

def draw_contourf(fig, ax_cross, xs,ys,da,):

    # xs = np.arange(0, da.shape[-1], 1)
    # ys = da.coords['vertical'].values
    # colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    colordict=['#191970','#005ffb','#5c9aff','#98ff98','#FFFFFF','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
    # colorlevel=[-120, -50, -30, -10, -5, -1, 1, 5, 10, 30, 50, 120]#雨量等级
    # colorlevel=[-120, -50, -30, -10,  -1, 1, 10, 30, 50, 120]#雨量等级
    # colorlevel=[-120, -50, -30, -10,  -1, 1, 10, 30, 50, 120]#雨量等级
    # colorlevel=[-80, -30, -10, -5, -1, 1, 5, 10, 30, 80]#雨量等级
    colorlevel=[-80, -8, -5, -3, -1, 1, 3, 5, 8,80]#雨量等级
    # colorlevel=[-280, -30, -20, -10, -5,  5, 10,20, 30, 280]#雨量等级
    colorticks = colorlevel[1:-1]
    dbz_contours = ax_cross.contourf(xs,
                                    ys,
                                    da.values,
                                    colors=colordict,
                                    levels=colorlevel
                                    )
    # ax_cross.set_ylim(1000, 200)
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks)
    cb_dbz.ax.tick_params(labelsize=10)


def set_ticks(fig,ax, x,y, pdic):
    """绘制图片的标签之类的
    """
    pass
    ## 标题之类的
    # ax.set_title('south', fontsize=10, loc='left')
    ax.set_xlabel("Time (Day/Hour)", fontsize=10)
    # ax.set_ylabel("Pressure (hPa)", fontsize=10)
    ax.set_ylabel("Height above ground level (km)", fontsize=10)

    ## 坐标标签
    x_ticks = x.values
    x_labels = x.time.dt.strftime('%d/%H')
    ax.set_xticks(x_ticks[::6])
    ax.set_xticklabels(x_labels[::6].values, rotation=0, fontsize=10)
    ax.tick_params(axis='both', labelsize=10, direction='out')

    ## 图例绘制



# %%
def draw(pdic):
    flnm = pdic['flnm']
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
    # dds = ds.interpolate_na(dim='height_agl',method='nearest')

    
    he = np.arange(0,20000, 100)
    dds = ds.interpolate_na(dim='height_agl')
    dds = dds.interp(height_agl=he, method='nearest', kwargs={'fill_value':'extrapolate'})
    
    

    dda = (dds['wind_speed']).T
    u = (dds['u']).T
    v = (dds['v']).T
    w = (dds['w']).T
    # q = (dds['q']).T*10**3
    div = dds['div'].T*10**5

    
    x = dda.time
    y = dda.height_agl
    y = y/1000
    wh = np.sqrt(u**2+v**2) # wind_horizontal

    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm,8*cm], dpi=300)
    ax = fig.add_axes([0.15,0.2,0.8,0.7])
    
    draw_contourf(fig,ax,x,y,w*10)
    draw_quiver(ax,x,y,u,v)
    crx = draw_contour(ax,x,y,div)
    ax.set_yticks(np.arange(0,15,1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_ylim(0,14)
    set_ticks(fig,ax,x,y,pdic)
    
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig_name = pdic['model']
    fig.savefig(fig_path+fig_name+'time_cross_left.png')
    
    
def draw_minus():
    # flnm = pdic['flnm']
    pdic = {'model':'minus'}
    flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3south.nc'
    flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0south.nc'
    def get_data(flnm):
        ds = xr.open_dataset(flnm)
        ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
        dds = ds.interpolate_na(dim='pressure',method='linear')
        dda = (dds['wind_speed']).T
        u = (dds['u']).T
        v = (dds['v']).T
        w = (dds['w']).T
        # q = (dds['q']).T*10**3
        # theta_e = (dds['theta_e']).T
        # theta = (dds['theta']).T
        # theta_e = (dds['theta_e']).T
        # theta_v = (dds['theta_v']).T
        x = dda.time
        y = dda.pressure
        wh = np.sqrt(u**2+v**2) # wind_horizontal
        return x,y, u,v,w,wh
    x3,y3,u3,v3,w3,wh3 = get_data(flnm3)
    x0,y0,u0,v0,w0,wh0 = get_data(flnm0)

    x = x0
    y = y0

    u = u3-u0    
    v = v3-v0
    w = w3-w0
    wh = wh3-wh0
    
    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm,8*cm], dpi=300)
    ax = fig.add_axes([0.15,0.15,0.8,0.8])
    draw_contourf(fig,ax,x,y,w*10)
    draw_quiver(ax,x,y,u,v, scale=20)
    crx = draw_contour(ax,x,y,wh)
    set_ticks(fig,ax,x,y,pdic)
    
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_time_cross/'
    fig_name = pdic['model']
    fig.savefig(fig_path+fig_name+'north_minus.png', bbox_inches = 'tight')


# %%

def main():

    model_list = ['gwd0', 'gwd3']
    # model_list = ['gwd3']
    flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'

    pdic = {}
    for model in model_list:    
        flpath1 = os.path.join(flpath,model)
        fname = 'time_height_'+model+'left.nc'

        flnm = os.path.join(flpath1, fname)
        pdic['flnm'] = flnm
        pdic['model'] = model
        draw(pdic)    

if __name__ == '__main__':
    main()
    # draw_minus()
    

