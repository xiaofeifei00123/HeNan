#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画时间高度廓线, 值是区域平均值(33-34N,111.5-113E)
多子图合并为一张图
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
def draw_contour(ax, x, y,da, **kw):
    """在地图上绘制等温线
    """
                        # levels=contour_levels,
    da = smooth2d(field=da, passes=8)                
    # levels = 
    levels=[-16,-8,-4, 4, 8, 16]
    if kw:
        if kw['levels']:
            levels = kw['levels']

    crx = ax.contour(x,
                        y,
                        da,
                        colors = 'black',

                        levels=levels,
                        linestyles = 'solid',
                        linewidths = 1,
                        alpha=1)

    ax.clabel(crx, inline=1, fontsize=12, colors='blue') # 等值线的标注
    return crx

def draw_quiver(ax,x,y,u,v, scale=60, ulength=10, xk=0.1, yk=0.1):
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
    # qk = ax.quiverkey(Q, X=0.69, Y=0.05, U=ulength, label=r'$(v, {} m/s)$'.format(ulength), labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢
    qk = ax.quiverkey(Q, X=xk, Y=yk, U=ulength, label=r'$(v, {} m/s)$'.format(ulength), labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢


def draw_contourf(ax_cross, xs,ys,da,):

    colordict=['#191970','#005ffb','#5c9aff','#98ff98','#FFFFFF','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
    colorlevel=[-80, -8, -5, -3, -1, 1, 3, 5, 8,80]#雨量等级
    colorticks = colorlevel[1:-1]
    dbz_contours = ax_cross.contourf(xs,
                                    ys,
                                    da.values,
                                    colors=colordict,
                                    levels=colorlevel
                                    )
    return dbz_contours
    # ax_cross.set_ylim(1000, 200)
    # cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks)
    # cb_dbz.ax.tick_params(labelsize=10)


def set_ticks(ax, x,y,):
    """绘制图片的标签之类的
    """
    pass
    ## 标题之类的
    # ax.set_title('south', fontsize=10, loc='left')
    # ax.set_xlabel("Time (Day/Hour)", fontsize=10)
    # ax.set_ylabel("Pressure (hPa)", fontsize=10)
    # ax.set_ylabel("Height above ground level (km)", fontsize=10)

    ## 坐标标签
    x_ticks = x.values
    x_labels = x.time.dt.strftime('%d/%H')
    ax.set_xticks(x_ticks[::6])
    ax.set_xticklabels(x_labels[::6].values, rotation=0, fontsize=10)
    ax.tick_params(axis='both', labelsize=10, direction='out')

    ## 图例绘制



def draw(ax,flnm):
    # flnm = pdic['flnm']

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

    # cm = 1/2.54
    # fig = plt.figure(figsize=[8*cm,8*cm], dpi=300)
    # ax = fig.add_axes([0.15,0.2,0.8,0.7])
    
    cf = draw_contourf(ax,x,y,w*10)
    draw_quiver(ax,x,y,u,v)
    crx = draw_contour(ax,x,y,div)
    ax.set_yticks(np.arange(0,15,1))  # 这个是选择哪几个坐标画上来的了,都有只是显不显示
    ax.set_ylim(0,14)
    set_ticks(ax,x,y)
    return cf
    
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    # fig_name = pdic['model']
    # fig.savefig(fig_path+fig_name+'time_cross_left.png')
    
    
def draw_minus(ax, loc='left'):
    # flnm = pdic['flnm']
    pdic = {'model':'minus'}
    # flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3left.nc'
    # flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0left.nc'
    flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3'+loc+'.nc'
    flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0'+loc+'.nc'
    def get_data(flnm):
        ds = xr.open_dataset(flnm)
        ds = ds.sel(time=slice('2021-07-20 00', '2021-07-21 00'))
        # dds = ds.interpolate_na(dim='pressure',method='linear')
        he = np.arange(0,20000, 100)
        dds = ds.interpolate_na(dim='height_agl')
        dds = dds.interp(height_agl=he, method='nearest', kwargs={'fill_value':'extrapolate'})

        dda = (dds['wind_speed']).T
        u = (dds['u']).T
        v = (dds['v']).T
        w = (dds['w']).T
        x = dda.time
        y = (dda.height_agl)/1000
        div = dds['div'].T*10**5
        wh = np.sqrt(u**2+v**2) # wind_horizontal
        return x,y, u,v,w,wh, div
    x3,y3,u3,v3,w3,wh3,div3 = get_data(flnm3)
    x0,y0,u0,v0,w0,wh0,div0 = get_data(flnm0)

    x = x0
    y = y0

    u = u3-u0    
    v = v3-v0
    w = w3-w0
    # wh = wh3-wh0
    div = div3-div0
    
    # cm = 1/2.54
    # fig = plt.figure(figsize=[8*cm,8*cm], dpi=300)
    # ax = fig.add_axes([0.15,0.2,0.8,0.7])
    # ax = fig.add_axes([0.15,0.15,0.8,0.8])
    draw_contourf(ax,x,y,w*10)
    draw_quiver(ax,x,y,u,v, scale=20, ulength=5, xk=0.9, yk=0.1)
    # crx = draw_contour(ax,x,y,wh)
    crx = draw_contour(ax,x,y,div, levels=[-16, -8,8, 16])
    ax.set_ylim(0,14)
    set_ticks(ax,x,y)
    ax.set_yticks(np.arange(0,15,1))
    

        
cm = 1/2.54
fig = plt.figure(figsize=(19*cm, 20*cm), dpi=600)
grid = plt.GridSpec(3,
                    3,
                    figure=fig,
                    left=0.1,
                    right=0.95,
                    bottom=0.15,
                    top=0.95,
                    wspace=0.2,
                    hspace=0.25)
                    # hspace=0.2)
                    # )
num = 9
axes = [None] * num  # 设置一个维度为8的空列表
for i in range(num):
    axes[i] = fig.add_subplot(grid[i])
    
draw_minus(axes[2],'all')
draw_minus(axes[5],'left')
draw_minus(axes[8],'right')

# model_list = ['gwd0', 'gwd3']
# model_list = ['gwd3']



flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/'
fname = 'time_height_'+'gwd0'+'all.nc'
flnm1 = os.path.join(flpath, fname)
draw(axes[0], flnm1)


fname = 'time_height_'+'gwd0'+'left.nc'
flnm1 = os.path.join(flpath, fname)
draw(axes[3], flnm1)

fname = 'time_height_'+'gwd0'+'right.nc'
flnm1 = os.path.join(flpath, fname)
draw(axes[6], flnm1)


flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/'
fname = 'time_height_'+'gwd3'+'all.nc'
flnm1 = os.path.join(flpath, fname)
draw(axes[1], flnm1)

fname = 'time_height_'+'gwd3'+'left.nc'
flnm1 = os.path.join(flpath, fname)
draw(axes[4], flnm1)

fname = 'time_height_'+'gwd3'+'right.nc'
flnm1 = os.path.join(flpath, fname)
cf = draw(axes[7], flnm1)


axes[3].set_ylabel("Height above ground level (km)", fontsize=12)
axes[7].set_xlabel("Time (Day/Hour)", fontsize=12)



ax = fig.add_axes([0.22,0.06,0.6,0.02])
# colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
colorlevel=[-80, -8, -5, -3, -1, 1, 3, 5, 8,80]#雨量等级
colorticks = colorlevel[1:-1]
cb = fig.colorbar(
    cf,
    cax=ax,
    orientation='horizontal',
    ticks=colorticks,
    # fraction = 0.05,  # 色标大小,相对于原图的大小
    # pad=0.1,  #  色标和子图间距离
)
fig.text(0.35,0.01, 'wind speed of vertical ($10^{-1} m/s$)')







title_list = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
for i in range(9):
    axes[i].set_title(title_list[i], y=0.98, loc='left', fontsize=10)
fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/time_height.png')






# %%
