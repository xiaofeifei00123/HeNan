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
matplotlib.rcParams['axes.unicode_minus']=False


# %%
def draw_contour(ax, x, y,da, **kw):
    """在地图上绘制等温线
    """
                        # levels=contour_levels,
    crx = ax.contour(x,
                        y,
                        da,
                        colors = 'blue',
                        # levels=contour_levels,
                        linestyles = 'solid',
                        # transform=ccrs.PlateCarree(),
                        # linewidth = 0.5,
                        alpha=0.8)

    ax.clabel(crx,inline=1, fontsize=20, colors='blue') # 等值线的标注
    return crx

def draw_quiver(ax,x,y,u,v):
    '''
    绘制风矢图
    '''
    u = u[::2,::2]
    v = v[::2,::2]
    x = x[::2]
    y = y[::2]
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=25,pivot='middle')  # 绘制风矢
    qk = ax.quiverkey(Q, X=0.67, Y=0.03, U=10, label=r'$(v, 10 m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢

def draw_contourf(fig, ax_cross, xs,ys,da,):

    # xs = np.arange(0, da.shape[-1], 1)
    # ys = da.coords['vertical'].values
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
    # ax_cross.set_ylim(0, 12000)
    # ax_cross.set_yticks(np.arange(0, 12000+1000, 1000))
    cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=colorticks)
    # ax_cross.tick_params(axis='both', labelsize=18, direction='out')
    cb_dbz.ax.tick_params(labelsize=12)

    # ## Set the x-ticks to use latitude and longitude labels
    # # coord_pairs = da.coords["xy_loc"].values
    # # x_ticks = np.arange(coord_pairs.shape[0])
    # x_ticks = da.time.values
    # x_labels = da.time.dt.strftime('%d/%H')
    # # x_labels = da.coords['xy_loc'].values
    # ## Set the desired number of x ticks below
    # # num_ticks = 6
    # # thin = int((len(x_ticks) / num_ticks) + .5)
    # ax_cross.set_xticks(x_ticks[::6])
    # ax_cross.set_xticklabels(x_labels[::6].values, rotation=30, fontsize=18)

    ## Set the x-axis and  y-axis labels
    # ax_cross.set_xlabel("Latitude, Longitude", fontsize=22)
    # ax_cross.set_ylabel("Height (m)", fontsize=22)
    ax_cross.invert_yaxis()

def set_ticks(fig,ax, x,y, pdic):
    """绘制图片的标签之类的
    """
    pass
    ## 标题之类的
    ax.set_title(pdic['model'], fontsize=26, loc='left')
    ax.set_xlabel("Time (Day/Hour)", fontsize=22)
    ax.set_ylabel("Pressure (hPa)", fontsize=22)

    ## 坐标标签
    x_ticks = x.values
    x_labels = x.time.dt.strftime('%d/%H')
    ax.set_xticks(x_ticks[::6])
    ax.set_xticklabels(x_labels[::6].values, rotation=0, fontsize=18)
    ax.tick_params(axis='both', labelsize=18, direction='out')

    ## 图例绘制



# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/time_height_gwd1.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3.nc'
def draw(pdic):
    flnm = pdic['flnm']
    ds = xr.open_dataset(flnm)
    dds = ds.interpolate_na(dim='pressure',method='linear')
    dda = (dds['wind_speed']).T
    u = (dds['u']).T
    v = (dds['v']).T
    w = (dds['w']).T
    # theta = (dds['theta']).T
    # theta_e = (dds['theta_e']).T
    # theta_v = (dds['theta_v']).T
    x = dda.time
    y = dda.pressure
    wh = np.sqrt(u**2+v**2) # wind_horizontal

    fig = plt.figure(figsize=[10,10])
    ax = fig.add_axes([0.11,0.1,0.88,0.8])
    # ax.set_title(pdic['model'], fontsize=26, loc='left')
    # ax.set_xlabel("Time (Day/Hour)", fontsize=22)
    # ax.set_ylabel("Pressure (hPa)", fontsize=22)
    # draw_quiver(ax,x,y,u,v)
    draw_contourf(fig,ax,x,y,w*10)
    draw_quiver(ax,x,y,u,v)
    crx = draw_contour(ax,x,y,wh)
    set_ticks(fig,ax,x,y,pdic)
    

    
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_time_cross/'
    fig_name = pdic['model']
    fig.savefig(fig_path+fig_name+'.png')
    
    


# %%

def main():

    model_list = ['gwd0', 'gwd1', 'gwd3']
    flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'

    pdic = {}
    for model in model_list:    
        flpath1 = os.path.join(flpath,model)
        fname = 'time_height_'+model+'.nc'
        flnm = os.path.join(flpath1, fname)
        pdic['flnm'] = flnm
        pdic['model'] = model

        draw(pdic)    

if __name__ == '__main__':
    main()
    


    # model = 'gwd0'
    # flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    # flpath1 = os.path.join(flpath,model)
    # fname = 'time_height_'+model+'.nc'
    # os.path.join(flpath1, fname)


    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/time_height_gwd1.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3.nc'




# %%
# ttt = w.time.dt.strftime('%d/%H')
# ttt[0].values

# model = 'gwd0'
# flpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
# flpath1 = os.path.join(flpath,model)
# fname = 'time_height_'+model+'.nc'
# flnm = os.path.join(flpath1, fname)
# flnm




