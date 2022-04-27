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
# from netCDF4 import Dataset
import matplotlib.pyplot as plt
# from matplotlib.cm import get_cmap
import cartopy.crs as crs
# from cartopy.feature import NaturalEarthFeature
from geopy.distance import distance  # 根据经纬度计算两点距离
plt.rcParams['axes.unicode_minus']=False 

from baobao.coord_transform import latlon2distance


# %%
class Draw():
    def __init__(self, ):
        pass
        self.colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        # self.colorlevel=[-16, -4, -3, -2,  -1, -0.5, 0.5, 1, 2,3, 4,16]  # 垂直速度
        self.colorlevel=[-16, -4, -2,  -1,-0.5,  -0.1, 0.1,0.5,  1, 2,4,16]  # 垂直速度
        # self.levels=np.arange(336, 372, 4)


    def draw_quiver(self, ax, u,v, scale=100, ulength=10):
        '''
        绘制风矢图
        '''
        x = u.cross_line_idx.values[::3]
        u = u[::3,::3]
        v = v[::3,::3]
        # x = u.cross_line_idx.values[::10]
        # u = u[::4,::10]
        # v = v[::4,::10]
        y = u.coords['vertical'].values
        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=scale,pivot='tip',minlength=0.001, width=0.015,zorder=2)  # 绘制风矢
        qk = ax.quiverkey(Q,
                        X=0.8, Y=1.04, 
                        U=ulength ,
                        label=r'${} m/s$'.format(ulength), 
                        labelpos='E',  # label在参考箭头的哪个方向
                        labelsep=0.05, # 箭头和标签之间的距离
                        fontproperties={'size':10}, 
                        coordinates='axes', # 是相对于ax的位置，还是相对于figure的位置
                        edgecolor='white',
                        facecolor='black'
                        )   # 设置参考风矢

    def draw_contour(self, ax, da, levels=np.arange(336, 372, 4)):
        pass

        xs = np.arange(0, da.shape[-1], 1)
        ys = da.coords['vertical'].values
        plt.contour
        cs = ax.contour(xs, ys, smooth2d(da.values, passes=16), levels=levels, colors='black', linewidths=0.5)
        ax.clabel(cs, inline=True, fontsize=10)

    def draw_contour2(self, ax,da):
        xs = np.arange(0, da.shape[-1], 1)
        ys = da.coords['vertical'].values
        # levels=np.arange(342, 372, 4)
        cs = ax.contour(xs, ys, smooth2d(da.values*10**4, passes=16), colors='red')
        ax.clabel(cs, inline=True, fontsize=10)

    def draw_contourf(self, ax_cross, da, ter_line):

        xs = np.arange(0, da.shape[-1], 1)
        ys = da.coords['vertical'].values
        cf = ax_cross.contourf(xs,
                                ys,
                                da.values,
                                colors=self.colordict,
                                levels=self.colorlevel
                            )

        ## 画山, 填充为灰色
        ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
                                        facecolor="#434343", zorder=3)
        return cf

    def draw_contourf_no_terrain(self, ax_cross, da, xs,ys):

        # xs = np.arange(0, da.shape[-1], 1)
        # ys = da.coords['vertical'].values
        cf = ax_cross.contourf(xs,
                                ys,
                                da.values,
                                colors=self.colordict,
                                levels=self.colorlevel
                            )

        return cf

    def set_ticks(self, ax_cross, da):
        pass
        ax_cross.set_ylim(0, 12000)
        ax_cross.set_yticks(np.arange(0, 12000+1000, 1000))
        # 可以设置为和纵坐标相关的数或者是啥东西    
        ax_cross.set_yticklabels(np.arange(0, 12000+1000, 1000)/1000)


        ax_cross.tick_params(axis='both', labelsize=10, direction='out')
        # cb_dbz = fig.colorbar(dbz_contours, ax=ax_cross, ticks=self.colorlevel[1:-1])
        # cb_dbz.ax.tick_params(labelsize=7)

        # Set the x-ticks to use latitude and longitude labels
        coord_pairs = da.coords["xy_loc"].values
        x_ticks = np.arange(coord_pairs.shape[0])

        ## set tick labels use distance
        da1 = latlon2distance(da)    
        x_labels = da1.distance.values.astype(int)
        ax_cross.set_xticks(x_ticks[::8])
        ax_cross.set_xticklabels(x_labels[::8], rotation=0, fontsize=10)

        # ax_cross.set_xlabel("Distance (km)", fontsize=10)
        ax_cross.set_ylabel("Height (km)", fontsize=10)



class GetData():
    def __init__(self, ):
        pass

    def drop_na(self, da):
        """处理数据, 这一步是必须要的，不然好像画不出来图
        好像是要把缺测的值去掉, 还要去掉白框
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

    def caculate_wind_hatch(self, u, v, w):
        ## 计算剖面风
        ldic = {
            'lat1':33.8,
            'lon1':111.7,
            'lat2':33.5,
            'lon2':113.2,
        }
        dy = (ldic['lat2']-ldic['lat1'])
        dx = (ldic['lon2']-ldic['lon1'])
        angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度
        hor = u*np.cos(angle)+v*np.sin(angle)
        ver = w
        return hor, ver

    def minus(self, da1, da2):
        """不同试验的剖面上，虽然垂直层数一样，但是其高度是有区别的, 最好的做法应该是插值，这里假设其高度是一样的

        Args:
            da1 (_type_): _description_
            da2 (_type_): _description_

        Returns:
            _type_: _description_
        """
        da = da2.values - da1.values    
        da = xr.DataArray(
            da,
            coords = da1.coords, 
            dims = da1.dims
        )
        return da

    def get_data(self, t='2021-07-20 00', flnm ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_zhengzhou.nc'):
        ds = xr.open_dataset(flnm)
        ds = ds.sel(time=t)

        w = ds['wa_cross']
        ter_line = ds['ter']
        u = ds['ua_cross']
        v = ds['va_cross']
        theta_e = ds['theta_e_cross']
        div = ds['div_cross']

        w = self.drop_na(w)
        u = self.drop_na(u)
        v = self.drop_na(v)
        theta_e = self.drop_na(theta_e)
        div = self.drop_na(div)

        hor, ver = self.caculate_wind_hatch(u,v,w)

        return {
            'u':u,
            'v':v,
            'w':w,
            'theta_e':theta_e,
            'div':div,
            'hor':hor,
            'ver':ver,
            'ter_line':ter_line
        }


def draw_one():
    pass
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm,7*cm), dpi=600)
    ax_cross = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    flnm ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_zhengzhou.nc'
    gd = GetData()
    dic = gd.get_data(t='2021-07-20 12', flnm=flnm)

    
    dr = Draw()
    dr.draw_contour(ax_cross, dic['theta_e'])
    dr.draw_contourf(ax_cross, dic['w'], dic['ter_line'])
    dr.draw_quiver(ax_cross,dic['hor'],dic['ver'], scale=60)
    dr.set_ticks(ax_cross, dic['theta_e'])
        
# cm = 1/2.54
# fig = plt.figure(figsize=(8*cm,7*cm), dpi=300)

def multi():
    cm = 1/2.54
    # fig = plt.figure(figsize=(19*cm, 8*cm), dpi=300)
    # ax_cross = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    fig = plt.figure(figsize=(8*cm, 24*cm), dpi=600)
    proj = crs.PlateCarree()  # 创建坐标系
    grid = plt.GridSpec(3,
                        1,
                        figure=fig,
                        left=0.15,
                        right=0.95,
                        bottom=0.05,
                        top=0.98,
                        wspace=0.2, # 两子图之间的宽
                        # hspace=0.25)
                        hspace=0.15)  # 两子图之间的高

    num = 3
    axes = [None]*num
    for i in range(num):
        axes[i] = fig.add_subplot(grid[i])


    flnm1 ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_zhengzhou.nc'
    flnm2 ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/cross_zhengzhou.nc'
    gd = GetData()
    dic1 = gd.get_data(t='2021-07-20 12', flnm=flnm1)
    dic2 = gd.get_data(t='2021-07-20 12', flnm=flnm2)


    dr0 = Draw()
    dr0.draw_contour(axes[0], dic1['theta_e'])
    cf0 = dr0.draw_contourf(axes[0], dic1['w'], dic1['ter_line'])
    dr0.draw_quiver(axes[0],dic1['hor'],dic1['ver'], scale=60)
    dr0.set_ticks(axes[0], dic1['theta_e'])

    dr1 = Draw()
    dr1.draw_contour(axes[1], dic2['theta_e'])
    cf1 = dr1.draw_contourf(axes[1], dic2['w'], dic2['ter_line'])
    dr1.draw_quiver(axes[1],dic2['hor'],dic2['ver'], scale=60)
    dr1.set_ticks(axes[1], dic1['theta_e'])


    dr2 = Draw()
    levels=np.arange(-12, 12+1, 4)
    dr2.draw_contour(axes[2], gd.minus(dic2['theta_e'],dic1['theta_e']), levels=levels)
    cf2 = dr2.draw_contourf(axes[2], gd.minus(dic2['w'],dic1['w']), dic1['ter_line'])
    dr2.draw_quiver(axes[2],gd.minus(dic2['hor'],dic1['hor']),gd.minus(dic2['ver'],dic1['ver']), scale=60)
    dr2.set_ticks(axes[2], gd.minus(dic2['theta_e'],dic1['theta_e']))



    cf_list = [cf0, cf1, cf2]
    laebel_list = ['(d)', '(e)', '(f)']
    for i in range(3):
        cb = fig.colorbar(
            # cf1,
            cf_list[i],
            ax=axes[i],
            orientation='horizontal',
            ticks=dr1.colorlevel[1:-1],
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小


    # axes[2].set_xlabel('Distance (km)')





    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/Figure8_zz.png')
# multi()



# %%
