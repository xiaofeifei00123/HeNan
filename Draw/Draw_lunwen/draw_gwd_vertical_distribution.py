#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制重力波拖曳应力
读取剖面数据
两类试验，一个是gwd1,一个是gwd3

绘图
-----------------------------------------
Time             :2021/11/23 10:29:43
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import numpy as np
import pandas as pd
import wrf
from wrf import getvar, CoordPair, vertcross, get_cartopy, smooth2d
from netCDF4 import Dataset
import os


import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cmaps
from baobao.coord_transform import latlon2distance
plt.rcParams['axes.unicode_minus']=False 

# %%



class CrossData():
    """获得垂直方向切向的数据
    提供剖面数据
    地形的剖面数据
    """
    def __init__(self, wrf_file) -> None:
        pass
        ## Create the start point and end point for the cross section
        # self.cross_start = CoordPair(lat=34, lon=110.5)
        # self.cross_end = CoordPair(lat=33.5, lon=113)
        self.cross_start = CoordPair(lat=33.8, lon=111.7)
        self.cross_end = CoordPair(lat=33.5, lon=113.2)
        # self.cross_start = CoordPair(lat=35.5, lon=113)
        # self.cross_end = CoordPair(lat=33.5, lon=113.5)
        self.ncfile = Dataset(wrf_file)
        ## 计算垂直坐标, 可以是离地高度、气压等
        # self.vert = getvar(self.ncfile, "height_agl")  # 离地高度坐标
        self.vert = getvar(self.ncfile, "z")  # 海拔高度
        # self.vert = getvar(self.ncfile, "pres")/100  # 气压坐标

    def get_vcross(self, var):
        """获得单个变量的切向数据, 竖着切

        Args:
            var ([type]): 变量名, 需要是wrf-python支持的

        Returns:
            [type]: [description]
        """

        var =  getvar(self.ncfile, var)
        var_vcross = vertcross(var, self.vert, wrfin=self.ncfile,
                                     start_point=self.cross_start,
                                        end_point=self.cross_end, 
                                        latlon=True, )
        ## 改变投影的attrs的格式
        pj = var_vcross.attrs['projection'].proj4()
        var_vcross = var_vcross.assign_attrs({'projection':pj})


        ## 改变xy_loc的coords的存储格式
        coord_pairs = var_vcross.coords["xy_loc"].values
        x_labels = [pair.latlon_str(fmt="{:.1f}, {:.1f}")
                    for pair in coord_pairs]
        var_vcross = var_vcross.assign_coords({'xy_loc':('cross_line_idx',x_labels)})
        return var_vcross

    def get_ter(self,):
        """获得地形高度
        """
        ter = wrf.getvar(self.ncfile, "ter", timeidx=-1)
        ter_line = wrf.interpline(ter, wrfin=self.ncfile, 
                            start_point=self.cross_start,
                            end_point=self.cross_end)
        ter_line = ter_line.assign_attrs({'projection':'lambert'})
        return ter_line

    def get_cross_data(self, var_list=['ua', 'va', 'wa', 'theta_e']):
        """获得垂直切一刀的数据

        Returns:
            [type]: [description]
        """
        da_cross_list = []
        for var in var_list:
            da = self.get_vcross(var)
            da_cross_list.append(da)
        
        ds = xr.merge(da_cross_list)
        return ds

    def get_proj(self):
        z = wrf.getvar(self.ncfile, "z")
        pj = get_cartopy(z)
        return pj




class DrawVertical():
    """绘制垂直剖面图
    """
    def draw_quiver(self, ax, u,v):
        '''
        绘制风矢图
        '''
        x = u.cross_line_idx.values[::5]
        u = u[::3,::5]
        v = v[::3,::5]
        y = u.coords['vertical'].values
        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=30,pivot='middle', zorder=1)  # 绘制风矢

    def draw_contour(self, ax, da):
        pass

        xs = np.arange(0, da.shape[-1], 1)
        ys = da.coords['vertical'].values
        levels=np.arange(342, 372, 2)
        cs = ax.contour(xs, ys, smooth2d(da.values, passes=16), levels=levels, colors='black')
        ax.clabel(cs, inline=True, fontsize=10)

    def draw_contourf(self, fig, ax_cross, da, ter_line):


        xs = np.arange(0, da.shape[-1], 1)
        ys = da.coords['vertical'].values
        # colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        # colordict=['white','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        # colorlevel=[0, 0.01, 0.02, 0.05, 0.1, 1,  5, 10, 100]
        colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
        # colorlevel=[-100,-0.1,-0.05,-0.03,-0.01, 0.01, 0.03, 0.05, 0.1, 100]
        # colorlevel=[-100,-0.1,-0.05,-0.01,-0.0001, 0.0001, 0.01, 0.05, 0.1, 100]
        # colorlevel=[-100,-10**(-1),-10**(-2),-10**(-3),-10**(-4), 0.0001, 0.01, 0.05, 0.1, 100]
        colorlevel=[-100,-0.1,-0.01,-0.001,-0.0001, 0.0001, 0.001, 0.01, 0.1, 100]
        # colorlevel=[-100, -0.35, -0.25, -0.15,-0.05, 0.05, 0.15, 0.25, 0.35, 100]  # 总共10个数
        
        colorticks = colorlevel[1:-1]
        dbz_contours = ax_cross.contourf(xs,
                                        ys,
                                        da.values,
                                        colors=colordict,
                                        # cmap=cmaps.precip3_16lev,
                                        # cmap=cmaps.WhiteBlueGreenYellowRed,
                                        levels=colorlevel
        )
        ax_cross.set_ylim(0, 12000)
        ax_cross.set_yticks(np.arange(0, 12000+1000, 1000))
        ax_cross.set_yticklabels((np.arange(0, 12000+1000, 1000)/1000).astype('int'))
        
        

        ax_cross.tick_params(axis='both', labelsize=10, direction='out')
        cb_dbz = fig.colorbar(dbz_contours, 
                              ax=ax_cross, 
                              ticks=colorticks,
                              orientation='horizontal',
                              fraction = 0.05,  # 色标大小,相对于原图的大小
                              pad=0.1,  #  色标和子图间距离
                              )
            
        # cb_dbz.ax.tick_params(labelsize=8)
        # cb_dbz.ax.set_xscale('log')
        # colorlabel = ['-0.1', '-0.05', '-0.01', '-0.001', '0.001', '0.01', '0.05', '0.1']
        colorlabel = ['$-10^{-1}$', '$-10^{-2}$', '$-10^{-3}$', '$-10^{-4}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$']
        cb_dbz.ax.set_xticklabels(colorlabel, fontsize=8)

        
        
        
        
        
        


        ## Set the x-ticks to use latitude and longitude labels
        coord_pairs = da.coords["xy_loc"].values
        x_ticks = np.arange(coord_pairs.shape[0])

        ax_cross.tick_params(axis='both', labelsize=10, direction='out')
        ## set tick labels use distance
        da1 = latlon2distance(da)    
        x_labels = da1.distance.values.astype(int)
        ax_cross.set_xticks(x_ticks[::8])
        ax_cross.set_xticklabels(x_labels[::8], fontsize=10)
        
        ## Set the x-axis and  y-axis labels
        # ax_cross.set_xlabel("Distance (km)", fontsize=10)
        ax_cross.set_ylabel("Height (km)", fontsize=10)


        # ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
        #                                 facecolor="saddlebrown", zorder=2)
        ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line),
                                        facecolor="#434343", zorder=2)


    def drop_na(self, da):
        """处理数据
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

# %%
def caculate_angle():
    """计算剖面的切角
    Returns:
        [type]: [description]
    """
    ldic = {
        'lat1':33.8,
        'lon1':111.7,
        'lat2':33.5,
        'lon2':113.2,
    }
    dy = (ldic['lat2']-ldic['lat1'])
    dx = (ldic['lon2']-ldic['lon1'])
    angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度
    return angle

def draw_cross_gwd1(flnm):
    # flnm  = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m_OGWD/wrfout_d04_2021-07-20_06:00:00'
    ds = xr.open_dataset(flnm)
    cd = CrossData(flnm)

    u = cd.get_vcross('DTAUX3D')
    v = cd.get_vcross('DTAUY3D')
    
    angle = caculate_angle()
    da = u*np.cos(angle)+v*np.sin(angle)   # 水平的拖曳力在剖面上的投影

    ter_line = cd.get_ter()
    dv = DrawVertical()
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm,8*cm), dpi=600)
    ax_cross = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    # title_t = ds.XTIME.dt.strftime('%d-%H').values[0]
    title_t = 'gwd1'
    # ax_cross.set_title(title_t, loc='left', fontsize=10)

    dv.draw_contourf(fig, ax_cross, dv.drop_na(da)*10**3, ter_line)
    # dv.draw_quiver(u,v)

    u = cd.get_vcross('ua')
    v = cd.get_vcross('va')
    dv.draw_quiver(ax_cross,u,v)
    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/cross+%s_right_d03.png'%title_t)
# %%

def draw_cross_gwd3(flnm):
    # flnm  = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m_OGWD/wrfout_d04_2021-07-20_06:00:00'
    ds = xr.open_dataset(flnm)
    cd = CrossData(flnm)

    ## 计算水平拖曳力
    u_list = [
        'DTAUX3D_LS',
        'DTAUX3D_SS',
        'DTAUX3D_BL',
        'DTAUX3D_FD',
    ]
    v_list = [
        'DTAUY3D_LS',
        'DTAUY3D_SS',
        'DTAUY3D_BL',
        'DTAUY3D_FD',
    ]
    dau = 0
    dav = 0
    for i,j in zip(u_list, v_list):
        u1 = cd.get_vcross(i)
        v1 = cd.get_vcross(j)
        dau = dau+u1
        dav = dav+v1
    ## 1. 计算切角
    angle = caculate_angle()
    print(angle)
    print(dau.mean())
    print(dav.mean())
    hor_drag = dau*np.cos(angle)+dav*np.sin(angle)   # 水平的拖曳力在剖面上的投影
    # hor_drag = dau*np.sin(angle)+dav*np.cos(angle)   # 水平的拖曳力在剖面上的投影
    ## 计算地形高度
    ter_line = cd.get_ter()

    ## 画图
    dv = DrawVertical()
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm,8*cm), dpi=600)
    ax_cross = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    # dv.draw_contourf(fig, ax_cross, dv.drop_na(hor_drag)*10**3, ter_line)
    dv.draw_contourf(fig, ax_cross, dv.drop_na(dau)*10**3, ter_line)

    u = cd.get_vcross('ua')
    v = cd.get_vcross('va')
    w = cd.get_vcross('wa')

    # self.cross_start = CoordPair(lat=35.5, lon=113)
    # self.cross_end = CoordPair(lat=33.5, lon=113.5)
    # ldic = {
    #     'lat1':35.5,
    #     'lon1':113,
    #     'lat2':33.5,
    #     'lon2':113.5,
    # }
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
    
    
    # dv.draw_quiver(ax_cross,u,ver*10)
    # dv.draw_quiver(ax_cross,dv.drop_na(dau)*10**5,dv.drop_na(dav)*10**5)
    # title_t = ds.XTIME.dt.strftime('%d-%H').values[0]
    title_t = 'gwd3'
    # ax_cross.set_title(title_t, loc='left', fontsize=10)
    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/cross+%s_right_d03.png'%title_t)



def main():
    # fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/wrfout_d03_' # 2021-07-20_05:00:00'
    # fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_' # 2021-07-20_05:00:00'
    # tt = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='1H')
    # tt = pd.date_range('2021-07-20 05', '2021-07-20 06', freq='1H')
    # tt = pd.date_range('2021-07-20 16', '2021-07-20 16', freq='1H')
    tt = pd.date_range('2021-07-20 12', '2021-07-20 12', freq='1H')
    # fl_path = /
    fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'

    # model_list = ['gwd1', 'gwd3']
    # model_list = ['gwd1', 'gwd3']
    model_list = ['gwd3']
    for model in model_list:
        for t in tt:
            fname = 'wrfout_d03_'+t.strftime('%Y-%m-%d_%H:%M:%S')
            fpath = os.path.join(fl_path, model)
            flnm = os.path.join(fpath, fname)
            print(fname)
            if model == 'gwd1':
                draw_cross_gwd1(flnm)
            elif model == 'gwd3':
                draw_cross_gwd3(flnm)
            # '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/'

    # flnm_gwd1 = ''

if __name__ == '__main__':
    main()
    