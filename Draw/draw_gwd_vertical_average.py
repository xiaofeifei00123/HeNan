#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
绘制重力波拖曳应力
读取剖面数据
绘图
-----------------------------------------
Time             :2021/11/23 10:29:43
Author           :Forxd
Version          :1.0
'''
# %%
from meteva.base.fun.statisticing import var_of_grd
import xarray as xr
import numpy as np
import pandas as pd
import wrf
from wrf import getvar, CoordPair, vertcross, get_cartopy, smooth2d
from netCDF4 import Dataset


import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import cmaps

# %%



class CrossData():
    """获得垂直方向切向的数据
    提供剖面数据
    地形的剖面数据
    """
    def __init__(self, wrf_file) -> None:
        pass
        ## Create the start point and end point for the cross section
        self.cross_start = CoordPair(lat=33.5, lon=110.5)
        self.cross_end = CoordPair(lat=33.5, lon=116)
        # self.cross_start = CoordPair(lat=35.5, lon=112.5)
        # self.cross_end = CoordPair(lat=33.5, lon=114)
        ## read the ncfile
        # wrf_file = '/mnt/zfm_18T/fengxiang/HeNan/Data/1900_90m/wrfout_d04_2021-07-20_08:00:00'
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
        x = u.cross_line_idx.values[::10]
        u = u[::3,::10]
        v = v[::3,::10]
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
        colordict=['white','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#颜色列表
        # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
        # colorlevel=np.array([-120, -50, -30, -10, -5, -1, 1, 5, 10, 30, 50, 120])#雨量等级
        # colorlevel=np.arange(0,1,0.01)
        colorlevel=[0, 0.01, 0.02, 0.05, 0.1, 1,  5, 10, 100]
        # colorlevel=np.arange(0,10,0.5)
        # colorticks=[-30, -20, -10, -5, -1, 1, 5, 10, 20, 30]#雨量等级
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
def draw_cross(flnm):
    # flnm  = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m_OGWD/wrfout_d04_2021-07-20_06:00:00'
    ds = xr.open_dataset(flnm)

    cd = CrossData(flnm)
    u = cd.get_vcross('DTAUX3D_LS')
    v = cd.get_vcross('DTAUY3D_LS')
    da = xr.ufuncs.sqrt(u**2+v**2)
    ter_line = cd.get_ter()
    dv = DrawVertical()
    fig = plt.figure(figsize=(10,8), dpi=400)
    ax_cross = fig.add_axes([0.2, 0.2, 0.75, 0.7])
    title_t = ds.XTIME.dt.strftime('%d-%H').values[0]
    ax_cross.set_title(title_t, loc='left', fontsize=18)

    dv.draw_contourf(fig, ax_cross, dv.drop_na(da)*10**3, ter_line)
    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_gwd3/gwd_cross+%s_33.5.png'%title_t)
# %%

def main():
    fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3-CTL/wrfout_d01_' # 2021-07-20_05:00:00'
    # fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3-BL/wrfout_d04_' # 2021-07-20_05:00:00'
    # tt = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='1H')
    # tt = pd.date_range('2021-07-20 05', '2021-07-20 06', freq='1H')
    tt = pd.date_range('2021-07-20 05', '2021-07-20 06', freq='1H')
    for t in tt:
        flnm = fl_path+t.strftime('%Y-%m-%d_%H:%M:%S')
        draw_cross(flnm)

# %%

fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3-CTL/wrfout_d01_' # 2021-07-20_05:00:00'
tt = pd.date_range('2021-07-20 05', '2021-07-20 06', freq='1H')
for t in tt:
    flnm = fl_path+t.strftime('%Y-%m-%d_%H:%M:%S')
ds = xr.open_dataset(flnm)
wrfin = Dataset(flnm)

var_list = ['DTAUX3D_LS', 'DTAUX3D_SS', 'DTAUX3D_FD', 'DTAUX3D_BL']


def get_da_vertical(var):
    # da = ds['DTAUX3D_LS'].squeeze()
    da = ds[var].squeeze()
    area = {
        'lon1':110.5,
        'lon2':116,
        'lat1':32,
        'lat2':36.5,
        'interval':0.125,
    }
    x1, y1 = wrf.ll_to_xy(wrfin,area['lat1'], area['lon1']).values
    x2, y2 = wrf.ll_to_xy(wrfin,area['lat2'], area['lon2']).values
    dda = da[:,y1:y2, x1:x2]
    ddda = dda.mean(dim={'south_north', 'west_east'})
    return ddda

dds = xr.Dataset()
for var in var_list:
    dds[var] = get_da_vertical(var)
dds
# %%
# da = dds['DTAUX3D_LS']
# da
# %%
import matplotlib.pyplot as plt
fig = plt.figure(figsize=[12,12])
ax = fig.add_axes([0.1,0.1, 0.8,0.8])
# LS = dds['DTAUX3D_LS']
# SS = dds['DTAUX3D_SS']
# BL = dds['DTAUX3D_BL']
# FD = dds['DTAUX3D_FD']
# ax.plot(LS.values, LS.bottom_top, label='LS')
# ax.plot(SS.values, SS.bottom_top, label='SS')
# ax.plot(BL.values, BL.bottom_top, label='BL')
# ax.plot(FD.values, FD.bottom_top, label='FD')

for var in var_list:
    dav = dds[var]
    ax.plot(dav.values, dav.bottom_top, label=var, lw=5)




# ax.set_xlim(0,0.00020)
ax.set_xscale('symlog')
ax.set_xlim(10**(-10),10**(-5))
ax.set_ylim(0,50)
# ax.set_xticks(np.linspace(10**(-10),10**(-5), 10))
# ax.set_yscale('symlog')
# ax.legend()

# ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
# ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
# ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
# ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
ax.legend(fontsize=20, edgecolor='white')
ax.set_xlabel('Drag(m/s^2)', fontsize=20)
ax.set_ylabel('Model Level', fontsize=20)
fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/aa.png')


    
    
    

# da.sel()
# da.sel(south_north=slice(y1,y2))
# da[:,:,1:2]
# x2


# index = ((da.lat<=area['lat2']) & (da.lat>=area['lat1']) & (da.lon>=area['lon1']) & (da.lon<=area['lon2']))
# da_obs = da.loc[:,index]  # 这里时间维度在前


# %%
if __name__ == '__main__':
    # main()
    pass
    