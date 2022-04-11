#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2022/04/04 19:07:43
Author           :Forxd
Version          :1.0
'''
# %%
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmaps
from pyproj import transform
from baobao.map import Map
import numpy as np

# import Figure7_multi.GetData as raingd
# from Figure7_multi import GetData as gd1
import draw_rain_wind as rw
import draw_distance_height_cross as dh


# %%
def set_colorbar(fig, cf, ax, dr):
    cb = fig.colorbar(
    cf,
    ax=ax,
    orientation='horizontal',
    ticks=dr.colorlevel[1:-1],
    fraction = 0.05,  # 色标大小,相对于原图的大小
    pad=0.08,  #  色标和子图间距离
            )
    cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小



cm = 1/2.54
# fig = plt.figure(figsize=(19*cm, 8*cm), dpi=300)
fig = plt.figure(figsize=(19*cm, 24*cm), dpi=600)
proj = ccrs.PlateCarree()  # 创建坐标系
grid = plt.GridSpec(3,
                    2,
                    figure=fig,
                    left=0.1,
                    right=0.95,
                    bottom=0.05,
                    top=0.98,
                    wspace=0.2, # 两子图之间的宽
                    # hspace=0.25)
                    hspace=0.2)  # 两子图之间的高

num = 6
axes = [None]*num

i = 0
while i < 6:
    axes[i] = fig.add_subplot(grid[i], projection=proj)
    axes[i+1] = fig.add_subplot(grid[i+1])
    i = i+2
    
    

gd = rw.GetData()
rain1, u1, v1 = gd.get_data_hourly('gwd0')
rain2, u2, v2 = gd.get_data_hourly('gwd3')
rain3, u3, v3 = gd.get_data_minus()
rain1.max()

ax1 = rw.add_map(axes[0])
ax2= rw.add_map(axes[2])
ax3= rw.add_map(axes[4])

ax1.set_title('(a)', loc='left', y=0.98, fontsize=9)
ax2.set_title('(b)', loc='left', y=0.98, fontsize=9)
ax3.set_title('(c)', loc='left', y=0.98, fontsize=9)

dr1 = rw.Draw()
cf1 = dr1.draw_contourf(rain1, ax1)
dr1.draw_quiver(u1,v1, ax1)

dr2 = rw.Draw()
cf2 = dr2.draw_contourf(rain2, ax2)
dr2.draw_quiver(u2,v2, ax2)

# dr.colorlevel = []
dr3 = rw.Draw()
dr3.colordict=['#0000fb','#3232fd','#6464fd','white','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
dr3.colorlevel= [-90,-20,-10,-3, 0, 3,10,20,90]  # 垂直速度的色标
cf3 = dr3.draw_contourf(rain3, ax3)
dr3.draw_quiver(u3,v3, ax3, scale=80, ulength=5)



set_colorbar(fig, cf1, axes[0], dr1)
set_colorbar(fig, cf2, axes[2], dr2)
set_colorbar(fig, cf3, axes[4], dr3)



## 距离高度图
flnm1 ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_zhengzhou.nc'
flnm2 ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/cross_zhengzhou.nc'
gd = dh.GetData()
dic1 = gd.get_data(t='2021-07-20 00', flnm=flnm1)
dic2 = gd.get_data(t='2021-07-20 00', flnm=flnm2)


dr0 = dh.Draw()
dr0.draw_contour(axes[1], dic1['theta_e'])
cf0 = dr0.draw_contourf(axes[1], dic1['w'], dic1['ter_line'])
dr0.draw_quiver(axes[1],dic1['hor'],dic1['ver'], scale=60)
dr0.set_ticks(axes[1], dic1['theta_e'])

dr1 = dh.Draw()
dr1.draw_contour(axes[3], dic2['theta_e'])
cf1 = dr1.draw_contourf(axes[3], dic2['w'], dic2['ter_line'])
dr1.draw_quiver(axes[3],dic2['hor'],dic2['ver'], scale=60)
dr1.set_ticks(axes[3], dic1['theta_e'])


dr2 = dh.Draw()
levels=np.arange(-12, 12+1, 4)
dr2.draw_contour(axes[5], gd.minus(dic1['theta_e'],dic2['theta_e']), levels=levels)
cf2 = dr2.draw_contourf(axes[5], gd.minus(dic1['w'],dic2['w']), dic1['ter_line'])
dr2.draw_quiver(axes[5],gd.minus(dic1['hor'],dic2['hor']),gd.minus(dic1['ver'],dic2['ver']), scale=60)
dr2.set_ticks(axes[5], gd.minus(dic1['theta_e'],dic2['theta_e']))




axes[1].set_title('(d)', loc='left', y=0.98, fontsize=9)
axes[3].set_title('(e)', loc='left', y=0.98, fontsize=9)
axes[5].set_title('(f)', loc='left', y=0.98, fontsize=9)
set_colorbar(fig, cf0, axes[1], dr0)
set_colorbar(fig, cf1, axes[3], dr1)
set_colorbar(fig, cf2, axes[5], dr2)




## 标注直线AB
# y = np.linspace(33.8, 33.5, 10)
# x = np.linspace(111.7, 113.2, 10)
# axes[4].plot(x,y, transform=ccrs.PlateCarree(), color='black', linewidth=2)
# axes[4].text(x[0], y[0], 'A')
# axes[4].text(x[-1], y[-1], 'B')




fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/Figure8_zz.png')

