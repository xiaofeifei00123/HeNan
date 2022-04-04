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
from baobao.map import Map

# import Figure7_multi.GetData as raingd
# from Figure7_multi import GetData as gd1
import rain_wind_multi as rw
import draw_distance_height_cross as dh


# %%



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
                    hspace=0.1)  # 两子图之间的高

num = 6
axes = [None]*num

i = 0
while i < 6:
    axes[i] = fig.add_subplot(grid[i], projection=proj)
    axes[i+1] = fig.add_subplot(grid[i])
    i = i+2
    
    
    
    
