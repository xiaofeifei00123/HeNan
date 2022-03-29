from email import header


#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
画多图, 多个子图， 合并图像
调用每张图的画图对象，调用每个图的数据
-----------------------------------------
Time             :2022/03/29 14:24:00
Author           :Forxd
Version          :1.0
'''

from cmath import pi
import sys,os
import xarray as xr
import numpy as np
import pandas as pd

# import salem  # 插值
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.io.shapereader import Reader, natural_earth
import matplotlib as mpl
from matplotlib.path import Path
import seaborn as sns
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import geopandas
import cmaps
from get_cmap import get_cmap_rain2
from multiprocessing import Pool
from baobao.map import Map

import draw_rain_distribution_24h as dd
import draw_rain_distribution_minus_24h as dm

# %%



if __name__ == '__main__':
    
    cm = round(1/2.54, 2)
    proj = ccrs.PlateCarree()  # 创建坐标系
    fig = plt.figure(figsize=(17*cm, 28*cm), dpi=600)
    # ax = fig.add_axes([0.1,0.08,0.85,0.85], projection=proj)
    # fig = plt.figure(figsize=(21, 20), dpi=400)  # 创建页面
    grid = plt.GridSpec(4,
                        2,
                        figure=fig,
                        left=0.05,
                        right=0.98,
                        bottom=0.12,
                        top=0.97,
                        wspace=0.1,
                        # hspace=0.25)
                        hspace=0.2)

    num = 8
    axes = [None] * num  # 设置一个维度为8的空列表
    for i in range(num):
        axes[i] = fig.add_subplot(grid[i], projection=proj)

    ### 画降水分布图
    gd1 = dd.GetData()
    ## 观测降水
    dr1 = dd.Draw(fig, axes[0])
    da = gd1.obs()
    cf = dr1.draw_tricontourf(da)    
    dr1.ax.set_title('(a)', loc='left', y=0.85)
    ## EC
    da = gd1.EC()
    dr1 = dd.Draw(fig, axes[2])
    cf = dr1.draw_single(da)    
    dr1.ax.set_title('(c)', loc='left', y=0.85)
    ## model
    da = gd1.onemodel('gwd0')
    dr1 = dd.Draw(fig, axes[4])
    cf = dr1.draw_single(da)    
    # axes[4].set_title('(e)')
    dr1.ax.set_title('(e)', loc='left', y=0.85)
    ## model
    da = gd1.onemodel('gwd3')
    dr1 = dd.Draw(fig, axes[6])
    cf = dr1.draw_single(da)    
    # axes[6].set_title('(f)')
    dr1.ax.set_title('(g)', loc='left', y=0.85)
    ## 用单独一个子图来画色标
    ax1 = fig.add_axes([0.07,0.06,0.40,0.02])
    colorlevel=[0, 0.1, 10, 25.0, 50, 100, 250,  700]#雨量等级
    colorticks = colorlevel[1:-1]
    cb = fig.colorbar(
        cf,
        cax=ax1,
        orientation='horizontal',
        ticks=colorticks,
    )
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小

    ### 画降水差值分布图
    gd2 = dm.GetData()

    ### gwd3-gwd0
    da = gd2.model_model()
    dr2 = dm.Draw(fig, axes[1])
    cf = dr2.draw_single(da)
    dr2.ax.set_title('(b)', loc='left', y=0.85)
    

    da = gd2.EC_OBS()
    dr2 = dm.Draw(fig, axes[3])
    cf = dr2.draw_single(da)
    dr2.ax.set_title('(d)', loc='left', y=0.85)
    
    ### gwd3-gwd0
    da = gd2.model_model()
    dr2 = dm.Draw(fig, axes[1])
    cf = dr2.draw_single(da)
    dr2.ax.set_title('(b)', loc='left', y=0.85)
    
    da = gd2.model_obs('gwd0')
    dr2 = dm.Draw(fig, axes[5])
    cf = dr2.draw_single(da)
    dr2.ax.set_title('(f)', loc='left', y=0.85)
    
    da = gd2.model_obs('gwd3')
    dr2 = dm.Draw(fig, axes[7])
    cf = dr2.draw_single(da)
    dr2.ax.set_title('(h)', loc='left', y=0.85)

    ax2 = fig.add_axes([0.55,0.06,0.40,0.02])
    colorlevel=[-700, -200, -100, -50, -20, 20, 50 , 100, 200,700 ]#雨量等级
    colorticks = colorlevel[1:-1]
    cb = fig.colorbar(
        cf,
        cax=ax2,
        orientation='horizontal',
        ticks=colorticks,
        # fraction = 0.05,  # 色标大小,相对于原图的大小
        # pad=0.1,  #  色标和子图间距离
    )
    
    cb.ax.tick_params(labelsize=10)  # 设置色标标注的大小
    
    fig_name = 'multi'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig.savefig(fig_path+fig_name)