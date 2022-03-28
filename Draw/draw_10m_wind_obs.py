#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
插值后观测的10m风场
-----------------------------------------
Time             :2021/12/30 23:37:35
Author           :Forxd
Version          :1.0
'''

# %%
from time import strftime
from cartopy.crs import Projection
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from pandas._libs.tslibs.timestamps import Timestamp
import xarray as xr
import meteva.base as meb
import numpy as np
import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
# import matplotlib as mpl
import cartopy.crs as ccrs
# import cmaps
# from wrf import constants, get_cartopy, smooth2d, getvar
# import metpy.calc as ca
# import metpy
# from metpy import calc as ca
# from metpy import calc as ca  # calc是一个文件夹
# from metpy.units import units  # units是units模块中的一个变量名(函数名？类名？)
# from metpy import constants  # constatns是一个模块名
# from pyproj import CRS
# from caculate_diag import Qv, QvDiv  ## 计算水汽通量散度的类
# from baobao.caculate import Qv, QvDiv


# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
ds = xr.open_dataset(flnm)
# %%
ds
# %%
t = '2021-07-20 12'
ds1 = ds.sel(time=t)
u = ds1['u']
v = ds1['v']
v

# %%
# %%

def draw_quiver(u,v, ax):
    '''
    绘制风矢图
    '''
    # u = u[::5,::5]
    # v = v[::5,::5]
    # y = u.coords['lat']
    y = u.lat.values
    x = u.lon.values
    # print(y)
    # print(type(u))
    # print(type(y))
    # x = u.coords['lon']
    # ax = self.ax
    # Q = ax.quiver(x,y,u,v,units='inches',scale=18,pivot='middle')  # 绘制风矢
    # Q = ax.quiver(x, y, u,v,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=18,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=25,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
    # qk = ax.quiverkey(Q, X=0.55, Y=0.15, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=0.55, Y=0.12, U=10, label=r'($10 ,g \cdot cm^{-1} \cdot hPa^{-1} \cdot s^{-1}$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    qk = ax.quiverkey(Q, X=0.75, Y=0.12, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':22})   # 设置参考风矢
    # qk = ax.quiverkey(Q, X=1.55, Y=0.05, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢

fig = plt.figure(figsize=[10,8])
ax = fig.add_axes([0.15,0.05,0.8,0.9], projection=ccrs.PlateCarree())
mb = mapview.BaseMap()
# mb.set_extent('中国陆地')
mb.drawcoastlines(linewidths=0.8, alpha=0.5)

# tt = (dic['time']).strftime('%Y-%m-%d %H%M')
# print(tt)
# ax.set_title(tt, loc='center', fontsize=25)
# ax.set_title(dic['model'], loc='left', fontsize=25)
# ax.set_title('OBS', loc='left', fontsize=25)
# ax.set_title(str(dic['level'])+'hPa', loc='right', fontsize=25)

mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
# mb.set_extent('中国陆地')
# mb.set_extent([107, 135, 20,40])
mb.set_extent([110, 116, 32, 36])
draw_quiver(u,v,ax)
# %%