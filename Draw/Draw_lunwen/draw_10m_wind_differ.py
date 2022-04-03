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
from baobao.map import Map
import os


from baobao.caculate import interp
import draw_10m_wind as d10

def get_wind_minus(t='2021-07-20 00'):
    # flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/Typhoon/weak_typhoon/10m_wind_station.nc'
    gd = d10.GetData(t)
    gd.get_wind_wrf()
    flnm_model1 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'
    u1,v1 = gd.get_wind_wrf(flnm_model1)  # gwd0的风速

    # model_list = ['gwd0', 'gwd1', 'gwd3','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    # model_list = ['gwd0',  'gwd3']
    model_list = ['gwd0',  'gwd3']
    flnm_model2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc'
    # fname = os.path.join(fpath, 'gwd3')
    # flnm = os.path.join(fname, '10m_wind_station.nc')
    u2,v2 = gd.get_wind_wrf(flnm_model2)
    ## 作差 
    u = u2-u1 
    v = v2-v1
    u = u.assign_coords({'lat':('sta',u1.lat.values), 'lon':('sta', u1.lon.values)})
    v = v.assign_coords({'lat':('sta',v1.lat.values), 'lon':('sta', v1.lon.values)})
    return u,v

def draw_wrf():
    u, v = get_wind_minus()
    pic_dic = {'model':'gwd3'}
    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm, 7*cm], dpi=300)
    ax = fig.add_axes([0.15,0.15,0.8,0.8], projection=ccrs.PlateCarree())
    dr = d10.DrawWind(fig, ax, scale=10, Ulength=2)
    # dw = d10.DrawWind()
    dr.draw(u,v,pic_dic)

    fig_name = '10mwind_minus'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save, pad_inches=0)



# %%
if __name__ == '__main__':
    # u,v = get_wind_wrf()
    # draw(u,v)
    # main()
    draw_wrf()
    