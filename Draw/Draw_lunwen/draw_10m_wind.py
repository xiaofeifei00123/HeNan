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

class Wind():
    pass

class DrawWind():

    def __init__(self, fig, ax, scale=40, Ulength=10):
        self.fig = fig 
        self.ax = ax
        # self.time = '2021-07-20 00'
        self.scale = scale   # 风矢的尺度
        self.Ulength =  Ulength  # 风矢参考的大小
        self.map_dic = {
                'proj':ccrs.PlateCarree(),
                'extent':[110.5, 116, 32, 36.5],
                'extent_interval_lat':1,
                'extent_interval_lon':1,
                }
        self.station = {
                'ZhengZhou': {
                    'abbreviation':'ZZ',
                    'lat': 34.76,
                    'lon': 113.65
                },
                'NanYang': {
                    'abbreviation':'NY',
                    'lat': 33.1,
                    'lon': 112.49,
                },
                'LuShi': {
                    'abbreviation':'LS',
                    'lat': 34.08,
                    'lon': 111.07,
                },
            }

    def draw_quiver(self, u,v, ax):
        '''
        绘制风矢图
        '''
        y = u.lat.values
        x = u.lon.values
        # Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=40,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        # qk = ax.quiverkey(Q, X=0.75, Y=0.08, U=10, label=r'($10 m/s$)', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢
        Q = ax.quiver(x, y, u.values,v.values,units='inches',scale=self.scale,pivot='middle', transform=ccrs.PlateCarree())  # 绘制风矢
        qk = ax.quiverkey(Q, X=0.75, Y=0.05, U=self.Ulength, label=r'(${} m/s$)'.format(self.Ulength), labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢
        # qk = ax.quiverkey(Q, x=0.5, y=0.08, U=self.Ulength, label=r'(${} m/s$)'.format(self.Ulength), labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢
        # qk = ax.quiverkey(Q, X=1.55, Y=0.05, U=10, label=r'$(\overrightarrow{qv_f}-\overrightarrow{qv_o}, 100\ g/kg \cdot m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':25})   # 设置参考风矢

    def draw(self, u,v,pic_dic={'model':'obs'}):

        fig = self.fig
        ax = self.ax
        mp = Map()
        ax = mp.create_map(ax, self.map_dic)
        ax.set_extent(self.map_dic['extent'])

        mp.add_station(ax, self.station, justice=True, delx=-0.1)
        if pic_dic['model'] == 'gwd0':
            pic_dic['model'] = 'no-gwd'
        # ax.set_title(pic_dic['model'], fontsize=10,loc='left')
        self.draw_quiver(u,v,ax)
        # fig_name = '10mwind_'+pic_dic['model']+pic_dic['time']+'typhoon'
        # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_10mwind/'
        # fig_save = os.path.join(fig_path, fig_name)
        # # fig.savefig(fig_save, bbox_inches='tight')
        # fig.savefig(fig_save,bbox_inces='tight', pad_inches=0)


class GetData():
    def __init__(self,time = '2021-07-20 0000'):
        self.time = time
        pass

    def get_wind_wrf(self,flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc',):
        """计算和筛选数据

        Args:
            flnm (str, optional): _description_. Defaults to '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc'.
            t (str, optional): _description_. Defaults to '2021-07-20 00'.

        Returns:
            _type_: _description_
        """
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3-test/10m_wind_station.nc'
        ds = xr.open_dataset(flnm)
        ds = ds.swap_dims({'time':'Time'})
        ds = ds.drop_vars('time')
        ds = ds.rename({'Time':'time'})
        # t = '2021-07-20 06'
        ds1 = ds.sel(time=self.time)
        da = ds1['uvmet10']
        u = da.sel(u_v='u')
        v = da.sel(u_v='v')
        return u,v

    def get_wind_obs(self,flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc', ):
        # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
        ds = xr.open_dataset(flnm)
        # t = '2021-07-20 06'
        ds1 = ds.sel(time=self.time)
        u = ds1['u']
        v = ds1['v']
        return u,v

def draw_wrf(t = '2021-07-20 12'):
    # model_list = ['gwd0', 'gwd1', 'gwd3','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    model_list = ['gwd0', 'gwd3']
    # model_list = ['weak_typhoon', 'strengthen_typhoon']
    fpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    gd = GetData(time=t)
    for model in model_list:
        fname = os.path.join(fpath, model)
        flnm = os.path.join(fname, '10m_wind_station.nc')
        u,v = gd.get_wind_wrf(flnm)
        pic_dic = {'model':model, 'time':t}
        # draw(u,v,pic_dic)
        cm = 1/2.54
        fig = plt.figure(figsize=[8*cm, 7*cm], dpi=300)
        ax = fig.add_axes([0.15,0.15,0.8,0.8], projection=ccrs.PlateCarree())
        dr = DrawWind(fig, ax)
        dr.draw(u,v, pic_dic)
        fig_name = '10mwind'+model
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
        fig_save = os.path.join(fig_path, fig_name)
        fig.savefig(fig_save, pad_inches=0)


def draw_obs(t = '2021-07-20 06'):
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
    t = '2021-07-20 06'
    gd = GetData(time=t)
    u,v = gd.get_wind_obs(flnm)
    pic_dic = {'model':'obs', 'time':t}
    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm, 7*cm], dpi=300)
    ax = fig.add_axes([0.15,0.15,0.8,0.8], projection=ccrs.PlateCarree())
    dr = DrawWind(fig, ax)
    dr.draw(u,v, pic_dic)

    fig_name = '10mwind'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save, pad_inches=0)


def main():
    t = '2021-07-20 12'
    draw_obs(t)
    draw_wrf(t)

# %%
if __name__ == '__main__':
    main()