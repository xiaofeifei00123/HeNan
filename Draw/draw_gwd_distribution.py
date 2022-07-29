#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
拖曳力的水平分布
-----------------------------------------
Time             :2021/09/27 15:45:32
Author           :Forxd
Version          :1.0
'''
# %%
from baobao.get_cmap import get_rgb

flnm_rgb = '/mnt/zfm_18T/fengxiang/HeNan/Draw/table/8colors.rgb'
rgb = get_rgb(flnm_rgb)

# %%
import sys,os
import xarray as xr
import numpy as np
import pandas as pd
import netCDF4 as nc
import wrf

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


from draw_10m_wind import GetData, DrawWind

# %%
def draw_wrf(t = '2021-07-20 12'):
    # model_list = ['gwd0', 'gwd1', 'gwd3','gwd3-FD', 'gwd3-BL','gwd3-SS', 'gwd3-LS']
    # model_list = ['gwd0', 'gwd3']
    model_list = ['gwd3']
    # model_list = ['weak_typhoon', 'strengthen_typhoon']
    fpath = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    gd = GetData(time=t)
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_12:00:00'

    ds = xr.open_dataset(flnm)
    # ds = ds.rename({'XTIME':'time'})
    # t = '2021-07-20 06'
    # ds1 = ds.sel(time=t)
    print(ds)
    # u = da.sel(u_v='u')
    # v = da.sel(u_v='v')
    
    
    
    

    # pic_dic = {'model':'tttt', 'time':t}
    # # draw(u,v,pic_dic)
    # cm = 1/2.54
    # fig = plt.figure(figsize=[8*cm, 7*cm], dpi=300)
    # ax = fig.add_axes([0.15,0.15,0.8,0.8], projection=ccrs.PlateCarree())
    # dr = DrawWind(fig, ax)
    # dr.draw(u,v, pic_dic)
    # # fig_name = '10mwind'
    # # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
    # # fig_save = os.path.join(fig_path, fig_name)
    # # fig.savefig(fig_save, pad_inches=0)
# draw_wrf()

# # %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_12:00:00'
# ds = xr.open_dataset(flnm)
# wrfnc = nc.Dataset(flnm)
# uv = wrf.getvar(wrfnc, 'uvmet10')
# u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
# # %%
# # uv1[0]
# v1
# cm = 1/2.54
# fig = plt.figure(figsize=[8*cm, 7*cm], dpi=300)
# ax = fig.add_axes([0.15,0.15,0.8,0.8], projection=ccrs.PlateCarree())

# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_12:00:00'

# ds = xr.open_dataset(flnm)
# wrfnc = nc.Dataset(flnm)
# # uv = wrf.getvar(wrfnc, 'uvmet10')
# # u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
# u = wrf.getvar(wrfnc, 'ua')
# v = wrf.getvar(wrfnc, 'va')
# u1 = u.isel(bottom_top=0).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
# v1 = v.isel(bottom_top=0).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
# dr = DrawWind(fig, ax)
# pic_dic = {'model':'tttt', 'time':'2021-07-20 00'}
# dr.draw(u1[::10, ::10],v1[::10, ::10], pic_dic)
# # %%
# fl_path ='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_12:00:00'
# ds = xr.open_dataset(fl_path)
# da = ds['DTAUX3D_FD'].sel(bottom_top=0)
# da
# da.max()
# ds['DUSFCG_SS'].max()

# %%
class Draw(object):

    def __init__(self) -> None:
        super().__init__()
        # self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/Province_9/Province_9.shp'
        self.path_province = '/mnt/zfm_18T/fengxiang/DATA/SHP/Province_shp/henan.shp'
        self.path_henan = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/henan.shp'
        # self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/Map/cn_shp/City_9/City_9.shp'
        self.path_city = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_henan/zhenzhou/zhenzhou_max.shp'
        self.path_tibet = '/mnt/zfm_18T/fengxiang/DATA/SHP/shp_tp/Tibet.shp'
        self.picture_path = '/mnt/zfm_18T/fengxiang/Asses_PBL/Rain/picture'


    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        ## 指定颜色
        # colorlevel=[0, 0.1, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 70]#雨量等级
        
        # colordict=['white','#A6F28F','#3DBA3D','#61BBFF','#0000FF','#FA00FA','#800040', '#EE0000']#降水颜色列表
        # colordict=['#0000fb','#3232fd','#6464fd','#a2a3fb','white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红

        
        # colorlevel=np.arange(-0.5, 0.6, 0.05)
        # colorlevel=[-0.45, -0.35, -0.25, -0.15,-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        # colorlevel=[-100, -0.35, -0.25, -0.15,-0.05, 0.05, 0.15, 0.25, 0.35, 100]
        # colorlevel=[-100,-0.1,-0.01,-0.001,-0.0001, 0.0001, 0.001, 0.01, 0.1, 100]
        colordict = rgb  # 8个颜色
        # colorlevel=[0.0001, 0.001, 0.01, 0.1,1, 10, 20,100]
        # colordict=['white','#fbbcbc', '#ff8383', '#fd4949', '#fd0000']#正负, 蓝-红
        colorlevel=[0, 0.0001, 0.001, 0.01, 0.1,100]
        x = data.lon
        y = data.lat
        
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                          levels=colorlevel,
                          colors = colordict,
                        #   levels=6,
                        #   cmap
                        #   cmap=cmaps.precip3_16lev,
                        #   cmap=cmaps.WhiteBlueGreenYellowRed,
                        #   cmap=cmaps.NCV_blue_red,
                          transform=ccrs.PlateCarree())
        return crx

    def draw_single(self, da, picture_dic, flnm):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        cm = 1/2.54
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.1,0.1,0.85,0.85], projection=proj)
        # ax.set_extent([])
        print("画{}时刻的图".format(str(picture_dic['date'])))
        date = picture_dic['date']
        dic = {'name':'HeNan',
               'cmap':cmaps.precip3_16lev,
               'time':str(date)}
               
        # ax = self.create_map(ax)
        mp = Map()
        map_dic = {
            'proj':ccrs.PlateCarree(),
            'extent':[110.5, 116, 32, 36.5],
            'extent_interval_lat':1,
            'extent_interval_lon':1,
        }

        ax = mp.create_map(ax, map_dic)
        ax.set_extent(map_dic['extent'])

        station = {
            'ZhengZhou': {
                'abbreviation':'ZZ',
                'lat': 34.76,
                'lon': 113.65
            },
        }
        mp.add_station(ax, station, justice=True)
        ax.set_title(date, fontsize=10,)
        ax.set_title(picture_dic['initial_time'], fontsize=10,loc='left')
        ax.set_title(picture_dic['type'], fontsize=10,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)
        # colorticks=[0.1,5,15,30.0,70,140]#雨量等级
        # colorlevel=[0, 0.1, 5, 10, 15.0, 20, 25, 30, 700]#雨量等级
        # colorlevel=[-0.45, -0.35, -0.25, -0.15,-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        colorlevel=[-100,-0.1,-0.01,-0.001,-0.0001, 0.0001, 0.001, 0.01, 0.1, 100]
        colorticks = colorlevel[1:-1]
        
        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小
        fig_name = picture_dic['type']+'_'+picture_dic['initial_time']+'_'+picture_dic['date']+'gai'

        
        ds = xr.open_dataset(flnm)
        wrfnc = nc.Dataset(flnm)
        # uv = wrf.getvar(wrfnc, 'uvmet10')
        # u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })

        u = wrf.getvar(wrfnc, 'ua')
        v = wrf.getvar(wrfnc, 'va')
        u1 = u.isel(bottom_top=0).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
        v1 = v.isel(bottom_top=0).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
        

        dr = DrawWind(fig, ax)
        pic_dic = {'model':'tttt', 'time':'2021-07-20 00'}
        dr.draw(u1[::10, ::10],v1[::10, ::10], pic_dic)
        

        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_gwd3/new_modify/'
        fig.savefig(fig_path+fig_name)


def get_gwd(ds, var_flag, bt=0):
    if var_flag == 'LS':
        da1 = ds['DTAUX3D_LS'].sel(bottom_top=bt)
        da2 = ds['DTAUY3D_LS'].sel(bottom_top=bt)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'SS':
        da1 = ds['DTAUX3D_SS'].sel(bottom_top=bt)
        da2 = ds['DTAUY3D_SS'].sel(bottom_top=bt)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'FD':
        da1 = ds['DTAUX3D_FD'].sel(bottom_top=bt)
        da2 = ds['DTAUY3D_FD'].sel(bottom_top=bt)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'BL':
        da1 = ds['DTAUX3D_BL'].sel(bottom_top=bt)
        da2 = ds['DTAUY3D_BL'].sel(bottom_top=bt)
        da = np.sqrt(da1**2+da2**2)
    else:
        da1 = ds['DTAUX3D_LS'].sel(bottom_top=bt)+ds['DTAUX3D_SS'].sel(bottom_top=bt)+ds['DTAUX3D_FD'].sel(bottom_top=bt)+ds['DTAUX3D_BL'].sel(bottom_top=bt)
        da2 = ds['DTAUY3D_LS'].sel(bottom_top=bt)+ds['DTAUY3D_SS'].sel(bottom_top=bt)+ds['DTAUY3D_FD'].sel(bottom_top=bt)+ds['DTAUY3D_BL'].sel(bottom_top=bt)
        da = np.sqrt(da1**2+da2**2)
    da = da*10**2
    return da



def draw_one(model):
    pass

    fl_path =os.path.join('/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify', model)
    # tt = pd.date_range('2021-07-20 00', '2021-07-21 00', freq='1H')
    # tt = pd.date_range('2021-07-20 12', '2021-07-20 12', freq='1H')
    tt = pd.date_range('2021-07-19 13', '2021-07-21 00', freq='1H')
    # tt = pd.date_range('2021-07-19 13', '2021-07-19 13', freq='1H')
    # tt = pd.date_range('2021-07-20 12', '2021-07-20 12', freq='1H')
    for t in tt:
        flnm = 'wrfout_d03_'+t.strftime('%Y-%m-%d_%H:%M:%S')
        flnm = os.path.join(fl_path, flnm)
    
        ds = xr.open_dataset(flnm)
        da = get_gwd(ds, 'all', 0)

        # da = da1
        # var = 'DTAUX'
        var = ''
        # print(da.max())
        gwd_sfc = da.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'}).squeeze()
        dr = Draw()
        picture_dic = {'date':gwd_sfc.time.dt.strftime("%Y%m%d-%H").values, 'type':model, 'initial_time':var}
        dr.draw_single(gwd_sfc, picture_dic, flnm)

def draw_dual():
    # model_list = ['gwd1', 'gwd3']
    # model_list = ['SS',]
    model_list = ['GWD3',]
    for model in model_list:
        draw_one(model)


if __name__ == '__main__':

    # draw_one()
    draw_dual()
    # draw_obs()
    # draw_forecast()
    # gd = GetData()
    # dd = gd.get_rain_obs()
    # gd = GetData()
    # da = gd.get_rain_wrf()
    # dr = Draw()
    # dr.draw_single(da, '2000_2012')
# %%
