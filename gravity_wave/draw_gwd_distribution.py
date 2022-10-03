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

# flnm_rgb = '/mnt/zfm_18T/fengxiang/HeNan/Draw/table/8colors.rgb' # rgb = get_rgb(flnm_rgb)

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
# import matplotlib.patches as patches
import matplotlib.pyplot as plt
import cmaps
# from get_cmap import get_cmap_rain2
from multiprocessing import Pool
from baobao.map import Map
from draw_10m_wind import GetData, DrawWind
from draw_rain_distribution_24h import Rain
# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify/GWD3/wrfout_d03_2021-07-19_12:00:00'
# # nc.Dataset(flnm)
# # ds = xr.open_dataset(flnm)
# wrfnc = nc.Dataset(flnm)
# u, v = wrf.getvar(wrfnc, 'uvmet10')
# # aa[0]
# # u
# v


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
        # self.colorlevel=[0, 0.1, 1, 5, 10,100]
        # self.colorlevel=[0, 10,  20,50,100, 1000]
        self.colorlevel=[0, 10,50,100,200, 2000]
        self.colorticks=self.colorlevel[1:-1]
        flnm_rgb = '/mnt/zfm_18T/fengxiang/HeNan/Draw/table/8colors.rgb'
        rgb = get_rgb(flnm_rgb)
        self.colordict = rgb  # 8个颜色


    def draw_contourf_single(self, data, ax, dic):
        """画填色图
        """

        ## 指定颜色
        x = data.lon
        y = data.lat
        
        
        crx = ax.contourf(x,
                          y,
                          data,
                          corner_mask=False,
                          levels=self.colorlevel,
                          colors = self.colordict,
                        #   levels=6,
                        #   cmap
                        #   cmap=cmaps.precip3_16lev,
                        #   cmap=cmaps.WhiteBlueGreenYellowRed,
                        #   cmap=cmaps.NCV_blue_red,
                          transform=ccrs.PlateCarree())
        return crx

    def draw_single(self, da, picture_dic, flnm, bt=0, level=500):
        """画单个的那种图

        Args:
            da (DataArray): 单个时次的降水
        """
        cm = 1/2.54
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
        proj = ccrs.PlateCarree()  # 创建坐标系
        ax = fig.add_axes([0.11,0.1,0.85,0.85], projection=proj)
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
        # mp.add_station(ax, station, justice=True)
        ax.set_title(date, fontsize=10,)
        ax.set_title(picture_dic['initial_time'], fontsize=10,loc='left')
        ax.set_title(picture_dic['type'], fontsize=10,loc='right')
        cf = self.draw_contourf_single(da, ax, dic)

        cb = fig.colorbar(
            cf,
            # cax=ax6,
            orientation='horizontal',
            ticks=self.colorticks,
            fraction = 0.05,  # 色标大小,相对于原图的大小
            pad=0.08,  #  色标和子图间距离
        )
        cb.ax.tick_params(labelsize=8)  # 设置色标标注的大小

        dr = Rain()
        ax.plot(np.linspace(dr.cross_start[0], dr.cross_end[0], 10), np.linspace(dr.cross_start[1], dr.cross_end[1], 10), color='black')

        # self.draw_wind_bt(fig, ax, bt=0)
        # self.draw_wind_level(fig, ax, level=level)
        
        

        # fig_name = picture_dic['type']+'_'+picture_dic['initial_time']+'_'+picture_dic['date']+'gai'
        # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_gwd3/new_modify/'
        # fig.savefig(fig_path+fig_name)
        return fig, ax

        
    def draw_wind_bt(self,fig, ax, flnm, bt):
        ds = xr.open_dataset(flnm)
        wrfnc = nc.Dataset(flnm)
        # uv = wrf.getvar(wrfnc, 'uvmet10')
        # u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })

        u = wrf.getvar(wrfnc, 'ua')
        v = wrf.getvar(wrfnc, 'va')
        u1 = u.isel(bottom_top=bt).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
        v1 = v.isel(bottom_top=bt).rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })

        dr = DrawWind(fig, ax)
        pic_dic = {'model':'tttt', 'time':'2021-07-20 00'}
        dr.draw(u1[::10, ::10],v1[::10, ::10], pic_dic)
        
    def draw_wind_sfc(self,fig, ax, flnm):
        # ds = xr.open_dataset(flnm)
        wrfnc = nc.Dataset(flnm)
        # uv = wrf.getvar(wrfnc, 'uvmet10')
        # u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })

        # u = wrf.getvar(wrfnc, 'ua')
        # v = wrf.getvar(wrfnc, 'va')
        u, v = wrf.getvar(wrfnc, 'uvmet10')
        u1 = u.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
        v1 = v.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })
        dr = DrawWind(fig, ax)
        pic_dic = {'model':'tttt', 'time':'2021-07-20 00'}
        dr.draw(u1[::10, ::10],v1[::10, ::10], pic_dic)

    def draw_wind_level(self,fig, ax, flnm, level):
        # ds = xr.open_dataset(flnm)
        wrfnc = nc.Dataset(flnm)
        # uv = wrf.getvar(wrfnc, 'uvmet10')
        # u1, v1 = uv.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time', })

        p = wrf.getvar(wrfnc, 'pressure')
        u = wrf.getvar(wrfnc, 'ua')
        v = wrf.getvar(wrfnc, 'va')

        u = wrf.interplevel(u, p, level).squeeze()
        v = wrf.interplevel(v, p, level).squeeze()

        u1 = u.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time', })
        v1 = v.rename({'XLAT':'lat', 'XLONG':'lon', 'Time':'time', })
        # print(u)

        dr = DrawWind(fig, ax)
        pic_dic = {'model':'tttt', 'time':'2021-07-20 00'}
        dr.draw(u1[::10, ::10],v1[::10, ::10], pic_dic)



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
    # da = da*10**2
    da = da*3600  # 一小时的加速度是多少(*3600s), 这里的单位变成m/s
    return da

# def get_gwd_level(flnm, var_flag, level=500):
def get_gwd_level(flnm, var_flag, level=500):

    wrfnc = nc.Dataset(flnm)
    p = wrf.getvar(wrfnc, 'pressure')
    ds = xr.open_dataset(flnm)
    

    if var_flag == 'LS':

        da1 = ds['DTAUX3D_LS']
        da2 = ds['DTAUY3D_LS']
        da1 = wrf.interplevel(da1, p, level)
        da2 = wrf.interplevel(da2, p, level)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'SS':
        da1 = ds['DTAUX3D_SS']
        da2 = ds['DTAUY3D_SS']
        da1 = wrf.interplevel(da1, p, level)
        da2 = wrf.interplevel(da2, p, level)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'FD':
        da1 = ds['DTAUX3D_FD']
        da2 = ds['DTAUY3D_FD']
        da1 = wrf.interplevel(da1, p, level)
        da2 = wrf.interplevel(da2, p, level)
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'BL':
        da1 = ds['DTAUX3D_BL']
        da2 = ds['DTAUY3D_BL']
        da1 = wrf.interplevel(da1, p, level)
        da2 = wrf.interplevel(da2, p, level)
        da = np.sqrt(da1**2+da2**2)
    else:
        da1 = ds['DTAUX3D_LS']+ds['DTAUX3D_SS']+ds['DTAUX3D_FD']+ds['DTAUX3D_BL']
        da2 = ds['DTAUY3D_LS']+ds['DTAUY3D_SS']+ds['DTAUY3D_FD']+ds['DTAUY3D_BL']
        da1 = wrf.interplevel(da1, p, level)
        da2 = wrf.interplevel(da2, p, level)
        da = np.sqrt(da1**2+da2**2)
    # da = da*10**2
    da = da*3600  # 一小时的加速度是多少(*3600s), 这里的单位变成m/s
    return da

def get_gwd_sfc(flnm, var_flag):

    wrfnc = nc.Dataset(flnm)
    p = wrf.getvar(wrfnc, 'pressure')
    ds = xr.open_dataset(flnm)
    

    if var_flag == 'LS':

        da1 = ds['DUSFCG_LS']
        da2 = ds['DVSFCG_LS']
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'SS':
        da1 = ds['DUSFCG_SS']
        da2 = ds['DVSFCG_SS']
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'FD':
        da1 = ds['DUSFCG_FD']
        da2 = ds['DVSFCG_FD']
        da = np.sqrt(da1**2+da2**2)
    elif var_flag == 'BL':
        da1 = ds['DUSFCG_BL']
        da2 = ds['DVSFCG_BL']
        da = np.sqrt(da1**2+da2**2)
    else:
        da1 = ds['DUSFCG_LS']+ds['DUSFCG_SS']+ds['DUSFCG_FD']+ds['DUSFCG_BL']
        da2 = ds['DVSFCG_LS']+ds['DVSFCG_SS']+ds['DVSFCG_FD']+ds['DVSFCG_BL']
        da = np.sqrt(da1**2+da2**2)
    # da = da*10**2
    da = da*3600  # 一小时的加速度是多少(*3600s), 这里的单位变成m/s
    return da

def draw_one(model):
    pass

    fl_path =os.path.join('/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/new_modify', model)
    # tt = pd.date_range('2021-07-20 00', '2021-07-20 00', freq='1H')
    tt = pd.date_range('2021-07-20 18', '2021-07-20 18', freq='1H')
    # tt = pd.date_range('2021-07-19 12', '2021-07-20 12', freq='1H')
    # tt = pd.date_range('2021-07-19 13', '2021-07-21 00', freq='1H')
    # tt = pd.date_range('2021-07-19 13', '2021-07-19 13', freq='1H')
    # tt = pd.date_range('2021-07-20 12', '2021-07-20 12', freq='1H')
    # force_flag = 'LS'   # 拖曳力种类
    # force_flag_list = ['LS', 'BL', 'SS', 'FD', 'ALL']
    force_flag_list = ['ALL']
    # force_flag_list = ['BL']
    for force_flag in force_flag_list:
        for t in tt:
            flnm = 'wrfout_d03_'+t.strftime('%Y-%m-%d_%H:%M:%S')
            flnm = os.path.join(fl_path, flnm)
        
            da = get_gwd_sfc(flnm, force_flag)
            print(da.max().values)
            var = ''
            gwd_sfc = da.rename({'XLAT':'lat', 'XLONG':'lon', 'XTIME':'time'}).squeeze()
            dr = Draw()
            picture_dic = {'date':gwd_sfc.time.dt.strftime("%Y%m%d-%H").values, 'type':model, 'initial_time':var}
            fig, ax = dr.draw_single(gwd_sfc, picture_dic, flnm)

            # dr.draw_wind_bt(fig, ax,flnm, bt=0)
            dr.draw_wind_sfc(fig, ax,flnm)
            ax.set_title(force_flag, fontsize=10,loc='left')

            fig_name = picture_dic['type']+'_'+picture_dic['initial_time']+'_'+picture_dic['date']+force_flag
            fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_gwd/'
            fig.savefig(fig_path+fig_name)

def draw_dual():
    # model_list = ['gwd1', 'gwd3']
    # model_list = ['SS',]
    # model_list = ['GWD3',]
    model_list = ['GWD3','SS','FD', 'CTRL']
    for model in model_list:
        print("**"*10)
        print(model)
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