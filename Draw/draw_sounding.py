#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
目的：
    绘制第三次科考Skewt探空曲线图，利用metpy库，
参考: 
    https://zhuanlan.zhihu.com/p/102032513
需要注意:
    这个是画单个探空站用的，所以画的要素比较全,

-----------------------------------------
Time             :2021/03/22 14:03:00
Author           :Forxd
Version          :1.0
'''

# %%
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # 这个还没用过
import numpy as np
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units
import pandas as pd
import xarray as xr
# %%

# %%

# %%
def get_data_micaps(ds, pic_dic):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    # ds = xr.open_dataset(flnm)
    # da = ds.sel(time='2021-07-20 00')
    t = pic_dic['time']
    da = ds.sel(time=t)
    p = da.pressure.values * units.hPa
    T = da.temp.values * units.degC
    Td = da.td.values * units.degC
    u = da.U.values * units('m/s')*2.5
    v = da.V.values * units('m/s')*2.5
    dic = {
        'p':p,
        'T':T,
        'Td':Td,
        'u':u,
        'v':v
    }
    return dic

def get_data_wrf(ds, pic_dic):

    ds = ds.rename({'ua':'u', 'va':'v',})
    t = pic_dic['time']
    # da = ds.sel(time=t)
    da = ds.sel(time=t, lat=34.71, lon=113.66, method='nearest')
    # da.dropna(dim='td', how='all')
    da = da.dropna(dim='pressure')
    p = da.pressure.values * units.hPa
    T = da.temp.values * units.K
    Td = da.td.values * units.degC
    u = da.u.values * units('m/s')*2.5
    v = da.v.values * units('m/s')*2.5
    dic = {
        'p':p,
        'T':T,
        'Td':Td,
        'u':u,
        'v':v
    }
    return dic

def draw_skewt(ds, pic_dic):
    # ds = get_data()
    print('draw %s时 %s类型的图'%(pic_dic['time'].strftime('%Y%M-%d %H%M'), pic_dic['type']))
    p = ds['p']
    T = ds['T']
    Td = ds['Td']
    u = ds['u']
    v = ds['v']
    
    # 计算抬升凝结高度上的温度和气压
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    # lcl_pressure, lcl_temperature = mpcalc.lcl(p, T, Td)
    # print(lcl_pressure, lcl_temperature)
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degK')

    fig = plt.figure(figsize=(10,10), dpi=300)
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])

    skew = SkewT(fig, rotation=30)
    skew.plot(p, T, 'r', linewidth=1.5)  # 画温度层节曲线
    skew.plot(p, Td, 'g', linewidth=1.5)  # 露点层节曲线
    # skew.plot_barbs(p[::80], u[::80], v[::80], length = 5)  # 风
    # skew.plot_barbs(p[::5], u[::5], v[::5], length = 5)  # 风
    # skew.plot_barbs(p[::5], u[::5], v[::5], length = 7)  # 风
    skew.plot_barbs(p[::5], u[::5], v[::5], length = 7)  # 风

    skew.ax.set_ylim(1000,100)
    skew.ax.set_xlim(-40,30)

    # Plot LCL temperature as black dot
    skew.plot(lcl_pressure, lcl_temperature, 'o', markerfacecolor='black',linewidth=20)

    # Plot the parcel profile as a black line
    skew.plot(p, parcel_prof, 'k', linewidth=1)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)

    # Plot a zero degree isotherm
    # skew.ax.axvline(0, color='c', linestyle='--', linewidth=1)

    # Add the relevant special lines
    skew.plot_dry_adiabats(alpha=0.1, colors='black')
    skew.plot_moist_adiabats(alpha=0.1, colors='black')
    skew.plot_mixing_lines(alpha=0.1)
    # plt.title("2021-07-20 12", fontsize=20)

    skew.ax.tick_params(axis='both', labelsize=20, direction='out')
    skew.ax.set_title(pic_dic['time'].strftime('%Y-%m-%d %H'), loc='right', size=20)
    skew.ax.set_title(pic_dic['type'], loc='left', size=20)
    skew.ax.set_ylabel('Pressure (hPa)', size=24)
    skew.ax.set_xlabel('Temperature (℃)', size=24)

    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_sounding/'
    fig_name = pic_dic['type']+'_'+pic_dic['time'].strftime('%Y%m%d_%H')
    fig.savefig(path+fig_name)


def draw_micaps():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    # t = pd.Timestamp('2021-07-20 00')
    tt = pd.date_range('2021-07-20 00', '2021-07-20 14', freq='6H')
    for t in tt:
        pic_dic = {
            'time':t,
            'type':'micaps'
        }
        dic = get_data_micaps(ds, pic_dic)
        # # print(da)
        draw_skewt(dic, pic_dic)

def draw_wrf(flnm, dic_model):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d03_latlon.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    # t = pd.Timestamp('2021-07-20 00')
    tt = pd.date_range('2021-07-20 00', '2021-07-20 14', freq='6H')
    
    for t in tt:
        pic_dic = {
            'time':t,
            # 'type':'ERA5_YSU_1800'
            'type':dic_model['file_type']+"_"+dic_model['initial_time']
        }
        print(pic_dic)
        dic = get_data_wrf(ds, pic_dic)

        draw_skewt(dic, pic_dic)

def draw_wrf_all():
    pass
    time_list = ['1800', '1812', '1900', '1912']
    initial_file_list = ['ERA5', 'GDAS']
    for f in initial_file_list:
        # time_list = ['1800']
        # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/'
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+f+'/'
        # gu = GetUpar()
        for t in time_list:
            # path_wrfout = path_main+'YSU_'+t+'/'
            # ds = gu.get_upar(path_wrfout)
            dic_model = {'initial_time':t, 'file_type':f}
            flnm = 'YSU_'+t
            path_in = path_main+flnm+'_upar_d03_latlon.nc'
            draw_wrf(path_in, dic_model)


# %%
if __name__ == '__main__':

    pass
    draw_wrf_all()
    # draw_micaps()
    # main()