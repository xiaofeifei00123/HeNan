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

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
ds = xr.open_dataset(flnm)
ds.time
ds.sel(time='2021-07-20 1200')

# %%
def get_data():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)
    da = ds.sel(time='2021-07-20 12')
    # da = ds.sel(time='2021-07-18 00')
    # da.temp
    # ds.time
    p = da.pressure.values * units.hPa
    # T = da.temp.values * units.K
    T = da.temp.values * units.degC
    Td = da.td.values * units.degC
    u = da.U.values * units('m/s')
    v = da.U.values * units('m/s')
    dic = {
        'p':p,
        'T':T,
        'Td':Td,
        'u':u,
        'v':v
    }
    return dic


# %%
if __name__ == '__main__':

    ## 读数据
    # flnm = './Data/gaize20140819.19.tsv'
    # col_names = ['T', 'v', 'u', 'P', 'TD']
    # df = pd.read_fwf(flnm,
    #                  skiprows=39,
    #                  usecols=[2, 4, 5, 7, 8],
    #                  names=col_names)

    # # 将缺省值设置-32678.00设置为None
    # df[df<=-32678] = None
    # # 将含有缺省值的行略去，how=any即只要该行含有Nan,略去整行
    # df = df.dropna(axis=0, subset=('T', 'v', 'u', 'P', 'TD'),
    #                how='any').reset_index(drop=True)

    # df = df[0:1600]
                
    
    # p = df['P'].values * units.hPa
    # T = df['T'].values * units.K
    # Td = df['TD'].values * units.K
    # u = df['u'].values * units('m/s')
    # v = df['v'].values * units('m/s')
    # print(type(p))
    ds = get_data()
    p = ds['p']
    T = ds['T']
    Td = ds['Td']
    u = ds['u']
    v = ds['v']
    # 计算抬升凝结高度上的温度和气压
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    # print(lcl_pressure, lcl_temperature)
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degK')

    fig = plt.figure(figsize=(10,10), dpi=300)
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])

    skew = SkewT(fig, rotation=30)
    skew.plot(p, T, 'r', linewidth=1.5)  # 画温度层节曲线
    skew.plot(p, Td, 'g', linewidth=1.5)  # 露点层节曲线
    # skew.plot_barbs(p[::80], u[::80], v[::80], length = 5)  # 风
    # skew.plot_barbs(p[::5], u[::5], v[::5], length = 5)  # 风
    skew.plot_barbs(p[::4], u[::4], v[::4], length = 5)  # 风

    skew.ax.set_ylim(1000,100)
    skew.ax.set_xlim(-50,40)

    # Plot LCL temperature as black dot
    skew.plot(lcl_pressure, lcl_temperature, 'o', markerfacecolor='black',linewidth=20)

    # Plot the parcel profile as a black line
    skew.plot(p, parcel_prof, 'k', linewidth=1)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)

    # Plot a zero degree isotherm
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=1)

    # Add the relevant special lines
    skew.plot_dry_adiabats(alpha=0.1, colors='black')
    skew.plot_moist_adiabats(alpha=0.1, colors='black')
    skew.plot_mixing_lines(alpha=0.1)
    plt.title("2021-07-20 12", fontsize=20)

    skew.ax.tick_params(axis='both', labelsize=20, direction='out')
    skew.ax.set_xlabel('')
    
    

    fig.savefig("test.png")