# %%
import xarray as xr
import numpy as np
import os,sys
import pandas as pd
from baobao.caculate import caculate_q_rh_thetav
from nmc_met_io.read_micaps import read_micaps_5


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes  # 这个还没用过
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import xarray as xr

# %%

def read_1sta_1time(flnm, id=57083):
    """单个时次单个站点的DataArry

    Args:
        flnm ([type]): [description]
    """
    df = read_micaps_5(flnm)
    ddf = df.loc[df.ID==str(id)]
    ## 这里是利用read_micaps_5, 只能是这些名称
    df1 = ddf.loc[:,['pressure', 'height', 'temperature', 'dewpoint', 'wind_angle', 'wind_speed']]
    df1 = df1.astype('float')

    df1['u'] = -1*np.sin(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    df1['v'] = -1*np.cos(ddf['wind_angle']/180*np.pi)*ddf['wind_speed']
    df1['height'] = df1['height']*10
    df1.rename(columns={'temperature':'temp', 'dewpoint':'td'}, inplace=True)


    df1.columns.name = 'vars'  # 还必须是这个名称, vars不能变
    df2 = df1.drop_duplicates('pressure', keep='last')
    df2 = df2.set_index('pressure')
    da = xr.DataArray(df2)
    da = da.assign_coords({'time':df.time[0]})
    da = da.assign_coords({'id':id})
    ## 对缺省的数据，进行插值
    dda = (da[::-1].interpolate_na(dim='pressure', method='nearest', fill_value="extrapolate"))[::-1]

    ds = dda.to_dataset(dim='vars')
    return ds

def draw_skewt(ds):
    p = ds['pressure'].values * units.hPa
    T = ds['temp'].values* units.degC
    Td = ds['td'].values* units.degC
    u = ds['u'].values* units('m/s')*2.5
    v = ds['v'].values* units('m/s')*2.5


    # 计算抬升凝结高度上的温度和气压
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    lcl_pressure

    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degK')  # 状态曲线
    cape, cin = mpcalc.cape_cin(p,T, Td, parcel_prof)

    ## 保存抬升凝结高度和cape值
    dic_diag = {
        'lcl_pressure':lcl_pressure.magnitude,
        'lcl_temperature':lcl_temperature.magnitude,
        'cape':cape.magnitude,
        'cin':cin.magnitude,
    }
    da = xr.DataArray(
        list(dic_diag.values()),
        coords={'var':list(dic_diag.keys())},
        dims=['var']
    )


    ## 画探空曲线
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm,7*cm), dpi=300)
    skew = SkewT(fig, rotation=30, rect=[0.2, 0.15, 0.70, 0.75])
    skew.plot(p, T, 'r', linewidth=1.5)  # 画温度层节曲线
    skew.plot(p, Td, 'g', linewidth=1.5)  # 露点层节曲线
    # skew.plot_barbs(p[::5], u[::5], v[::5], length = 7)  # 画风
    skew.plot_barbs(p, u, v, length = 5)  # 画风


    ## 设置范围
    skew.ax.set_ylim(1000,200)
    # skew.ax.set_ylim(990,200)
    skew.ax.set_xlim(-30,30)

    # Plot LCL temperature as black dot, 画LCL
    skew.plot(lcl_pressure, lcl_temperature, 'o', markerfacecolor='black',linewidth=20)


    # Plot the parcel profile as a black line, 画状态曲线
    skew.plot(p, parcel_prof, 'k', linewidth=1)

    # Shade areas of CAPE and CIN, 画cape和cin的填色
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)

    # skew.ax.text(30, 300, 'AAAAA', zorder=10)
    skew.ax.set_title('Cape: {}'.format(cape.magnitude.round(1)), loc='right', y=0.8, x=0.85)

    # Plot a zero degree isotherm
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=1)

    # Add the relevant special lines, 画感绝热和湿绝热线
    skew.plot_dry_adiabats(alpha=0.1, colors='black')
    skew.plot_moist_adiabats(alpha=0.1, colors='black')
    skew.plot_mixing_lines(alpha=0.1)

    ## 画风玫瑰子图, 有点看不大懂
    # ax_hod = inset_axes(skew.ax, '30%', '%30', loc=1)
    # ax_hod = fig.add_axes([0.6, 0.65, 0.2, 0.2])
    # h = Hodograph(ax_hod, component_range=80.)
    # h.add_grid(increment=10)
    # h.plot(u,v)
    # h.plot_colormapped(u,v,)


    skew.ax.tick_params(axis='both', labelsize=10, direction='out')
    # skew.ax.set_title(pic_dic['time'].strftime('%Y-%m-%d %H'), loc='right', size=20)
    # skew.ax.set_title(pic_dic['type'], loc='left', size=20)
    skew.ax.set_ylabel('Pressure (hPa)', size=10)
    skew.ax.set_xlabel('Temperature (℃)', size=10)
    skew.ax.set_yticks([990, 925,  850, 800, 700, 600, 500, 400, 300, 200])

    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_sounding/'
    # fig_name = pic_dic['type']+'_'+pic_dic['time'].strftime('%Y%m%d_%H')
    fig_name = 'aaaaaaaaaaaaaaa'
    fig.savefig(path+fig_name)
#     return da
flnm = '/mnt/zfm_18T/fengxiang/aa/2014_07_11-12/sounding/data/micaps/14071120.000'
ds = read_1sta_1time(flnm, id='58203')
draw_skewt(ds)