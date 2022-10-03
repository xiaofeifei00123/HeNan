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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import xarray as xr
# %%

# %%

# %%
def get_data_micaps(ds, pic_dic):
    """_summary_

    Args:
        ds (_type_): _description_
        pic_dic (_type_): _description_

    Returns:
        _type_: _description_
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    # ds = xr.open_dataset(flnm)
    # da = ds.sel(time='2021-07-20 00')
    t = pic_dic['time']
    da = ds.sel(time=t)
    p = da.pressure.values * units.hPa
    T = da.temp.values * units.degC
    Td = da.td.values * units.degC
    # u = da.U.values * units('m/s')*2.5
    # v = da.V.values * units('m/s')*2.5
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

def get_data_wrf(ds, pic_dic):
    """从聚合的wrf, nc文件中读取需要的站点数据

    Args:
        ds ([type]): [description]
        pic_dic ([type]): [description]

    Returns:
        [type]: [description]
    """

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
    """计算需要的诊断变量

    Args:
        ds ([type]): [description]
        pic_dic ([type]): [description]

    Returns:
        [type]: [description]
    """
    # ds = get_data()
    print('draw %s时 %s类型的图'%(pic_dic['time'].strftime('%Y%M-%d %H%M'), pic_dic['type']))
    p = ds['p']
    T = ds['T']
    Td = ds['Td']
    u = ds['u']
    v = ds['v']
    
    # 计算抬升凝结高度上的温度和气压
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degK')  # 状态曲线
    cape, cin = mpcalc.cape_cin(p,T, Td, parcel_prof)
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
    # return da

    # fig = plt.figure(figsize=(10,10), dpi=400)
    # skew = SkewT(fig, rotation=30)
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm,7*cm), dpi=300)
    skew = SkewT(fig, rotation=30, rect=[0.2, 0.15, 0.70, 0.75])

    skew.plot(p, T, 'r', linewidth=1.5)  # 画温度层节曲线
    skew.plot(p, Td, 'g', linewidth=1.5)  # 露点层节曲线
    skew.plot_barbs(p[::5], u[::5], v[::5], length = 5)  # 画风
    # skew.plot_barbs(p, u, v, length = 5)  # 画风

    ## 设置范围
    # skew.ax.set_ylim(1000,200)
    skew.ax.set_ylim(990,200)
    skew.ax.set_xlim(-30,30)

    # Plot LCL temperature as black dot, 画LCL
    skew.plot(lcl_pressure, lcl_temperature, 'o', markerfacecolor='black',linewidth=20)

    # Plot the parcel profile as a black line, 画状态曲线
    skew.plot(p, parcel_prof, 'k', linewidth=1)

    # Shade areas of CAPE and CIN, 画cape和cin的填色
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)

    # Plot a zero degree isotherm
    # skew.ax.axvline(0, color='c', linestyle='--', linewidth=1)

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

    tt = pic_dic['time']+pd.Timedelta('8H')
    skew.ax.set_title(tt.strftime('%Y-%m-%d %H'), loc='left', size=10, y=0.97)
    # skew.ax.set_title(pic_dic['time'].strftime('%Y-%m-%d %H'), loc='right', size=10)
    # skew.ax.set_title(pic_dic['type'], loc='left', size=10)
    skew.ax.set_title('CAPE: {} \n CIN:{}'.format(cape.magnitude.round(1), cin.magnitude.round(1)),  y=0.8, x=0.7)
    # skew.ax.set_title('CIN: {}'.format(cin.magnitude.round(1)),  y=0.9, x=0.7)
    skew.ax.set_ylabel('Pressure (hPa)', size=10)
    skew.ax.set_xlabel('Temperature (℃)', size=10)
    skew.ax.set_yticks([990, 925,  850, 800, 700, 600, 500, 400, 300, 200])

    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_sounding/micaps/'
    fig_name = pic_dic['type']+'_'+pic_dic['time'].strftime('%Y%m%d_%H')
    fig.savefig(path+fig_name)
    return da


def draw_micaps():
    pass
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    # t = pd.Timestamp('2021-07-20 00')
    tt = pd.date_range('2021-07-19 00', '2021-07-19 14', freq='6H')
    # tt = pd.date_range('2021-07-20 00', '2021-07-20 14', freq='6H')
    tt_list = []
    for t in tt:
        pic_dic = {
            'time':t,
            'type':'micaps'
        }
        dic = get_data_micaps(ds, pic_dic)
        # # print(da)
        da = draw_skewt(dic, pic_dic)
        tt_list.append(da)
    ds1 = xr.concat(tt_list, pd.Index(tt, name='time'))
    return ds1

def draw_wrf(flnm, dic_model):
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d03_latlon.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    tt = pd.date_range('2021-07-20 00', '2021-07-20 14', freq='6H')
    
    tt_list = []
    for t in tt:
        pic_dic = {
            'time':t,
            # 'type':'ERA5_YSU_1800'
            'type':dic_model['file_type']+"_"+dic_model['initial_time']
        }
        # print(pic_dic)
        dic = get_data_wrf(ds, pic_dic)
        da = draw_skewt(dic, pic_dic)
        tt_list.append(da)
        # ds1[t] = da
    # print(ds1.var)
    ds1 = xr.concat(tt_list, pd.Index(tt, name='time'))
    return ds1
    

def draw_wrf_all():
    pass
    # time_list = ['1800', '1812', '1900', '1912']
    # initial_file_list = ['ERA5', 'GDAS']

    model_list = ['1900_90m', '1900_900m','1912_900m', '1912_90m']
    for model in model_list:
        path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'+model+'/'
        dic_model = {'initial_time':'', 'file_type':model}
        # flnm = 'YSU_'+t
        flnm = 'upar.nc'
        path_in = path_main+flnm
        diag = draw_wrf(path_in, dic_model)
    #     # ds2[f+t] = diag
    #     # dic2[f+t] = diag
    # return ds2

def draw_1km():
    pass
    ds2 = xr.Dataset()
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/'
    fl_list = ['1km_1912_hgt_upar_d04_latlon.nc', '1km_1912_upar_d04_latlon.nc']
    # gu = GetUpar()
    type_list = ['1km_hgt', '1km']
    for fl, type in zip(fl_list, type_list):
        dic_model = {'initial_time':'1912', 'file_type':type}
        path_in = path_main+fl
        diag = draw_wrf(path_in, dic_model)


def test():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800_upar_d03_latlon.nc'
    ds = xr.open_dataset(flnm)  # 所有时次探空的集合
    t = pd.Timestamp('2021-07-20 00')
    pic_dic = {
        'time':t,
        'type':'ERA5_YSU_1800',
    }
    dic = get_data_wrf(ds, pic_dic)
    diag = draw_skewt(dic, pic_dic)


def get_diag():
    ds2 = draw_wrf_all()
    ds3 = ds2.to_array(dim='model').to_dataset(dim='var')
    ds3['lcl_temperature'] = ds3['lcl_temperature']-273.15
    ds4 = ds3.to_array(dim='var').to_dataset(dim='model')
    da_micaps = draw_micaps()
    ds4['micaps'] = da_micaps
    ds4
    # %%
    t = ds4.time
    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/'
    flnm = path+'./diagnostic.xlsx'
    writer = pd.ExcelWriter(flnm)
    for i in t:
        title = str(i.dt.strftime('%Y-%m-%d_%H%M').values)
        print(title)
        ds = ds4.sel(time=i).drop_vars(['time'])
        df = ds.to_dataframe().T
        df.to_excel(writer, sheet_name=title)
    writer.save()
    writer.close()




# %%
if __name__ == '__main__':
    pass
    #### 
    # draw_1km()
    draw_micaps()
    # draw_wrf_all()
    ## 需要决定是否注释117行的return, 不注释了则运行快一些，否则会慢一些
    # get_diag() ## 如果要计算的话，draw里面的return 位置要放在上面