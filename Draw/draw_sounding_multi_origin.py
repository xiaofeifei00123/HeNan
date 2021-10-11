#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
目的：
    绘制第三次科考Skewt探空曲线图，利用metpy库，
    将多个模式的数据画在一起进行对比
    需要解决的一个问题是怎么把不同时次的有数据的模式挑出来
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
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou_all.nc'
ds = xr.open_dataset(flnm)
t = pd.Timestamp('2021-07-19 00')
dds = ds.sel(time=t)
# %%
dds.dropna(dim='model', how='all')
# dds.sel(model='ERA51900')

# %%
def get_data(da):
    """获得一个试验所需要的资料
    将数据转换成画图的需要单位

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou_all.nc'
    # ds = xr.open_dataset(flnm)
    # da = ds.sel(model=model_dic['type'], time=model_dic['time'])
    da = da.dropna(dim='pressure')

    # pressure_levels = [980, 900, 800, 700, 600, 500, 400, 300, 200]
    pressure_levels = [980, 900, 800, 700, 600, 500]
    # pressure_levels = [980, 900, 800, 700, 600, 500, 400, 300]
    da = da.sel(pressure=pressure_levels)
    ## 反转一下纵坐标， 地面数据在最前面
    da = da.reindex(pressure=da.pressure[::-1])
    p = da.pressure.values * units.hPa
    T = da.temp.values * units.degC
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

def draw_skewt(skew, dic, pic_dic):
    """画探空曲线

    Args:
        skew ([type]): [description]
        dic ([type]): [description]
        pic_dic ([type]): [description]
    """
    print('draw %s时 %s类型的图'%(pic_dic['time'].strftime('%Y%M-%d %H%M'), pic_dic['model']))
    # p = dic['p'][::10]
    # T = dic['T'][::10]
    # Td = dic['Td'][::10]
    
    p = dic['p']
    T = dic['T']
    Td = dic['Td']
    # 计算抬升凝结高度上的温度和气压
    if pic_dic['model'] == 'micaps':
        pic_dic['model'] = 'Micaps' 
    skew.plot(p, T,  linewidth=1.5, color = pic_dic['color'],linestyle='solid', marker=pic_dic['marker'], label=pic_dic['model'], markersize=5)  # 画温度层节曲线
    skew.plot(p, Td,  linewidth=1.5, color = pic_dic['color'],linestyle='dashed', marker=pic_dic['marker'], markersize=5)  # 画温度层节曲线
    

def get_pic_dic_list(ds):
    """控制数据和图形的字典的列表
    把所有的循环放在这个里面来做
    设置不同模式的线条颜色等
    parameters:
        ds: 包含所有数据的Dataset
        t: 需要那个时次的

    Returns:
        [type]: [description]
    """
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou_all.nc'
    # ds = xr.open_dataset(flnm)
    # t = pd.Timestamp('2021-07-20 00')
    model_list = list(ds.model.values)
    color = None
    marker = None
    pic_dic_list = []
    for model in model_list:
        if model == 'micaps':
            color = 'black'
            marker = 'o'
        else:
            f = model[0:4]
            t = model[4:]
            if t == '1800':
                color = 'red'
            elif t == '1812':
                color = 'green'
            elif t == '1900':
                color = 'blue'
            elif t == '1912':
                color = 'm'

            if f == 'ERA5':
                marker = 's'
            elif f == 'GDAS':
                marker = '^'

        pic_dic = {
            'model':model,
            'color':color,
            'marker':marker,
        }
        # pic_dic.update({'Micaps':pic_dic.pop('micaps')})
        # print(pic_dic)
        pic_dic_list.append(pic_dic)
    return pic_dic_list



def main(t):
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou_all.nc'
    ds = xr.open_dataset(flnm)
    ## 对某一个时次来看
    # t = pd.Timestamp('2021-07-20 06')
    fig = plt.figure(figsize=(10,10), dpi=300)
    skew = SkewT(fig, rotation=45, rect=[0.05, 0.08, 0.94, 0.88], aspect=70)
    # skew = SkewT(fig, rotation=45, rect=[0.1, 0.1, 0.85, 0.6],aspect=1)
    # skew.ax.set_extent([10,30, 200, 1000])
    ## 筛选时间，将不合要求的时间删除
    ds = ds.sel(time=t)
    ds = ds.dropna(dim='model', how='all')
    pic_dic_list = get_pic_dic_list(ds)
    skew.ax.set_ylim(1000,500)
    skew.ax.set_xlim(15,30)
    for pic_dic in pic_dic_list:
        ## 获取变量
        # da = ds.sel(model=pic_dic['model'], time=t)
        ## 筛选模式
        da = ds.sel(model=pic_dic['model'])
        dic = get_data(da)
        pic_dic.update({'time':t})
        # print(dic)
        draw_skewt(skew, dic, pic_dic)
    ## 设置范围
    skew.ax.legend(edgecolor='white', fontsize=20,loc='lower left')
    # skew.ax.legend(edgecolor='white')
    # Add the relevant special lines, 画感绝热和湿绝热线
    skew.plot_dry_adiabats(alpha=0.1, colors='black')
    skew.plot_moist_adiabats(alpha=0.1, colors='black')
    skew.plot_mixing_lines(alpha=0.1)

    skew.ax.tick_params(axis='both', labelsize=20, direction='out')
    skew.ax.set_title(pic_dic['time'].strftime('%Y-%m-%d %H'), loc='right', size=20)
    # skew.ax.set_title(pic_dic['model'], loc='left', size=20)
    skew.ax.set_ylabel('Pressure (hPa)', size=24)
    skew.ax.set_xlabel('Temperature (℃)', size=24)

    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_sounding/'
    fig_name = 'diff1'+'_'+pic_dic['time'].strftime('%Y%m%d_%H')
    fig.savefig(path+fig_name)
# main()


# %%
if __name__ == '__main__':
    pass
    t1 = pd.date_range('2021-07-18 00', '2021-07-21 12', freq='6H')
    t2 = pd.date_range('2021-07-18 18', '2021-07-21 12', freq='24H')
    ttt = t1.drop(t2)
    ttt
    for t in ttt:
        main(t)
    #### 