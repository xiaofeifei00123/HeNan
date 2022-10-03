# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from draw_upar_850_d03_div import draw
import metpy
from metpy.units import units
from baobao.caculate import caculate_q_rh_thetaev

# %%
def get_data():
    def cal_ri(ds1):
        pre = ds1.pressure*units.hPa
        temp = ds1['temp']*units.degC
        # theta = metpy.calc.potential_temperature(pre, temp)
        ds2 = ds1.reset_coords()
        ds3 = caculate_q_rh_thetaev(ds1)
        theta = ds3['theta_v']*units.K
        height = ds1['geopt']*units.m
        u = ds1['u']
        v = ds1['v']
        u = u*units('m/s')
        v = v*units('m/s')
        ri = metpy.calc.gradient_richardson_number(height, theta, u, v, vertical_dim=0)
        # bv = metpy.calc.brunt_vaisala_frequency(height, theta)
        return ri

        

    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)

    ri_list = []
    for t in ds.time:
        dss = ds.sel(time=t)
        ri = cal_ri(dss)
        ri_list.append(ri)
    rii = xr.concat(ri_list, dim='time')
    ri2 = rii.sel(time=slice('2021-07-19 00', '2021-07-20 20'))
    return ri2

# %%

def draw(ri2):
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    line_list = ['--', '--', '--', '-', '-', '-']
    color_list = ['red', 'blue', 'green', 'red', 'blue', 'green']
    i = 0
    for t in ri2.time:
        ri1 = ri2.sel(time=t)
        x = ri1.values
        y = ri1.pressure.values
        label = (t+pd.Timedelta('8H')).dt.strftime('%d/%H').values
        ax.plot(x,y, label=label, linestyle=line_list[i], color=color_list[i])
        i += 1
    ax.set_xlim(-0.5, 1)
    # ax.set_xlim(-2, -0.5)
    ax.invert_yaxis()
    # ax.vlines(1)
    ax.axvline(x=0.25, color='black')
    ax.set_ylim(1000, 850)
    ax.set_yticks(np.arange(1000, 849, -50))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.legend(edgecolor='white')
    ax.set_ylabel('Pressure (hPa)')
    ax.set_xlabel('Richardson number')
    path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_sounding/micaps/'
    fig_name = path+'richard'
    fig.savefig(fig_name)

if __name__ == '__main__':
    
    ri2 = get_data()
    draw(ri2)