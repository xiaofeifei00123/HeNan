# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from draw_upar_850_d03_div import draw
import metpy
import metpy.calc as ca
from metpy.units import units
from baobao.caculate import caculate_q_rh_thetaev

# %%

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
ds = xr.open_dataset(flnm)
ds
# %%
# (ds['geopt']/9.86).max()
ds['geopt']
# %%
# ca.
# ds.sel(pressure=200, method='nearest')['geopt']
geo = ds['geopt'].isel(time=0)
u = ds['u'].isel(time=0)
# %%
# u.shape
delta = geo[1:].values-geo[0:-1].values
# delta.shape
# geo.shape
ca.gradient(f=u)
# %%
# u
# ds.reset_coords()
ht = ds['geopt'].mean(dim='time').values
ds2 = ds.assign_coords({'height':('pressure',ht)})
ds3 =ds2.swap_dims({'pressure':'height'})
du = ca.gradient(ds3['u'].isel(time=0))[0]
dv = ca.gradient(ds3['v'].isel(time=0))[0]
# db
# %%
# db[0].plot()
# db[0]
np.sqrt(du**2+dv**2)
# dv**2
# db.plot()
# u.values.shape
# delta.shape

# %%
def get_data():
    """一个时次的
    """

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
        bv = ca.brunt_vaisala_frequency_squared(height, theta)

        
        ds3 =ds1.swap_dims({'pressure':'height'})
        du = ca.gradient(ds3['u'])[0]
        dv = ca.gradient(ds3['v'])[0]
        # duvz = np.sqrt(du**2+dv**2)
        duvz = du**2+dv**2
        duvz.time.values = u.time.values
        
        return ri, bv, duvz


    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
    ds = xr.open_dataset(flnm)
    ht = ds['geopt'].mean(dim='time').values
    ds = ds.assign_coords({'height':('pressure',ht)})



    ri_list = []
    bv_list = []
    duvz_list = []
    for t in ds.time:
        dss = ds.sel(time=t)
        ri, bv, duvz = cal_ri(dss)
        ri_list.append(ri)
        bv_list.append(bv)
        duvz_list.append(duvz)
    rii = xr.concat(ri_list, dim='time')
    bvv = xr.concat(bv_list, dim='time')
    duvzz = xr.concat(duvz_list, dim='time')

    ri2 = rii.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    bv2 = bvv.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    duvz2 = duvzz.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    duvz2 = duvz2.swap_dims({'height':'pressure'})
    return ri2, bv2, duvz2

ri, bv, duvz = get_data()
# %%
# duvz.plot()


# ri.shape
# ri.sel(press)
# duvz
r = ri.sel(pressure=slice(1000, 950)).mean(dim='pressure')
b = bv.sel(pressure=slice(1000, 950)).mean(dim='pressure')
du = duvz.sel(pressure=slice(1000, 950)).mean(dim='pressure')

# duvz.plot()
# duvz.plot()
cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 5*cm), dpi=300)
ax = fig.add_axes([0.1, 0.27, 0.85, 0.7])
# ax.plot()
# x = r.time.dt.strftime('%d\n%H')
x = r.time.dt.strftime('%H\n%d')
y = r.values
y1 = b.values
y2 = du.values
y3 = b/du
# y4 = b/ri
ax.plot(x,y,label='$R_i$', color='black', marker='.')
ax.plot(x,y1*10000, label='$N^2*10^4$', color='red', marker='.')
# ax.plot(x,y2*1000, label='$shear^2$', color='blue', marker='.')
# ax.plot(x,y2*100, label='$shear^2$')
# ax.plot(x,y1/y2, label='$ N^2/shear^2$', marker='.', color='black')
# ax.plot(x,y4, label='$ N^2/shear^2$')
# ax.set_xticks(x[2::2])
# ax.set_xticks(x[4::1])
ax.set_ylim(-2, 5)
ax.axhline(y=0.25, color='black')
ax.legend(edgecolor='white')
ax.set_xlabel('Time (Hour/Date)')
# ax.set_ylabel('Time (Hour/Date)')
figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture/'
fig.savefig(figpath+'richard')



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

# if __name__ == '__main__':
    
#     ri2 = get_data()
#     draw(ri2)