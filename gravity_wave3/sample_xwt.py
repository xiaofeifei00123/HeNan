"""
Python sample script for cross-wavelet analysis and the statistical approach
suggested by Torrence and Compo (1998) using the wavelet module. To run
this script successfully, the `matplotlib` and `progressbar` modules have to
be installed.


Disclaimer
----------
This module is based on routines provided by C. Torrence and G. P. Compo
available at <http://paos.colorado.edu/research/wavelets/>, on routines
provided by A. Grinsted, J. Moore and S. Jevrejeva available at
<http://noc.ac.uk/using-science/crosswavelet-wavelet-coherence>, and
on routines provided by A. Brazhe available at
<http://cell.biophys.msu.ru/static/swan/>.

This software is released under a BSD-style open source license. Please read
the license file for furter information. This routine is provided as is
without any express or implied warranties whatsoever.

Authors
-------
Nabil Freij, Sebastian Krieger

"""
# %%
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
import pycwt as wavelet

from pycwt.helpers import find
from matplotlib.image import NonUniformImage
from baobao import caculate as ca
import xarray as xr
from geopy.distance import distance  # 根据经纬度计算两点距离
import cmaps
import matplotlib.ticker as ticker
# %%

def get_data_div_vor():
    def drop_na(da):
        """处理数据, 这一步是必须要的，不然好像画不出来图
        """
        for i in range(da.shape[-1]):
            column_vals = da[:,i].values
            # Let's find the lowest index that isn't filled. The nonzero function
            # finds all unmasked values greater than 0. Since 0 is a valid value
            # for dBZ, let's change that threshold to be -200 dBZ instead.
            first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
            da[0:first_idx, i] = da[first_idx, i]
        da = da.dropna(dim='vertical')
        return da

        
    def latlon2distance(da2):
        """将剖面数据的经纬度横坐标变为距离坐标

        Args:
            da2 (_type_): _description_

        Returns:
            _type_: _description_
        """
        # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
        # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross1.nc'
        # ds = xr.open_dataset(flnm)
        # ds1 = ds.sel(time='2021-07-20 08')
        # da = ds1['wa_cross']
        # da1 = da.interpolate_na(dim='vertical', method='linear',  fill_value="extrapolate")
        # da2 = da1.sel(vertical=2000, method='nearest')
        dd = da2.xy_loc
        def str_latlon(string):
            # d1 = dd.values[0]
            lat = float(string.split(',')[0])
            lon = float(string.split(',')[1])
            return lat, lon

        d2 = dd.values
        lat_list = []
        lon_list = []
        for i in d2:
            # print(i)
            lat, lon = str_latlon(i)
            lat_list.append(lat)
            lon_list.append(lon)

        dis_list = [0]
        di = 0
        for i in range(len(lat_list)-1):
            # print(i)
            lat1 = lat_list[i]
            lon1 = lon_list[i]
            loc1 = (lat1, lon1)
            lat2 = lat_list[i+1]
            lon2 = lon_list[i+1]
            loc2 = (lat2, lon2)
            dist = distance(loc1,loc2).km
            di = di+dist
            dis_list.append(di)
        dis_list
        dis_array = (np.array(dis_list)).round(1)
        dis_array
        da2 = da2.assign_coords({'distance':('cross_line_idx',dis_array)})
        da3 = da2.swap_dims({'cross_line_idx':'distance'})
        return da3
        
    flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/cross4_1time.nc'
    # flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/cross4_1time.nc'
    # flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/SS/wrfout/cross4_1time.nc'
    # flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/FD/wrfout/cross4_1time.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.squeeze()


    def select_vertical(data):
        # data = data.interp(vertical=2000, method='slinear')
        data = data.interp(vertical=5200, method='slinear')
        data = latlon2distance(data)
        data = data*10**5
        return data

    da = drop_na(ds['div_cross'])
    db = drop_na(ds['vor_cross'])
    div = select_vertical(da)
    vor = select_vertical(db)
    # print(div.max(), vor.max())
    # return div, vor
    return vor, div

def get_data_cross_wrf(var='div'):
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_C.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
    hh = ds['height'].mean(dim='time').values
    ds2 = ds.assign_coords({'z':('pressure',hh)})
    ds3 = ds2.swap_dims({'pressure':'z'})
    da = ds3[var]#.sel(z=1000, method='nearest')
    db = da.sel(z = np.sort(da.z))
    # dc = db.sel(z=1500, method='nearest')
    dc = db.sel(z=1000, method='nearest')
    if var == 'div':
        dc = dc*10**5
    elif var == 'wa':
        dc = dc*10**2
    elif var == 'vor':
        dc = dc*10**5
    else:
        print('出现错误')
    da = dc
    time = da.time.values
    sst = da.values
    return sst, time

# d1, t1 = get_data_cross_wrf('vor')
# d2, t2 = get_data_cross_wrf('wa')
# t1 = np.arange(len(t1))
# t2 = np.arange(len(t2))
# t1, s1 = t1, d1
# t2, s2 = t2, d2

# da1, da2 = get_data_div_vor()
da1, da2 = get_data_div_vor()
# print(da1.mean().values, da2.mean().values)
# da1.shape

# cm = 1/2.54
# fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# x = np.arange(da1.shape[0])
# y = da1.values
# y2 = da2.values
# ax.plot(x,y)
# ax.plot(x,y2)

# np.where()
# da1.where(da1==0., da1, 1)
# xr.where(da1==0., 0, 2).plot()




t1 = da1.distance.values
t2 = da2.distance.values
s1 = da1.values
s2 = da2.values
mother = 'morlet'
# s1 = s1+50

dt = np.diff(t1)[0]
n1 = t1.size
n2 = t2.size
n = min(n1, n2)

# Change the probablity density function (PDF) of the data. The time series
# of Baltic Sea ice extent is highly bi-modal and we therefore transform the
# timeseries into a series of percentiles. The transformed series probably
# reacts 'more linearly' to climate.
s2, _, _ = wavelet.helpers.boxpdf(s2)

# Calculates the standard deviation of each time series for later
# normalization.
std1 = s1.std()
std2 = s2.std()


# I. Continuous wavelet transform
# ===============================

# Calculate the CWT of both normalized time series. The function wavelet.cwt
# returns a a list with containing [wave, scales, freqs, coi, fft, fftfreqs]
# variables.
# mother = wavelet.Morlet(6)          # Morlet mother wavelet with m=6
mother = wavelet.Morlet(6)          # Morlet mother wavelet with m=6
slevel = 0.95                       # Significance level
dj = 1/12                           # Twelve sub-octaves per octaves
s0 = -1  # 2 * dt                  # Starting scale, here 6 months
J = -1  # 7 / dj                    # Seven powers of two with dj sub-octaves
if True:
    alpha1, _, _ = wavelet.ar1(s1)  # Lag-1 autocorrelation for red noise
    alpha2, _, _ = wavelet.ar1(s2)  # Lag-1 autocorrelation for red noise
else:
    alpha1 = alpha2 = 0.0           # Lag-1 autocorrelation for white noise
# The following routines perform the wavelet transform and siginificance
# analysis for two data sets.
W1, scales1, freqs1, coi1, _, _ = wavelet.cwt(s1/std1, dt, dj, s0, J, mother)

signif1, fft_theor1 = wavelet.significance(1.0, dt, scales1, 0, alpha1,
                                           significance_level=slevel,
                                           wavelet=mother)
W2, scales2, freqs2, coi2, _, _ = wavelet.cwt(s2/std2, dt, dj, s0, J, mother)
signif2, fft_theor2 = wavelet.significance(1.0, dt, scales2, 0, alpha2,
                                           significance_level=slevel,
                                           wavelet=mother)

power1 = (np.abs(W1)) ** 2             # Normalized wavelet power spectrum
power2 = (np.abs(W2)) ** 2             # Normalized wavelet power spectrum
period1 = 1/freqs1
period2 = 1/freqs2
sig95_1 = np.ones([1, n1]) * signif1[:, None]
sig95_1 = power1 / sig95_1             # Where ratio > 1, power is significant
sig95_2 = np.ones([1, n2]) * signif2[:, None]
sig95_2 = power2 / sig95_2             # Where ratio > 1, power is significant

s2 = s2[np.argwhere((t2 >= min(t1)) & (t2 <= max(t1))).flatten()]

# Calculate the cross wavelet transform (XWT). The XWT finds regions in time
# frequency space where the time series show high common power. Torrence and
# Compo (1998) state that the percent point function -- PPF (inverse of the
# cumulative distribution function) -- of a chi-square distribution at 95%
# confidence and two degrees of freedom is Z2(95%)=3.999. However, calculating
# the PPF using chi2.ppf gives Z2(95%)=5.991. To ensure similar significance
# intervals as in Grinsted et al. (2004), one has to use confidence of 86.46%.
W12, cross_coi, freq, signif = wavelet.xwt(s1, s2, dt, dj=1/12, s0=-1, J=-1,
                                           significance_level=0.8646,
                                           wavelet='morlet', normalize=True)

cross_power = np.abs(W12)**2
cross_sig = np.ones([1, n]) * signif[:, None]
cross_sig = cross_power / cross_sig  # Power is significant where ratio > 1
cross_period = 1/freq

# s1 = s1**2
WCT, aWCT, corr_coi, freq, sig = wavelet.wct(s1, s2, dt, dj=1/12, s0=-1, J=-1,
                                             significance_level=0.8646,
                                             wavelet='morlet', normalize=True,
                                             cache=True)
# print(WCT.mean())
# WCT, aWCT, corr_coi, freq, sig = wavelet.wct(s1, s2, dt, dj=1/12, s0=-1, J=-1,
#                                              significance_level=0.8646,
#                                              wavelet='morlet', normalize=True,
#                                              cache=False)

cor_sig = np.ones([1, n]) * sig[:, None]
cor_sig = np.abs(WCT) / cor_sig  # Power is significant where ratio > 1
cor_period = 1 / freq

# left (W).
angle = 0.5 * np.pi - aWCT
u, v = np.cos(angle), np.sin(angle)

cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# cr = ax2.contourf(t1,cor_period, WCT, cmap=cmaps.NCV_blue_red)
cr = ax2.contourf(t1,cor_period, WCT, cmap=cmaps.WhiteBlueGreenYellowRed)

fig.colorbar(cr, orientation='horizontal')



ax2.fill(np.concatenate([t1, t1[-1:]+dt, t1[-1:]+dt, t1[:1]-dt, t1[:1]-dt]),
         np.concatenate([corr_coi, [1e-9], cor_period[-1:], cor_period[-1:],
                         [1e-9]]),
         'k', alpha=0.3, hatch='x')
ax2.set_title('Cross-Correlation')
# ax2.quiver(t1[::3], cor_period[::3], u[::3, ::3], v[::3, ::3], units='height',
#            angles='uv', pivot='mid', linewidth=1, edgecolor='k',
#            headwidth=10, headlength=10, headaxislength=5, minshaft=2,
#            minlength=5)
ax2.quiver(t1[::5], cor_period[::5], u[::5, ::5], v[::5, ::5], units='height',
           angles='uv', pivot='mid', linewidth=0.8, edgecolor='k',
           headwidth=5, headlength=5, headaxislength=4, minshaft=1,
           minlength=5, width=0.002)
ax2.set_ylim(6, 256)
ax2.set_yscale('log', base=2, subs=None)
ax2.set_xlim(max(t1.min(), t2.min()), min(t1.max(), t2.max()))
ax2.set_yticks([8, 16, 32, 64, 128])
axx = plt.gca().yaxis
axx.set_major_formatter(ticker.ScalarFormatter())
ax2.ticklabel_format(axis='y', style='plain')
fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_wave'
fig.savefig(fig_path+'/xwt_GWD3.png')
# %%

# flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/cross4_1time.nc'
# ds = xr.open_dataset(flnm)