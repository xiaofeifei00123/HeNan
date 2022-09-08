# %%
import netCDF4 as nc
from netCDF4 import Dataset
import wrf
from wrf import getvar
# from baobao.caculate import caculate_div3d
from baobao.caculate import get_div_wrfout
import xarray as xr
import pandas as pd
import numpy as np

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
import pycwt

# %%
"""
In this example we will load the NINO3 sea surface temperature anomaly dataset
between 1871 and 1996. This and other sample data files are kindly provided by
C. Torrence and G. Compo at
<http://paos.colorado.edu/research/wavelets/software.html>.
"""
# We begin by importing the relevant libraries. Please make sure that PyCWT is
# properly installed in your system.
from __future__ import division
import numpy
from matplotlib import pyplot

import pycwt as wavelet
from pycwt.helpers import find

# Then, we load the dataset and define some data related parameters. In this
# case, the first 19 lines of the data file contain meta-data, that we ignore,
# since we set them manually (*i.e.* title, units).
url = 'http://paos.colorado.edu/research/wavelets/wave_idl/nino3sst.txt'
dat = numpy.genfromtxt(url, skip_header=19)
title = 'NINO3 Sea Surface Temperature'
label = 'NINO3 SST'
units = 'degC'
t0 = 1871.0
dt = 0.25  # In years

# We also create a time array in years.
N = dat.size
t = numpy.arange(0, N) * dt + t0

# We write the following code to detrend and normalize the input data by its
# standard deviation. Sometimes detrending is not necessary and simply
# removing the mean value is good enough. However, if your dataset has a well
# defined trend, such as the Mauna Loa CO\ :sub:`2` dataset available in the
# above mentioned website, it is strongly advised to perform detrending.
# Here, we fit a one-degree polynomial function and then subtract it from the
# original data.
p = numpy.polyfit(t - t0, dat, 1)
dat_notrend = dat - numpy.polyval(p, t - t0)
std = dat_notrend.std()  # Standard deviation
var = std ** 2  # Variance
dat_norm = dat_notrend / std  # Normalized dataset

# The next step is to define some parameters of our wavelet analysis. We
# select the mother wavelet, in this case the Morlet wavelet with
# :math:`\omega_0=6`.
mother = wavelet.Morlet(6)
s0 = 2 * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
dj = 1 / 12  # Twelve sub-octaves per octaves
J = 7 / dj  # Seven powers of two with dj sub-octaves
alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise

# The following routines perform the wavelet transform and inverse wavelet
# transform using the parameters defined above. Since we have normalized our
# input time-series, we multiply the inverse transform by the standard
# deviation.
wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J,
                                                      mother)
iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

# We calculate the normalized wavelet and Fourier power spectra, as well as
# the Fourier equivalent periods for each wavelet scale.
power = (numpy.abs(wave)) ** 2
fft_power = numpy.abs(fft) ** 2
period = 1 / freqs

# We could stop at this point and plot our results. However we are also
# interested in the power spectra significance test. The power is significant
# where the ratio ``power / sig95 > 1``.
signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                         significance_level=0.95,
                                         wavelet=mother)
sig95 = numpy.ones([1, N]) * signif[:, None]
sig95 = power / sig95

# Then, we calculate the global wavelet spectrum and determine its
# significance level.
glbl_power = power.mean(axis=1)
dof = N - scales  # Correction for padding at edges
glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha,
                                        significance_level=0.95, dof=dof,
                                        wavelet=mother)

# We also calculate the scale average between 2 years and 8 years, and its
# significance level.
sel = find((period >= 2) & (period < 8))
Cdelta = mother.cdelta
scale_avg = (scales * numpy.ones((N, 1))).transpose()
scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha,
                                             significance_level=0.95,
                                             dof=[scales[sel[0]],
                                                  scales[sel[-1]]],
                                             wavelet=mother)

# Finally, we plot our results in four different subplots containing the
# (i) original series anomaly and the inverse wavelet transform; (ii) the
# wavelet power spectrum (iii) the global wavelet and Fourier spectra ; and
# (iv) the range averaged wavelet spectrum. In all sub-plots the significance
# levels are either included as dotted lines or as filled contour lines.

# Prepare the figure
pyplot.close('all')
pyplot.ioff()
figprops = dict(figsize=(11, 8), dpi=72)
fig = pyplot.figure(**figprops)

# First sub-plot, the original time series anomaly and inverse wavelet
# transform.
ax = pyplot.axes([0.1, 0.75, 0.65, 0.2])
ax.plot(t, iwave, '-', linewidth=1, color=[0.5, 0.5, 0.5])
ax.plot(t, dat, 'k', linewidth=1.5)
ax.set_title('a) {}'.format(title))
ax.set_ylabel(r'{} [{}]'.format(label, units))

# Second sub-plot, the normalized wavelet power spectrum and significance
# level contour lines and cone of influece hatched area. Note that period
# scale is logarithmic.
bx = pyplot.axes([0.1, 0.37, 0.65, 0.28], sharex=ax)
levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
bx.contourf(t, numpy.log2(period), numpy.log2(power), numpy.log2(levels),
            extend='both', cmap=pyplot.cm.viridis)
extent = [t.min(), t.max(), 0, max(period)]
bx.contour(t, numpy.log2(period), sig95, [-99, 1], colors='k', linewidths=2,
           extent=extent)
bx.fill(numpy.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                           t[:1] - dt, t[:1] - dt]),
        numpy.concatenate([numpy.log2(coi), [1e-9], numpy.log2(period[-1:]),
                           numpy.log2(period[-1:]), [1e-9]]),
        'k', alpha=0.3, hatch='x')
bx.set_title('b) {} Wavelet Power Spectrum ({})'.format(label, mother.name))
bx.set_ylabel('Period (years)')
#
Yticks = 2 ** numpy.arange(numpy.ceil(numpy.log2(period.min())),
                           numpy.ceil(numpy.log2(period.max())))
bx.set_yticks(numpy.log2(Yticks))
bx.set_yticklabels(Yticks)

# Third sub-plot, the global wavelet and Fourier power spectra and theoretical
# noise spectra. Note that period scale is logarithmic.
cx = pyplot.axes([0.77, 0.37, 0.2, 0.28], sharey=bx)
cx.plot(glbl_signif, numpy.log2(period), 'k--')
cx.plot(var * fft_theor, numpy.log2(period), '--', color='#cccccc')
cx.plot(var * fft_power, numpy.log2(1./fftfreqs), '-', color='#cccccc',
        linewidth=1.)
cx.plot(var * glbl_power, numpy.log2(period), 'k-', linewidth=1.5)
cx.set_title('c) Global Wavelet Spectrum')
cx.set_xlabel(r'Power [({})^2]'.format(units))
cx.set_xlim([0, glbl_power.max() + var])
cx.set_ylim(numpy.log2([period.min(), period.max()]))
cx.set_yticks(numpy.log2(Yticks))
cx.set_yticklabels(Yticks)
pyplot.setp(cx.get_yticklabels(), visible=False)

# Fourth sub-plot, the scale averaged wavelet spectrum.
dx = pyplot.axes([0.1, 0.07, 0.65, 0.2], sharex=ax)
dx.axhline(scale_avg_signif, color='k', linestyle='--', linewidth=1.)
dx.plot(t, scale_avg, 'k-', linewidth=1.5)
dx.set_title('d) {}--{} year scale-averaged power'.format(2, 8))
dx.set_xlabel('Time (year)')
dx.set_ylabel(r'Average variance [{}]'.format(units))
ax.set_xlim([t.min(), t.max()])

pyplot.show()


# %%
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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
import pycwt as wavelet

from pycwt.helpers import find
from matplotlib.image import NonUniformImage

data1 = dict(name='Arctic Oscillation', nick='AO', file='jao.dat')
data2 = dict(name='Baltic Sea ice extent', nick='BMI', file='jbaltic.dat')
mother = 'morlet'

# Loads the data to be analysed.
t1, s1 = np.loadtxt(data1['file'], unpack=True)
t2, s2 = np.loadtxt(data2['file'], unpack=True)
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
mother = wavelet.Morlet(6)          # Morlet mother wavelet with m=6
slevel = 0.95                       # Significance level
dj = 1/12                           # Twelve sub-octaves per octaves
s0 = -1  # 2 * dt                   # Starting scale, here 6 months
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

# First plot is of both CWT
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)

extent1 = [t1.min(), t1.max(), 0, max(period1)]
extent2 = [t2.min(), t2.max(), 0, max(period2)]
im1 = NonUniformImage(ax1, interpolation='bilinear', extent=extent1)
im1.set_data(t1, period1, power1)
ax1.images.append(im1)
ax1.contour(t1, period1, sig95_1, [-99, 1], colors='k', linewidths=2,
            extent=extent1)
ax1.fill(np.concatenate([t1, t1[-1:]+dt, t1[-1:]+dt, t1[:1]-dt, t1[:1]-dt]),
         np.concatenate([coi1, [1e-9], period1[-1:], period1[-1:], [1e-9]]),
         'k', alpha=0.3, hatch='x')
ax1.set_title('{} Wavelet Power Spectrum ({})'.format(data1['nick'],
                                                      mother.name))

im2 = NonUniformImage(ax2, interpolation='bilinear', extent=extent2)
im2.set_data(t2, period2, power2)
ax2.images.append(im2)
ax2.contour(t2, period2, sig95_2, [-99, 1], colors='k', linewidths=2,
            extent=extent2)
ax2.fill(np.concatenate([t2, t2[-1:]+dt, t2[-1:]+dt, t2[:1]-dt, t2[:1]-dt]),
         np.concatenate([coi2, [1e-9], period2[-1:], period2[-1:], [1e-9]]),
         'k', alpha=0.3, hatch='x')
ax2.set_xlim(max(t1.min(), t2.min()), min(t1.max(), t2.max()))
ax2.set_title('{} Wavelet Power Spectrum ({})'.format(data2['nick'],
                                                      mother.name))


# II. Cross-wavelet transform
# ===========================

# Due to the difference in the time series, the second signal has to be
# trimmed for the XWT process.
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

# Calculate the wavelet coherence (WTC). The WTC finds regions in time
# frequency space where the two time seris co-vary, but do not necessarily have
# high power.
WCT, aWCT, corr_coi, freq, sig = wavelet.wct(s1, s2, dt, dj=1/12, s0=-1, J=-1,
                                             significance_level=0.8646,
                                             wavelet='morlet', normalize=True,
                                             cache=True)

cor_sig = np.ones([1, n]) * sig[:, None]
cor_sig = np.abs(WCT) / cor_sig  # Power is significant where ratio > 1
cor_period = 1 / freq

# Calculates the phase between both time series. The phase arrows in the
# cross wavelet power spectrum rotate clockwise with 'north' origin.
# The relative phase relationship convention is the same as adopted
# by Torrence and Webster (1999), where in phase signals point
# upwards (N), anti-phase signals point downwards (S). If X leads Y,
# arrows point to the right (E) and if X lags Y, arrow points to the
# left (W).
angle = 0.5 * np.pi - aWCT
u, v = np.cos(angle), np.sin(angle)

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.55, 0.05, 0.35])
cbar_ax_1 = fig.add_axes([0.85, 0.05, 0.05, 0.35])

extent_cross = [t1.min(), t1.max(), 0, max(cross_period)]
extent_corr = [t1.min(), t1.max(), 0, max(cor_period)]
im1 = NonUniformImage(ax1, interpolation='bilinear', extent=extent_cross)
im1.set_data(t1, cross_period, cross_power)
ax1.images.append(im1)
ax1.contour(t1, cross_period, cross_sig, [-99, 1], colors='k', linewidths=2,
            extent=extent_cross)
ax1.fill(np.concatenate([t1, t1[-1:]+dt, t1[-1:]+dt, t1[:1]-dt, t1[:1]-dt]),
         np.concatenate([cross_coi, [1e-9], cross_period[-1:],
                         cross_period[-1:], [1e-9]]),
         'k', alpha=0.3, hatch='x')
ax1.set_title('Cross-Wavelet')
ax1.quiver(t1[::3], cross_period[::3], u[::3, ::3], v[::3, ::3],
           units='width', angles='uv', pivot='mid', linewidth=1,
           edgecolor='k', headwidth=10, headlength=10, headaxislength=5,
           minshaft=2, minlength=5)
fig.colorbar(im1, cax=cbar_ax)

im2 = NonUniformImage(ax2, interpolation='bilinear', extent=extent_corr)
im2.set_data(t1, cor_period, WCT)
ax2.images.append(im2)
ax2.contour(t1, cor_period, cor_sig, [-99, 1], colors='k', linewidths=2,
            extent=extent_corr)
ax2.fill(np.concatenate([t1, t1[-1:]+dt, t1[-1:]+dt, t1[:1]-dt, t1[:1]-dt]),
         np.concatenate([corr_coi, [1e-9], cor_period[-1:], cor_period[-1:],
                         [1e-9]]),
         'k', alpha=0.3, hatch='x')
ax2.set_title('Cross-Correlation')
ax2.quiver(t1[::3], cor_period[::3], u[::3, ::3], v[::3, ::3], units='height',
           angles='uv', pivot='mid', linewidth=1, edgecolor='k',
           headwidth=10, headlength=10, headaxislength=5, minshaft=2,
           minlength=5)
ax2.set_ylim(2, 35)
ax2.set_xlim(max(t1.min(), t2.min()), min(t1.max(), t2.max()))
fig.colorbar(im2, cax=cbar_ax_1)

plt.draw()
plt.show()





# %%

wrf_file = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d01_2021-07-19_18:00:00'


div = get_div_wrfout(wrf_file)
div
# %%
flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
# ds_model = xr.open_dataset(flnm_model)
ds_obs = xr.open_dataset(flnm_obs)
ds_obs.time
# ds_obs.time.values = ds_obs.time.values+pd.Timedelta('12H')
# %%
tt = ds_obs.time.values+pd.Timedelta('8H')
ds_obs = ds_obs.assign_coords({'time':tt})
# %%
ds_obs.time
# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_A.nc'
ds = xr.open_dataset(flnm)
ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
hh = ds['height'].mean(dim='time').values
ds2 = ds.assign_coords({'z':('pressure',hh)})
ds3 = ds2.swap_dims({'pressure':'z'})
da = ds3['div']#.sel(z=1000, method='nearest')
db = da.sel(z = np.sort(da.z))
dc = db.sel(z=1000, method='nearest')
dc = dc*10**5
da = dc
time = da.time.values

sst = da.values
# %%
# dc.plot()
a = da.values
a
# %%
import numpy as np

# np.fft(a)
# np.fft.fft2(a)
np.fft.ifftshift(a)

# np.fft(a)
# %%
# %%
"""
1. 原始数据是2021年7月17日00时~7月23日00时的降水数据
2. 每隔1小时记录一次， 采样频率f_s = 1/(60*60)
3. 总的数据点个数N是145

有关 [离散傅里叶变化] 的量
0. 最大周期， N*T_s, N小时， N*60*60 s
1. 最小频率(最大周期对应的, 最慢的波) = 1/(N*T_s) = f_s/N = 1/(60*60)/145, 其他频率都是这个最小频率的n倍
2. 对原始数据了解的频率范围: 0 -- 1/2*f_s  ?

"""
cm = 1/2.54
fig = plt.figure(figsize=(16*cm, 8*cm), dpi=600)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
x = da.time.dt.strftime('%d%H')
y = a
ax.plot(x,y)

# %%

X = fft(da.values[1:])   # 变为圆频率的函数
N = len(X)   # 采样点数
n = np.arange(N)
sr = 1/(60*60)   # 采样频率， 每秒采样多少次
T = N/sr   # 最大周期, 整个区间作为一个完整的波对应的周期
freq = n/T   # 不同波对应的频率, 1/T是最小频率，所有的频率是1/T的整数倍

n_oneside = int(N/2)  # 有个共轭的复数, 这里最好保证这个N是偶数
f_oneside = freq[0:n_oneside]

t_h = 1/f_oneside/(60*60)   # 把横坐标变为小时







fig = plt.figure(figsize=(16*cm, 4*cm), dpi=600)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# ax.plot(f_oneside, np.abs(X[:n_oneside]), 'black')
ax.plot(f_oneside, np.abs(X[:n_oneside]**2/n_oneside), 'black')
# ax.plot(t_h, np.abs(X[:n_oneside]**2/n_oneside), 'black')
# ax.plot(t_h, np.abs(X[:n_oneside]), 'red')
# ax.stem(t_h, np.abs(X[:n_oneside]), 'b', markerfmt=' ', basefmt='-b')


# ax.stem(f_oneside, np.abs(X[:n_oneside]), 'b', markerfmt=' ', basefmt='-b')
# ax.stem(t_h, np.abs(X[:n_oneside]), 'b', markerfmt=' ', basefmt='-b')

# ax.set_xlim(2, 150)
# ax.set_xlim(2, 20)
# ax.set_ylim(0, 200)
# %%
# 1/145
# freq[1]-freq[0]
fig = plt.figure(figsize=(8*cm, 4*cm), dpi=600)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])


## 去除噪声
A = np.abs(X)
# idx = A>200
# idx = A>150 
# idx = A<120 
idx = (A>150) & (A<180)
# idx = (A>100) & (A<120)
# idx = (A<100)
freq_clean = idx[:n_oneside] * X[:n_oneside]
t_h = 1/f_oneside/(60*60)   # 把横坐标变为小时
# ax.stem(f_oneside, np.abs(freq_clean), 'b', markerfmt=' ', basefmt='-b')
ax.stem(t_h, np.abs(freq_clean), 'b', markerfmt=' ', basefmt='-b')
ax.set_xlim(0, 24)
# ax.set_ylim(0, 200)

# %%
## 逆傅里叶变换
iX = ifft(freq_clean)
fig = plt.figure(figsize=(8*cm, 4*cm), dpi=600)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
y = da.values[0:int(N/2)]
# ax.plot(t_h, iX)
ax.plot(f_oneside, y, color='blue')
ax.plot(f_oneside, iX, color='red')

# ax.plot(t_h, y, color='blue')
# ax.plot(t_h, iX, color='red')
# ax.set_xlim(0, 24)
# %%
# da.values[:N/2]
# %%
# type(N)
# int(N/2)

# %%
import numpy as np
from scipy.stats.distributions import chi2
'''
功率谱分析
输入：
x：需要分析的时间序列(原始序列，未标准化或距平处理)
m：最大滞后相关长度，m取值范围最好在(n/10)~(n/3)之间，n为样本数，可以多次调整m获得最佳效果，通常取m=n/3
alpha1：红噪音检验信度
alpha2：白噪音检验信度
输出：
l：功率谱图的X坐标，对应的周期为2m/l，使用时自己调整tick labels
Sl：功率谱估计值
Sr：红噪音
Sw：白噪音
r1：落后一个时刻的自相关函数，用于查看使用哪种噪音检验
'''
def specx_anal(x,m,alpha1,alpha2):
    n = x.shape[0]
    x = (x - np.mean(x))/np.std(x)
    r1 = np.zeros((n-6))
    r2 = np.zeros((n-7))
    for i in np.arange(0,n-6):
        r1[i]=np.sum(x[:n-i]*x[i:])/x[:n-i].shape[0]
    for i in np.arange(1,n-6):
        r2[i-1]=np.sum(x[:n-i]*x[i:])/x[:n-i].shape[0]
    r2 = r2[::-1]
    r = np.hstack((r2,r1))
    l = np.arange(0,m+1,1)
    tao = np.arange(1,m,1)
    Sl  = np.zeros((m+1))
    Tl  = np.zeros((m+1))
    S0l = np.zeros((m+1))
    a = np.array((r.shape[0]+1)/2).astype('int32')
    r = r[a-1:a+m]
    a=r[1:-1]*(1+np.cos(np.pi*tao/m))
    for i in np.arange(2,m+1,1):
        Sl[i-1]=(r[0]+np.sum(a*np.cos(l[i-1]*np.pi*tao/m)))/m 
    Sl[0]=(r[0]+np.sum(a*np.cos(l[0]*np.pi*tao/m)))/(2*m)
    Sl[-1]=(r[0]+np.sum(a*np.cos(l[-1]*np.pi*tao/m)))/(2*m)
    for i in range(l.shape[0]):
        Tl[i]=2*m/l[i]
    f=(2*n-m/2)/m
    S=np.mean(Sl)
    for i in range(l.shape[0]):
        S0l[i]=S*(1-r[1]*r[1])/(1+r[1]*r[1]-2*r[1]*np.cos(l[i]*np.pi/m))
    x2r = chi2.ppf(1-alpha1,df = f)
    Sr=S0l*x2r/f
    x2w = chi2.ppf(1-alpha2,df = f)
    Sw=S*x2w/f;
    r1=r[1]
    return l,Sl,Sr,Sw,r1

# %%
x = da.values
# l,Sl,Sr,Sw,r1 = specx_anal(x, 15, 0.1, 0.1)
l,Sl,Sr,Sw,r1 = specx_anal(x, 15, 0.1, 0.1)
plt.plot(l,Sl,'-b',label='Real')
plt.plot(l,Sr,'--r',label='red noise')
plt.plot(l,np.linspace(Sw,Sw,l.shape[0]),'--m',label='white noise')
plt.legend()
plt.show()
print(r1)
