#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2022/04/21 22:02:25
Author           :Forxd
Version          :1.0
'''


# %%
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import xarray as xr
import pandas as pd
from geopy.distance import distance
import numpy as np
from waveletFunctions import wave_signif, wavelet


# %%
def  get_data(model='CTRL'):
    if model == 'CTRL':
        flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/cross_rain.nc'
    elif model == 'GWD3':
        flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_rain.nc'
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_rain.nc'
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross2.nc'
    ds = xr.open_dataset(flnm)
    ds1 = ds.sel(time='2021-07-20 12')
    da = ds1['wa_cross']
    da1 = da.interpolate_na(dim='vertical', method='linear',  fill_value="extrapolate")
    da2 = da1.sel(vertical=2000, method='nearest')
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
    da3
    return da3


def caculate(sst):
    # sst = get_data()
    time = sst.distance.values
    sst = sst.values

    # sst = sst - np.mean(sst)
    variance = np.std(sst, ddof=1) ** 2
    print("variance = ", variance)
    if 0:
        variance = 1.0
        sst = sst / np.std(sst, ddof=1)
    n = len(sst)
    dt = 0.5
    # dt = 3   # 对应着小波的平移
    pad = 0  # pad the time series with zeroes (recommended)
    dj = 0.25  # this will do 4 sub-octaves per octave
    s0 = 2 * dt  # this says start at a scale of 6 months
    j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.72  # lag-1 autocorrelation for red noise background
    # print("lag1 = ", lag1)
    mother = 'MORLET'

    # Wavelet transform:
    wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
    power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
    global_ws = (np.sum(power, axis=1) / n)  # time-average over all times

    # Significance levels:
    n = len(sst)
    signif = wave_signif(([variance]), dt=dt, sigtest=0, scale=scale,
        lag1=lag1, mother=mother)
    # expand signif --> (J+1)x(N) array
    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    # power.shape
    # sig95.shapek
    sig95 = power / sig95  # where ratio > 1, power is significant

    # Global wavelet spectrum & significance levels:
    dof = n - scale  # the -scale corrects for padding at edges
    global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1,
        lag1=lag1, dof=dof, mother=mother)

    # Scale-average between El Nino periods of 2--8 years
    avg = np.logical_and(scale >= 2, scale < 8)

    Cdelta = 0.776  # this is for the MORLET wavelet
    # expand scale --> (J+1)x(N) array
    scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    scale_avg = power / scale_avg  # [Eqn(24)]
    scale_avg = dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
    scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2,
        lag1=lag1, dof=([2, 7.9]), mother=mother)
    # return time, sst, period, power
    return time, period, power, sst, sig95, coi, global_signif, global_ws
# ------------------------------------------------------ Plotting

def draw(time, period, power, sst, sig95, coi, ax3, fig):

    # --- Contour plot wavelet power spectrum
    # levels = [0, 0.05, 0.1, 0.2, 0.4, 999]
    # # levels = [-0.2, -0.1,0, 0.1,  0.2, 999]
    # # *** or use 'contour'
    # CS = ax2.contourf(time, period, power, len(levels))
    # im = ax2.contourf(CS, levels=levels,
    #     colors=['white', 'bisque', 'orange', 'orangered', 'darkred'])
    # ticks = levels[1:-1]
    # print(sst.max())

    # ax2.set_xlabel('距离 (km)', labelpad=1)
    # ax2.set_ylabel('波长 (km)')
    # # plt.title('b) Wavelet Power Spectrum (contours at 0.5,1,2,4\u00B0C$^2$)')
    # # plt.xlim(xlim[:])
    # # 95# significance contour, levels at -99 (fake) and 1 (95# signif)
    # ax2.contour(time, period, sig95, [-99, 1], colors='k')
    # # cone-of-influence, anything "below" is dubious
    # ax2.fill_between(time, coi * 0 + period[-1], coi, facecolor="none",
    #     edgecolor="#00000040", hatch='x')
    # # plt.plot(time, coi, 'k')
    # ax2.plot(time, coi, 'k')
    # # format y-scale
    # ax2.set_yscale('log', base=2, subs=None)
    # ax2.set_ylim([np.min(period), np.max(period)])
    # # ax2.set_ylim(0, 16)
    # # ax2.set_yticks([0, 16])
    # ax = plt.gca().yaxis
    # ax.set_major_formatter(ticker.ScalarFormatter())
    # ax2.ticklabel_format(axis='y', style='plain')
    # plt.colorbar(im, cax=CS, orientation='horizontal')

    #   , fraction=0.05, pad=0.5)

    # plt.subplots_adjust(right=0.7, top=0.9)

    # --- Plot global wavelet spectrum
    # ax3 = plt.subplot(gs[1, -1])
    # ax3 = ax2
    ax3.plot(global_ws, period)
    ax3.plot(global_signif, period, '--')
    ax3.set_xlabel('Power $(m^2/s^2)$')
    # ax3.set_title('c) Global Wavelet Spectrum')
    plt.xlim([0, 1.25 * np.max(global_ws)])
    # format y-scale
    ax3.set_yscale('log', base=2, subs=None)
    plt.ylim([np.min(period), np.max(period)])
    ax = plt.gca().yaxis
    ax.set_major_formatter(ticker.ScalarFormatter())
    ax3.ticklabel_format(axis='y', style='plain')
    ax3.set_xlim(0,3)

    '''
    # --- Plot 2--8 yr scale-average time series
    ax4 = plt.subplot(gs[2, 0:3])
    ax4.plot(time, scale_avg, 'k')
    # plt.xlim(xlim[:])
    ax4.set_xlabel('距离(km)')
    ax4.set_ylim(0,0.7)
    ax4.set_ylabel('Avg variance $(m/s)$')
    '''
    # ax1.set_title()
    # ax3 = fig.add_axes([0.2,0.08,0.7,0.03])
    # cb = fig.colorbar(
    #     im,
    #     cax=ax3,
    #     orientation='horizontal',
    #     ticks=ticks,
    #     fraction = 0.1,  # 色标大小,相对于原图的大小
    #     pad=0.5,  #  色标和子图间距离
    #     )
    # # ax1.set_xticks([0,50,100,150])
    # ax2.set_xticks([0,50,100,150])
    # # ax1.set_xlim(0,150)
    # ax2.set_xlim(0,150)
    # ax2.set_yticks([1, 2,4,8,16, 32, 64])
    # ax2.set_ylim(0,64)
    # # ax1.set_yticks([-2,-1,0,1,2])
    # # ax1.set_title('(a)', loc='left', y=0.97)
    # # ax2.set_title('(b)', loc='left', y=0.97)
    # return im, ticks

if __name__ == '__main__':
    
    cm = 1/2.4
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=600)
    ax1 = fig.add_axes([0.2,0.6,0.6,0.3])

    sst = get_data('CTRL')
    time, period, power, sst, sig95, coi, global_signif, global_ws = caculate(sst)
    draw(time, period, power, sst, sig95, coi, ax1, fig)

    ax2 = fig.add_axes([0.2,0.15,0.6,0.3])
    sst = get_data('GWD3')
    time, period, power, sst, sig95, coi, global_signif, global_ws = caculate(sst)
    draw(time, period, power, sst, sig95, coi, ax2, fig)
    # time, period, power, sst, sig95, coi,global_signif, global_ws = caculate(sst)
    # # im, ticks = draw(time, period, power, sst, sig95, coi, ax3,ax4, fig)
    # draw(time, period, power, sst, sig95, coi, global_signif, global_ws, ax1,ax3, fig)

    # ax1.set_title()
    # ax3 = fig.add_axes([0.2,0.08,0.7,0.03])
    # ax3 = fig.add_axes([0.2,0.08,0.7,0.03])
    # ax5 = fig.add_axes([0.3,0.06,0.5,0.03])
    # ax5 = fig.add_axes([0.84,0.3,0.03,0.5])
    # # # ax3.set_title('')
    # cb = fig.colorbar(
    #     im,
    #     cax=ax5,
    #     # orientation='horizontal',
    #     orientation='vertical',
    #     ticks=ticks,
    #     # fraction = 0.05,  # 色标大小,相对于原图的大小
    #     # pad=0.5,  #  色标和子图间距离
    #     label='aaaaa'
    #     )
    
    # ax1.set_title('(a)', loc='left', y=0.97)
    # ax2.set_title('(b)', loc='left', y=0.97)

    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/wavelet_power.png')