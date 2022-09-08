# %%
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import xarray as xr
import pandas as pd
from geopy.distance import distance
import numpy as np
from waveletFunctions import wave_signif, wavelet
import wrf
import netCDF4 as nc
# %%

def get_data_rain():
    area = {
        'lat1':33.5,
        'lat2':36.0,
        'lon1':112,
        'lon2':115,
        }        
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_model = xr.open_dataset(flnm_model)
    ds_obs = xr.open_dataset(flnm_obs)
    # ds_obs.time.values = ds_obs.time.values+pd.Timedelta('12H')

    # tt = ds_obs.time.values+pd.Timedelta('8H')
    # ds_obs = ds_obs.assign_coords({'time':tt})

    # cm = Common()
    gd = GetData()
    ds_model_mean = gd.caculate_area_mean(ds_model, area)
    ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)
    dsall = xr.merge([ds_obs_mean, ds_model_mean])
    ds = dsall.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    # ds = dsall.sel(time=slice('2021-07-17 00', '2021-07-19 00'))
    # ds = ds.resample(time='12H').sum()
    return ds


def get_data_div():
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_D.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_D.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
    hh = ds['height'].mean(dim='time').values
    ds2 = ds.assign_coords({'z':('pressure',hh)})
    ds3 = ds2.swap_dims({'pressure':'z'})
    da = ds3['div']#.sel(z=1000, method='nearest')
    db = da.sel(z = np.sort(da.z))
    dc = db.sel(z=1500, method='nearest')
    dc = dc*10**5
    da = dc
    time = da.time.values
    sst = da.values
    return sst, time, da

def  get_data_vs():
    """vertical wind speed
    Returns:
        _type_: _description_
    """
    # flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/cross_rain.nc'
    flnm='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross_rain.nc'
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



sst, time, da = get_data_div()
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

# ------------------------------------------------------ Plotting
cm = 1/2.54
fig = plt.figure(figsize=(8*cm, 6*cm), dpi=600)
ax = fig.add_axes([0.15, 0.2, 0.7, 0.7])


# --- Contour plot wavelet power spectrum
levels = [-40, -20,0, 20,  40, 999]
# *** or use 'contour'
# CS = ax.contourf(time, period, power, len(levels))
# im = ax.contourf(CS, levels=levels,
#     colors=['white', 'bisque', 'orange', 'orangered', 'darkred'])

# colorlevel=[0, 1, 10, 25, 50, 100, 250, 400,600,800,1000, 2000]#雨量等级
# colorlevel=[0, 10, 20, 30, 50, 70, 100, 150, 200, 400,800,1000]
colorlevel=[0, 1, 2, 3, 5, 7, 10, 15, 20, 40,80,100]
# colorlevel=[-60, -40,-30,  -20 , -10, -5, 5,10,20 ,30, 40, 60]
# colorticks=[-40, -20 , 0, 20 , 40]
colorticks = colorlevel[1:-1]
rgbtxt = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/rain_6d/11colors.txt'
def get_rgb(fn):
    """
    fn: rgb_txt文件存储路径
    """
    # fn = './11colors.txt'
    df = pd.read_csv(fn, skiprows=4, sep='\s+',encoding='gbk',header=None, names=['r','g','b'])
    rgb = []
    for ind, row in df.iterrows():
        rgb.append(row.tolist())
    rgb = np.array(rgb)/255.
    return rgb
rgb = get_rgb(rgbtxt)
colordict = rgb  

crx = ax.contourf(time,
                    period,
                    power,
                    corner_mask=False,
                    levels=colorlevel,
                    colors = colordict,
                    )


cb = fig.colorbar(
    crx,
    # cax=ax6,
    orientation='vertical',
    # orientation='horizontal',
    ticks=colorticks,
    fraction = 0.06,  # 色标大小,相对于原图的大小
    pad=0.05,  #  色标和子图间距离
    )



# ax.set_xlabel('时间 （日期/小时）')
ax.set_xlabel('时间 （日期）')
ax.set_ylabel('周期 （小时）')
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
ax.contour(time, period, sig95, [-99, 1], colors='k')
# cone-of-influence, anything "below" is dubious
ax.fill_between(time, coi * 0 + period[-1], coi, facecolor="none",
    edgecolor="#00000040", hatch='x')
# plt.plot(time, coi, 'k')
ax.plot(time, coi, 'k')
# format y-scale
ax.set_yscale('log', base=2, subs=None)
ax.set_yticks([0, 1, 2, 3, 4, 6,10, 16])
ax.set_ylim(1, 32)
x = da.time.dt.strftime('%d')
ax.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))

## 不用指数形式标注纵坐标
axx = plt.gca().yaxis
axx.set_major_formatter(ticker.ScalarFormatter())
ax.ticklabel_format(axis='y', style='plain')
# %%
figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/'
fig.savefig(figpath+'wave')
