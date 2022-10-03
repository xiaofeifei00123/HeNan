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
from read_rain_wrf import GetData
from common import Common
# %%


# %%
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
    da3
    return da3

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

# %%
def get_div_distance():
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/cross2.nc' 
    ds = xr.open_dataset(flnm)
    div = ds['div_cross'].sel(time='2021-07-20 00')#.sel(vertical=2000, method='nearest')
    div = drop_na(div)
    # div2 = div.sel(vertical=1000, method='nearest')
    div2 = div.sel(vertical=4200, method='nearest')
    sst = div2.values*10**5
    tt = latlon2distance(div2).distance.values
    da = latlon2distance(div2)
    # da = div2
    return sst, tt, da

def get_vs_distance():
    """垂直速度随时间变化

    Returns:
        _type_: _description_
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/cross2.nc' 
    ds = xr.open_dataset(flnm)
    div = ds['wa_cross'].sel(time='2021-07-20 00')#.sel(vertical=2000, method='nearest')
    div = drop_na(div)
    div2 = div.sel(vertical=4200, method='nearest')
    div2
    sst = div2.values*10
    tt = latlon2distance(div2).distance.values
    da = latlon2distance(div2)
    # da = div2
    return sst, tt, da

def get_data_rain_obs():
    def caculate_area_mean_obs(da,area):
        mask = (
            (da.coords['lat']>area['lat1'])
            &(da.coords['lat']<area['lat2'])
            &(da.coords['lon']<area['lon2'])
            &(da.coords['lon']>area['lon1'])
        )
        aa = xr.where(mask, 1, np.nan)
        db = da*aa
        dsr = db.mean(dim=['lat', 'lon'])
        return dsr
    # area = {
    #     'lat1':33.5,
    #     'lat2':36.0,
    #     'lon1':112.2,
    #     'lon2':114.8,
    #     }        
    # area = {
    #     'lat1':32,
    #     'lat2':36.5,
    #     'lon1':110.5,
    #     'lon2':116,
    #     }        
    area = {
        'lat1':34.4,
        'lat2':34.9,
        'lon1':113.0,
        'lon2':113.7,
    }        
    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_obs = xr.open_dataset(flnm_obs)
    gd = GetData()
    ds_obs_mean  = caculate_area_mean_obs(ds_obs, area)
    sst = ds_obs_mean['PRCP'].values
    tt = ds_obs_mean.time.values
    da = ds_obs_mean['PRCP']
    return sst, tt, da

def get_rain_wrf():
    # flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/rain_model_da.nc'
    ds_model = xr.open_dataset(flnm_model)
    def caculate_area_mean(da, area,):
        lon = da['lon'].values
        lat = da['lat'].values
        #     ## 构建掩膜, 范围内的是1， 不在范围的是nan值
            
        clon = xr.where((lon<area['lon2']) & (lon>area['lon1']), 1, np.nan)
        clat = xr.where((lat<area['lat2']) & (lat>area['lat1']), 1, np.nan)
        da = da*clon*clat
        # if 'south_north' in list(da.dims):
        da_mean = da.mean(dim=['south_north', 'west_east'])
        da_mean
        return da_mean

    # area = {
    #     'lat1':33.5,
    #     'lat2':36.0,
    #     'lon1':112,
    #     'lon2':115,
    #     }        
    # area = {
    #     'lat1':32,
    #     'lat2':36.5,
    #     'lon1':110.5,
    #     'lon2':116,
    #     }        
    area = {
        'lat1':34.4,
        'lat2':34.9,
        'lon1':113.0,
        'lon2':113.7,
    }        
    ds_model_mean = caculate_area_mean(ds_model, area)
    ds_model_mean = ds_model_mean['GWD3']
    da = ds_model_mean.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    sst = da.values
    tt = da.time.values
    return sst, tt, da






def get_data_div():
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_D.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_D.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_C.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_D.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/FD/wrfout/time_cross_D.nc'
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
    return sst, time, da

def get_data_vertical_speed():
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_C.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_A.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_D.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_D.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_E.nc'
    
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_B.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_B.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/FD/wrfout/time_cross_B.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/SS/wrfout/time_cross_B.nc'

    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_A.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
    hh = ds['height'].mean(dim='time').values
    ds2 = ds.assign_coords({'z':('pressure',hh)})
    ds3 = ds2.swap_dims({'pressure':'z'})
    # da = ds3['div']#.sel(z=1000, method='nearest')
    da = ds3['wa']#.sel(z=1000, method='nearest')
    db = da.sel(z = np.sort(da.z))
    # dc = db.sel(z=10000, method='nearest')
    # dc = db.sel(z=[8000, 9000, 10000,11000,  12000], method='nearest').mean(dim='z')

    # dcc = xr.where((db.z>3000)&(db.z<12000), db, np.nan)
    # dcc = xr.where((db.z>5000)&(db.z<12000), db, np.nan)

    # dcc = xr.where((db.z>0)&(db.z<3000), db, np.nan)
    # dc = dcc.mean(dim='z')

    # dc = db.sel(z=3000, method='nearest')
    # dc = db.sel(z=4200, method='nearest')
    dc = db.sel(z=1500, method='nearest')
    # dc = dc*10
    dc = dc*10**2
    da = dc
    time = da.time.values
    sst = da.values
    return sst, time, da

def  get_data_vs_distance():
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
    da2 = da1.sel(vertical=1000, method='nearest')
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



# sst, time, da = get_data_vertical_speed()
# sst, time, da = get_div_distance()
# sst, time, da = get_vs_distance()
# sst, time, da = get_data_div()

sst, time, da = get_data_rain_obs()
# sst, time, da = get_rain_wrf()
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
fig = plt.figure(figsize=(8*cm, 5*cm), dpi=600)
ax = fig.add_axes([0.15, 0.2, 0.7, 0.7])


# --- Contour plot wavelet power spectrum
levels = [-40, -20,0, 20,  40, 999]
# *** or use 'contour'
# CS = ax.contourf(time, period, power, len(levels))
# im = ax.contourf(CS, levels=levels,
#     colors=['white', 'bisque', 'orange', 'orangered', 'darkred'])

# colorlevel=[0, 1, 10, 25, 50, 100, 250, 400,600,800,1000, 2000]#雨量等级
colorlevel=[0, 10, 20, 30, 50, 70, 100, 150, 200, 400,800,50000]  # 垂直速度
# colorlevel=[0, 100, 150, 200, 400,800,1000, 1500, 2000, 3000, 4000, 50000]
# colorlevel=[0,  50, 100, 150, 200,250, 400,800,1200, 1600, 2000, 50000]
# colorlevel=[0, 1, 2, 3, 5, 7, 10, 15, 20, 40,80,100]
# colorlevel=[0, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2.0, 4.0,8.0,100.0]
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
# ax.set_xlabel('距离（km）')
# ax.set_ylabel('波长（km）')
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
ax.contour(time, period, sig95, [-99, 1], colors='k')
# cone-of-influence, anything "below" is dubious
ax.fill_between(time, coi * 0 + period[-1], coi, facecolor="none",
    edgecolor="#00000040", hatch='x')
# plt.plot(time, coi, 'k')
ax.plot(time, coi, 'k')
# format y-scale
ax.set_yscale('log', base=2, subs=None)
ax.set_yticks([0, 1, 2, 3, 4, 6,8, 10,13,16, 24, 32, 64])
ax.set_ylim(2, 32)
# ax.set_ylim([np.min(period), np.max(period)])
## 横坐标为距离
# x = da.distance
# ax.set_xticks(time[::24])
# ax.set_xticklabels(x.values[::24], rotation=30, fontsize=10)


## 横坐标为时间
x = da.time.dt.strftime('%d')
ax.set_xticks(time[::24])
ax.set_xticklabels(x.values[::24], fontsize=10)

## 不用指数形式标注纵坐标
axx = plt.gca().yaxis
axx.set_major_formatter(ticker.ScalarFormatter())
ax.ticklabel_format(axis='y', style='plain')
# %%
figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_wave/'
fig.savefig(figpath+'wave_time_rain_obs')
# fig.savefig(figpath+'wave_time_4km_SS')
# fig.savefig(figpath+'wave_time_4km_gwd3')
# fig.savefig(figpath+'tttt')
# fig.savefig(figpath+'wave_rain')
# fig.savefig(figpath+'wave_rain_obs')

# %%
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
# ds_obs = xr.open_dataset(flnm_obs)
# ds_obs
