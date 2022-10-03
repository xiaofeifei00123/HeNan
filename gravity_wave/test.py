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
from netCDF4 import Dataset
from wrf import getvar, interpline, CoordPair, vertcross
from geopy.distance import distance  # 根据经纬度计算两点距离
# %%
flnm1 = '/home/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d01_2021-07-20_00:00:00'
flnm2 = '/home/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d02_2021-07-20_00:00:00'
flnm3 = '/home/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d03_2021-07-20_00:00:00'
ds1 = xr.open_dataset(flnm1)
ds2 = xr.open_dataset(flnm2)
ds3 = xr.open_dataset(flnm3)
# %%
ds1['DUSFCG_BL'].max()
# ds2['DUSFCG_BL'].max()
# ds3['DUSFCG_BL'].max()

ds1['DTAUX3D_BL'].max()
ds2['DTAUX3D_BL'].max()
ds3['DTAUX3D_BL'].max()
# %%
for i in list(ds3.data_vars):
    print(i)
# ds1['DTAUX3D_BL'].max()

# %%
# flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
# ds_obs = xr.open_dataset(flnm_obs)
# ds_obs.time
# tt = ds_obs.time.values+pd.Timedelta('8H')
# ds_obs = ds_obs.assign_coords({'time':tt})
# # ds_obs.time
# ds_obs['PRCP']

# %%
def get_rain_wrf():
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
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
    area = {
        'lat1':32,
        'lat2':36.5,
        'lon1':110.5,
        'lon2':116,
        }        
    ds_model_mean = caculate_area_mean(ds_model, area)
    ds_model_mean = ds_model_mean['GWD3']
    da = ds_model_mean.sel(time=slice('2021-07-17 00', '2021-07-23 00'))
    sst = da.values
    tt = da.time.values
    return sst, tt, da
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

a = sst
# np.fft(a)
# np.fft.fft2(a)
# np.fft.ifftshift(a)

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
# ax.plot(f_oneside, np.abs(X[:n_oneside]**2/n_oneside), 'black')
ax.set_xlim(0, 20)
ax.set_ylim(0, 4)
ax.plot(t_h, np.abs(X[:n_oneside]**2/n_oneside), 'black')
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
# idx = (A>150) & (A<180)
# idx = (A>100) & (A<120)
idx = (A<200)
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
sst, tt, da = get_rain_wrf()
sst

# %%
sst.shape

# %%
# x = da.values
x = sst
# l,Sl,Sr,Sw,r1 = specx_anal(x, 15, 0.1, 0.1)
l,Sl,Sr,Sw,r1 = specx_anal(x, 14, 0.1, 0.1)
plt.plot(l,Sl,'-b',label='Real')
plt.plot(l,Sr,'--r',label='red noise')
plt.plot(l,np.linspace(Sw,Sw,l.shape[0]),'--m',label='white noise')
plt.legend()
plt.show()
print(r1)
