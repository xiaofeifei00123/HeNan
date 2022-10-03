#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算水平散度的剖面数据
-----------------------------------------
Time             :2022/01/12 17:01:24
Author           :Forxd
Version          :1.0
'''
# %%

    
# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from draw_upar_850_d03_div import draw
import metpy
from metpy.units import units

# %%
def get_data():
    def cal_ri(ds1):
        pre = ds1.pressure*units.hPa
        temp = ds1['temp']*units.degC
        theta = metpy.calc.potential_temperature(pre, temp)
        height = ds1['geopt']*units.m
        u = ds1['u']
        v = ds1['v']
        u = u*units('m/s')
        v = v*units('m/s')
        ri = metpy.calc.gradient_richardson_number(height, theta, u, v, vertical_dim=0)
        # bv = metpy.calc.brunt_vaisala_frequency(height, theta)
        return ri

        

# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/upar_zhenzhou.nc'
ds = xr.open_dataset(flnm)
ds1 = ds.sel(time='2021-07-20 00')
ds1
# %%
from baobao.caculate import caculate_q_rh_thetaev
ds2 = ds1.reset_coords()
ds3 = caculate_q_rh_thetaev(ds1)
# %%
ds3['theta_v']


# %%
pre = ds1.pressure*units.hPa
temp = ds1['temp']*units.degC
theta = metpy.calc.potential_temperature(pre, temp)
height = ds1['geopt']*units.m
u = ds1['u']
v = ds1['v']
u = u*units('m/s')
v = v*units('m/s')
ri = metpy.calc.gradient_richardson_number(height, theta, u, v, vertical_dim=0)
ri
# %%
import metpy.calc as ca
duu = ca.gradient(u, deltas=height.values)
duu
duu[0].plot()
# %%
du = xr.DataArray((u[1:].values - u[0:-1].values), coords=u[0:-1].coords)
dv = xr.DataArray((v[1:].values - v[0:-1].values), coords=u[0:-1].coords)
dz = xr.DataArray((height[1:].values - height[0:-1].values), coords=u[0:-1].coords)
# du = xr.DataArray((u[1:].values - u[0:-1].values), coords=u[0:-1].coords)
# du
# dz
# (du/dz).plot()
ws = ((du/dz)**2+(dv/dz)**2)
# %%
dth = xr.DataArray((theta[1:].values - theta[0:-1].values), coords=u[0:-1].coords)
dth
the = theta[0:-1]
g = 9.86
bn = g/the*(dth/dz)
# %%
bn.plot()
# theta
# height.values
# %%
bnn = ca.brunt_vaisala_frequency_squared(height, theta)
bnn.plot()
# u.values
# duu = ca.gradient(u.values, coordinates=height.values)
# %%
(bnn/ws).plot(label='fri')
(bnn*10**4).plot(label='bn')
(ws*10**4).plot(label='ws')
plt.ylim(-0.5, 1)
plt.legend()


# u.shape

# %%

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

# %%

# h[0:-1]


# %%

def highPassFiltering(img,size):#传递参数为傅里叶变换后的频谱图和滤波尺寸
    h, w = img.shape[0:2]#获取图像属性
    h1,w1 = int(h/2), int(w/2)#找到傅里叶频谱图的中心点
    img[h1-int(size/2):h1+int(size/2), w1-int(size/2):w1+int(size/2)] = 0#中心点加减滤波尺寸的一半，刚好形成一个定义尺寸的滤波大小，然后设置为0
    return img

def lowPassFiltering(img,size):#传递参数为傅里叶变换后的频谱图和滤波尺寸
    h, w = img.shape[0:2]#获取图像属性
    h1,w1 = int(h/2), int(w/2)#找到傅里叶频谱图的中心点
    img2 = np.zeros((h, w), np.uint8)#定义空白黑色图像，和傅里叶变换传递的图尺寸一致
    img2[h1-int(size/2):h1+int(size/2), w1-int(size/2):w1+int(size/2)] = 1#中心点加减滤波尺寸的一半，刚好形成一个定义尺寸的滤波大小，然后设置为1，保留低频部分
    img3=img2*img #将定义的低通滤波与传入的傅里叶频谱图一一对应相乘，得到低通滤波
    return img3

def bandpass(img, w, radius):    
    """_summary_

    Args:
        img (_type_): 二维的DataArray
        w (_type_): 带宽, 
        radius (_type_): 带中心到频率平面原点的距离
    """
    pass
    rows, cols = img.shape
    crow,ccol = int(rows/2), int(cols/2) #中心位置
    w = w
    radius = radius
    mask = np.ones((rows, cols, 2), np.uint8)
    for i in range(0, rows):
        for j in range(0, cols):
            # 计算(i, j)到中心点的距离
            d = np.sqrt(pow(i - crow, 2) + pow(j - ccol, 2))
            if radius - w / 2 < d < radius + w / 2:
                mask[i, j, 0] = mask[i, j, 1] = 0
            else:
                mask[i, j, 0] = mask[i, j, 1] = 1
    #掩膜图像和频谱图像乘积

    f = img * mask[:,:,0]
    return f
    


# img_path = './picture_upar/850/div/gwd3_850_2021072017w_speed.png'
# img = cv.imread(img_path, 2)


flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/upar_latlon.nc'
ds = xr.open_dataset(flnm)
da = ds.sel(time='2021-07-20 12').sel(pressure=850)['w']
# da.max()
bb = da.interpolate_na(dim='lat', method='linear',  fill_value="extrapolate")
img = bb


f = np.fft.fft2(img)
#将左上角低频部分移动到中间
fshift = np.fft.fftshift(f)
#调用高通滤波函数
# img1=highPassFiltering(fshift,10)
# img1=lowPassFiltering(fshift,30)
img1 = bandpass(fshift, w=20, radius=0)

#复数化整，方便观察频谱图
res = np.log(np.abs(img1))
#傅里叶逆变换
ishift = np.fft.ifftshift(img1)
iimg = np.fft.ifft2(ishift)
iimg = np.real(iimg)
dda = xr.DataArray(iimg,
                    coords=da.coords,
                    dims=da.dims,
                    )
  


# dda.max()

# bb-iimg
# plt.subplot(131), plt.imshow(img,'gray'), plt.title('原图像')
# plt.axis('off')
# plt.subplot(132), plt.imshow(res,'gray'), plt.title('高通滤波')
# plt.axis('off')
# plt.subplot(133), plt.imshow(iimg,'gray'), plt.title('高通滤波结果')
# plt.axis('off')
# plt.show()

# print(flnm_wrf)
dic= {
    'model':'gwd3_once',
    'level':850,
    'time':pd.Timestamp('2021-07-20 16'),
    'flnm':flnm,
}    
ds_wrf = xr.open_dataset(flnm)
# ds_wrf = ds_wrf.rename({'ua':'u', 'va':'v'})
t = dic['time']
level = dic['level']
ds2 = ds_wrf.sel(time=t, pressure=level)
ds2
draw(dda*15, ds2['u'], ds2['v'], dic)


# %%
# iimg.values



# %%
# %%
flnm = '/home/fengxiang/HeNan/Data/OBS/rain_ec.nc'
da2 = xr.open_dataarray(flnm)
da2 = da2.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
da2 = da2.sel(lat=slice(32,36.5)).sel(lon=slice(110.5,116))
da2 = da2.sum(dim='time')


# %%
# da1
area = {
    'lon1':110.5,
    'lon2':116,
    'lat1':32,
    'lat2':36.5,
    'interval':0.125,
}
da = xr.open_dataarray('/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain.nc')
da1 = rain_station2grid(da, area)
da1 = da1.sel(time=slice('2021-07-20 01', '2021-07-21 00'))  # 24小时逐小时降水
da1 = da1.sum(dim='time')
# %%
# ddc.lat
# da2.lat
# dd
da2-da1



# da1-da_rain
# type(da_rain.dims)
# da_rain.plot()
# dd = da_rain.sum(dim='time').plot()
# dd
# da
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'+model+'/'+'rain.nc'
# flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain.nc'
# flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/rain.nc'

# %%
# flnm3 = 
flnm3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/rain_latlon.nc'
flnm0 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/rain_latlon.nc'
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d04/'+model+'/'+'rain.nc'
da3 = xr.open_dataarray(flnm3)
da3 = da3.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
da3 = da3.sum(dim='time') 
da0 = xr.open_dataarray(flnm0)
da0 = da0.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
da0 = da0.sum(dim='time') 
da0
# %%
# da0
# da3-da0
# da3.lat.values
da = da3.values-da0.values
# da = xr.DataArray
da
# %%
dda = xr.DataArray(
    da,
    coords={
        'lat':da3.lat.values,
        'lon':da3.lon.values,
    },
    dims=['lat', 'lon']
    )
dda




# da0
# k
# da0










# %%
import xarray as xr
import numpy as np
import pandas as pd
from wrf import getvar, CoordPair, vertcross, get_cartopy
import wrf
from netCDF4 import Dataset
from multiprocessing import Pool
from baobao.caculate import caculate_div
from metpy.units import units  # units是units模块中的一个变量名(函数名？类名？)
from metpy import calc as ca  # calc是一个文件夹
import matplotlib.pyplot as plt

# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/cross.nc'
ds = xr.open_dataset(flnm)
da = ds['div_cross']
da = da.isel(time=0)
# %%

# self.cross_start = CoordPair(lat=33, lon=111)
# self.cross_end = CoordPair(lat=35.5, lon=114.5)

ldic = {
    'lat1':33,
    'lon1':111,
    'lat2':35.5,
    'lon2':114.5,
}
dy = (ldic['lat2']-ldic['lat1'])
dx = (ldic['lon2']-ldic['lon1'])
angle = np.arctan2(dy,dx)  # 对边和直角边, 弧度

# np.cos(angle)
# deg = 180.0/np.pi # 角度和弧度之间的转换
# rad = np.pi/180.0


# angle
# np.arctan(dy/dx)




# %%
def draw_contour2(ax,da):
    xs = np.arange(0, da.shape[-1], 1)
    ys = da.coords['vertical'].values
    # levels=np.arange(342, 372, 4)
    cs = ax.contour(xs, ys, da.values*10**4, colors='red', zorder=3)
    ax.clabel(cs, inline=True, fontsize=10)

    

fig = plt.figure(figsize=[10,8])
ax = fig.add_axes([0.1,0.1,0.8,0.8])



# %%

def concat_dxy(var, index):
    """将二维数据在垂直方向上累加
    变成三维的

    Args:
        var ([type]): 二维变量
        index ([type]): 垂直方向的索引

    Returns:
        var_concat: 累加到一块的数据
    """
    var_list = []
    # index = u.bottom_top.values
    for i in index:
        # print(i)
        aa = var.magnitude  # 转为numpy
        bb = xr.DataArray(aa, dims=['south_north', 'west_east']) # 转为DataArray
        var_list.append(bb)
    var_concat = xr.concat(var_list, dim=pd.Index(index, name='bottom_top'))
    return var_concat

def caculate_div3d(u, v, lon, lat):
    """求wrfout数据中三维的u,v数据对应的散度
    就先原始的wrfout数据吧

    Args:
        u ([type]): 三维
        v ([type]): 三维
        lon ([type]): 二维
        lat ([type]): 二维
    """
    pass
    # u =  getvar(ncfile, 'ua')
    # v =  getvar(ncfile, 'va')
    # lon = u.XLONG
    # lat = u.XLAT
    u = u*units('m/s')
    v = v*units('m/s')
    dx, dy = ca.lat_lon_grid_deltas(lon.values, lat.values)
    ## 重组dx和dy, 其实就是把dx和dy的垂直维度加上，虽然每个垂直层上数据一样, 这是由于metpy计算时的问题导致的
    index = u.bottom_top.values
    ddx = concat_dxy(dx, index)
    ddy = concat_dxy(dy, index)
    dddx = ddx.values*units('m')
    dddy = ddy.values*units('m')
    ### 因为这个函数的问题，所以dx必须是和u维度相对应的
    div = ca.divergence(u=u, v=v, dx=dddx, dy=dddy)
    div
    return div
    # div.min()*10**3*10

# if __name__ == '__main__':
#     # main()
    
wrf_file = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d01_2021-07-19_18:00:00'
ncfile = Dataset(wrf_file)
u =  getvar(ncfile, 'ua')
v =  getvar(ncfile, 'va')
lon = u.XLONG
lat = u.XLAT
div = caculate_div3d(u,v, lon, lat)
div = div.rename('div')
div
# %%
# xr.merge([u,div])
u.projection
div.attrs = u.attrs
# %%
div.to_netcdf('aa.nc')
