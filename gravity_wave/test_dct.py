#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

图像滤波, 高通、低通、带通
参考:
fft变换(高斯) https://zhuanlan.zhihu.com/p/332427422
dct 变换  https://blog.csdn.net/weixin_45355387/article/details/123994314
-----------------------------------------
Time             :2022/07/20 12:12:15
Author           :Forxd
Version          :1.0
'''


# %%
import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import scipy
import netCDF4 as nc
import wrf
import matplotlib.pyplot as plt
import cmaps
import scipy

# %%

def high_pass_filter(img, radius=10):
    r = radius

    rows, cols = img.shape
    # center = int(rows / 2), int(cols / 2)
    center = 0,0
    mask = np.zeros((rows, cols), np.uint8)
    x, y = np.ogrid[:rows, :cols]
    mask_area = (x - center[0]) ** 2 + (y - center[1]) ** 2 <= r * r
    mask[mask_area] = 1
    mask = 1-mask
    return mask

def low_pass_filter(img, radius=10):
    r = radius

    rows, cols = img.shape
    # center = int(rows / 2), int(cols / 2)
    center = 0,0

    mask = np.zeros((rows, cols), np.uint8)
    x, y = np.ogrid[:rows, :cols]
    mask_area = (x - center[0]) ** 2 + (y - center[1]) ** 2 <= r * r
    mask[mask_area] = 1
    return mask

def band_pass_filters(img, r_in=10, r_out=30):
    rows, cols = img.shape
    crow, ccol = int(rows / 2), int(cols / 2)

    # radius_out = r_out
    # radius_in = r_in

    mask = np.zeros((rows, cols), np.uint8)
    center = [0, 0]
    x, y = np.ogrid[:rows, :cols]
    mask_area = np.logical_and(((x - center[0]) ** 2 + (y - center[1]) ** 2 >= r_in ** 2),
                               ((x - center[0]) ** 2 + (y - center[1]) ** 2 <= r_out ** 2))
    mask[mask_area] = 1  # 带通
    # mask = 1 - mask    # 带阻
    return mask



def filter_lambda_high(img, dx, lam):
    """根据波长，进行高通滤波，得到小于这个波长的短波
    """
    pass
    L = min(img.shape)*dx  # 单位km
    r = L/lam  # 波数
    mask = high_pass_filter(img,radius=r)   # 左下角为1， 其它为0
    Y = cv2.dct(img)
    Y = Y*mask  # 波数大于r_h的短波
    y2 = cv2.idct(Y)
    return y2


def filter_lambda_low(img, dx, lam):
    """进行低通,波长大于lam的长波

    Args:
        img (_type_): _description_
        dx (_type_): _description_
        lam (_type_): _description_
    """

    L = min(img.shape)*dx  # 单位km
    r = L/lam  # 波数
    mask = low_pass_filter(img,radius=r)   # 左下角为1， 其它为0
    Y = cv2.dct(img)
    Y = Y*mask  # 波数大于r_h的短波
    y2 = cv2.idct(Y)
    return y2

def filter_lambda_band(img, dx, lam1, lam2):
    """带通滤波， lam1<lam< lam2的波

    Args:
        img (_type_): _description_
        dx (_type_): _description_
        lam1 (_type_): _description_
        lam2 (_type_): _description_
    """
    L = min(img.shape)*dx  # 单位km
    r2 = L/lam1  # 波数, 波长越小，波数越大
    r1 = L/lam2  # 波数
    mask = band_pass_filters(img,r_in=r1, r_out=r2)   # 左下角为1， 其它为0
    Y = cv2.dct(img)
    Y = Y*mask
    y2 = cv2.idct(Y)
    return y2


def draw_img(img_origin, img_low, img_band, img_high,lon,lat):
    colordict=['#0000fb','#3232fd','#6464fd','white','white','#fbbcbc', '#fd4949', '#fd0000']#正负, 蓝-红
    colorlevel= [-5,-0.5,-0.1,-0.01, 0, 0.01,0.1,0.5,5]  # 垂直速度的色标
    # colorlevel = [8700,8800,8900,9000,9100,9200,9300,9400,9500]
    # colorlevel = [9050, 9100,9150,9200,9250,9300,9350,9400,9450]
    # colorlevel = [-200, -100, -50, -10, 0, 10, 50, 100, 200]
    # colorlevel = [0, 1, 10, 50, 100, 200,300,400,500]
    # colorlevel=[0, 0.1, 5, 15.0, 30, 70, 140, 250.0, 700]#雨量等级
    # colorlevel = [3000, 3050, 3100, 3150, 3200, 3250, 3300, 3350, 3400]
    # colorlevel = np.linspace(0, 4000, 9)

    cm = 1/2.54
    fig = plt.figure(figsize=[17*cm, 16*cm], dpi=300)

    ax1 = fig.add_axes([0.1, 0.5, 0.4, 0.35])
    ax2 = fig.add_axes([0.55, 0.5, 0.4, 0.35])

    ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.35])
    ax4 = fig.add_axes([0.55, 0.1, 0.4, 0.35])
    ax5 = fig.add_axes([0.1, 0.01, 0.85, 0.04])
    colorticks = colorlevel[1:-1]

    # cmap = cmaps.BlAqGrWh2YeOrReVi22
    # cmap = cmaps.prcp_1
    # crx = ax1.contourf(lon, lat,img_origin, cmap=cmap, levels=colorlevel)
    # crx = ax2.contourf(lon, lat,img_filter, cmap=cmap, levels=colorlevel)

    crx = ax1.contourf(lon, lat,img_origin, colors=colordict, levels=colorlevel)
    crx = ax2.contourf(lon, lat,img_low, colors=colordict, levels=colorlevel)
    crx = ax3.contourf(lon, lat,img_band, colors=colordict, levels=colorlevel)
    crx = ax4.contourf(lon, lat,img_high, colors=colordict, levels=colorlevel)

    ax1.set_title('a) Orig.', loc='left', y=0.85, x=0.1, bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
    ax2.set_title(r'$b) \quad \lambda>100km$', loc='left', y=0.85, x=0.1, bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
    ax3.set_title(r'$c) \quad 50<\lambda<100km$', loc='left', y=0.85, x=0.1, bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
    ax4.set_title(r'$d) \quad \lambda<50km$', loc='left', y=0.85, x=0.1, bbox=dict(facecolor='white', alpha=1, edgecolor='black'))

    cb = fig.colorbar(
        crx,
        cax=ax5,
        orientation='horizontal',
        ticks=colorticks,
        # fraction = 0.05,  # 色标大小,相对于原图的大小
        # pad=0.1,  #  色标和子图间距离
    )


# draw_img(img,imgr_low, imgr_band, imgr_high,lon,lat)

# %%

if __name__ == '__main__':
# Read image
# img = cv2.imread("imori.jpg").astype(np.float32)


    # flnm = '/mnt/zfm_18T/fengxiang/aa/2014_07_11-12/analysis/test2_20140711_1200-1500-cycle1h/data/exp1_noda/wrfout_d02_2014-07-11_15:00:00'
    # wrfnc = nc.Dataset(flnm)
    # w = wrf.getvar(wrfnc, 'wa')
    # z = wrf.getvar(wrfnc, 'z')
    # zz = w[10,:,:]
    # # zz = z[10,:,:]

    # lon = zz.XLONG.values
    # lat = zz.XLAT.values
    # img = zz.values

    
    flnm = '/mnt/zfm_18T/fengxiang/aa/2014_07_11-12/analysis/test3_20140711_1200-1500-cycle30min/data/interplevel_vars_multitime_multimodel_old.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(level=600).sel(time='2014-07-11 1530').sel(model='radar_surface')
    w = ds['wa']
    lon = w.lon.values
    lat = w.lat.values
    img = w.values
    


    
    
    

    # imgr_high = filter_lambda_high(img, dx=3, lam=50)
    imgr_high = filter_lambda_high(img, dx=3, lam=50)
    imgr_low = filter_lambda_low(img, dx=3, lam=100)
    imgr_band = filter_lambda_band(img, dx=3, lam1=50, lam2=100)
    
    draw_img(img,imgr_low, imgr_band, imgr_high,lon,lat)

# %%
# flnm = '/mnt/zfm_18T/fengxiang/aa/2014_07_11-12/analysis/test3_20140711_1200-1500-cycle30min/data/interplevel_vars_multitime_multimodel_old.nc'
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/wrfout_d01_2021-07-19_12:00:00'
ds = xr.open_dataset(flnm)
# ds = ds.sel(level=600).sel(time='2014-07-11 1530').sel(model='radar_surface')
# w = ds['wa']
# lon = zz.XLONG.values
# lat = zz.XLAT.values
# img = w.values
