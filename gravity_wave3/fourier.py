# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
# %%


# flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/SS/wrfout/cross4_1time.nc'
# flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/SS/wrfout/cross4_times.nc'
flnm = '/home/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/cross4_times.nc'
ds = xr.open_dataset(flnm)
da = ds['wa_cross']
# da = ds['ua_cross']
# da = ds['div_cross']
# da1 = da.interpolate_na(dim='vertical', method='linear',  fill_value="extrapolate")
da1 = da.interpolate_na(dim='vertical', method='nearest',  fill_value="extrapolate")
da1
# %%
ds
# %%
da2 = da1.sel(vertical=12000, method='nearest')
# da2 = da1.interp(vertical=12000, method='linear')
da3 = da2[6:, 1:]
# %%
# da3 = da2
da3 = da2[0:-6,0:-1]

# da3.plot()
# da3.mean()
da3.shape
# %%
db = da3
grat = db.values
cm = 1/2.54
# fig = plt.figure(figsize=(10*cm, 3*cm), dpi=300)
fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# ax1 = plt.subplot(1,3,1)
# ax2 = plt.subplot(1,3,2)
# ax3 = plt.subplot(1,3,3)
# ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])

# crx1 = ax1.contourf(db)
# fig.colorbar(crx1,orientation='horizontal')

ft = np.fft.ifftshift(grat)
ft = np.fft.fft2(ft)   # 傅里叶便换
ft = np.fft.fftshift(ft)

# ft = np.fft.fft2(grat)
# ft = ft[69:,69:]

x = grat.shape[0]
x2 = np.arange(0,x,1)-(x-1)/2
y2 = x2

# ax2.set_xlim(0, 0.1)
# ax2.set_ylim(0, 0.1)

crx2 = ax2.contourf(x2, y2, abs(ft)**2, levels=100, cmap='rainbow')
# crx2 = ax2.contourf(x2, y2, abs(ft))
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
# ax2.set_xlim(-2, 2)
# ax2.set_ylim(-2, 2)
# ax2.set_xlim(-10, 10)
# ax2.set_ylim(-10, 10)
# # ax2.set_ylim(0, 5)
# ax2.set_xticks(np.arange(-5, 6, 1))
# ax2.set_yticks(np.arange(-5, 6, 1))
# ax2.set_xlim(-5, 5)
# ax2.set_ylim(-5, 5)
# ax2.set_ylim(65, 75)
# %%
# ax2.set_xticks(np.arange(200, 221, 10))

# %%

# plt.show()
# fig.colorbar(crx2)
# fig.colorbar(crx2,orientation='horizontal')
# ax2.contourf(np.log(abs(ft)))
# ax2.imshow(np.log(abs(ft)))

## 滤波
def Ideal(src, d0, ftype):
    template = np.zeros(src.shape, dtype=np.float32)  # 构建滤波器
    r, c = src.shape
    for i in range(r):
        for j in range(c):
            distance = np.sqrt((i - r/2)**2 + (j - c/2)**2)
            if distance < d0:
                template[i, j] = 1
            else:
                template[i, j] = 0

    if ftype == 'high':
        template = 1 - template
    return template

# filter = Ideal(ft, 20,'high')  # 低通滤波
filter = Ideal(ft, 10,'low')  # 低通滤波
ft = ft*filter


ift = np.fft.ifftshift(ft)
ift = np.fft.ifft2(ift)  # 傅里叶逆变换
ift = np.fft.fftshift(ift)
ift = ift.real  # Take only the real part
crx3 = ax3.contourf(ift)
# fig.colorbar(crx3)
