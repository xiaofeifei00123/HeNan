#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
存放公共的绘图程序
-----------------------------------------
Time             :2021/10/12 09:19:46
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.crs as ccrs
import meteva.base as meb
from nmc_met_graphics.plot import mapview
from wrf import get_cartopy, smooth2d, getvar
import cmaps

def draw_contourf_quick(da):
    """对二维的含经纬度的数据，
    快速的绘制填色图
    判断数据正确性的

    Args:

        da ([type]): 二维数据，lat, lon

    Returns:
        [type]: [description]
    """
    pass
    def draw_contourf(ax, da):
        """在地图上绘制填色图
        """
        x = da.lon
        y = da.lat
        contour_levels = np.arange(0, 660, 20)
        # colormap = cmaps.ViBlGrWhYeOrRe
        colormap = cmaps.precip3_16lev  # 反转色标
        crx = ax.contourf(x,
                            y,
                            da,
                            cmap=colormap,
                        #   norm=norm,
                        #   extend='both',
                        #   extend='max',
                            levels=contour_levels,
                            transform=ccrs.PlateCarree())
        return crx

    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0.12,0.03,0.8,0.97], projection=ccrs.PlateCarree())
    mb = mapview.BaseMap()
    # mb.set_extent('中国陆地')
    mb.set_extent('河南')
    mb.drawcoastlines(linewidths=0.8, alpha=0.5)

    mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
    # mb.set_extent([107, 135, 20,40])
    cs = draw_contourf(ax,da)

if __name__ == '__main__':
    # main()
    pass
    