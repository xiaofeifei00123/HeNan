# %%
# import metdig
# import datetime
# from metdig.io.cassandra import get_obs_stations
import pandas as pd
import numpy as np
import xesmf as xe
from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units
import xarray as xr
import datetime
import meteva.base as meb


import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import xarray as xr
import meteva.base as meb
import numpy as np
# import os
import pandas as pd
from nmc_met_io.read_micaps import read_micaps_1, read_micaps_2, read_micaps_14
import meteva.base as meb
from nmc_met_graphics.plot import mapview
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmaps
from wrf import get_cartopy, smooth2d, getvar


# %%
def draw_contourf(ax, da):
    """在地图上绘制填色图
    """
    x = da.lon
    y = da.lat
    contour_levels = np.arange(0, 200, 10)
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


flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc'
ds = xr.open_dataset(flnm)
ds = ds.drop_vars('EC')
dds = ds.sel(time=slice('2021-07-20 0000', '2021-07-20 2000')).sum(dim='time')
# %%
da = dds['ERA51900']
fig = plt.figure(figsize=[10,8])
ax = fig.add_axes([0.12,0.03,0.8,0.97], projection=ccrs.PlateCarree())
mb = mapview.BaseMap()
# mb.set_extent('中国陆地')
mb.set_extent('河南')
mb.drawcoastlines(linewidths=0.8, alpha=0.5)

mb.drawstates(linewidths=0.8, alpha=0.5) # 省界
# mb.set_extent([107, 135, 20,40])
cs = draw_contourf(ax,da)
