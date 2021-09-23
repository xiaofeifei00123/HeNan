# %%
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim)
# %%

flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/ERA5/YSU_1800/wrfout_d03_2021-07-18_00:00:00'
ncfile = Dataset(flnm)
z = getvar(ncfile, 'z')
u = getvar(ncfile, 'ua')

# %%
get_cartopy(u)