# %%
import netCDF4 as nc
from netCDF4 import Dataset
import wrf
from wrf import getvar
# from baobao.caculate import caculate_div3d
from baobao.caculate import get_div_wrfout

# %%

wrf_file = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/wrfout_d01_2021-07-19_18:00:00'


div = get_div_wrfout(wrf_file)
div