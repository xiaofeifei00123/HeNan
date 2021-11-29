# %%
import xarray as xr
import os
# %%
# path = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m_gwd3/wrfout_d01*'
path = '/mnt/zfm_18T/fengxiang/HeNan/Data/1912_90m_OGWD/wrfout_d04*'
fl_list = os.popen('ls {}'.format(path)) # 打开一个管道
fl_list = fl_list.read().split()
fl_list

# %%
def read_dusfcg(flnm):
    ds = xr.open_dataset(flnm)
    da = ds['DUSFCG']
    # da = ds['OA1LS']
    print(da.max().values)
for fl in fl_list:
    read_dusfcg(fl)

# %%
