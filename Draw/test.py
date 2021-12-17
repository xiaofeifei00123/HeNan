# %%
import xarray as xr
import os

# %%
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3/wrfout_d03_2021-07-20_05:00:00'
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd1/wrfout_d03_2021-07-20_06:00:00'
ds = xr.open_dataset(flnm)
# ds['DUSFCG'].plot()


flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/gwd3/wrfout_d03_2021-07-20_06:00:00'
ds = xr.open_dataset(flnm)
da = ds['DVSFCG_LS']+ds['DVSFCG_SS']+ds['DVSFCG_FD']+ds['DVSFCG_BL']
# da = ds['DUSFCG_LS']+ds['DUSFCG_SS']+ds['DUSFCG_FD']+ds['DUSFCG_BL']
da.min()
# da
# ds['DUSFCG_LS']

# ds1['DVSFCG_FD']+.mean()
# var_list = [
#         'DTAUX3D_LS', 
#         'DTAUX3D_SS',
#         'DTAUX3D_BL',
#         'DTAUX3D_FD',
#         'DTAUY3D_LS', 
#         'DTAUY3D_SS',
#         'DTAUY3D_BL',
#         'DTAUY3D_FD',
#         'DUSFCG_LS', 
#         'DUSFCG_SS',
#         'DUSFCG_BL',
#         'DUSFCG_FD',
#         'DVSFCG_LS', 
#         'DVSFCG_SS',
#         'DVSFCG_BL',
#         'DVSFCG_FD',
#         ]

# for var in var_list:
#     da = ds[var]
#     print(da.max().values)



