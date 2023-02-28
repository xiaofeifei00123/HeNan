#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
区域平均的降水变化
-----------------------------------------
Time             :2022/08/23 20:07:20
Author           :Forxd
Version          :1.0
'''

# %%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from read_rain_wrf import GetData
from common import Common
from draw_rain_distribution_minus import Rain


import meteva.method as mem
import numpy as np


# %%


def replace_rain():
    """将19日18-20日00时的降水，换成1912起报的结果
    """
    import xarray as xr
    flnm_19_gwd3 = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/DA/GWD3/2021-07-19-12__2021-07-21-00/all.nc'
    flnm_19_ctrl = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/DA/GWD3/2021-07-19-12__2021-07-21-00/all.nc'
    flnm_all = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model_da.nc'
    flnm_new = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model_da_new.nc'
    da19_gwd3 = xr.open_dataarray(flnm_19_gwd3)
    da19_ctrl = xr.open_dataarray(flnm_19_ctrl)
    ds = xr.open_dataset(flnm_all)
    tt = pd.date_range('2021-07-19 18', '2021-07-20 00', freq='1H')
    ds['GWD3'].loc[tt, :, ::] = da19_gwd3.sel(time=tt)
    ds['CTRL'].loc[tt, :, ::] = da19_ctrl.sel(time=tt)
    ds.to_netcdf(flnm_new)





# %%
def caculate_area_mean_obs(da,area):
    mask = (
        (da.coords['lat']>area['lat1'])
        &(da.coords['lat']<area['lat2'])
        &(da.coords['lon']<area['lon2'])
        &(da.coords['lon']>area['lon1'])
    )
    aa = xr.where(mask, 1, np.nan)
    db = da*aa
    dsr = db.mean(dim=['lat', 'lon'])
    return dsr

def get_data():
    """降水核心区域平均降水随时间变化曲线

    Returns:
        _type_: _description_
    """
    flnm_model= '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model_da_new.nc'
    ds_model = xr.open_dataset(flnm_model)
    gd = GetData()
    com = Common()
    ds_model_mean = gd.caculate_area_mean(ds_model, com.areaD)
    ds_model_mean


    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_obs = xr.open_dataset(flnm_obs)
    da_obs = ds_obs['PRCP']
    da_obs_mean  = caculate_area_mean_obs(da_obs, com.areaE)
    da_obs_mean

    ds_model_mean['OBS'] = da_obs_mean
    ds_mean = ds_model_mean
    ds_mean
    return ds_mean

def draw(ds, fig, ax, *args, **kw):
    # color_list = ['black', 'green', 'blue', 'red', 'orange']
    # color_list = ['black', 'blue', 'red','green', 'orange']
    # color_list = ['black', 'red', 'blue','orange', 'green']
    color_list = ds.color_list
    linestyle_list = ds.line_list
    var_list = list(ds.data_vars)
    i = 0
    for var in var_list:
        da = ds[var]
        x = da.time.dt.strftime('%d/%H')
        y = da.values
        if var == 'PRCP':
            var = 'OBS'
        # ax.plot(x,y, label=var, color=color_list[i], **kw)
        ax.plot(x,y, label=var, color=color_list[i],linestyle=linestyle_list[i], **kw)
        i+=1
    ax.legend(edgecolor='white')

    # ax.set_xticks(x[::24])
    # ax.set_xticklabels(x[::24].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(6))

    ax.set_xticks(x[12::12])
    ax.set_xticklabels(x[12::12].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(6))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ax.set_xlim(47, 124)
    # ax.set_ylim(0, 20)
    # ax.set_xticks(x[::2])
    # ax.set_xticklabels(x[::2].values, rotation=0, fontsize=10)
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(1))

# %%
def draw_skill(ds):
    fo = ds[['GWD3', 'CTRL']].to_array().values
    fo_name_list = ["GWD3","CTRL"]
    ob = ds['OBS'].values
    cm = 1/2.54
    width = 8*cm
    height = 8*cm
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
    save_path = fig_path+'rain_time_skill'
    mem.taylor_diagram(ob,fo,member_list = fo_name_list, width=width, height=height,save_path=save_path)
# # %%

def main(ds):
    cm = 1/2.54
    fig = plt.figure(figsize=(16*cm, 8*cm), dpi=300)
    ax  = fig.add_axes([0.1, 0.15, 0.85, 0.8])
    ax.set_ylabel('Precipitation (mm)')
    ax.set_xlabel('Time (Date/Hour)')
    # ds.attrs['color_list'] = ['red', 'green', 'black', 'red', 'green', 'black']
    # ds.attrs['color_list'] = ['green', 'red', 'black', 'red', 'green', 'black']
    ds.attrs['color_list'] = ['red', 'black', 'red', 'green', 'black']
    ds.attrs['line_list'] = ['-', '-', '-', '--', '--', '--']
    draw(ds, fig, ax)
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
    fig.savefig(fig_path+'model_A')

if __name__ == "__main__":
    ds = get_data()
    main(ds[['GWD3', 'OBS']])
    # draw_skill(ds)

# # %%

# # import xarray as xr
# # # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/2021-07-19-12__2021-07-21-00/all.nc'
# # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/2021-07-20-12__2021-07-22-00/all.nc'
# # da = xr.open_dataarray(flnm)
# # da
# # # %%
# # gd = GetData()
# # cm = Common()
# # area = {
# #     'lat1':33.5,
# #     'lat2':36.0,
# #     'lon1':112.5,
# #     'lon2':114.5,
# #     }        
# # daa = gd.caculate_area_mean(da, area)
# # daa
# # # %%
# # daa.plot()
# # # da.mean(dim=['south'])
