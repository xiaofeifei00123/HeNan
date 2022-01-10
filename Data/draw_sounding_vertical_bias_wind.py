#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算bias偏差和均方根误差
不同时次和不同站点就是不同的试验罢了
-----------------------------------------
Time             :2021/12/17 10:40:47
Author           :Forxd
Version          :1.0
'''
# %%
import xarray as xr
import pandas as pd
import numpy as np
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt
import os

# import meteva as meb
# import meteva.base as meb
# meb.xarray_to_griddata()


# %%
class Data():
    
    def __init__(self, station='zhengzhou', data_time='2021-07-20 00') -> None:
        flnm1 = '/home/fengxiang/HeNan/Data/GWD/d03/sounding_all.nc'
        flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_all.nc'
        self.ds1 = xr.open_dataset(flnm1)
        self.ds2 = xr.open_dataset(flnm2)
        self.model_list = ['OBS', 'gwd0', 'gwd1', 'gwd3']
        self.color_list = ['black',  'blue','green', 'red']
        self.station = station
        self.time = data_time

    def get_data(self, var='u', station='zhengzhou', select_time='2021-07-20 00'):
        """筛选变量、站点和时间

        Args:
            var (str, optional): [description]. Defaults to 'u'.
            station (str, optional): [description]. Defaults to 'zhengzhou'.
            select_time (str, optional): [description]. Defaults to '2021-07-20 00'.

        Returns:
            [type]: [description]
        """
        pass
        select_time = self.time
        var_model = self.ds1[var].sel(station=station).sel(time=select_time).dropna(dim='pressure', how='all')
        var_obs = self.ds2[var].sel(station=station).sel(time=select_time).dropna(dim='pressure', how='all')
        var_obs = var_obs.expand_dims(dim='model')
        var_obs = var_obs.assign_coords({'model':['OBS']})
        vvv = xr.concat([var_model, var_obs], dim='model')
        if var=='q':
            vvv = vvv*10**3
        return vvv
        # return var_obs

    def get_data_one_var(self, var='wind_speed'):
        aa = self.get_data(var, station=self.station)
        data_dic = {}
        for model in self.model_list:
            data_dic[model] = aa.sel(model=model).dropna(dim='pressure') # 对不同模式的数据进行筛选
        return data_dic
        

    def get_var(self,):
        var_dic = {}
        # var_list = ['temp', 'theta_v','rh', 'q']
        var_list = ['u', 'v', 'wind_speed', 'wind_angle']
        for var in var_list:
            var_dic[var] = self.get_data_one_var(var)
        return var_dic


# %%

def interp_onetime(var_dic, model='OBS'):
    # var_list = ['temp', 'theta_v','rh', 'q']
    var_list = ['u', 'v', 'wind_speed', 'wind_angle']
    ds1 = xr.Dataset()
    for var in var_list:
        da_obs = var_dic[var][model].interp(pressure=np.arange(990, 190, -10)).round(1)
        ds1[var] = da_obs
    return ds1

station_list = ['zhengzhou', 'nanyang','lushi']
def get_bias_data():
    tt = pd.date_range('2021-07-20 00', '2021-07-20 12', freq='6H')
    model_list = ['OBS', 'gwd0', 'gwd1', 'gwd3']
    ds2_list = []
    for t in tt:

        ds_sta_list = []
        for sta in station_list:
            gd = Data(data_time=t, station=sta)
            var_dic = gd.get_var()
            ds1_list = []
            for model in model_list:
                ds1 = interp_onetime(var_dic,model)
                ds1_list.append(ds1)
            dds_sta = xr.concat(ds1_list, dim='model')
            ds_sta_list.append(dds_sta)

        dds =  xr.concat(ds_sta_list,dim='station')
        ds2_list.append(dds)
    ds3 = xr.concat(ds2_list, dim='time')
    return ds3
ds3 = get_bias_data()
# %%
# ds3
## 计算偏差

model_list_new = ['gwd0', 'gwd1', 'gwd3']

ds_bias_list = []
for model in model_list_new:
    ds_b = (ds3.sel(model=model)-ds3.sel(model='OBS')).mean(dim=['time','station'])
    cc = ds_b.expand_dims(dim='model')
    ds_c = cc.assign_coords({'model':[model]}) 
    ds_bias_list.append(ds_c)
ds_bias = xr.concat(ds_bias_list, dim='model')

## 计算均方根误差
ds_rmse_list = []
for model in model_list_new:
    ds_r = xr.ufuncs.sqrt((((ds3.sel(model=model)-ds3.sel(model='OBS'))**2).sum(dim=['time', 'station'])/6))
    cc = ds_r.expand_dims(dim='model')
    ds_c = cc.assign_coords({'model':[model]}) 
    ds_rmse_list.append(ds_c)
ds_rmse = xr.concat(ds_rmse_list, dim='model')
# ds_rmse

# %%
# ds_rmse

    



# %%
# ds_bias_gwd0 = (ds3.sel(model='gwd0')-ds3.sel(model='OBS')).mean(dim='time')
# ds_bias_gwd1 = (ds3.sel(model='gwd1')-ds3.sel(model='OBS')).mean(dim='time')
# ds_bias_gwd3 = (ds3.sel(model='gwd3')-ds3.sel(model='OBS')).mean(dim='time')
# ds_rmse_gwd0 = xr.ufuncs.sqrt((((ds3.sel(model='gwd0')-ds3.sel(model='OBS'))**2).sum(dim='time')/3))
# ds_rmse_gwd1 = xr.ufuncs.sqrt((((ds3.sel(model='gwd1')-ds3.sel(model='OBS'))**2).sum(dim='time')/3))
# ds_rmse_gwd3 = xr.ufuncs.sqrt((((ds3.sel(model='gwd3')-ds3.sel(model='OBS'))**2).sum(dim='time')/3))

# ds_rmse_gwd3






# %%
def draw(var_dic, station='zhengzhou'):

    # model_list = ['OBS', 'gwd0', 'gwd1', 'gwd3']
    model_list = ['gwd0', 'gwd1', 'gwd3']
    color_list = ['blue','green', 'red']
    var_list = ['u', 'v', 'wind_speed', 'wind_angle']
    # var_list = ['temp', 'theta_v', 'rh', 'q']
    

    fig = plt.figure(figsize=[12,8])
    axes = [None]*4

    axes[0] = fig.add_subplot(141)
    axes[1] = fig.add_subplot(142, sharey=axes[0], )
    axes[2] = fig.add_subplot(143, sharey=axes[0])
    axes[3] = fig.add_subplot(144, sharey=axes[0])
    fig.subplots_adjust(hspace=0.4, wspace=0.3)

    axes[1].tick_params('y', labelleft=False)
    axes[2].tick_params('y', labelleft=False)
    axes[3].tick_params('y', labelleft=False)




    for model,color in zip(model_list,color_list):
        # axes[0].plot(var_dic['temp'].sel(model=model).values, var_dic['temp'].sel(model=model).pressure, color=color, label=model)
        # axes[1].plot(var_dic['theta_v'].sel(model=model).values, var_dic['theta_v'].sel(model=model).pressure, color=color, label=model)
        # axes[2].plot(var_dic['rh'].sel(model=model).values, var_dic['rh'].sel(model=model).pressure, color=color, label=model)
        # axes[3].plot(var_dic['q'].sel(model=model).values, var_dic['q'].sel(model=model).pressure, color=color, label=model)
        axes[0].plot(var_dic['u'].sel(model=model).values, var_dic['u'].sel(model=model).pressure, color=color, label=model)
        axes[1].plot(var_dic['v'].sel(model=model).values, var_dic['v'].sel(model=model).pressure, color=color, label=model)
        axes[2].plot(var_dic['wind_speed'].sel(model=model).values, var_dic['wind_speed'].sel(model=model).pressure, color=color, label=model)
        axes[3].plot(var_dic['wind_angle'].sel(model=model).values, var_dic['wind_angle'].sel(model=model).pressure, color=color, label=model)

        
    axes[0].invert_yaxis()
    axes[0].legend(loc='upper center', bbox_to_anchor=(2.4,1.0,0.5,0.15), ncol=4, edgecolor='white', fontsize=18)
    # axes[0].legend(loc='best')
    axes[0].set_ylim(1000,200)
    fts = 12
    axes[0].set_ylabel('pressure (hPa)', fontsize=fts*1.5)

    # axes[0].set_xlim(-5,5)
    # axes[1].set_xlim(-10,10)
    # axes[3].set_xlim(0,20)

    def set_ticks(ax, var):
        pass
        fts = 12
        if var=='wind_angle':
            var = 'wind_direction'
        ax.set_xlabel(var, fontsize=fts*1.8)
        ax.xaxis.set_tick_params(labelsize=fts*1.3)
        ax.yaxis.set_tick_params(labelsize=fts*1.3)
        ax.tick_params(which='major',length=8,width=1.0) # 控制标签大小 
        ax.tick_params(which='minor',length=4,width=0.5)  #,colors='b')
        if station=='bias':
            ax.axvline(x=0, color='black')
        return ax

    i = 0
    for var in var_list:
        set_ticks(axes[i], var)
        i+=1
    # # ## 添加高度坐标
    # da = var_dic['temp'].sel(model='gwd0')
    # ax5 = axes[3].twinx()
    # ax5.plot(da.values, da.pressure, alpha=0)
    # ax5.invert_yaxis()
    # dda = da.swap_dims({'pressure':'height'})
    # ddda = dda.interp(height=[500, 1000, 1500, 2000, 3000, 5000, 10000, 20000], method='linear', kwargs={'fill_value':'extrapolate'})
    # ax5.set_yticks(ddda.pressure.values)
    # ax5.set_yticklabels(ddda.height.values.round(1))
    # ax5.set_ylabel('Height (m)', fontsize=fts*1.5)
    # set_ticks(ax5, 'wind')
    # fig_name= 'nanyang_wind'
    fig_name=station+"_wind"
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/sounding_upar/'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save)

if __name__ == '__main__':

    station_list = ['zhengzhou', 'nanyang']    

    # for sta in station_list:
        # gd = Data(station=sta) 
        # var_dic = gd.get_var()
        # print(var_dic['q']['gwd3'].max())
    draw(ds_rmse, station='rmse')
    draw(ds_bias, station='bias')
    


