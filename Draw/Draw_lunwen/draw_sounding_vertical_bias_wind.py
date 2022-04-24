#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算探空廓线的bias偏差和均方根误差
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
# plt.rcParams['axes.unicode_minus']=False 


# %%
class Data():
    
    def __init__(self, station='zhengzhou', data_time='2021-07-20 00'):
        flnm1 = '/home/fengxiang/HeNan/Data/GWD/d03/sounding_all.nc'
        flnm2 = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/micaps_sounding_station_all.nc'
        self.ds1 = xr.open_dataset(flnm1)
        self.ds2 = xr.open_dataset(flnm2)
        # self.model_list = ['OBS', 'gwd0', 'gwd1', 'gwd3']
        # self.color_list = ['black',  'blue','green', 'red']
        self.model_list = ['OBS', 'gwd0',  'gwd3']
        self.color_list = ['black',  'blue', 'red']
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
    # tt = pd.date_range('2021-07-20 00', '2021-07-20 12', freq='6H')
    tt = pd.DatetimeIndex([
        '2021-07-20 00',
        '2021-07-20 06',
        '2021-07-20 12',
        '2021-07-21 00',
        ])
    model_list = ['OBS', 'gwd0', 'gwd3']
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
ds3
# %%
# ds3
## 计算偏差

# model_list_new = ['gwd0', 'gwd1', 'gwd3']
model_list_new = ['gwd0',  'gwd3']

ds_bias_list = []
for model in model_list_new:
    # ds_b = (ds3.sel(model=model)-ds3.sel(model='OBS')).mean(dim=['time','station'])
    ds_b = (ds3.sel(model=model)-ds3.sel(model='OBS')).mean(dim='station').mean(dim='time')
    # ds_b = abs(ds3.sel(model=model)-ds3.sel(model='OBS')).mean(dim='station').mean(dim='time')
    cc = ds_b.expand_dims(dim='model')
    ds_c = cc.assign_coords({'model':[model]}) 
    ds_bias_list.append(ds_c)
ds_bias = xr.concat(ds_bias_list, dim='model')

## 计算均方根误差的时间平均值
ds_rmse_list = []
for model in model_list_new:
    # ds_r = xr.ufuncs.sqrt((((ds3.sel(model=model)-ds3.sel(model='OBS'))**2).sum(dim=['time', 'station'])/6))
    ds_r = (np.sqrt(((ds3.sel(model=model)-ds3.sel(model='OBS'))**2).mean(dim='station')).mean(dim='time'))
    cc = ds_r.expand_dims(dim='model')
    ds_c = cc.assign_coords({'model':[model]}) 
    ds_rmse_list.append(ds_c)
ds_rmse = xr.concat(ds_rmse_list, dim='model')
ds_rmse




def draw_big(var_dic, station='zhengzhou'):
    """这个只有风速和风向两个变量了

    Args:
        var_dic (_type_): _description_
        station (str, optional): _description_. Defaults to 'zhengzhou'.

    Returns:
        _type_: _description_
    """

    # model_list = ['gwd0', 'gwd1', 'gwd3']
    # color_list = ['black','blue', 'red']
    model_list = ['gwd0',  'gwd3']
    color_list = ['red', 'blue']
    # var_list = ['wind_speed', 'wind_angle']
    var_list = ['风速', '风向']
    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm,9*cm], dpi=600)
    axes = [None]*2

    # axes[0] = fig.add_axes([0.13,0.15, 0.4,0.8])
    # axes[1] = fig.add_axes([0.58,0.15, 0.4,0.8], sharey=axes[0])
    axes[0] = fig.add_axes([0.20,0.15, 0.36,0.78])
    axes[1] = fig.add_axes([0.6,0.15, 0.36,0.78], sharey=axes[0])
    # axes[0].set_ylabel('Pressure (hPa)', fontsize=10)

    axes[1].tick_params('y', labelleft=False)
    for model,color in zip(model_list,color_list):
        label = model
        if model == 'gwd0':
            label = 'CTRL'
        elif model == 'gwd3':
            label = 'GWD3'
        axes[0].plot(var_dic['wind_speed'].sel(model=model).values, var_dic['wind_speed'].sel(model=model).pressure, color=color, label=label, linewidth=1)
        axes[1].plot(var_dic['wind_angle'].sel(model=model).values, var_dic['wind_angle'].sel(model=model).pressure, color=color, label=label, linewidth=1)

        
    axes[0].invert_yaxis()
    # axes[0].legend(loc='upper center', bbox_to_anchor=(2.4,1.0,0.5,0.15), ncol=4, edgecolor='white', fontsize=18)
    # axes[0].legend(loc='upper center', bbox_to_anchor=(1.0, 1.0,0.5,0.1),ncol=3)
    # axes[0].legend(loc='upper right', fontsize=10, edgecolor='white')
    axes[0].legend(loc='upper left', fontsize=10, edgecolor='white')
    fts = 10
    # axes[0].set_ylabel('pressure (hPa)', fontsize=fts)

    axes[0].set_xlim(0,6)
    axes[1].set_xlim(20,70)
    axes[0].set_ylim(1000,150)
    axes[1].set_ylim(1000,150)

    def set_ticks(ax, var):
        pass
        fts = 10
        if var=='wind_angle':
            var = 'wind_direction'
        ax.set_xlabel(var, fontsize=fts)
        ax.xaxis.set_tick_params(labelsize=fts)
        ax.yaxis.set_tick_params(labelsize=fts)
        ax.tick_params(which='major',length=4,width=0.6) # 控制标签大小 
        ax.tick_params(which='minor',length=2,width=0.3)  #,colors='b')
        if station=='bias':
            ax.axvline(x=0, color='black')
        return ax

    i = 0
    for var in var_list:
        set_ticks(axes[i], var)
        i+=1
    fig_name=station+"_wind"
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/sounding_upar/'
    axes[0].set_ylabel('Pressure (hPa)', fontsize=12)
    axes[0].set_xlabel('RMSE (风速，m/s)', fontsize=10)
    # axes[0].set_title('风速', loc='right',fontsize=10)
    axes[1].set_xlabel('RMSE (风向，$^{\circ}$)', fontsize=10)
    # axes[1].set_title('风向', loc='right',fontsize=12)
    axes[0].set_title('(b)', loc='left',y=0.98, fontsize=10)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save)
    
def draw_bias(var_dic, station='zhengzhou'):
    """这个只有风速和风向两个变量了

    Args:
        var_dic (_type_): _description_
        station (str, optional): _description_. Defaults to 'zhengzhou'.

    Returns:
        _type_: _description_
    """

    # model_list = ['gwd0', 'gwd1', 'gwd3']
    # color_list = ['black','blue', 'red']
    model_list = ['gwd0',  'gwd3']
    color_list = ['red', 'blue']
    # var_list = ['wind_speed', 'wind_angle']
    var_list = ['风速', '风向']
    cm = 1/2.54
    fig = plt.figure(figsize=[8*cm,9*cm], dpi=600)
    axes = [None]*2

    # axes[0] = fig.add_axes([0.13,0.15, 0.4,0.8])
    # axes[1] = fig.add_axes([0.58,0.15, 0.4,0.8], sharey=axes[0])
    axes[0] = fig.add_axes([0.20,0.15, 0.36,0.78])
    axes[1] = fig.add_axes([0.6,0.15, 0.36,0.78], sharey=axes[0])
    # axes[0].set_ylabel('Pressure (hPa)', fontsize=10)

    axes[1].tick_params('y', labelleft=False)
    for model,color in zip(model_list,color_list):
        label = model
        if model == 'gwd0':
            label = 'CTRL'
        if model == 'gwd3':
            label = 'GWD3'
        axes[0].plot(var_dic['wind_speed'].sel(model=model).values, var_dic['wind_speed'].sel(model=model).pressure, color=color, label=label, linewidth=1)
        axes[1].plot(var_dic['wind_angle'].sel(model=model).values, var_dic['wind_angle'].sel(model=model).pressure, color=color, label=label, linewidth=1)

        
    axes[0].invert_yaxis()
    # axes[0].legend(loc='upper center', bbox_to_anchor=(2.4,1.0,0.5,0.15), ncol=4, edgecolor='white', fontsize=18)
    # axes[0].legend(loc='upper center', bbox_to_anchor=(1.0, 1.0,0.5,0.1),ncol=3)
    # axes[0].legend(loc='upper right', fontsize=10, edgecolor='white')
    axes[0].legend(loc='upper right', fontsize=10, edgecolor='white')
    fts = 10
    # axes[0].set_ylabel('pressure (hPa)', fontsize=fts)

    # axes[0].set_xlim(0,6)
    # axes[1].set_xlim(20,70)
    axes[0].set_xlim(-5,5)
    axes[1].set_xlim(-40,40)
    axes[0].set_ylim(1000,150)
    axes[1].set_ylim(1000,150)

    def set_ticks(ax, var):
        pass
        fts = 10
        if var=='wind_angle':
            var = 'wind_direction'
        ax.set_xlabel(var, fontsize=fts)
        ax.xaxis.set_tick_params(labelsize=fts)
        ax.yaxis.set_tick_params(labelsize=fts)
        ax.tick_params(which='major',length=4,width=0.6) # 控制标签大小 
        ax.tick_params(which='minor',length=2,width=0.3)  #,colors='b')
        if station=='bias':
            ax.axvline(x=0, color='black')
        return ax

    i = 0
    for var in var_list:
        set_ticks(axes[i], var)
        i+=1
    fig_name=station+"_wind"
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/sounding_upar/'
    axes[0].set_ylabel('Pressure (hPa)', fontsize=12)
    axes[0].set_xlabel('BIAS (风速，m/s)', fontsize=10)
    axes[1].set_xlabel('BIAS(风向，$^{\circ}$)', fontsize=10)
    axes[0].set_title('(a)', loc='left',y=0.98, fontsize=10)

    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save)
    
    

if __name__ == '__main__':

    draw_big(ds_rmse, station='rmse')
    draw_bias(ds_bias, station='bias')
    



# %%
