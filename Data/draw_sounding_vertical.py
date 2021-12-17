# %%
import xarray as xr
import wrf
import netCDF4 as nc
import matplotlib.pyplot as plt
import os



# %%
class Data():
    
    def __init__(self, station='zhengzhou', data_time='2021-07-20 00') -> None:
        flnm1 = '/home/fengxiang/HeNan/Data/GWD/sounding_all.nc'
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
        # var_list = ['theta_v', 'temp','q', 'rh']
        var_list = ['u', 'v', 'wind_speed', 'wind_angle']
        for var in var_list:
            var_dic[var] = self.get_data_one_var(var)
        return var_dic




# data_dic = get_data_one_var('wind_speed')


# %%
def draw(var_dic, station='zhengzhou'):

    model_list = ['OBS', 'gwd0', 'gwd1', 'gwd3']
    color_list = ['black',  'blue','green', 'red']
    var_list = ['u', 'v', 'wind_speed', 'wind_angle']
    

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
        axes[0].plot(var_dic['u'][model].values, var_dic['u'][model].pressure, color=color, label=model)
        axes[1].plot(var_dic['v'][model].values, var_dic['v'][model].pressure, color=color, label=model)
        axes[2].plot(var_dic['wind_speed'][model].values, var_dic['wind_speed'][model].pressure, color=color, label=model)
        axes[3].scatter(var_dic['wind_angle'][model].values, var_dic['wind_angle'][model].pressure, color=color, label=model)

        
    axes[0].invert_yaxis()
    axes[0].legend(loc='upper center', bbox_to_anchor=(2.4,1.0,0.5,0.15), ncol=4, edgecolor='white', fontsize=18)
    # axes[0].legend(loc='best')
    axes[0].set_ylim(1000,200)
    fts = 12
    axes[0].set_ylabel('pressure (hPa)', fontsize=fts*1.5)

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
        return ax

    i = 0
    for var in var_list:
        set_ticks(axes[i], var)
        i+=1
    # # ## 添加高度坐标
    da = var_dic['wind_angle']['gwd0']
    ax5 = axes[3].twinx()
    ax5.plot(da.values, da.pressure, alpha=0)
    ax5.invert_yaxis()
    dda = da.swap_dims({'pressure':'height'})
    ddda = dda.interp(height=[100,500, 1000, 1500, 2000, 3000, 5000, 10000, 20000], method='linear', kwargs={'fill_value':'extrapolate'})
    ax5.set_yticks(ddda.pressure.values)
    ax5.set_yticklabels(ddda.height.values.round(1))
    ax5.set_ylabel('Height (m)', fontsize=fts*1.5)
    set_ticks(ax5, 'wind')
    # fig_name= 'nanyang_wind'
    fig_name=station+'_wind'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_upar/'
    fig_save = os.path.join(fig_path, fig_name)
    fig.savefig(fig_save)

if __name__ == '__main__':
    # main()

    station_list = ['zhengzhou', 'nanyang']    

    for sta in station_list:
        gd = Data(station=sta) 
        var_dic = gd.get_var()
        draw(var_dic, station=sta)
    


