# %%
import xarray as xr
import numpy as np
import pandas as pd
import os
from baobao.caculate import uv2wind


class Bias():
    """计算偏差, 均方根误差这些误差值
    """
    def __init__(self, wind_diag='wind_speed'):
        self.wind_diag = wind_diag   # 判断是风向还是风速
        pass

    def get_ds(self, flnm_wrf= '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'):
        """获得观测和预测的数据

        Args:
            flnm_wrf (str, optional): [description]. Defaults to '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'.

        Returns:
            [type]: [description]
        """
        flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
        # flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'

        ds_obs = xr.open_dataset(flnm_obs)
        ds_obs = ds_obs.rename({'idx':'sta'})

        ds = xr.open_dataset(flnm_wrf)
        ds = ds.swap_dims({'time':'Time'})
        ds = ds.drop_vars('time')
        ds_wrf = ds.rename({'Time':'time'})

        ## 筛选相同时间的值
        t = (ds_obs+ds_wrf).time
        ds_obs = ds_obs.sel(time=t)
        ds_wrf = ds_wrf.sel(time=t)
        return ds_wrf, ds_obs


    def get_bias(self, flnm_wrf = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'):

        ds_wrf, ds_obs = self.get_ds(flnm_wrf)
        u2 = ds_wrf['uvmet10'].sel(u_v='u')
        v2 = ds_wrf['uvmet10'].sel(u_v='v')
        ws2,wa2 = uv2wind(u2,v2)  # wind_speed, wind_angle
        ws1 = ds_obs['wind_speed']
        wa1 = ds_obs['wind_angle']


        ## 所有时次的平均值(对时间做平均, 只保留一个值)
        # bias_ws = (ws2-ws1).mean(dim=['time', 'sta'])
        bias_ws = (ws2-ws1).mean(dim='time').mean(dim='sta')
        mae_ws = np.absolute((ws2-ws1)).mean(dim='time').mean(dim='sta')
        rmse_ws = np.sqrt(((ws2-ws1)**2).mean(dim='time')).mean(dim='sta')
        # bias_wa = (wa2-wa1).mean(dim=['time', 'sta'])
        bias_wa = (wa2-wa1).mean(dim='time').mean(dim='sta')
        mae_wa = np.absolute((wa2-wa1)).mean(dim='time').mean(dim='sta')
        rmse_wa = np.sqrt(((wa2-wa1)**2).mean(dim='sta')).mean(dim='time')
        ## 风素的, 二选一
        if self.wind_diag == 'wind_speed':
            dic = {
                'BIAS':float('%.1f'%bias_ws),
                'RMSE':float('%.1f'%rmse_ws),
                'MAE':float('%.1f'%mae_ws),
            }
        # elif 
        elif self.wind_diag == 'wind_direction':
        ## 风向的
            dic = {
                'BIAS':float('%.1f'%bias_wa),
                'RMSE':float('%.1f'%rmse_wa),
                'MAE':float('%.1f'%mae_wa),
            }
        else:
            print("输入的判断风向还是风速的变量有问题")
        bb = pd.Series(dic)
        bb.index.name = 'model'
        cc = xr.DataArray(bb)
        return cc

    def get_bias_dual(self, ):
        """获得多个模式的偏差值

        Returns:
            [type]: [description]
        """
        model_list = ['gwd0', 'gwd1', 'gwd3']
        fl_path = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
        # fname = model'
        bias_list = []
        for model in model_list:
            fname = model+'/10m_wind_station.nc'
            flnm = os.path.join(fl_path,fname)
            aa = self.get_bias(flnm)
            print(aa)
            aa.rename(model)
            bias_list.append(aa)
        cc = xr.concat(bias_list, pd.Index(model_list, name='modle'))

        ## 处理成DataFrame
        dd = cc.to_dataset(dim='modle')
        d2 = dd.to_dataframe()
        ddf = d2.T
        df = ddf
        return df

    ## 计算偏差的改进率
    def caculate_delta(self, df):
        """计算偏差的改进率

        Args:
            df ([type]): [description]
        """
        def reshape_df(df2,df3):
            """改变原来的数组形状

            Args:
                df2 ([type]): [description]
                df3 ([type]): [description]

            Returns:
                [type]: [description]
            """
            df2['dd1'] = df3.iloc[:,0].values
            df2['dd2'] = df3.iloc[:,1].values
            df2['dd3'] = df3.iloc[:,2].values
            order = ['BIAS', 'dd1','RMSE', 'dd2','MAE', 'dd3']
            df3 = df2[order]
            df3['dd1'] = df3['dd1'].apply(lambda x: format(x, '.2%'))
            df3['dd2'] = df3['dd2'].apply(lambda x: format(x, '.2%'))
            df3['dd3'] = df3['dd3'].apply(lambda x: format(x, '.2%'))
            return df3
        a = (df.loc['gwd1']-df.loc['gwd0'])#/ddf.loc['gwd0']
        b = (df.loc['gwd3']-df.loc['gwd0'])#/ddf.loc['gwd0']
        c = df.loc['gwd0']
        delta1 = -1*a/c
        delta3 = -1*b/c
        df.loc['delta0'] = 0
        df.loc['delta1'] = delta1
        df.loc['delta3'] = delta3
        df1 = df.round(4)

        df2 = df1.iloc[0:3, 0:3]
        df3 = df1.iloc[3:, 0:3]
        df4 = reshape_df(df2,df3)
        return df4

# bias = Bias('wind_direction')
bias = Bias('wind_speed')
df = bias.get_bias_dual()
ddf = bias.caculate_delta(df)
ddf
# %%