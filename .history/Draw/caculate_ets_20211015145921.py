#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:
计算TS评分，并绘制bar图
-----------------------------------------
Time             :2021/10/12 10:50:16
Author           :Forxd
Version          :1.0
'''
# %%
# import metdig
# import datetime
# from metdig.io.cassandra import get_obs_stations
import pandas as pd
import numpy as np
import xarray as xr
from draw_global import draw_contourf_quick
import meteva.method as mem
import meteva.base as meb


import matplotlib.pyplot as plt
from cycler import cycler
from global_variable import station_dic
import seaborn as sns
# plt.rcParams['font.sans-serif'] = 'STHupo'
# plt.rcParams['font.family'] = ['sans-serif']
# plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
# plt.rcParams['font.sans-serif'] = ['SimSun']  # 指定默认字体为黑体
# plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 指定默认字体为黑体
# plt.rcParams['font.sans-serif'] = ['Helvetica']  # 指定默认字体为黑体
# %%
# from matplotlib import font_manager
 
# for font in font_manager.fontManager.ttflist:
#     # 查看字体名以及对应的字体文件名
#     # print(font.name, '-', font.fname)
#     print(font.name,)



# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc'
ds = xr.open_dataset(flnm)
ds





# %%
def contourf():
    """验证数据的正确性
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.drop_vars('EC')
    dds = ds.sel(time=slice('2021-07-19 1200', '2021-07-20 1200')).sum(dim='time')
    da = dds['ERA51912']
    da
    draw_contourf_quick(da)

# %%

class Caculate():

    def __init__(self, rain,) -> None:
        # self.flag = flag
        # self.area = area
        self.rain = rain
        # self.threshold = threshold
        self.model_list = ['ERA51800','ERA51812','ERA51900','ERA51912','GDAS1800','GDAS1812','GDAS1900','GDAS1912',]

    def get_two_scale(self, threshold):
        # flag = 'all'
        """计算ETS评分这些值
            这个阈值是要变化的
        """
        rain = self.rain  # 空间分布的总降水
        # model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        model_list = self.model_list

        ### 计算二分类评分
        hfmc_dic = {}
        Accuracy = {}
        ETS = {}
        TS = {}
        POFD = {}
        HK = {}
        ODR = {}
        BIAS = {}
        for model in model_list:
            # hfmc = mem.hfmc(rain['obs'].values/30, rain[model].values/30, grade_list=[0.25])

            ## mask高原外的数据
            rain_obs = rain['OBS'].values
            rain_model = rain[model].values
            # print(rain_model.max())
            # print("计算 %s" %model)
            hfmc = mem.hfmc(rain_obs, rain_model, grade_list=[threshold])   # 这里算的都是平均值的评分, 我想要的是评分的平均值, 也就是要算每个时次的评分
            ETS[model] = mem.ets_hfmc(hfmc) 
            # TS[model] = mem.ts_hfmc(hfmc) 
            Accuracy[model] = mem.pc_hfmc(hfmc)
            # POFD[model] = mem.pofd_hfmc(hfmc)
            # POFD[model] = mem.far_hfmc(hfmc)
            # HK[model] = mem.hk_yesorno_hfmc(hfmc)
            # ODR[model] = mem.odds_ratio_hfmc(hfmc)
            BIAS[model] = mem.bias_hfmc(hfmc)

        # grade_list = [Accuracy, ETS, TS, POFD, HK, ODR, BIAS]
        # grade_list_name = ['Accuracy', 'ETS', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']
        # grade_list = [Accuracy, TS,  BIAS]
        grade_list = [Accuracy, ETS,  BIAS]
        # grade_list_name = ['Accuracy', 'TS', 'BIAS']
        grade_list_name = ['Accuracy', 'ETS', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        # print(df)
        # df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/tt.csv')
        return df



    def get_space_scale(self, rain_threshold):
        """计算空间评分
        """
        rain = self.rain
        """空间降水"""
        # model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']
        model_list = self.model_list

        sal_list = []
        for model in model_list:
            pass
            print("计算 %s" %model)
            # tttime = pd.date_range()
            rain_model = rain[model]
            rain_obs  = rain['OBS']
            grid_obs = meb.xarray_to_griddata(rain_obs)
            grid_model = meb.xarray_to_griddata(rain_model)

            look_ff = mem.mode.feature_finder(grid_obs, grid_model, smooth=10, threshold=rain_threshold/15, minsize=5)
            look_match = mem.mode.centmatch(look_ff)
            look_merge = mem.mode.merge_force(look_match)
            sal = mem.sal(look_ff)
            ssal = pd.Series(sal)
            sal_list.append(ssal)
        
        df = pd.concat(sal_list, axis=1)
        df.columns = model_list
        return df

    def get_time_scale(self, rain):
        """计算时间序列的评分
        rain是时间序列的降水, 就是一个一维的数组
        """
        # rd = rd(self.flag, self.area)
        # rain = rd.get_rain_total()
        # rain = rd.get_rain_times()
        # print(rain)
        model_list = ['ACM2', 'YSU', 'QNSE', 'QNSE_EDMF', 'TEMF']

        MAE = {} # 平均绝对误差
        RMSE = {} # 均方根误差
        SD = {}  # 预报标准差/观测标准差
        CORR = {}
        for model in model_list:
            pass
            rain_model = rain[model].values
            rain_obs  = rain['obs'].values
            # MAE[model] = mem.mae(rain_obs, rain_model)
            sd = mem.ob_fo_std(rain_obs, rain_model)
            SD[model] = sd[1]/sd[0]
            RMSE[model] = mem.rmse(rain_obs, rain_model)
            CORR[model] = mem.corr(rain_obs, rain_model)
        # grade_list = [MAE, RMSE, SD, CORR]
        # grade_list_name = ['MAE', 'RMSE', 'SD', 'CORR']
        # grade_list = [RMSE, SD, CORR]
        grade_list = [RMSE, SD]
        # grade_list_name = ['RMSE', 'SD[Model]/SD[OBS]']
        grade_list_name = ['RMSE', 'SD[Fcst]/SD[Obs]']
        # grade_list_name = ['RMSE', 'SD', 'CORR']
        # grade_list_name = ['Accuracy', 'TS', 'POFD', 'HK', 'ODR', 'BIAS']

        b = []
        for i in grade_list:
            a = pd.Series(i)
            b.append(a)
            # print(a)
        df = pd.concat(b, axis=1)
        df.columns = grade_list_name
        df = df.T
        print(df)
        # df.to_csv('/home/fengxiang/Project/Asses_PBL/Data/time_score.csv')
        return df


def draw_taylor(dds):
    """绘制泰勒图

    Args:
        dds ([type]): [description]
    """
    da = dds.to_array(dim='model')
    ob = dds['OBS'].values
    fo = da.drop_sel(model=['OBS', 'GFS'])
    fo = fo.values
    model_list = ['ERA51800','ERA51812','ERA51900','ERA51912','GDAS1800','GDAS1812','GDAS1900','GDAS1912',]
    grade_list=[50, ]
    fo_name_list = model_list
    mem.performance(ob,fo,grade_list = grade_list,member_list = fo_name_list,)




# %%
def get_ts():
    """获取不同阈值的ts评分
    """
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/rain_all.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.drop_vars('EC')
    dds = ds.sel(time=slice('2021-07-20 0000', '2021-07-20 1200')).sum(dim='time')
    ca = Caculate(dds)
    # threshold_list = [0.1, 5, 15, 30, 70, 140]
    # rain_list = ['小雨', '中雨', '大雨', '暴雨', '大暴雨', '特大暴雨']
    threshold_list = [15, 30, 70, 140]
    rain_list = ['大雨', '暴雨', '大暴雨', '特大暴雨']
    ts_list = []
    for threshold in threshold_list:
        df = ca.get_two_scale(threshold).T
        print(df)
        ddf = df['ETS']
        ddf.index.name = 'model'
        da = xr.DataArray.from_series(ddf)
        # da = get_ts(threshold)
        ts_list.append(da)
    dda = xr.concat(ts_list, pd.Index(threshold_list, name='threshold'))
    return dda
    # dda.sel(threshold=0.1)
# ts = get_ts()
# ts
# %%
##### 画图
def draw_bar(ts):
    # threshold_list = ts.threshold
    # labels = list(threshold_list.astype('str').values)
    # rain_list = ['小雨', '中雨', '大雨', '暴雨', '大暴雨', '特大暴雨']
    # threshold_list = [0.1, 5, 15, 30, 70, 140]
    # rain_list = ['小雨≥0.1', '中雨≥5', '大雨≥15', '暴雨≥30', '大暴雨≥70', '特大暴雨≥140']
    rain_list = ['大雨≥15', '暴雨≥30', '大暴雨≥70', '特大暴雨≥140']
    labels = rain_list
    x = np.arange(len(labels))
    ds = ts
    ds = ts.to_dataset(dim='model')

    model_list = ['ERA51800','ERA51812','ERA51900','ERA51912','GDAS1800','GDAS1812','GDAS1900','GDAS1912',]
    # color_list = ['green', 'blue','orange', 'red',  'green', 'blue','orange', 'red', ]
    # color_list = ['green', 'blue','orange', 'red',  'darkgreen', 'darkblue', 'darkorange', 'darkred',]
    color_list = ['#7cfc00', '#3366cc','#eee', '#ff0066',  '#32cd32', '#0000cc', '#ff4500', '#ff0000',]
    hatch_list = ['.','.','.','.','x','x','x','x',]
    width = 0.10
    fig = plt.figure(figsize=(12, 6), dpi=300)  # 创建页面
    ax = fig.add_axes([0.05,0.1, 0.93,0.8])
    i = -3
    j = 0
    for model in model_list:
        # rects = ax.bar(x+width*i, ds[model], width, label=model, facecolor='white', edgecolor=color_list[j], hatch=hatch_list[j])
        rects = ax.bar(x+width*i, ds[model], width, label=model, color=color_list[j])
        i += 1
        j += 1
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=18)
    ax.set_yticks(np.arange(0,1.1,0.1))
    ax.set_ylim(-0.0,0.5)
    # ax.set_yticklabels(fontsize=24)
    ax.tick_params(axis='both', labelsize=19, direction='out')
    # ax.legend(fontsize=20, edgecolor='white', loc='upper right', ncol=2)
    # ax.legend(fontsize=20, edgecolor='white', loc='upper right', ncol=1, bbox_to_anchor=(0.7, 0.2, 0.5, 0.5))
    ax.legend(fontsize=18, edgecolor='white', loc='upper right', ncol=1)
    ax.set_title('ETS', fontsize=24)
    fig.savefig('/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_rain/ets_score.png')
if __name__ == '__main__':
    pass
    ts = get_ts()
    draw_bar(ts)