import draw_10m_wind as dw
import draw_10m_wind_differ as dwd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os



cm = round(1/2.54, 2)
proj = ccrs.PlateCarree()  # 创建坐标系
fig = plt.figure(figsize=(17*cm, 27*cm), dpi=600)
grid = plt.GridSpec(4,
                    2,
                    figure=fig,
                    left=0.05,
                    right=0.98,
                    bottom=0.1,
                    top=0.97,
                    wspace=0.1,
                    # hspace=0.25)
                    hspace=0.2)
num = 8
axes = [None] * num  # 设置一个维度为8的空列表
for i in range(num):
    axes[i] = fig.add_subplot(grid[i], projection=proj)

    
## 画风场的分布图
t = '2021-07-20 00'
flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
gd = dw.GetData(time=t)
## 观测风场
u,v = gd.get_wind_obs(flnm_obs)
pic_dic = {'model':'obs', 'time':t}
dr = dw.DrawWind(fig, axes[0])
dr.ax.set_title('(a)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
dr.draw(u,v, pic_dic)


## gwd0风场
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'
u,v = gd.get_wind_wrf(flnm)
# pic_dic = {'model':'gwd0', 'time':t}
pic_dic = {'model':'gwd0'}
dr = dw.DrawWind(fig, axes[2])
dr.ax.set_title('(c)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
dr.draw(u,v, pic_dic)

## gwd3风场
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc'
u,v = gd.get_wind_wrf(flnm)
pic_dic = {'model':'gwd3'}
dr = dw.DrawWind(fig, axes[4])
dr.draw(u,v, pic_dic)
dr.ax.set_title('(e)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))

## gwd3-gwd0风场
u, v = dwd.get_wind_minus(t)
pic_dic = {'model':'gwd3-gwd0'}
# dr = dw.DrawWind(fig, axes[6], scale=10, Ulength=1)
dr = dw.DrawWind(fig, axes[6])
dr.draw(u,v,pic_dic)
dr.ax.set_title('(g)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))



t = '2021-07-20 12'
flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/Data/OBS/10m_wind_station.nc'
gd = dw.GetData(time=t)
## 观测风场
u,v = gd.get_wind_obs(flnm_obs)
pic_dic = {'model':'obs', 'time':t}
dr = dw.DrawWind(fig, axes[1])
dr.draw(u,v, pic_dic)
dr.ax.set_title('(b)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))


## gwd0风场
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/10m_wind_station.nc'
u,v = gd.get_wind_wrf(flnm)
# pic_dic = {'model':'gwd0', 'time':t}
pic_dic = {'model':'gwd0'}
dr = dw.DrawWind(fig, axes[3])
dr.draw(u,v, pic_dic)
dr.ax.set_title('(d)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))

## gwd3风场
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/10m_wind_station.nc'
u,v = gd.get_wind_wrf(flnm)
pic_dic = {'model':'gwd3'}
dr = dw.DrawWind(fig, axes[5])
dr.draw(u,v, pic_dic)
dr.ax.set_title('(f)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))

## gwd3-gwd0风场
u, v = dwd.get_wind_minus(t)
pic_dic = {'model':'gwd3-gwd0'}
# dr = dw.DrawWind(fig, axes[6], scale=10, Ulength=1)
dr = dw.DrawWind(fig, axes[7])
dr.draw(u,v,pic_dic)
dr.ax.set_title('(h)', loc='left', y=0.845,x=0.04 ,bbox=dict(facecolor='white', alpha=1, edgecolor='white'))



fig_name = '10mwind'
fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_lunwen/'
fig_save = os.path.join(fig_path, fig_name)
fig.savefig(fig_save, pad_inches=0)
