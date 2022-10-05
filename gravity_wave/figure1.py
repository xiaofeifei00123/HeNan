# %%
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
# from draw_rain_distribution_24h import set_colorbar
import draw_rain_distribution_24h as drd
import draw_rain_time_obs  as drt
import matplotlib.colorbar as cb


# %%
cm = 1/2.54
proj = ccrs.PlateCarree()  # 创建坐标系
fig = plt.figure(figsize=(17*cm, 8*cm), dpi=300)
ax1 = fig.add_axes([0.05,0.1,0.3,0.82], projection=proj)
ax2  = fig.add_axes([0.5, 0.2, 0.4, 0.6])
ax_colorbar  = fig.add_axes([0.12, 0.05, 0.3, 0.05])

cf = drd.draw_obs(ax1, fig)
# drd.set_colorbar(cf, ax_colorbar, fig)
levels = cf.levels
colorticks = levels[1:-1]
cb = fig.colorbar(
    cf,
    cax=ax_colorbar,
    # ax,
    orientation='horizontal',
    ticks=colorticks,
    fraction = 0.06,  # 色标大小,相对于原图的大小
    pad=0.05,  #  色标和子图间距离
    )
ax_colorbar.tick_params(labelsize=10)  # 设置色标标注的大小
tic = cb.get_ticks()
labels = list(map(lambda x: str(x) if x<1 else str(int(x)), tic))  # 将colorbar的标签变为字符串
cb.set_ticklabels(labels)






drt.draw_obs(ax2)
fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture/'
fig.savefig(fig_path+'ttt.png')
# %%


# def draw_obs_main():
#     cm = 1/2.54
#     fig = plt.figure(figsize=(9*cm, 5*cm), dpi=300)
#     ax2  = fig.add_axes([0.12, 0.2, 0.73, 0.7])
#     draw_obs(ax2, fig)

#     fig.legend(loc='upper left', edgecolor='white', bbox_to_anchor=(0.13, 0.7, 0.2, 0.2))
#     figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_time/'
#     fig.savefig(figpath+'rain_time')



