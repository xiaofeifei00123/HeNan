# %%
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from draw_rain_distribution import Draw as drs
# %%




def set_colorbar(cf, ax, fig):    
    ## 设置colorbar
    levels = cf.levels
    colorticks = levels[1:-1]
    cb = fig.colorbar(
        cf,
        # cax=ax6,
        orientation='horizontal',
        ticks=colorticks,
        fraction = 0.06,  # 色标大小,相对于原图的大小
        pad=0.1,  #  色标和子图间距离
        )
    ax.tick_params(labelsize=10)  # 设置色标标注的大小
    tic = cb.get_ticks()
    labels = list(map(lambda x: str(x) if x<1 else str(int(x)), tic))  # 将colorbar的标签变为字符串
    cb.set_ticklabels(labels)

def draw_obs(tl):

    flnm_obs = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_obs.nc'
    ds_obs = xr.open_dataset(flnm_obs)
    # ds_obs = ds_obs.sel(time=slice('2021-07-20 01', '2021-07-21 00'))
    ds_obs = ds_obs.sel(time=tl)
    da_obs = ds_obs.sum(dim='time')['PRCP']
    

    cm = round(1/2.54, 2)
    proj = ccrs.PlateCarree()  # 创建坐标系
    fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
    ax = fig.add_axes([0.13,0.1,0.82,0.8], projection=proj)
    # dr = Draw(fig, ax)
    dr = drs()
    
    cf = dr.draw_rain_24h(da_obs, ax)    

    ## 设置colorbar
    set_colorbar(cf, ax, fig)

    ax.set_title('OBS', fontsize=10,loc='left')
    t1 = str(ds_obs.time.dt.strftime('%d/%H')[0].values)
    t2 = str(ds_obs.time.dt.strftime('%d/%H')[-1].values)
    tfig = str(ds_obs.time.dt.strftime('%d%H')[-1].values)
    tt = t1+'-'+t2
    ax.set_title(tt, fontsize=10,loc='center')
    
    fig_name = 'obs'+tfig+'aa'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_new/'
    fig.savefig(fig_path+fig_name)



def draw_model(tl):
    # flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model.nc'
    flnm_model = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/data/'+'rain_model_da.nc'
    ds_model = xr.open_dataset(flnm_model)
    # for model in ['gwd3','gwd1', 'gwd0']:
    # model_list = ['SS2']
    # model_list = ['CTRL','FD','GWD3','SS']
    # model_list = ['CTRL', 'GWD3']
    model_list = ['GWD3']
    for model in model_list:

        ds = ds_model.sel(time=tl)
        # da = ds[model].sum(dim='time')
        db = ds_model[model].isel(time=0)
        da = ds_model[model].sum(dim='time').assign_coords(db.coords)
        # print(da)

        cm = 1/2.54
        proj = ccrs.PlateCarree()  # 创建坐标系
        fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
        ax = fig.add_axes([0.13,0.1,0.82,0.8], projection=proj)
        dr = drs()


        cf = dr.draw_rain_24h(da, ax)    
        set_colorbar(cf, ax, fig)

        ax.set_title(model, fontsize=10,loc='left')
        t1 = str(ds.time.dt.strftime('%d/%H')[0].values)
        t2 = str(ds.time.dt.strftime('%d/%H')[-1].values)
        tt = t1+'-'+t2
        tfig = str(ds.time.dt.strftime('%d%H')[-1].values)
        ax.set_title(tt, fontsize=10,loc='center')
        
        fig_name = model+tfig
        fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_rain/rain_da/'
        fig.savefig(fig_path+fig_name)


if __name__ == '__main__':
    
    time1=slice('2021-07-17 01', '2021-07-18 00')
    time2=slice('2021-07-18 01', '2021-07-19 00')
    time3=slice('2021-07-19 01', '2021-07-20 00')
    time4=slice('2021-07-20 01', '2021-07-21 00')
    time5=slice('2021-07-21 01', '2021-07-22 00')
    # time6=slice('2021-07-22 01', '2021-07-23 00')
    time7=slice('2021-07-17 01', '2021-07-23 00')
    # time8=slice('2021-07-19 13', '2021-07-20 00')
    # time_list = [time1, time2, time3, time4, time5, time6]
    # time_list = [time4,time7]
    time_list = [time7]
    # for tl in time_list[0:1]:
    for tl in time_list:
        # draw_model(tl)
        draw_obs(tl)
    # %%
