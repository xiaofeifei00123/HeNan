# %%
from src import read_time_height_wrf 
from src import draw_time_height
import xarray as xr
import matplotlib.pyplot as plt
# %%
# sd = read_time_height_wrf.Sounding()

# %%
# sd
# flnm = 
def save_data_one(path_wrfout='/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/'):
    # south
    areaA = {
        'lat1':33.5,
        'lat2':33.8,
        'lon1':112.3,
        'lon2':112.8,
        }        
    
    areaB = {
        'lat1':34,
        'lat2':35,
        'lon1':112.4,
        'lon2':113.4,
        }        

    areaC = {
        'lat1':35.3,
        'lat2':35.8,
        'lon1':113.4,
        'lon2':114.2,
        }        
    areaD = {
        'lat1':33.5,
        'lat2':36.0,
        'lon1':112.2,
        'lon2':114.8,
        }        
    areaE = {
        'lat1':32,
        'lat2':36.5,
        'lon1':110.5,
        'lon2':116,
        }        
    areaF = {
        'lat1':33.5,
        'lat2':34.5,
        'lon1':112.0,
        'lon2':113,
        }        
    areaG = {
        'lat1':33.8,
        'lat2':34.4,
        'lon1':113.2,
        'lon2':114.2,
        }        
    ### 背风坡
    areaH = {
        'lat1':33.4,
        'lat2':34.2,
        'lon1':111.4,
        'lon2':112.2,
    }        
    ### 迎风坡
    areaI = {
        'lat1':33.4,
        'lat2':34.2,
        'lon1':112.6,
        'lon2':113.4,
    }        

    # area = areaA
    # area_list = [areaA, areaB,areaC] 
    # area_list = [areaD, areaE]
    # area_list = [areaB,areaA]
    area_list = [areaH, areaI]
    # arname_list = ['A']

    # arlist = ['A', 'B', 'C']
    # arlist = ['D', 'E']
    # arlist = ['D', 'C']
    # arlist = ['B', 'A']
    arlist = ['H', 'I']

    i = 0
    for area in area_list:
        # area = 'area'+ar

        sd = read_time_height_wrf.Sounding(area=area)
        ds = sd.sounding_main(path_wrfout)
        # path_save = '/mnt/zfm_18T/fengxiang/HeNan/draw_gravity/data/'
        # flnm_save = path_wrfout+'time_cross_south.nc'
        # ar = area[-1]
        flnm_save = path_wrfout+'time_cross_'+arlist[i]+'.nc'
        ds.to_netcdf(flnm_save)
        i+=1

def save_data_all():
    model_list = ['GWD3','CTRL', 'SS', 'FD']
    # model_list = ['GWD3', 'CTRL']
# model_list = ['GWD3',]
    # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'
    for model in model_list:
        path = path_main+model+'/wrfout/'
        save_data_one(path)
        # save_one_model_mp(path)


# %%
def draw_one(flnm, fig, ax):
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=slice('2021-07-19 12', '2021-07-21 12'))
    dr = draw_time_height.Draw(fig, ax)
    dr.draw(ds)

def draw_all():
    cm = 1/2.54


    model_list = ['CTRL', 'SS', 'FD', 'GWD3']
    # model_list = ['GWD3',]
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'
    # fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_cross/newall/time_cross/new/'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_cross/time_vs/'
    # area_list = ['south', 'middle', 'north']
    # area_list = ['A', 'B', 'C']
    # area_list = ['B','D', 'E']
    area_list = ['H','I']
    for area in area_list:
        for model in model_list:
            # fig = plt.figure(figsize=(16*cm, 6*cm), dpi=300)
            # ax = fig.add_axes([0.1,0.15, 0.85, 0.75])
            # ax= fig.add_axes([0.1, 0.2, 0.88, 0.7])

            cm = 1/2.54
            fig = plt.figure(figsize=(8*cm, 5*cm), dpi=600)
            ax = fig.add_axes([0.15,0.22, 0.72, 0.67])
            
            path = path_main+model+'/wrfout/'
            # flnm = path+'time_cross_south.nc'
            flnm = path+'time_cross_'+area+'.nc'
            draw_one(flnm,fig, ax)
            fig_name = fig_path+model+area
            print(model+area)
            ax.set_title(model, loc='left')
            ax.set_title(area, loc='right')
            fig.savefig(fig_name)

if __name__ == '__main__':
    
    # save_data_all()
    draw_all()
# %%

# import xarray as xr
# flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_H.nc'
# ds = xr.open_dataset(flnm)
# ds
# # %%
# ds.sel(time=slice('2021-07-19 12', '2021-07-21 12'))