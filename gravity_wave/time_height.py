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
        'lat1':34.4,
        'lat2':34.8,
        'lon1':113,
        'lon2':113.8,
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

    # area = areaA
    # area_list = [areaA, areaB,areaC] 
    # area_list = [areaD, areaE]
    area_list = [areaD]
    # arname_list = ['A']

    # arlist = ['A', 'B', 'C']
    # arlist = ['D', 'E']
    arlist = ['D']

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
    model_list = ['CTRL', 'SS', 'FD', 'GWD3']
    # path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/'
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'
    for model in model_list[1:]:
        path = path_main+model+'/wrfout/'
        save_data_one(path)
        # save_one_model_mp(path)


# %%
def draw_one(flnm, fig, ax):
    ds = xr.open_dataset(flnm)
    dr = draw_time_height.Draw(fig, ax)
    dr.draw(ds)

def draw_all():
    cm = 1/2.54


    model_list = ['CTRL', 'SS', 'FD', 'GWD3']
    # model_list = ['FD',]
    path_main = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/'
    fig_path = '/mnt/zfm_18T/fengxiang/HeNan/Draw/picture_cross/newall/time_cross/new/'
    # area_list = ['south', 'middle', 'north']
    # area_list = ['A', 'B', 'C']
    area_list = ['D', 'E']
    for area in area_list:
        for model in model_list:
            fig = plt.figure(figsize=(17*cm, 8*cm), dpi=300)
            ax = fig.add_axes([0.1,0.15, 0.85, 0.75])
            path = path_main+model+'/wrfout/'
            # flnm = path+'time_cross_south.nc'
            flnm = path+'time_cross_'+area+'.nc'
            draw_one(flnm,fig, ax)
            fig_name = fig_path+model+area
            print(model+area)
            ax.set_title(model, loc='left')
            ax.set_title(area, loc='right')
            fig.savefig(fig_name)
# save_data_all()
draw_all()
# %%