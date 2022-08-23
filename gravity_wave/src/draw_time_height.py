# %%
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt


class Draw():
    """画单张的图
    """
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax

    def draw_contour(self, x, y,da, **kw):
        """在地图上绘制等温线
        """
        levels=np.arange(336, 372, 4)
                        
        crx = self.ax.contour(x,
                            y,
                            da,
                            colors = 'black',
                            levels=levels,
                            linestyles = 'solid',
                            # transform=ccrs.PlateCarree(),
                            # linewidth = 0.1,
                            alpha=0.8)

        self.ax.clabel(crx,inline=1, fontsize=10, colors='blue', zorder=2) # 等值线的标注
        return crx

    def draw_quiver(self, x,y,u,v, scale=50):
        '''
        绘制风矢图
        '''

        u = u[::24,::6]
        v = v[::24,::6]
        x = x[::6]
        y = y[::24]
        Q = self.ax.quiver(x, y, u.values,v.values,units='inches',scale=scale,pivot='middle')  # 绘制风矢
        qk = self.ax.quiverkey(Q, X=0.67, Y=0.05, U=10, label=r'$(v, 10 m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢

    def draw_contourf(self, xs,ys,da,):

        # xs = np.arange(0, da.shape[-1], 1)
        # ys = da.coords['vertical'].values
        # colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        colordict=['#191970','#005ffb','#5c9aff','#98ff98','#FFFFFF','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
        # colorlevel=[-120, -50, -30, -10, -5, -1, 1, 5, 10, 30, 50, 120]#雨量等级
        # colorlevel=[-120, -50, -30, -10,  -1, 1, 10, 30, 50, 120]#雨量等级
        # colorlevel=[-120, -5, -3, -1,  -0.1, 0.1, 1, 3, 5, 120]#雨量等级
        colorlevel=[-120, -10,-5, -3, -1, 1, 3, 5, 10, 120]#雨量等级
        # colorticks=[-30, -20, -10, -5, -1, 1, 5, 10, 20, 30]#雨量等级
        colorticks = colorlevel[1:-1]
        dbz_contours = self.ax.contourf(xs,
                                        ys,
                                        da.values,
                                        colors=colordict,
                                        levels=colorlevel
                                        )
        # ax_cross.set_ylim(0, 12000)
        # ax_cross.set_yticks(np.arange(0, 12000+1000, 1000))
        # cb_dbz = self.fig.colorbar(dbz_contours, ticks=colorticks)
        cb_dbz = self.fig.colorbar(dbz_contours, ax=self.ax, ticks=colorticks, orientation='vertical', fraction=0.06, pad=0.02)
        # ax_cross.tick_params(axis='both', labelsize=18, direction='out')
        # cb_dbz.ax.tick_params(labelsize=10)

        # self.ax.invert_yaxis()

    def set_ticks(self,x,y, ):
        """绘制图片的标签之类的
        """
        pass
        ## 标题之类的
        # self.ax.set_title('south', fontsize=10, loc='left')
        self.ax.set_xlabel("Time (Day/Hour)", fontsize=10)
        self.ax.set_ylabel("Hieght (km)", fontsize=10)

        ## 坐标标签
        x_ticks = x.values
        x_labels = x.time.dt.strftime('%d/%H')
        self.ax.set_xticks(x_ticks[::24])
        self.ax.set_xticklabels(x_labels[::24].values, rotation=0, fontsize=10)
        self.ax.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
        self.ax.set_yticks(np.arange(0, 10000.1, 1000).astype('int'))
        self.ax.set_yticklabels(np.arange(0, 10.1, 1), rotation=0, fontsize=10)
        self.ax.tick_params(axis='both', labelsize=10, direction='out')

        ## 图例绘制



    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/time_height_gwd1.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3.nc'
    def draw(self, ds):

        ds = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
        hh = ds['height'].mean(dim='time').values
        ds2 = ds.assign_coords({'z':('pressure',hh)})
        ds3 = ds2.swap_dims({'pressure':'z'})
        da = ds3['ua'].T
        x = da.time
        y = da.z
        w = (ds3['wa']).T
        
        u = ds3['ua'].T
        v = ds3['va'].T
        w = ds3['wa'].T
        div = ds3['div'].T*10**5
        thetae = ds3['theta_e'].T+273.15
        print(thetae.max())
        
        self.ax.set_xlabel("Time (Day/Hour)", fontsize=22)
        self.ax.set_ylabel("Pressure (hPa)", fontsize=22)
        print(x.shape, y.shape, w.shape)
        self.draw_contourf(x,y,div)
        self.draw_contour(x,y, thetae)
        self.ax.set_ylim(0, 10000)
        # self.draw_quiver(x,y,u,v)
        self.set_ticks(x,y)

if __name__ == '__main__':

    # path_save = '/mnt/zfm_18T/fengxiang/HeNan/draw_gravity/data/'
    # flnm_save = path_save+'gwd3.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_south.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_middle.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/CTRL/wrfout/time_cross_A.nc'
    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_A.nc'
    ds = xr.open_dataset(flnm)
    ds
    # cm = 1/2.54
    # fig = plt.figure(figsize=(17*cm, 8*cm), dpi=600)
    # ax = fig.add_axes([0.1,0.15, 0.85, 0.8])
    # dr = Draw(fig, ax)
    # dr.draw(ds)
    # fig.savefig('aa.png')


# %%
flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_A.nc'
ds = xr.open_dataset(flnm)
ds
    