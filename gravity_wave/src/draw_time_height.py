# %%
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from wrf import smooth2d
from matplotlib.ticker import AutoMinorLocator

class Draw():
    """画单张的图
    """
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax

    def draw_contour(self, x, y,da, **kw):
        """在地图上绘制等温线
        """
        levels=np.arange(348, 356+1, 4)
                        
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

        # u = u[::24,::6]
        # v = v[::24,::6]
        # x = x[::6]
        # y = y[::24]
        # u = u[::400,::8]
        # v = v[::400,::8]
        # x = x[::8]
        # y = y[::400]
        u = u[::20,::8]
        v = v[::20,::8]
        x = x[::8]
        y = y[::20]
        Q = self.ax.quiver(x, y, u.values,v.values,units='inches',scale=scale,pivot='middle')  # 绘制风矢
        # qk = self.ax.quiverkey(Q, X=0.67, Y=0.05, U=10, label=r'$(v, 10 m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢
        qk = self.ax.quiverkey(Q, X=0.7, Y=0.05, U=10, label=r'$(v, 10 m/s)$', labelpos='E',coordinates='figure',  fontproperties={'size':10})   # 设置参考风矢

    def draw_contourf(self, xs,ys,da,):

        # xs = np.arange(0, da.shape[-1], 1)
        # ys = da.coords['vertical'].values
        # colordict=['#191970','#005ffb','#5c9aff','#98ff98','#ddfddd','#FFFFFF','#fffde1','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        colordict=['#191970','#005ffb','#5c9aff','#98ff98','#FFFFFF','#ffff9e', '#ffc874','#ffa927', '#ff0000']#颜色列表
        # colorlevel=[-80, -30, -20, -10, -5, -1, 1, 5, 10, 20, 30, 80]#雨量等级
        # colorlevel=[-120, -50, -30, -10, -5, -1, 1, 5, 10, 30, 50, 120]#雨量等级
        # colorlevel=[-120, -50, -30, -10,  -1, 1, 10, 30, 50, 120]#雨量等级
        # colorlevel=[-120, -5, -3, -1,  -0.1, 0.1, 1, 3, 5, 120]#雨量等级
        # colorlevel=[-120, -4, -2, -1,  -0.2, 0.2, 1, 2, 4, 120]#散度等级
        # colorlevel=[-120, -3, -2, -1,  -0.5, 0.5, 1, 2, 3, 120]#散度等级
        colorlevel=[-120,  -2,-1, -0.5,  -0.1, 0.1, 0.5,1, 2, 120]#散度等级
        # colorlevel=[-120, -10,-5, -3, -1, 1, 3, 5, 10, 120]#雨量等级
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

    # def set_ticks(self,x,y, ):
    #     """绘制图片的标签之类的
    #     """
    #     pass
    #     ## 标题之类的
    #     # self.ax.set_title('south', fontsize=10, loc='left')
    #     self.ax.set_xlabel("Time (Day/Hour)", fontsize=10)
    #     self.ax.set_ylabel("Hieght (km)", fontsize=10)

    #     ## 坐标标签
    #     x_ticks = x.values
    #     x_labels = x.time.dt.strftime('%d/%H')
    #     self.ax.set_xticks(x_ticks[::24])
    #     self.ax.set_xticklabels(x_labels[::24].values, rotation=0, fontsize=10)
    #     self.ax.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
    #     self.ax.set_yticks(np.arange(0, 10000.1, 1000).astype('int'))
    #     self.ax.set_yticklabels(np.arange(0, 10.1, 1), rotation=0, fontsize=10)
    #     self.ax.tick_params(axis='both', labelsize=10, direction='out')

    #     ## 图例绘制



    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd0/time_height_gwd0.nc'
    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd1/time_height_gwd1.nc'
    

    # flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/gwd3/time_height_gwd3.nc'

    def switch_p2z(self, ds):    
        """将气压坐标转换为高度坐标, ds中必须含有height变量

        Args:
            ds (_type_): _description_

        Returns:
            _type_: _description_
        """
        ds = ds.reset_coords()
        ds1 = ds.interpolate_na(dim='pressure',method='linear',fill_value="extrapolate")
        ds2 = ds1.interp(pressure=np.arange(0, 1050, 5), method='linear')
        ds3 = ds2.dropna(dim='pressure')
        hh = ds3['height'].mean(dim='time').values
        ds4 = ds3.assign_coords({'z':('pressure',hh)})
        ds5 = ds4.swap_dims({'pressure':'z'})
        return ds5
    
    

    def draw(self, ds):
        

        ## 将气压坐标转为高度坐标
        ds3 = self.switch_p2z(ds)

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
        # self.draw_contourf(x,y,div)

        self.draw_contourf(x,y,w*10)
        # self.ax.contourf(x,y,w*10)
        self.draw_contour(x,y, thetae)
        self.draw_quiver(x,y,u,v)

        self.ax.set_ylim(0, 16000)
        self.ax.set_yticks(np.arange(0, 16000+1, 2000))
        y_labels = np.arange(0, 16+1, 2).astype(int)

        self.ax.set_yticklabels(y_labels)
        self.ax.set_xlabel("Time (Day/Hour)", fontsize=10)
        self.ax.set_ylabel("Hieght (km)", fontsize=10)
        # self.set_ticks(x,y)
        x_ticks = x.values
        x_labels = x.time.dt.strftime('%d/%H')
        self.ax.set_xticks(x_ticks[::12])
        self.ax.set_xticklabels(x_labels[::12].values, rotation=0, fontsize=10)
        
        self.ax.xaxis.set_minor_locator(plt.MultipleLocator(0.125))
        # self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.tick_params(axis='both', labelsize=10, direction='out')

# %%

if __name__ == '__main__':

    flnm = '/mnt/zfm_18T/fengxiang/HeNan/Data/GWD/d03/newall/GWD3/wrfout/time_cross_I.nc'
    ds = xr.open_dataset(flnm)
    ds = ds.sel(time=slice('2021-07-19 12', '2021-07-21 12'))
    cm = 1/2.54
    fig = plt.figure(figsize=(8*cm, 5*cm), dpi=600)
    ax = fig.add_axes([0.15,0.22, 0.72, 0.7])
    dr = Draw(fig, ax)
    dr.draw(ds)
    figpath = '/mnt/zfm_18T/fengxiang/HeNan/gravity_wave/figure/picture_cross/time_vs/'
    figname = 'aa'
    fig.savefig(figpath+figname)
# %%

