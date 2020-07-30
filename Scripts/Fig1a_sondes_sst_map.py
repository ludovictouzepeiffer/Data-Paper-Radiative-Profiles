# -*- coding: utf-8 -*-
"""
Drawing Figure 1a.

B. Fildier
"""

# load modules

import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable


#%% Main script

if __name__ == "__main__":

#%% 
    
    # get the mean 
    
    print('import sonde data')

    # REPLACE PATHS BELOW
    inputdir_sondes = '/Users/bfildier/Data/EUREC4A/merged/sondes'
    inputdir_SST = '/Users/bfildier/Data/EUREC4A/ERA'
    
    radiosondes = xr.open_dataset(os.path.join(inputdir_sondes,'all_radiosondes.nc'))
    dropsondes = xr.open_dataset(os.path.join(inputdir_sondes,'all_dropsondes.nc'))
    
    
    print('import SST data')    
    sst_source = 'ERA'
    if sst_source == 'ERA':
        sst_all = xr.open_dataset(os.path.join(inputdir_SST,'SST_2020_01_02_ERA5_hourly_Barbados.nc'))
        sst = sst_all.mean(dim='time')
    else:
        sst = xr.open_dataset('../Input/SST_mean.nc')
        
    print('define coordinates for map')

    if sst_source == 'ERA':
        x_sst,y_sst = np.meshgrid(sst.longitude.values,sst.latitude.values)
    else:
        x_sst,y_sst = np.meshgrid(sst.lon.values, sst.lat.values)
        
    x_ds,y_ds = dropsondes.lon.values.flatten(),dropsondes.lat.values.flatten()
    x_rs,y_rs = radiosondes.lon.values.flatten(),radiosondes.lat.values.flatten()    
    
    print('draw figure')
    
    fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree()},figsize=(10,10))
    
    ax.coastlines(resolution='50m')
    ax.set_extent([-60,-49,6,17])
    gl = ax.gridlines(color='Grey',draw_labels=True)
    fs = 12
    gl.xlabel_style = {'size': fs}
    gl.ylabel_style = {'size': fs}

    # plot
    sst_cmap = plt.cm.viridis_r
    alpha = 1
    
    if sst_source == 'ERA':
        sst_values = sst.sstk.values
    else:
        sst_values = sst.SST.values
    
    ssthm = ax.pcolormesh(x_sst,y_sst,sst_values,cmap=sst_cmap, alpha=alpha,edgecolors='face', antialiased=True, transform=ccrs.PlateCarree())
    
    dm = ax.scatter(x_ds,y_ds,marker='o',color='white',alpha=0.1,s=1,linewidths=0.01,label='Dropsondes')
    rm = ax.scatter(x_rs,y_rs,marker='o',color='tomato',alpha=0.03,s=1,linewidths=0.01,label='Radiosondes')
        
    # ax.legend(handles=[dm,rm],fontsize=18)
    
    # add colorbar
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="3%", pad=0.7, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(ssthm, cax=ax_cb)
    cbar.solids.set(alpha=1)
        
    # Change alpha
    fig.canvas.draw()
    colors = ssthm.get_facecolor()
    def alpha_to_white(color):
        white = np.array([1,1,1])
        alpha = color[-1]
        color = color[:-1]
        return alpha*color + (1 - alpha)*white
    colors = np.array([alpha_to_white(color) for color in colors])
    ssthm.set_facecolor(colors)

    # labels
    cbar.ax.tick_params(labelsize=fs)
    cbar.ax.set_ylabel('SST (K)',fontsize=fs)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.text(0.03,0.94,'(a)',transform = ax.transAxes,fontsize=18,color='white')
    
    # save
    plt.savefig('../Figures/Figure1a.png',bbox_inches='tight')
    plt.show()