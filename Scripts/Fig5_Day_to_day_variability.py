import numpy as np
from netCDF4 import Dataset, num2date # to work with NetCDF files
from os.path import expanduser
import matplotlib.pyplot as plt
home = expanduser("~") # Get users home directory
import statsmodels.api as sm
from scipy import stats

import xarray as xr
import pytz
import glob, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt     
import matplotlib.dates as mdates
import datetime
from dateutil import tz
import metpy.calc as mpcalc
from metpy.units import units

import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import matplotlib.cm as cmx

matplotlib.rcParams.update({'font.size': 24})

#Directory where sondes are stored
dir_profile = "/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Products/"
path_to_sonde_profiles = os.path.join(dir_profile,"rad_profiles_all_sondes_ERA.nc")

sonde_profiles = xr.open_dataset(path_to_sonde_profiles)

sonde_BCO = sonde_profiles.where(sonde_profiles.platform=="BCO", drop=True)

def get_variables_day_to_day(profiles):

    data = profiles["q_rad"]
    data["q_rad_sw"] = profiles["q_rad_lw"]
    data["q_rad_lw"] = profiles["q_rad_sw"]

    data["time"] = profiles["launch_time"]
    data = data.drop_vars(["lay","col"])
    
    data = data.to_dataframe()
    data["time"] = data["time"].dt.tz_localize(pytz.UTC).dt.tz_convert('America/Barbados').dt.strftime("%Y%m%dT%H")
    
    data["time"] = pd.to_datetime(data["time"], format="%Y%m%dT%H")

    data = data.reset_index()
    data = data.set_index(["time","zlay"])
    data = data.groupby([pd.Grouper(freq='1H', level='time', label="right"), 
                                 pd.Grouper(level='zlay')]).mean()
    
    #come back to xarray and get q_rad
    data = data.to_xarray()
    
    time = data.time.values
    zlay = data.zlay.values
    q_rad = np.transpose(data.q_rad.values)
    q_rad_sw = np.transpose(data.q_rad_lw.values)
    q_rad_lw = np.transpose(data.q_rad_sw.values)

    return time, zlay, q_rad, q_rad_lw, q_rad_sw

def plot_day_to_day(profiles):
    
    time, zlay, q_rad, q_rad_lw, q_rad_sw = get_variables_day_to_day(profiles)   
 
    dates_list = [date for date in time]    
        
    fig, ax = plt.subplots(3,1,figsize=(20,30))
    
    fig.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=0.2, hspace=0.2)

    pad=10
    fs=24
    loc="left"
    fw="bold"
    ax[0].set_title(r'a) Shortwave', loc=loc, pad=pad,fontsize=fs, fontweight=fw)
    ax[1].set_title(r'b) Longwave', loc=loc, pad=pad, fontsize=fs, fontweight=fw)
    ax[2].set_title(r'c) Net', loc=loc, pad=pad, fontsize=fs, fontweight=fw)

    ax[1].set_ylabel('Altitude (km)')
    ax[2].set_xlabel('Date')

    ymin=0
    ymax=10
    
    colormap = matplotlib.cm.get_cmap("RdBu_r")
    val_min = -4
    val_max = 4

    zlay=zlay/1000
    
    ax[0].pcolormesh(dates_list, zlay, q_rad_sw, cmap=colormap,vmin=val_min, vmax=val_max)
    ax[1].pcolormesh(dates_list, zlay, q_rad_lw, cmap=colormap,vmin=val_min, vmax=val_max)
    im = ax[2].pcolormesh(dates_list, zlay, q_rad, cmap=colormap,vmin=val_min, vmax=val_max)

    myFmt = mdates.DateFormatter('%m-%d')
    
    ini = np.datetime64('2020-01-19 00:00:00')
    end = np.datetime64('2020-02-17 00:00:00')
    
    for k in range(3):
        ax[k].xaxis.set_major_formatter(myFmt)
        ax[k].set_ylim([0,ymax])
        ax[k].set_xlim([ini,end])

    for k in range(3):
        ticks = ax[k].get_xticks()
        ax[k].set_xticks(np.linspace(ticks[0], ticks[-1], 10))

    ax[0].tick_params(labelbottom=False)    
    ax[1].tick_params(labelbottom=False)    

    x,y,w,h = ax[2].get_position().bounds
    c_map_ax = fig.add_axes([x, y-0.25*h, 1*w, 0.06*h])
    cbar = fig.colorbar(im,cax=c_map_ax, orientation="horizontal", extend="both")
    cbar.ax.set_xlabel('Heating Rate (K/day)',color='k') # cbar legend   
    
    fig.savefig('../Figures/Fig5_Day_to_day_variability.jpg')  
    
plot_day_to_day(sonde_BCO)
