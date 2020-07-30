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

def get_variables_to_plot(profiles):
    
    #get only some coordinates of the original xarray
    data = profiles["q_rad"]  
    data["q_rad_lw"] = profiles["q_rad_lw"]
    data["q_rad_sw"] = profiles["q_rad_sw"]
    data["time"] = profiles["launch_time"]
    data = data.drop_vars(["lay","col"])
    
    data = data.to_dataframe()
    data["time"] = data["time"].dt.tz_localize(pytz.UTC).dt.tz_convert('America/Barbados').dt.strftime("%H:%M")
    
    data["time"] = pd.to_datetime(data["time"], format="%H:%M")

    data = data.reset_index()
    data = data.set_index(["time","zlay"])
    data = data.groupby([pd.Grouper(freq='10min', level='time'), 
                                 pd.Grouper(level='zlay')]).mean()
        
    #come back to xarray and get q_rad
    data = data.to_xarray()
    time = data.time.values
    zlay = data.zlay.values
        
    #fill values with 0 in a new array
    ini = np.datetime64('1900-01-01 00:00:00')
    end = ini + np.timedelta64(24,'h')
    count_time = np.arange(ini, end, np.timedelta64(10, 'm'))

    q_rad = np.zeros((len(count_time), len(zlay)))
    q_rad_sw = np.zeros((len(count_time), len(zlay)))
    q_rad_lw = np.zeros((len(count_time), len(zlay)))

    ds =  xr.Dataset({'q_rad': (['count_time', 'zlay'],  q_rad),
                       'q_rad_sw': (['count_time', 'zlay'],  q_rad_sw),
                        'q_rad_lw': (['count_time', 'zlay'],  q_rad_lw)},
                    coords={"count_time": count_time, "zlay": zlay})

    array = ds.to_dataframe()

    for itime in time:
        for izlay in zlay:
            array.q_rad.loc[itime, izlay] = data["q_rad"].sel(time=itime).sel(zlay=izlay).values
            array.q_rad_lw.loc[itime, izlay] = data["q_rad_lw"].sel(time=itime).sel(zlay=izlay).values
            array.q_rad_sw.loc[itime, izlay] = data["q_rad_sw"].sel(time=itime).sel(zlay=izlay).values
    
    data = array.to_xarray()
    
    q_rad = np.transpose(data.q_rad.values)
    q_rad_lw = np.transpose(data.q_rad_lw.values)
    q_rad_sw = np.transpose(data.q_rad_sw.values)
    zlay = data.zlay.values
    time = data.count_time.values
        
    return time, zlay, q_rad, q_rad_lw, q_rad_sw

def plot_diurnal_cycle(profiles):
    
    time, zlay, q_rad, q_rad_lw, q_rad_sw = get_variables_to_plot(profiles)   
 
    dates_list = [date for date in time]    
        
    fig, ax = plt.subplots(1,3,figsize=(20,10))
    
    ax[0].set_title(r'Shortwave')
    ax[0].set_ylabel('Altitude (km)')
    ax[1].set_title('Longwave')
    ax[1].set_xlabel('Time (Local)')
    ax[2].set_title(r'Net')

    ymin=0
    ymax=10
    
    colormap = matplotlib.cm.get_cmap("RdBu_r")
    val_min = -4
    val_max = 4

    zlay=zlay/1000
    
    ax[0].pcolormesh(dates_list, zlay, q_rad_sw, cmap=colormap,vmin=val_min, vmax=val_max)
    ax[1].pcolormesh(dates_list, zlay, q_rad_lw, cmap=colormap,vmin=val_min,vmax=val_max)
    im = ax[2].pcolormesh(dates_list, zlay, q_rad, cmap=colormap,vmin=val_min, vmax=val_max)
    
    myFmt = mdates.DateFormatter('%-H')
    
    for k in range(3):
        ax[k].xaxis.set_major_formatter(myFmt)
        ax[k].set_ylim([0,ymax])
    

    for k in range(3):
        ticks = ax[k].get_xticks()
        ax[k].set_xticks(np.linspace(ticks[0], mdates.date2num(mdates.num2date(ticks[-1])), 4))

    ax[1].tick_params(labelleft=False)    
    ax[2].tick_params(labelleft=False)    

    cb = fig.colorbar(im, ax=ax[2], extend="both")
    cb.ax.set_ylabel('Heating Rate (K/day)')
        
    fig.savefig('../Figures/Fig4_Diurnal_composite.jpg')    
    
plot_diurnal_cycle(sonde_profiles)
