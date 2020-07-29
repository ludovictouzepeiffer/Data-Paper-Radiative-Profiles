import numpy as np
from netCDF4 import Dataset, num2date # to work with NetCDF files
from os.path import expanduser
import matplotlib.pyplot as plt
from scipy import stats

import xarray as xr
import pytz
import glob, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt     
import matplotlib.dates as mdates
import datetime
import metpy.calc as mpcalc
from metpy.units import units

import pandas as pd

matplotlib.rcParams.update({'font.size': 24})


#Directory where sondes are stored
dir_profile = "/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Products/"
path_to_sonde_profiles = os.path.join(dir_profile,"rad_profiles_all_sondes_ERA.nc")

sonde_profiles = xr.open_dataset(path_to_sonde_profiles)

def get_relative_humidity(profiles):
    
    profiles["play"].attrs['units'] = 'hPa'
    profiles["tlay"].attrs['units'] = 'kelvin'
    
    rh = mpcalc.relative_humidity_from_mixing_ratio(profiles["mr"], profiles["tlay"], profiles["play"])
    profiles["rh"] = (["launch_time","zlay"], rh.magnitude)  
    
    return profiles

def get_specific_humidity(profiles):
        
    qv = mpcalc.specific_humidity_from_mixing_ratio(profiles["mr"])
    profiles["qv"] = (["launch_time","zlay"], qv.magnitude)  
    
    return profiles

def plot_average_profiles(profiles):
    
    profiles = get_relative_humidity(profiles)
    profiles = get_specific_humidity(profiles)
       
    tlay_mean = profiles["tlay"].quantile(0.5, dim="launch_time")
    tlay_fq = (profiles["tlay"]).quantile(0.25, dim="launch_time")
    tlay_lq = (profiles["tlay"]).quantile(0.75, dim="launch_time")

    qv_mean = profiles["qv"].quantile(0.5, dim="launch_time")*1000
    qv_fq = (profiles["qv"]).quantile(0.25, dim="launch_time")*1000
    qv_lq = (profiles["qv"]).quantile(0.75, dim="launch_time")*1000
    
    rh_mean = profiles["rh"].quantile(0.5, dim="launch_time")
    rh_fq = (profiles["rh"]).quantile(0.25, dim="launch_time")
    rh_lq = (profiles["rh"]).quantile(0.75, dim="launch_time")

    q_rad_mean = profiles["q_rad"].quantile(0.5, dim="launch_time")
    q_rad_fq = profiles["q_rad"].quantile(0.25, dim="launch_time")
    q_rad_lq = profiles["q_rad"].quantile(0.75, dim="launch_time")
  
    q_rad_lw_mean = profiles["q_rad_lw"].quantile(0.5, dim="launch_time")
    q_rad_lw_fq = profiles["q_rad_lw"].quantile(0.25, dim="launch_time")
    q_rad_lw_lq = profiles["q_rad_lw"].quantile(0.75, dim="launch_time")
    
    q_rad_sw = profiles["q_rad_sw"].where(profiles["q_rad_sw"].mean(dim="zlay") > 0, drop=True)
    
    q_rad_sw_mean = q_rad_sw.quantile(0.5, dim="launch_time", skipna=True)
    q_rad_sw_fq = q_rad_sw.quantile(0.25, dim="launch_time", skipna=True)
    q_rad_sw_lq = q_rad_sw.quantile(0.75, dim="launch_time", skipna=True)
    
    zlay = profiles["zlay"]/1000
    
    fig, ax = plt.subplots(2,3,figsize=(20,20))
    
    ax[0,0].set_xlabel('Temperature (K)')
    ax[0,1].set_xlabel('Specific humidity (g/kg)')
    ax[0,2].set_xlabel('Relative humidity (%)')
    
    ax[0,0].set_ylabel('Altitude (km)')
    ax[1,0].set_ylabel('Altitude (km)')

    fs=24
    ax[0,1].set_title('Environmental means', fontsize=fs)
    ax[1,0].set_title('Shortwave', fontsize=fs)
    ax[1,1].set_title('Longwave', fontsize=fs)
    ax[1,2].set_title('Net', fontsize=fs)

    ax[1,1].set_xlabel('Heating rates (K/day)')
    
    ymin=0
    ymax=10
    for k in range(3):
        ax[1,k].set_xlim([-4.9,4.9])
        for i in range (2):
            ax[i,k].grid(color='k', linestyle='--', linewidth=0.8)
            ax[i,k].set_ylim([ymin,ymax]) 
            ax[i,k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
            for axis in ['top','bottom','left','right']:
                ax[i,k].spines[axis].set_linewidth(1.3)
            ax[i,k].spines['right'].set_visible(False)
            ax[i,k].spines['top'].set_visible(False)
        
    cl= "k"
    alpha=0.20
    
    ax[0,0].fill_betweenx(zlay,tlay_fq, tlay_lq, alpha=alpha, color=cl)    
    ax[0,0].plot(tlay_mean, zlay, color=cl)

    ax[0,1].fill_betweenx(zlay,qv_fq, qv_lq, alpha=alpha, color=cl)    
    ax[0,1].plot(qv_mean, zlay, color=cl)
    
    ax[0,2].fill_betweenx(zlay,rh_fq, rh_lq, alpha=alpha, color=cl)    
    ax[0,2].plot(rh_mean, zlay, color=cl)
    
    ax[1,0].fill_betweenx(zlay,q_rad_sw_fq, q_rad_sw_lq, alpha=alpha, color=cl)    
    ax[1,0].plot(q_rad_sw_mean, zlay, color=cl)

    ax[1,1].fill_betweenx(zlay,q_rad_lw_fq, q_rad_lw_lq, alpha=alpha, color=cl)    
    ax[1,1].plot(q_rad_lw_mean, zlay, color=cl)
    
    ax[1,2].fill_betweenx(zlay,q_rad_fq, q_rad_lq, alpha=alpha, color=cl)    
    ax[1,2].plot(q_rad_mean, zlay, color=cl)
    
    x_text=0.8
    y_text=0.9
    ax[0,0].text(x_text,y_text,'(a)',transform = ax[0,0].transAxes,fontsize=fs)
    ax[0,1].text(x_text,y_text,'(b)',transform = ax[0,1].transAxes,fontsize=fs)
    ax[0,2].text(x_text,y_text,'(c)',transform = ax[0,2].transAxes,fontsize=fs)
    ax[1,0].text(x_text,y_text,'(d)',transform = ax[1,0].transAxes,fontsize=fs)
    ax[1,1].text(x_text,y_text,'(e)',transform = ax[1,1].transAxes,fontsize=fs)
    ax[1,2].text(x_text,y_text,'(f)',transform = ax[1,2].transAxes,fontsize=fs)
    
    fig.savefig('../Figures/Fig3_Average_profiles.jpg')
    
plot_average_profiles(sonde_profiles)