#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from netCDF4 import Dataset, num2date # to work with NetCDF files
from os.path import expanduser
import matplotlib.pyplot as plt
import xarray as xr
import glob, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt     
import matplotlib.dates as mdates
import datetime
import metpy.calc as mpcalc
from metpy.units import units
import seaborn as sns
import pandas as pd

sns.set(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=3.2, color_codes=True, rc=None)

#Directory where sondes are stored
#dir_profile = "/media/ludo/DATA/google-drive/Thèse/EUREC4a/github/Input/Products/"
#path_to_sonde_profiles = os.path.join(dir_profile,"rad_profiles_all_sondes_ERA.nc")

dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Data/"
fp_rad_profiles = os.path.join(dir_profile, "rad_profiles_all_sondes_ERA.nc")
sonde_profiles = xr.open_dataset(fp_rad_profiles)

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
       
    tlay_median = profiles["tlay"].quantile(0.5, dim="launch_time")
    tlay_25 = (profiles["tlay"]).quantile(0.25, dim="launch_time")
    tlay_75 = (profiles["tlay"]).quantile(0.75, dim="launch_time")
    tlay_05 = (profiles["tlay"]).quantile(0.05, dim="launch_time")
    tlay_95 = (profiles["tlay"]).quantile(0.95, dim="launch_time")

    qv_median = profiles["qv"].quantile(0.5, dim="launch_time")*1000
    qv_25 = (profiles["qv"]).quantile(0.25, dim="launch_time")*1000
    qv_75 = (profiles["qv"]).quantile(0.75, dim="launch_time")*1000
    qv_05 = (profiles["qv"]).quantile(0.05, dim="launch_time")*1000
    qv_95 = (profiles["qv"]).quantile(0.95, dim="launch_time")*1000

    rh_median = profiles["rh"].quantile(0.5, dim="launch_time")
    rh_25 = (profiles["rh"]).quantile(0.25, dim="launch_time")
    rh_75 = (profiles["rh"]).quantile(0.75, dim="launch_time")
    rh_05 = (profiles["rh"]).quantile(0.05, dim="launch_time")
    rh_95 = (profiles["rh"]).quantile(0.95, dim="launch_time")

    q_rad_median = profiles["q_rad"].quantile(0.5, dim="launch_time")
    q_rad_25 = (profiles["q_rad"]).quantile(0.25, dim="launch_time")
    q_rad_75 = (profiles["q_rad"]).quantile(0.75, dim="launch_time")
    q_rad_05 = (profiles["q_rad"]).quantile(0.05, dim="launch_time")
    q_rad_95 = (profiles["q_rad"]).quantile(0.95, dim="launch_time")

    q_rad_lw_median = profiles["q_rad_lw"].quantile(0.5, dim="launch_time")
    q_rad_lw_25 = (profiles["q_rad_lw"]).quantile(0.25, dim="launch_time")
    q_rad_lw_75 = (profiles["q_rad_lw"]).quantile(0.75, dim="launch_time")
    q_rad_lw_05 = (profiles["q_rad_lw"]).quantile(0.05, dim="launch_time")
    q_rad_lw_95 = (profiles["q_rad_lw"]).quantile(0.95, dim="launch_time")

    q_rad_sw = profiles["q_rad_sw"].where(profiles["q_rad_sw"].mean(dim="zlay") > 0, drop=True)
    q_rad_sw_median = q_rad_sw.quantile(0.5, dim="launch_time")
    q_rad_sw_25 = q_rad_sw.quantile(0.25, dim="launch_time")
    q_rad_sw_75 = q_rad_sw.quantile(0.75, dim="launch_time")
    q_rad_sw_05 = q_rad_sw.quantile(0.05, dim="launch_time")
    q_rad_sw_95 = q_rad_sw.quantile(0.95, dim="launch_time")

    zlay = profiles["zlay"]/1000
    
    fig, ax = plt.subplots(2,3,figsize=(30,30))
    
    ax[0,0].set_xlabel('Temperature (K)')
    ax[0,1].set_xlabel('Specific humidity (g/kg)')
    ax[0,2].set_xlabel('Relative humidity (%)')
    
    ax[0,0].set_ylabel('Altitude (km)')
    ax[1,0].set_ylabel('Altitude (km)')

    fs=43
    ax[0,1].set_title('Environmental means', fontsize=fs)
    ax[1,0].set_title('Shortwave', fontsize=fs)
    ax[1,1].set_title('Longwave', fontsize=fs)
    ax[1,2].set_title('Net', fontsize=fs)

    ax[1,1].set_xlabel('Heating rates (K/day)')
    
    ymin=0.03
    ymax=10
    for k in range(3):
        ax[1,k].set_xlim([-6.5,6.5])
        for i in range (2):
            ax[i,k].grid(color='k', linestyle='--', linewidth=0.8)
            ax[i,k].set_ylim([ymin,ymax]) 
            ax[i,k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
            for axis in ['top','bottom','left','right']:
                ax[i,k].spines[axis].set_linewidth(1.3)
            ax[i,k].spines['right'].set_visible(False)
            ax[i,k].spines['top'].set_visible(False)
        
    cl= "k"
    alpha=0.30
    alpha1=0.10

    ax[0,0].plot(tlay_median, zlay, color=cl, linewidth=3, label="median")
    ax[0,0].fill_betweenx(zlay,tlay_25, tlay_75, alpha=alpha, color=cl, label="25-75%")    
    ax[0,0].fill_betweenx(zlay,tlay_05, tlay_25, alpha=alpha1, color=cl, label="5-95%")    
    ax[0,0].fill_betweenx(zlay,tlay_75, tlay_95, alpha=alpha1, color=cl)    
    ax[0,0].legend(loc="lower left")

    ax[0,1].fill_betweenx(zlay, qv_25, qv_75, alpha=alpha, color=cl)    
    ax[0,1].fill_betweenx(zlay,qv_05, qv_25, alpha=alpha1, color=cl)    
    ax[0,1].fill_betweenx(zlay,qv_75, qv_95, alpha=alpha1, color=cl)    
    ax[0,1].plot(qv_median, zlay, color=cl, linewidth=3)
    
    ax[0,2].fill_betweenx(zlay, rh_25, rh_75, alpha=alpha, color=cl)    
    ax[0,2].fill_betweenx(zlay,rh_05, rh_25, alpha=alpha1, color=cl)    
    ax[0,2].fill_betweenx(zlay,rh_75, rh_95, alpha=alpha1, color=cl)    
    ax[0,2].plot(rh_median, zlay, color=cl, linewidth=3)

    ax[1,0].fill_betweenx(zlay, q_rad_sw_25, q_rad_sw_75, alpha=alpha, color=cl)    
    ax[1,0].fill_betweenx(zlay,q_rad_sw_05, q_rad_sw_25, alpha=alpha1, color=cl)    
    ax[1,0].fill_betweenx(zlay,q_rad_sw_75, q_rad_sw_95, alpha=alpha1, color=cl)    
    ax[1,0].plot(q_rad_sw_median, zlay, color=cl, linewidth=3)
    
    ax[1,1].fill_betweenx(zlay, q_rad_lw_25, q_rad_lw_75, alpha=alpha, color=cl)    
    ax[1,1].fill_betweenx(zlay,q_rad_lw_05, q_rad_lw_25, alpha=alpha1, color=cl)    
    ax[1,1].fill_betweenx(zlay,q_rad_lw_75, q_rad_lw_95, alpha=alpha1, color=cl)    
    ax[1,1].plot(q_rad_lw_median, zlay, color=cl, linewidth=3)
    
    ax[1,2].fill_betweenx(zlay, q_rad_25, q_rad_75, alpha=alpha, color=cl)    
    ax[1,2].fill_betweenx(zlay,q_rad_05, q_rad_25, alpha=alpha1, color=cl)    
    ax[1,2].fill_betweenx(zlay,q_rad_75, q_rad_95, alpha=alpha1, color=cl)    
    ax[1,2].plot(q_rad_median, zlay, color=cl, linewidth=3)
    
    x_text=0.9
    y_text=0.9
    ax[0,0].text(x_text,y_text,'(a)',transform = ax[0,0].transAxes,fontsize=fs)
    ax[0,1].text(x_text,y_text,'(b)',transform = ax[0,1].transAxes,fontsize=fs)
    ax[0,2].text(x_text,y_text,'(c)',transform = ax[0,2].transAxes,fontsize=fs)
    ax[1,0].text(x_text,y_text,'(d)',transform = ax[1,0].transAxes,fontsize=fs)
    ax[1,1].text(x_text,y_text,'(e)',transform = ax[1,1].transAxes,fontsize=fs)
    ax[1,2].text(x_text,y_text,'(f)',transform = ax[1,2].transAxes,fontsize=fs)
    fig.tight_layout() 
    
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    
    fig.savefig('/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Figures/Paper_figures/Fig2_Average_profiles_edit.png')
    
plot_average_profiles(sonde_profiles)
