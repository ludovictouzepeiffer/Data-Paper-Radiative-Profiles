#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bar plot of number of sondes (dropsondes and radiosondes) launched per 10-minute bin

@author: annaleaalbright
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pytz
import os
import datetime
import pandas as pd
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.dates as mdates
plt.rcParams.update({'font.size': 18})


#%% 
# =============================================================================
#             open files
# =============================================================================  
dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/Dropsondes/Data/"
path_to_dropsonde_profiles = os.path.join(dir_profile,"rad_profiles_all_dropsondes.nc")
path_to_radiosonde_profiles = os.path.join(dir_profile, "rad_profiles_all_radiosondes.nc")

dropsonde_profiles = xr.open_dataset(path_to_dropsonde_profiles)
radiosonde_profiles = xr.open_dataset(path_to_radiosonde_profiles)
all_profiles = xr.concat([dropsonde_profiles, radiosonde_profiles], dim="launch_time")

#%%

# =============================================================================
#             PDF of sondes per hour, by platform
# =============================================================================  

def calc_pdf_sondes(profiles):
    
    # select coordinates
    data = profiles["q_rad"]  
    data["time"] = profiles["launch_time"]
    data = data.drop_vars(["lay","col", "play"])
    data = data.to_dataframe()
    data["time"] = data["time"].dt.tz_localize(pytz.UTC).dt.tz_convert('America/Barbados').dt.strftime("%H:%M")
    data["time"] = pd.to_datetime(data["time"], format="%H:%M")

    data = data.reset_index()
    data = data.set_index(["time","zlay"])
    data = data.groupby([pd.Grouper(freq='10min', level='time'), pd.Grouper(level='zlay')]).count()
    data = data.groupby(["time"]).mean()

    # get time
    data = data.to_xarray()
    time = data.time.values
    
    ini = np.datetime64('1900-01-01 00:00:00')
    end = ini + np.timedelta64(24,'h')
    count_time = np.arange(ini, end, np.timedelta64(10, 'm'))
    count = np.zeros(len(count_time))

    array = {'count_time': count_time, 'count': count}
    array = pd.DataFrame(data=array)
    array = array.set_index(["count_time"])

    for itime in time:
        array.loc[itime, "count"] = data["q_rad"].sel(time=itime).values
       
    return count_time, array["count"].values
   
time, count_vec_dropsondes = calc_pdf_sondes(dropsonde_profiles)
time, count_vec_radiosondes = calc_pdf_sondes(radiosonde_profiles)


#%% try out different plots

# ordinary line plot (not ideal)

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1, 1, 1) 
ax.plot(time, count_vec_dropsondes, linewidth=4, color='crimson', alpha =1, label="dropsondes (1140 total)")  
ax.plot(time, count_vec_radiosondes, linewidth=2, color='royalblue', alpha = 0.8, label="radiosondes (1177 total)")  
myFmt = mdates.DateFormatter("%-H")
ax.xaxis.set_major_formatter(myFmt)
ticks = ax.get_xticks()
ax.set_xticks(np.linspace(ticks[0], mdates.date2num(mdates.num2date(ticks[-2])), 8))
ax.set_ylabel("sondes per interval")
ax.set_xlabel("time (10-minute intervals)")
ax.grid(True, alpha=0.5)
ax.set_ylim([0, 120])
ax.legend(loc='best', frameon=False, fontsize=20)

# bar plot, but without dates (finish)
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1, 1, 1) 
test_x_radiosondes = np.arange(len(count_vec_radiosondes))
test_x_dropsondes = np.arange(len(count_vec_dropsondes))
ax.bar(test_x_radiosondes, count_vec_radiosondes, color='royalblue', label="radiosondes (1177 total)")
ax.bar(test_x_dropsondes, count_vec_dropsondes, color="crimson", alpha = 1, label="dropsondes (1140 total)")
ax.set_ylabel("sondes per interval")
ax.set_xlabel("bin count")
ax.grid(True, alpha=0.5)
ax.set_ylim([0, 120])
ax.legend(loc='best', frameon=False, fontsize=20)

# not working
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1, 1, 1) 
ax.bar(dates_list_radiosondes, count_vec_radiosondes, color='royalblue', label="radiosondes (1177 total)")
ax.bar(dates_list_dropsondes, count_vec_dropsondes, color="crimson", alpha = 1, label="dropsondes (1140 total)")
ax.set_ylabel("sondes per interval")
ax.set_xlabel("time (10-minute intervals)")
myFmt = mdates.DateFormatter("%-H")
ax.xaxis.set_major_formatter(myFmt)
ticks = ax.get_xticks()
ax.set_xticks(np.linspace(ticks[0], mdates.date2num(mdates.num2date(ticks[-1])), 8))
ax.grid(True, alpha=0.5)
ax.set_ylim([0, 120])
ax.legend(loc='best', frameon=False, fontsize=20)
