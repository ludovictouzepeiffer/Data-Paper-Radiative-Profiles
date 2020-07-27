#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 6

# plot mixing ratio, theta, and LW / SW / net radiative heating rates for four example days (one for each pattern)
# and metric of spatial variance in radiative heating rates

"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
plt.rcParams.update({'font.size': 24})
import os
import seaborn as sns
sns.set(context='notebook', style='white', palette='deep', font='sans-serif', font_scale=3, color_codes=True, rc=None)

#%% 
# =============================================================================
#            pre-processing functions
# =============================================================================  

def choose_HALO_circle(sondes): 
    HALO_sondes = sondes.where(sondes.Platform=='HALO')
    # HALO circle, source: http://eurec4a.eu/fileadmin/user_upload/eurec4a/documents/Flight-Planning-Coordinates.pdf
    lon_center, lat_center = -57.717,13.3
    lon_pt_circle, lat_pt_circle = -57.245,14.1903
    r_circle = np.sqrt((lon_pt_circle-lon_center)**2+(lat_pt_circle-lat_center)**2)
    lat_N = lat_center + r_circle 
    lat_S = lat_center - r_circle 
    lon_W = lon_center - r_circle 
    lon_E = lon_center + r_circle    
    sondes_circle = HALO_sondes.where((HALO_sondes['lat']<lat_N) & (HALO_sondes['lat']>lat_S)\
                                      & (HALO_sondes['lon']<lon_E) &(HALO_sondes['lon']>lon_W), drop=True)
    nsondes_HALO_circle = len(sondes_circle.launch_time)
    print(nsondes_HALO_circle, "total sondes launched in HALO circle")
    return sondes_circle

def get_data_one_day(fp_dropsondes, day_str):
    all_sondes = xr.open_dataset(fp_dropsondes).swap_dims({"sounding": "launch_time"})
    sondes_oneday = all_sondes.sel(launch_time=day_str) 
    Sondes_Circle = choose_HALO_circle(sondes_oneday)
    
    # remove sondes with large nunber of NaNs
    # threshold : require this many non-NA values
    sondes_circle = Sondes_Circle.dropna(dim="launch_time", \
                                         subset=['q'], \
                                         how='any', thresh=300) 
    
    nsondes_qc = len(sondes_circle['launch_time'])
    print(nsondes_qc, "sondes after quality control on " + day_str)

    # =============================================================================
    #          load environmental data
    # ============================================================================= 

    kgtog = 1000
    CtoK = 273.15
    
    # load specific humidity, convert from kg/kg to g/kg
    if (sondes_circle['q'].max().values > 10): 
        specific_humidity_g_kg  = sondes_circle['q']  * kgtog
    else:
        specific_humidity_g_kg  = sondes_circle['q']
            
    # convert temperature from C to K
    if (sondes_circle['T'].max().values < 100): 
        temp_K  = sondes_circle['T'] + CtoK
    else:
        temp_K  = sondes_circle['T']
            
    # convert RH to % units
    RH  = sondes_circle['rh']
    if (sondes_circle['rh'].max().values < 1): 
        RH  = sondes_circle['rh'] * 100
    else:
        RH  = sondes_circle['rh']
    
    launch_time = sondes_circle['launch_time']
    alt_vec = sondes_circle['height']

    # save as xarray Dataset    
    xr_day = xr.Dataset(
            data_vars={'specific_humidity': (('launch_time', 'height'), specific_humidity_g_kg),
                       'temp_K':    (('launch_time', 'height'), temp_K),
                       'RH':    (('launch_time', 'height'), RH),},
            coords={'launch_time': launch_time,
                    'height': alt_vec})
    
    return xr_day

#%%
# =============================================================================
# Load preprocessed environmental data in HALO circle 
# =============================================================================  

# load JOANNE dropsondes
input_dir = '/Users/annaleaalbright/Dropbox/EUREC4A/Dropsondes/Data/'
fp_dropsondes = os.path.join(input_dir,'EUREC4A_JOANNE_Dropsonde-RD41_Level_3_v0.5.7-alpha+0.g45fe69d.dirty.nc')

day_str_fish = '2020-01-22' 
xr_fish = get_data_one_day(fp_dropsondes, day_str_fish)

day_str_flower = '2020-02-02' 
xr_flower = get_data_one_day(fp_dropsondes, day_str_flower)

day_str_gravel = '2020-02-05' 
xr_gravel = get_data_one_day(fp_dropsondes, day_str_gravel)

day_str_sugar = '2020-02-09' 
xr_sugar = get_data_one_day(fp_dropsondes, day_str_sugar)

#%%

def plot_mean_environmental_profiles(xr_fish, xr_flower, xr_gravel, xr_sugar):
    """ 
    Plot mean temperature, specific humidity, and relative humidity for each day
    
    """
    fig, ax = plt.subplots(1,3,figsize=(20,8))

    ax[0].set_ylabel('Altitude (km)')
    ax[0].set_xlabel('Temperature / K ')
    ax[1].set_xlabel('Specific humidity / g/kg')
    ax[2].set_xlabel('Relative humidity / %')
    ax[1].set_title('Environmental means')

    Dates = ['fish: 2020-01-22', 'flower: 2020-02-02', 'gravel: 2020-02-05', 'sugar: 2020-02-09']
    var_vec = ["temp_K","specific_humidity", "RH"]

    ymin=0
    ymax=10
    for k in range(3):
        #ax[k].grid(color='k', linestyle='-', linewidth=0.8)
        ax[k].grid(True, alpha=0.7)
        ax[k].set_ylim([ymin,ymax]) 
        ax[k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
        ax[k].spines['right'].set_visible(False)
        ax[k].spines['top'].set_visible(False)
        for axis in ['top','bottom','left','right']:
              ax[k].spines[axis].set_linewidth(1.3)
    
    # =============================================================================
    #      temperature
    # =============================================================================
    var = var_vec[0]
    var1 = xr_fish[var]
    var2 = xr_flower[var]
    var3 = xr_gravel[var]
    var4 = xr_sugar[var]
    km_vec = var4.height.values / 1000
      
    ax[0].plot(var1.mean(dim='launch_time'), km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[0].plot(var2.mean(dim='launch_time'), km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[0].plot(var3.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[0].plot(var4.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])


    # =============================================================================
    #     specific humidity
    # =============================================================================
    var = var_vec[1]
    var1 = xr_fish[var]
    var2 = xr_flower[var]
    var3 = xr_gravel[var]
    var4 = xr_sugar[var]
    km_vec = var4.height.values / 1000

    ax[1].plot(var1.mean(dim='launch_time'), km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[1].plot(var2.mean(dim='launch_time'), km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[1].plot(var3.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[1].plot(var4.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])
    
    ax[1].tick_params(labelleft=False)
    ax[0].set_ylim([0,9])
    ax[1].set_ylim([0,9])
    
    # =============================================================================
    #      relative humidity
    # =============================================================================
    var = var_vec[2]
    var1 = xr_fish[var]
    var2 = xr_flower[var]
    var3 = xr_gravel[var]
    var4 = xr_sugar[var]
    km_vec = var4.height.values / 1000

    ax[2].plot(var1.mean(dim='launch_time'), km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[2].plot(var2.mean(dim='launch_time'), km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[2].plot(var3.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[2].plot(var4.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])    
    ax[2].tick_params(labelleft=False)
    ax[2].set_ylim([0,9])

    fig.tight_layout() 
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    
#%% call function to plot

plot_mean_environmental_profiles(xr_fish, xr_flower, xr_gravel, xr_sugar)

#%% 
# =============================================================================
#          radiative profiles
# ============================================================================= 
    
#%% pre-process

def get_rad_data_one_day(day_str, all_profiles):

    # only choose HALO sondes 
    profiles_oneday = all_profiles.sel(launch_time=day_str) 
    HALO_profiles = profiles_oneday.where(profiles_oneday.platform=='HALO', drop=True)
    
    # =============================================================================
    #          LOAD data
    # ============================================================================= 
    q_rad = HALO_profiles['q_rad']
    q_rad_lw = HALO_profiles['q_rad_lw']
    q_rad_sw = HALO_profiles['q_rad_sw']
    launch_time = HALO_profiles['launch_time']
    alt_vec = HALO_profiles['zlev'][:-1].values

    # compress into xarray object    
    rad_day = xr.Dataset(
            data_vars={'Q_rad':    (('launch_time', 'height'), q_rad),
                       'Q_rad_lw':    (('launch_time', 'height'), q_rad_lw),
                       'Q_rad_sw':    (('launch_time', 'height'), q_rad_sw)},
            coords={'launch_time': launch_time,
                    'height': alt_vec})
    
    return rad_day

#%% load data

# =============================================================================
#             open files
# =============================================================================  
dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Data/"
path_to_rad_profiles = os.path.join(dir_profile,"rad_profiles_all_sondes_ERA.nc")
all_profiles = xr.open_dataset(path_to_rad_profiles)

#%% call function

day_str_fish = '2020-01-22' 
rad_fish = get_rad_data_one_day(day_str_fish, all_profiles)

day_str_flower = '2020-02-02' 
rad_flower = get_rad_data_one_day(day_str_flower, all_profiles)

day_str_gravel = '2020-02-05' 
rad_gravel = get_rad_data_one_day(day_str_gravel, all_profiles)

day_str_sugar = '2020-02-09' 
rad_sugar = get_rad_data_one_day(day_str_sugar, all_profiles)

#%% subplots

def plot_rad_mean_profiles(rad_fish, rad_flower, rad_gravel, rad_sugar):
    
    fig, ax = plt.subplots(1,3,figsize=(20,8))
    
    ax[0].set_ylabel('Altitude (km)')  
    ax[0].set_title('Shortwave')
    ax[1].set_title('Longwave')
    ax[2].set_title('Net')
    ax[1].set_xlabel('Mean heating rates / K/day')
    
    xmin=-6
    xmax=5
    ymin=0.1
    ymax=10
    for k in range(3):
        ax[k].grid(True, alpha=0.6)
        ax[k].set_xlim([xmin,xmax]) 
        ax[k].set_ylim([ymin,ymax]) 
        ax[k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
        ax[k].spines['right'].set_visible(False)
        ax[k].spines['top'].set_visible(False)
        for axis in ['top','bottom','left','right']:
              ax[k].spines[axis].set_linewidth(1.3)
        
    Dates = ['fish', 'flower', 'gravel', 'sugar']
    var_vec = ["Q_rad_sw", "Q_rad_lw", "Q_rad"]
    
    # =============================================================================
    #     SW 
    # =============================================================================
    var = var_vec[0]
    
    fish_var = rad_fish[var].where(rad_fish["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    flower_var = rad_flower[var].where(rad_flower["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    gravel_var = rad_gravel[var].where(rad_gravel["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    sugar_var = rad_sugar[var].where(rad_sugar["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    km_vec = sugar_var.height.values / 1000

    ax[0].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[0].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[0].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[0].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    # =============================================================================
    #      LW
    # =============================================================================
    var = var_vec[1]
    fish_var = rad_fish[var]
    flower_var = rad_flower[var]
    gravel_var = rad_gravel[var]
    sugar_var = rad_sugar[var]
    
    ax[1].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[1].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[1].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[1].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    # =============================================================================
    #      Net
    # =============================================================================
    var = var_vec[2]
    fish_var = rad_fish[var]
    flower_var = rad_flower[var]
    gravel_var = rad_gravel[var]
    sugar_var = rad_sugar[var]
    
    ax[2].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[2].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[2].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[2].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    
    ax[1].tick_params(labelleft=False)  
    ax[2].tick_params(labelleft=False)
    
    fig.tight_layout() 
    
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    

#%% call function to plot
    
plot_rad_mean_profiles(rad_fish, rad_flower, rad_gravel, rad_sugar)

#%% plot standard deviation

def plot_radiative_sigma(rad_fish, rad_flower, rad_gravel, rad_sugar):
    

    fig, ax = plt.subplots(1,3,figsize=(20,8))
    ax[0].set_ylabel('Altitude (km)')  
    
    ax[0].set_title('Shortwave')
    ax[1].set_title('Longwave')
    ax[2].set_title('Net')
    
    ax[1].set_xlabel('Standard deviation / K/day')
    
    xmin=0
    xmax=5
    ymin=0.1
    ymax=10
    for k in range(3):
        #ax[k].grid(color='k', linestyle='-', linewidth=0.8)
        ax[k].grid(True, alpha=0.6)
        ax[k].set_xlim([xmin,xmax]) 
        ax[k].set_ylim([ymin,ymax]) 
        ax[k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
        ax[k].spines['right'].set_visible(False)
        ax[k].spines['top'].set_visible(False)
        for axis in ['top','bottom','left','right']:
              ax[k].spines[axis].set_linewidth(1.3)
        
    Dates = ['fish', 'flower', 'gravel', 'sugar']
    var_vec = ["Q_rad_sw", "Q_rad_lw", "Q_rad"]
    
    # =============================================================================
    #      SW
    # =============================================================================
    var = var_vec[0]
    km_vec = rad_fish[var].height.values / 1000

    fish_var = rad_fish[var].where(rad_fish["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    flower_var = rad_flower[var].where(rad_flower["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    gravel_var = rad_gravel[var].where(rad_gravel["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    sugar_var = rad_sugar[var].where(rad_sugar["Q_rad_sw"].mean(dim="height") > 0, drop=True)
    
    fish_spatial_variance = fish_var.std(dim="launch_time")
    flower_spatial_variance = flower_var.std(dim="launch_time")
    gravel_spatial_variance = gravel_var.std(dim="launch_time")
    sugar_spatial_variance = sugar_var.std(dim="launch_time")

    ax[0].plot(fish_spatial_variance,km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[0].plot(flower_spatial_variance,km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[0].plot(gravel_spatial_variance, km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[0].plot(sugar_spatial_variance, km_vec,linewidth=5, color='lightgrey', label= Dates[3])
    #ax[0].legend(loc='best')
    
    # =============================================================================
    #      LW
    # =============================================================================
    var = var_vec[1]
    fish_spatial_variance = rad_fish[var].std(dim="launch_time")
    flower_spatial_variance = rad_flower[var].std(dim="launch_time")
    gravel_spatial_variance = rad_gravel[var].std(dim="launch_time")
    sugar_spatial_variance = rad_sugar[var].std(dim="launch_time")
    
    ax[1].plot(fish_spatial_variance,km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[1].plot(flower_spatial_variance,km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[1].plot(gravel_spatial_variance, km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[1].plot(sugar_spatial_variance, km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    # =============================================================================
    #      net
    # =============================================================================
    var = var_vec[2]
    
    fish_spatial_variance = rad_fish[var].std(dim="launch_time")
    flower_spatial_variance = rad_flower[var].std(dim="launch_time")
    gravel_spatial_variance = rad_gravel[var].std(dim="launch_time")
    sugar_spatial_variance = rad_sugar[var].std(dim="launch_time")
    
    ax[2].plot(fish_spatial_variance,km_vec,linewidth=5, color='blue', label= Dates[0])
    ax[2].plot(flower_spatial_variance,km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    ax[2].plot(gravel_spatial_variance, km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    ax[2].plot(sugar_spatial_variance, km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    
    ax[1].tick_params(labelleft=False)  
    ax[2].tick_params(labelleft=False)
    
    fig.tight_layout() 
    
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    

#%%
plot_radiative_sigma(rad_fish, rad_flower, rad_gravel, rad_sugar)
