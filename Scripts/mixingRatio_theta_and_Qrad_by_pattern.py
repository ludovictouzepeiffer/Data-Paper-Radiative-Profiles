# -*- coding: utf-8 -*-

#!/usr/bin/env python
"""
# plot mixing ratio, theta, and LW / SW / net heating rates for four example days (one for each pattern)
@author: annaleaalbright
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
plt.rcParams.update({'font.size': 24})
import os

#%% 
# =============================================================================
#            Functions
# =============================================================================  

# choose only dropsondes in HALO circle
def choose_HALO_circle(sondes): 
    HALO_sondes = sondes.where(sondes.Platform=='HALO')
    # HALO circle, source: http://eurec4a.eu/fileadmin/user_upload/eurec4a/documents/Flight-Planning-Coordinates.pdf
    lon_center, lat_center = -57.717,13.3
    lon_pt_circle, lat_pt_circle = -57.245,14.1903
    r_circle = np.sqrt((lon_pt_circle-lon_center)**2+(lat_pt_circle-lat_center)**2)
    buffer = 0.01
    lat_N = lat_center + r_circle + buffer*r_circle
    lat_S = lat_center - r_circle - buffer*r_circle
    lon_W = lon_center - r_circle - buffer*r_circle
    lon_E = lon_center + r_circle + buffer*r_circle     
    sondes_circle = HALO_sondes.where((HALO_sondes['lat']<lat_N) & (HALO_sondes['lat']>lat_S)\
                                      & (HALO_sondes['lon']<lon_E) &(HALO_sondes['lon']>lon_W), drop=True)
    nsondes_HALO_circle = len(sondes_circle.launch_time)
    print(nsondes_HALO_circle, "total sondes launched in HALO circle")
    return sondes_circle

# select dropsonde data for one day with date string
def get_data_one_day(fp_dropsondes, day_str):
    all_sondes = xr.open_dataset(fp_dropsondes)
    sondes_oneday = all_sondes.sel(launch_time=day_str) 
    Sondes_Circle = choose_HALO_circle(sondes_oneday)
    
    # remove sondes with large nunber of NaNs
    # threshold : require this many non-NA values
    sondes_circle = Sondes_Circle.dropna(dim="launch_time", \
                                         subset=['alt','pres','u_wind','v_wind', 'wspd','lat','lon','mr', 
                                         'theta', 'theta_e', 'theta_v', 'rh', 'T', 'dz', 'q', 'dp'], \
                                         how='any', thresh=8_000) 
    
    nsondes_qc = len(sondes_circle.launch_time)
    print(nsondes_qc, "sondes after quality control on " + day_str)

    # =============================================================================
    #          LOAD environmental data
    # ============================================================================= 
    mixing_ratio = sondes_circle['mr']
    theta = sondes_circle['theta']
    temp = sondes_circle['T']
    RH = sondes_circle['rh']
    pressure = sondes_circle['pres']
    launch_time = sondes_circle['launch_time']
    alt_vec = sondes_circle['alt']

    # compress into xarray object    
    xr_day = xr.Dataset(
            data_vars={'mixing_ratio':    (('launch_time', 'alt'), mixing_ratio),
                       'theta':    (('launch_time', 'alt'), theta),
                       'temp':    (('launch_time', 'alt'), temp),
                       'RH':    (('launch_time', 'alt'), RH),
                       'pressure':    (('launch_time', 'alt'), pressure)},
            coords={'launch_time': launch_time,
                    'alt': alt_vec})
    
    return xr_day

#%%
    
# =============================================================================
# Load sondes in HALO circle + pre-process
# =============================================================================  

fp_dropsondes = '/Users/annaleaalbright/Dropbox/EUREC4A/Dropsondes/Data/all_sondes_w_cloud_flag.nc'

# 22 Jan (fish), 2 Feb (flowers), 5 Feb (gravel) and 9 Feb (sugar). 
day_str_fish = '2020-01-22' 
xr_fish = get_data_one_day(fp_dropsondes, day_str_fish)

day_str_flower = '2020-02-02' 
xr_flower = get_data_one_day(fp_dropsondes, day_str_flower)

day_str_gravel = '2020-02-05' 
xr_gravel = get_data_one_day(fp_dropsondes, day_str_gravel)

day_str_sugar = '2020-02-09' 
xr_sugar = get_data_one_day(fp_dropsondes, day_str_sugar)

#%%

def plot_vars_together(xr_fish, xr_flower, xr_gravel, xr_sugar):
    
    plt.rcParams.update({'font.size': 30})

    fig, ax = plt.subplots(1,2,figsize=(15,10))
    ax[0].set_ylabel('Altitude (km)')
    ax[0].set_xlabel('Mixing ratio / g/kg')
    ax[1].set_xlabel('Theta / K ')
    Dates = ['fish: 2020-01-22', 'flower: 2020-02-02', 'gravel: 2020-02-05', 'sugar: 2020-02-09']
    alpha_all = 0.05
    var_vec = ["mixing_ratio", "theta"]


    #xmin=-5
    #xmax=5
    ymin=0
    ymax=10
    for k in range(2):
        #ax[k].grid(color='k', linestyle='-', linewidth=0.8)
        ax[k].grid(True, alpha=0.7)
        #ax[k].set_xlim([xmin,xmax]) 
        ax[k].set_ylim([ymin,ymax]) 
        ax[k].tick_params(direction='in', bottom=True, top=True, left=True, right=True,grid_alpha=0.6)
        ax[k].spines['right'].set_visible(False)
        ax[k].spines['top'].set_visible(False)
        for axis in ['top','bottom','left','right']:
              ax[k].spines[axis].set_linewidth(1.3)
    

    # mixing ratio
    var = var_vec[0]
    var1 = xr_fish[var]
    var2 = xr_flower[var]
    var3 = xr_gravel[var]
    var4 = xr_sugar[var]
    km_vec = var4.alt.values / 1000

        
    for i in range(len(var1.launch_time)):
        ax[0].plot(var1.isel(launch_time=i), km_vec,color="blue", linewidth=2,alpha=alpha_all)
    ax[0].plot(var1.mean(dim='launch_time'), km_vec,linewidth=5, color='blue', label= Dates[0])
    for i in range(len(var2.launch_time)):
        ax[0].plot(var2.isel(launch_time=i), km_vec,color="palevioletred", linewidth=2,alpha=alpha_all + 0.1)
    ax[0].plot(var2.mean(dim='launch_time'), km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    for i in range(len(var3.launch_time)):
        ax[0].plot(var3.isel(launch_time=i), km_vec,color="deepskyblue", linewidth=2,alpha=alpha_all)
    ax[0].plot(var3.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    for i in range(len(var4.launch_time)):
        ax[0].plot(var4.isel(launch_time=i), km_vec,color="lightgrey", linewidth=2,alpha=alpha_all+0.05)
    ax[0].plot(var4.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])
    #ax[0].legend(loc='best')

    
    #ax[0].gca().spines['right'].set_visible(False)
    #ax[0].gca().spines['top'].set_visible(False)

    # theta
    var = var_vec[1]
    var1 = xr_fish[var]
    var2 = xr_flower[var]
    var3 = xr_gravel[var]
    var4 = xr_sugar[var]
    km_vec = var4.alt.values / 1000

        
    for i in range(len(var1.launch_time)):
        ax[1].plot(var1.isel(launch_time=i), km_vec,color="blue", linewidth=2,alpha=alpha_all)
    ax[1].plot(var1.mean(dim='launch_time'), km_vec,linewidth=5, color='blue', label= Dates[0])
    for i in range(len(var2.launch_time)):
        ax[1].plot(var2.isel(launch_time=i), km_vec,color="palevioletred", linewidth=2,alpha=alpha_all + 0.1)
    ax[1].plot(var2.mean(dim='launch_time'), km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    for i in range(len(var3.launch_time)):
        ax[1].plot(var3.isel(launch_time=i), km_vec,color="deepskyblue", linewidth=2,alpha=alpha_all)
    ax[1].plot(var3.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    for i in range(len(var4.launch_time)):
        ax[1].plot(var4.isel(launch_time=i), km_vec,color="lightgrey", linewidth=2,alpha=alpha_all+0.05)
    ax[1].plot(var4.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])
    #ax[1].legend(loc='best')
    
    ax[1].tick_params(labelleft=False)
    
    fig.tight_layout() 
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    
#%% call function

plot_vars_together(xr_fish, xr_flower, xr_gravel, xr_sugar)

#%% 
# =============================================================================
#          radiative profiles
# ============================================================================= 
    
#%% 
def get_rad_data_one_day(day_str, all_profiles):

    profiles_oneday = all_profiles.sel(launch_time=day_str) 

    
    # =============================================================================
    #          LOAD data
    # ============================================================================= 
    q_rad = profiles_oneday['q_rad']
    q_rad_lw = profiles_oneday['q_rad_lw']
    q_rad_sw = profiles_oneday['q_rad_sw']
    launch_time = profiles_oneday['launch_time']
    alt_vec = profiles_oneday['zlev'][:-1].values

    # compress into xarray object    
    rad_day = xr.Dataset(
            data_vars={'Q_rad':    (('launch_time', 'alt'), q_rad),
                       'Q_rad_lw':    (('launch_time', 'alt'), q_rad_lw),
                       'Q_rad_sw':    (('launch_time', 'alt'), q_rad_sw)},
            coords={'launch_time': launch_time,
                    'alt': alt_vec})
    
    return rad_day

#%% load data

dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/Dropsondes/Data/"
path_to_dropsonde_profiles = os.path.join(dir_profile,"rad_profiles_all_dropsondes.nc")
path_to_radiosonde_profiles = os.path.join(dir_profile, "rad_profiles_all_radiosondes.nc")

dropsonde_profiles = xr.open_dataset(path_to_dropsonde_profiles)
radiosonde_profiles = xr.open_dataset(path_to_radiosonde_profiles)
all_profiles = xr.concat([dropsonde_profiles, radiosonde_profiles], dim="launch_time")

#%% call function

# 22 Jan (fish), 2 Feb (flower), 5 Feb (gravel) and 9 Feb (sugar). 
day_str_fish = '2020-01-22' 
rad_fish = get_rad_data_one_day(day_str_fish, all_profiles)

day_str_flower = '2020-02-02' 
rad_flower = get_rad_data_one_day(day_str_flower, all_profiles)

day_str_gravel = '2020-02-05' 
rad_gravel = get_rad_data_one_day(day_str_gravel, all_profiles)

day_str_sugar = '2020-02-09' 
rad_sugar = get_rad_data_one_day(day_str_sugar, all_profiles)


#%% subplots

def plot_rad_rates(rad_fish, rad_flower, rad_gravel, rad_sugar):
    
    plt.rcParams.update({'font.size': 30})

    fig, ax = plt.subplots(1,3,figsize=(20,10))
    ax[0].set_ylabel('Altitude (km)')  
    
    ax[0].set_title('Shortwave')
    ax[1].set_title('Longwave')
    ax[2].set_title('Net')
    
    ax[1].set_xlabel('Heating rate (K/day)')
    
    xmin=-5
    xmax=5
    ymin=0
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
        
    #Dates = ['fish: 2020-01-22', 'flower: 2020-02-02', 'gravel: 2020-02-02', 'sugar: 2020-02-09']
    Dates = ['fish', 'flower', 'gravel', 'sugar']

    var_vec = ["Q_rad_sw", "Q_rad_lw", "Q_rad"]
    alpha_all = 0.02
    
    # SW 
    var = var_vec[0]
    fish_var = rad_fish[var]
    flower_var = rad_flower[var]
    gravel_var = rad_gravel[var]
    sugar_var = rad_sugar[var]
    km_vec = sugar_var.alt.values / 1000

    for i in range(0, len(fish_var.launch_time)):
        ax[0].plot(fish_var.isel(launch_time=i), km_vec,color="blue", linewidth=2,alpha=alpha_all)
    ax[0].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    for i in range(0, len(flower_var.launch_time)):
        ax[0].plot(flower_var.isel(launch_time=i), km_vec,color="palevioletred", linewidth=2,alpha=alpha_all)
    ax[0].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    for i in range(0, len(gravel_var.launch_time)):
        ax[0].plot(gravel_var.isel(launch_time=i), km_vec,color="deepskyblue", linewidth=2,alpha=alpha_all)
    ax[0].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    for i in range(0, len(sugar_var.launch_time)):
        ax[0].plot(sugar_var.isel(launch_time=i), km_vec,color="lightgrey", linewidth=2,alpha=alpha_all)
    ax[0].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])
    #ax[0].legend(loc='best')
    # LW
    var = var_vec[1]
    fish_var = rad_fish[var]
    flower_var = rad_flower[var]
    gravel_var = rad_gravel[var]
    sugar_var = rad_sugar[var]
    
    # fish mediumslateblue
    for i in range(len(fish_var.launch_time)):
        ax[1].plot(fish_var.isel(launch_time=i), km_vec,color="blue", linewidth=2,alpha=alpha_all)
    ax[1].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    for i in range(len(flower_var.launch_time)):
        ax[1].plot(flower_var.isel(launch_time=i), km_vec,color="palevioletred", linewidth=2,alpha=alpha_all)
    ax[1].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    for i in range(len(gravel_var.launch_time)):
        ax[1].plot(gravel_var.isel(launch_time=i), km_vec,color="deepskyblue", linewidth=2,alpha=alpha_all)
    ax[1].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    for i in range(len(sugar_var.launch_time)):
        ax[1].plot(sugar_var.isel(launch_time=i), km_vec,color="lightgrey", linewidth=2,alpha=alpha_all)
    ax[1].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    # Net
    var = var_vec[2]
    fish_var = rad_fish[var]
    flower_var = rad_flower[var]
    gravel_var = rad_gravel[var]
    sugar_var = rad_sugar[var]
    
    for i in range(len(fish_var.launch_time)):
        ax[2].plot(fish_var.isel(launch_time=i), km_vec,color="blue", linewidth=2,alpha=alpha_all)
    ax[2].plot(fish_var.mean(dim='launch_time'),km_vec,linewidth=5, color='blue', label= Dates[0])
    for i in range(len(flower_var.launch_time)):
        ax[2].plot(flower_var.isel(launch_time=i), km_vec,color="palevioletred", linewidth=2,alpha=alpha_all)
    ax[2].plot(flower_var.mean(dim='launch_time'),km_vec,linewidth=5, color='palevioletred', label= Dates[1])
    for i in range(len(gravel_var.launch_time)):
        ax[2].plot(gravel_var.isel(launch_time=i), km_vec,color="deepskyblue", linewidth=2,alpha=alpha_all)
    ax[2].plot(gravel_var.mean(dim='launch_time'), km_vec,linewidth=5, color='deepskyblue', label= Dates[2])
    for i in range(len(sugar_var.launch_time)):
        ax[2].plot(sugar_var.isel(launch_time=i), km_vec,color="lightgrey", linewidth=2,alpha=alpha_all)
    ax[2].plot(sugar_var.mean(dim='launch_time'), km_vec,linewidth=5, color='lightgrey', label= Dates[3])

    
    ax[1].tick_params(labelleft=False)  
    ax[2].tick_params(labelleft=False)
    
    fig.tight_layout() 
    
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    

#%%
plot_rad_rates(rad_fish, rad_flower, rad_gravel, rad_sugar)
