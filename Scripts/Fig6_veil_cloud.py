#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A.L. Albright

Figure 6. Plot environmental variables and radiative heating rates for a persistent 'veil cloud' 
occurring on 24 January 2020

"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

plt.rcParams.update({"font.size": 24})
import os
import seaborn as sns

sns.set(
    context="notebook",
    style="white",
    palette="deep",
    font="sans-serif",
    font_scale=3,
    color_codes=True,
    rc=None,
)

#%%
# =============================================================================
#            pre-processing functions
# =============================================================================


def choose_HALO_circle(sondes):
    HALO_sondes = sondes.where(sondes.Platform == "HALO")
    # HALO circle, source: http://eurec4a.eu/fileadmin/user_upload/eurec4a/documents/Flight-Planning-Coordinates.pdf
    lon_center, lat_center = -57.717, 13.3
    lon_pt_circle, lat_pt_circle = -57.245, 14.1903
    r_circle = np.sqrt(
        (lon_pt_circle - lon_center) ** 2 + (lat_pt_circle - lat_center) ** 2
    )
    lat_N = lat_center + r_circle
    lat_S = lat_center - r_circle
    lon_W = lon_center - r_circle
    lon_E = lon_center + r_circle
    sondes_circle = HALO_sondes.where(
        (HALO_sondes["lat"] < lat_N)
        & (HALO_sondes["lat"] > lat_S)
        & (HALO_sondes["lon"] < lon_E)
        & (HALO_sondes["lon"] > lon_W),
        drop=True,
    )
    nsondes_HALO_circle = len(sondes_circle.launch_time)
    print(nsondes_HALO_circle, "total sondes launched in HALO circle")
    return sondes_circle


# =============================================================================
#          load environmental data
# =============================================================================

def get_env_data_one_day(fp_dropsondes, day_str):
    all_sondes = xr.open_dataset(fp_dropsondes).swap_dims({"sounding": "launch_time"})
    sondes_oneday = all_sondes.sel(launch_time=day_str)
    Sondes_Circle = choose_HALO_circle(sondes_oneday)

    # remove sondes with large nunber of NaNs
    sondes_circle = Sondes_Circle.dropna(
        dim="launch_time", subset=["q"], how="any", thresh=300
    )

    nsondes_qc = len(sondes_circle["launch_time"])
    print(nsondes_qc, "sondes after quality control on " + day_str)

    kgtog = 1000
    CtoK = 273.15

    # load specific humidity, convert from kg/kg to g/kg
    if sondes_circle["q"].max().values > 10:
        specific_humidity_g_kg = sondes_circle["q"] * kgtog
    else:
        specific_humidity_g_kg = sondes_circle["q"]

    # convert temperature from C to K
    if sondes_circle["T"].max().values < 100:
        temp_K = sondes_circle["T"] + CtoK
    else:
        temp_K = sondes_circle["T"]

    # convert RH to % units
    RH = sondes_circle["rh"]
    if sondes_circle["rh"].max().values < 1:
        RH = sondes_circle["rh"] * 100
    else:
        RH = sondes_circle["rh"]

    launch_time = sondes_circle["launch_time"]
    alt_vec = sondes_circle["height"]

    # save as xarray Dataset
    xr_day = xr.Dataset(
        data_vars={
            "specific_humidity": (("launch_time", "height"), specific_humidity_g_kg),
            "temp_K": (("launch_time", "height"), temp_K),
            "RH": (("launch_time", "height"), RH),
        },
        coords={"launch_time": launch_time, "height": alt_vec},
    )

    return xr_day


def get_rad_data_one_day(fp_rad_profiles, day_str):

    all_rad_profiles = xr.open_dataset(fp_rad_profiles)
    # only choose HALO sondes
    profiles_oneday = all_rad_profiles.sel(launch_time=day_str)
    HALO_profiles = profiles_oneday.where(profiles_oneday.platform == "HALO", drop=True)
    
    q_rad = HALO_profiles["q_rad"]
    q_rad_lw = HALO_profiles["q_rad_lw"]
    q_rad_sw = HALO_profiles["q_rad_sw"]
    launch_time = HALO_profiles["launch_time"]
    alt_vec = HALO_profiles["zlev"][:-1].values

    # compress into xarray object
    rad_day = xr.Dataset(
        data_vars={
            "Q_rad": (("launch_time", "height"), q_rad),
            "Q_rad_lw": (("launch_time", "height"), q_rad_lw),
            "Q_rad_sw": (("launch_time", "height"), q_rad_sw),
        },
        coords={"launch_time": launch_time, "height": alt_vec},
    )

    return rad_day


#%%
# =============================================================================
# load Joanne dropsondes
# =============================================================================

# REPLACE PATHS BELOW

input_dir = "/Users/annaleaalbright/Dropbox/EUREC4A/Dropsondes/Data/"
fp_dropsondes = os.path.join(
    input_dir, "EUREC4A_JOANNE_Dropsonde-RD41_Level_3_v0.5.7-alpha+0.g45fe69d.dirty.nc"
)

# =============================================================================
#             radiative profiles
# =============================================================================
dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Data/"
fp_rad_profiles = os.path.join(dir_profile, "rad_profiles_all_sondes_ERA.nc")

#%%

day_veil_clouds = "2020-01-24"

# environmental conditions
env_veil = get_env_data_one_day(fp_dropsondes, day_veil_clouds)
env_veil_0 = env_veil.sel(launch_time=slice("2020-01-24T12:53", "2020-01-24T12:56"))
# env_veil_1 = env_veil.sel(launch_time = slice('2020-01-24T13:53', '2020-01-24T13:56'))
# env_veil_2 = env_veil.sel(launch_time = slice('2020-01-24T14:53', '2020-01-24T14:56'))
env_veil_3 = env_veil.sel(launch_time=slice("2020-01-24T15:53", "2020-01-24T15:56"))
# env_veil_concat = xr.concat([env_veil_0, env_veil_1, env_veil_2, env_veil_3], dim="launch_time")
env_veil_concat = xr.concat([env_veil_0, env_veil_3], dim="launch_time")

rad_veil = get_rad_data_one_day(fp_rad_profiles, day_veil_clouds)
rad_veil_0 = rad_veil.sel(launch_time=slice("2020-01-24T12:53", "2020-01-24T12:56"))
# rad_veil_1 = rad_veil.sel(launch_time = slice('2020-01-24T13:53', '2020-01-24T13:56'))
# rad_veil_2 = rad_veil.sel(launch_time = slice('2020-01-24T14:53', '2020-01-24T14:56'))
rad_veil_3 = rad_veil.sel(launch_time=slice("2020-01-24T15:53", "2020-01-24T15:56"))
# rad_veil_concat = xr.concat([rad_veil_0, rad_veil_1, rad_veil_2, rad_veil_3], dim="launch_time")
rad_veil_concat = xr.concat([rad_veil_0, rad_veil_3], dim="launch_time")
km_vec = rad_veil_concat.height.values / 1000


#%%
def plot_mean_environmental_veil(xr_day, xr_veil):
    """ 
    Plot mean temperature, specific humidity, and relative humidity for each day
    
    """
    fig, ax = plt.subplots(1, 3, figsize=(20, 8))

    ax[0].set_ylabel("Altitude (km)")
    ax[0].set_xlabel("Temperature (K)")
    ax[1].set_xlabel("Specific humidity (g/kg)")
    ax[2].set_xlabel("Relative humidity (%)")
    ax[1].set_title("Environmental profiles")

    var_vec = ["temp_K", "specific_humidity", "RH"]

    ymin = 0
    ymax = 10
    for k in range(3):
        ax[k].grid(True, alpha=0.7)
        ax[k].set_ylim([ymin, ymax])
        ax[k].tick_params(
            direction="in", bottom=True, top=True, left=True, right=True, grid_alpha=0.7
        )
        ax[k].spines["right"].set_visible(False)
        ax[k].spines["top"].set_visible(False)
        for axis in ["top", "bottom", "left", "right"]:
            ax[k].spines[axis].set_linewidth(1.3)

    # =============================================================================
    #      temperature
    # =============================================================================
    var = var_vec[0]
    var1 = xr_veil[var]
    var1_interp = var1.interpolate_na(dim="height", method="linear")

    km_vec = var1.height.values / 1000

    colors1 = iter(["royalblue", "navy"])

    for i in range(0, len(var1_interp.launch_time)):
        ax[0].plot(
            var1_interp.isel(launch_time=i),
            km_vec,
            color=next(colors1),
            linewidth=5,
            alpha=0.7,
        )

    # =============================================================================
    #     specific humidity
    # =============================================================================
    var = var_vec[1]
    var2 = xr_veil[var]
    var2_interp = var2.interpolate_na(dim="height", method="linear")

    colors2 = iter(["royalblue", "navy"])

    for i in range(0, len(var2_interp.launch_time)):
        ax[1].plot(
            var2_interp.isel(launch_time=i),
            km_vec,
            color=next(colors2),
            linewidth=5,
            alpha=0.7,
        )

    ax[1].tick_params(labelleft=False)
    ax[0].set_ylim([0.1, 9])
    ax[1].set_ylim([0.1, 9])

    # =============================================================================
    #      relative humidity
    # =============================================================================
    var = var_vec[2]
    var3 = xr_veil[var]
    var3_interp = var3.interpolate_na(dim="height", method="linear")

    colors3 = iter(["royalblue", "navy"])

    for i in range(0, len(var3_interp.launch_time)):
        ax[2].plot(
            var3_interp.isel(launch_time=i),
            km_vec,
            color=next(colors3),
            linewidth=5,
            alpha=0.7,
        )

    ax[2].set_ylim([0.1, 9])

    fig.tight_layout()
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)


#%%

plot_mean_environmental_veil(env_veil, env_veil_concat)

#%%


def plot_rad_profiles_veil_cloud(xr_veil):
    """ 
    Plot mean temperature, specific humidity, and relative humidity for each day
    
    """
    fig, ax = plt.subplots(1, 3, figsize=(20, 8))

    ax[0].set_ylabel("Altitude (km)")
    ax[0].set_ylabel("Altitude (km)")
    ax[0].set_title("Shortwave")
    ax[1].set_title("Longwave")
    ax[2].set_title("Net")
    ax[1].set_xlabel("Heating rates (K/day)")

    var_vec = ["Q_rad_sw", "Q_rad_lw", "Q_rad"]

    ymin = 0
    ymax = 10
    xmin = -20
    xmax = 15
    for k in range(3):
        ax[k].grid(True, alpha=0.7)
        ax[k].set_ylim([ymin, ymax])
        ax[k].set_xlim([xmin, xmax])
        ax[k].tick_params(
            direction="in", bottom=True, top=True, left=True, right=True, grid_alpha=0.7
        )
        ax[k].spines["right"].set_visible(False)
        ax[k].spines["top"].set_visible(False)
        for axis in ["top", "bottom", "left", "right"]:
            ax[k].spines[axis].set_linewidth(1.3)

    # =============================================================================
    #      SW
    # =============================================================================
    var = var_vec[0]
    var1 = xr_veil[var]
    km_vec = var1.height.values / 1000

    colors1 = iter(["royalblue", "navy"])

    for i in range(0, len(var1.launch_time)):
        ax[0].plot(
            var1.isel(launch_time=i),
            km_vec,
            color=next(colors1),
            linewidth=5,
            alpha=0.7,
        )

    # =============================================================================
    #     LW
    # =============================================================================
    var = var_vec[1]
    var2 = xr_veil[var]

    colors2 = iter(["royalblue", "navy"])

    for i in range(0, len(var2.launch_time)):
        ax[1].plot(
            var2.isel(launch_time=i),
            km_vec,
            color=next(colors2),
            linewidth=5,
            alpha=0.7,
        )
    ax[1].tick_params(labelleft=False)
    ax[0].set_ylim([0.1, 9])
    ax[1].set_ylim([0.1, 9])

    # =============================================================================
    #      net
    # =============================================================================
    var = var_vec[2]
    var3 = xr_veil[var]

    colors3 = iter(["royalblue", "navy"])

    for i in range(0, len(var3.launch_time)):
        ax[2].plot(
            var3.isel(launch_time=i),
            km_vec,
            color=next(colors3),
            linewidth=5,
            alpha=0.7,
        )
    ax[2].set_ylim([0.1, 9])

    fig.tight_layout()
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)


#%% call function to plot

plot_rad_profiles_veil_cloud(rad_veil_concat)
