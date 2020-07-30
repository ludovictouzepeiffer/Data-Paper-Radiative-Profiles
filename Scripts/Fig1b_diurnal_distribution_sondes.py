#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pytz
import os
import seaborn as sns
import pandas as pd
import matplotlib.dates as mdates

plt.rcParams.update({"font.size": 18})


#%%
# =============================================================================
#             open files
# =============================================================================

dir_profile = "/Users/annaleaalbright/Dropbox/EUREC4A/RadiativeProfiles/Data/"
path_to_rad_profiles = os.path.join(dir_profile, "rad_profiles_all_sondes_ERA.nc")
all_profiles = xr.open_dataset(path_to_rad_profiles)

#%%
# =============================================================================
# count number of sondes per platform / dropsondes vs. radiosondes
# =============================================================================

# dropsondes
HALO_sondes = all_profiles.where(all_profiles.platform == "HALO", drop=True)
P3_sondes = all_profiles.where(all_profiles.platform == "P3", drop=True)
all_dropsondes = xr.concat([HALO_sondes, P3_sondes], dim="launch_time")

# radiosondes
ATL_radiosondes = all_profiles.where(all_profiles.platform == "ATL", drop=True)
BCO_radiosondes = all_profiles.where(all_profiles.platform == "BCO", drop=True)
RHB_radiosondes = all_profiles.where(all_profiles.platform == "RHB", drop=True)
MER_radiosondes = all_profiles.where(all_profiles.platform == "MER", drop=True)
MET_radiosondes = all_profiles.where(all_profiles.platform == "MET", drop=True)
all_radiosondes = xr.concat(
    [
        ATL_radiosondes,
        BCO_radiosondes,
        RHB_radiosondes,
        MER_radiosondes,
        MET_radiosondes,
    ],
    dim="launch_time",
)

num_dropsondes = int(len(HALO_sondes.launch_time) + len(P3_sondes.launch_time))
num_radiosondes = int(
    len(BCO_radiosondes.launch_time)
    + len(RHB_radiosondes.launch_time)
    + len(MER_radiosondes.launch_time)
    + len(MET_radiosondes.launch_time)
    + len(ATL_radiosondes.launch_time)
)
print("total number of dropsondes:", num_dropsondes)
print("total number of radiosondes:", num_radiosondes)

#%%

# =============================================================================
#             PDF of sondes per hour, by platform
# =============================================================================


def calc_pdf_sondes(profiles):

    # select coordinates
    data = profiles["q_rad"]
    data["time"] = profiles["launch_time"]
    data = data.drop_vars(["lay", "col", "play"])
    data = data.to_dataframe()
    data["time"] = (
        data["time"]
        .dt.tz_localize(pytz.UTC)
        .dt.tz_convert("America/Barbados")
        .dt.strftime("%H:%M")
    )
    data["time"] = pd.to_datetime(data["time"], format="%H:%M")

    data = data.reset_index()
    data = data.set_index(["time", "zlay"])
    data = data.groupby(
        [pd.Grouper(freq="10min", level="time"), pd.Grouper(level="zlay")]
    ).count()
    data = data.groupby(["time"]).mean()

    # get time
    data = data.to_xarray()
    time = data.time.values

    ini = np.datetime64("1900-01-01 00:00:00")
    end = ini + np.timedelta64(24, "h")
    count_time = np.arange(ini, end, np.timedelta64(10, "m"))
    count = np.zeros(len(count_time))

    array = {"count_time": count_time, "count": count}
    array = pd.DataFrame(data=array)
    array = array.set_index(["count_time"])

    for itime in time:
        array.loc[itime, "count"] = data["q_rad"].sel(time=itime).values

    return count_time, array["count"].values


#%% call function

time, count_vec_dropsondes = calc_pdf_sondes(all_dropsondes)
time, count_vec_radiosondes = calc_pdf_sondes(all_radiosondes)

#%% plot
sns.set(
    context="notebook",
    style="whitegrid",
    palette="deep",
    font="sans-serif",
    font_scale=2,
    color_codes=True,
    rc=None,
)

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(1, 1, 1)
ax.plot(
    time,
    count_vec_radiosondes,
    linewidth=5,
    color="salmon",
    alpha=1,
    label=f"{num_radiosondes} radiosondes",
)
ax.plot(
    time,
    count_vec_dropsondes,
    linewidth=5,
    color="navy",
    alpha=1,
    label=f"{num_dropsondes} dropsondes",
)
myFmt = mdates.DateFormatter("%-H")
ax.xaxis.set_major_formatter(myFmt)
ticks = ax.get_xticks()
ax.set_xticks(np.linspace(ticks[0], mdates.date2num(mdates.num2date(ticks[-2])), 8))
ax.set_ylabel("number of sondes")
ax.set_xlabel("local time")
ax.grid(True, alpha=0.5)
ax.set_ylim([0, 120])
ax.legend(loc="best", frameon=False)
plt.minorticks_on()
plt.gca().tick_params(axis="y", which="minor", left=False)
plt.tight_layout()
sns.despine()
ax.autoscale_view()
