import numpy as np
import matplotlib.pyplot as plt


# I/O directories
satdir ='../output/TERRA_Modis/'
satid="MODIS_Terra_CorrectedReflectance_TrueColor"

# HALO circle
lon_center, lat_center = -57.717,13.3
lon_pt_circle, lat_pt_circle = -57.245,14.1903
r_circle = np.sqrt((lon_pt_circle-lon_center)**2+(lat_pt_circle-lat_center)**2)

# Image box
lonmin,lonmax = -60,-55
dlon = lonmin-lonmax
latmin,latmax = 11.5,15
dlat = latmin-latmax
width = 10000
height = int(width*dlat/dlon)

date_str="2020-02-02"
start_time="00:00"
