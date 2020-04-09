# by Ludovic Touze-Peiffer
# thanks to Hauke Schulz for showing his initial method

import urllib.request
import datetime
from calendar import monthrange
import os
import argparse

if __name__ == "__main__":

    ##-- import movie parameters

    from movie_params import *

    ##-- load arguments if present

    # Arguments to be used if want to change options while executing script
    parser = argparse.ArgumentParser(description="Downloads MODIS images")
    parser.add_argument("--date", type=str, default=date_str,help="Flight date, YYY-MM-DD")
    parser.add_argument("--sat_id",type=str,default=satid,
        help='SAT variable ID')
    args = parser.parse_args()

    # Define ouput path
    path_dir = satdir
    os.makedirs(path_dir,exist_ok=True)

    # Define dates
    start_date = datetime.datetime.strptime(args.date+start_time,"%Y-%m-%d%H:%M")


    str_date = start_date.strftime("%Y-%m-%dT%H:%M:%SZ")
    lon_lat = str(latmin)+','+str(lonmin)+ ',' + str(latmax)+','+str(lonmax)

    # define download url
    url = ('https://wvs.earthdata.nasa.gov/api/v1/snapshot?'+
    'REQUEST=GetSnapshot&TIME='+
    str_date+
    '&BBOX='+
    lon_lat+
    '&CRS=EPSG:4326&LAYERS='+
    args.sat_id+
    '&WRAP=day&FORMAT=image/jpeg&WIDTH='+
    str(width)+
    '&HEIGHT='+
    str(height)+
    '&ts=1585928448079')

    # save
    save_str = os.path.join(path_dir,'%s.jpg'%str_date)
    urllib.request.urlretrieve(url, save_str)


    
