# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import pandas as pd
from numpy import *
import os
import sys
import settings

## Function that returns index of closest latitude in the tracking array to the one provided by the tracking zone
def get_closest_index(lats,lat):
        import operator
        lat_index, min_value = min(enumerate(abs(lats-lat)), key=operator.itemgetter(1))
        return lat_index

## Get monthly moisture footprints (forwards & backwards)
def track_footprints(month, latitude, longitude, evap, precip):
    
    ## Create lat,lon arrays equal to dimensions of tracking dataset
    lats=np.arange(90,-90,-0.5)
    lons=np.arange(0,360,0.5)
    
    latidx=get_closest_index(lats,latitude)
    lonidx=get_closest_index(lons,longitude)
    
    ## UTrack atmospheric moisture trajectory dataset
    UTrack = xr.open_dataset(os.path.join(settings.PATH_UTRACK, 'utrack_climatology_0.5_'+ str(month).zfill(2)+'.nc')).moisture_flow

    ## Forward tracking footprint
    fp_fw = UTrack[latidx, lonidx,:,:].values
    fp_fw = fp_fw * -0.1
    fp_fw = e**fp_fw
    fp_fw[fp_fw == 1] = 0
    forward_fp = fp_fw / np.nansum(fp_fw)
    
    ET_source = (evap[month-1])[latidx, lonidx].values
    forward_fp = forward_fp * ET_source
    
    ## Backward tracking footprint
    fp_bw = UTrack[:,:,latidx, lonidx].values
    fp_bw = fp_bw * -0.1
    fp_bw = e**fp_bw
    fp_bw[fp_bw == 1] = 0
    ET = (evap[month-1]).values
    # backward_fp = fp_bw * ET
        
    fp_bw = ET * fp_bw
    backward_fp = fp_bw / np.nansum(fp_bw)
    
    P_sink = (precip[month-1])[latidx,lonidx].values
    backward_fp = backward_fp * P_sink

    return forward_fp, backward_fp

## Runner Function: Takes DataFrame with lons,lats and aggregates footprints for each sink/source cell given by (lon,lat) for each month (1-12)
## output: --> Datasets containing monthly moisture footprints
def moisture_tracking_runner(Screened_data):
    
    evap = xr.open_dataset(os.path.join(settings.PATH_DATA_ERA5_E)).e.sel(time=slice('2008-01-01','2017-12-01')) 
    ## Convert to multi-year mean
    evap = evap.groupby('time.month').mean(dim = 'time')
    evap = (xr.where(evap > 0, evap, 0))

    precip = xr.open_dataset(os.path.join(settings.PATH_DATA_ERA5_P)).tp.sel(time=slice('2008-01-01','2017-12-01')) 
    ## Convert to multi-year mean
    precip = precip.groupby('time.month').mean(dim = 'time')
    precip = (xr.where(precip > 0, precip, 0))

    Forward_footprint_monthly_final  = np.zeros(shape=(12,360,720))
    Backward_footprint_monthly_final  = np.zeros(shape=(12,360,720))
    
    for j in range(Screened_data.shape[0]): # iterate over source/sink cells given by df(lats,lons)
        latitude, longitude = np.array(Screened_data.loc[j])
        
        Forward_footprint_monthly = np.zeros(shape=(12,360,720))
        Backward_footprint_monthly = np.zeros(shape=(12,360,720))

        for i in range(12): # iterate over months (1:12)
            Forward_footprint_monthly[i,:,:], Backward_footprint_monthly[i,:,:] = track_footprints(i+1,latitude,longitude,evap,precip)
            
        Forward_footprint_monthly = np.where(np.isnan(Forward_footprint_monthly),0,Forward_footprint_monthly)
        Backward_footprint_monthly = np.where(np.isnan(Backward_footprint_monthly),0,Backward_footprint_monthly)
        
        Forward_footprint_monthly_final = Forward_footprint_monthly_final + Forward_footprint_monthly
        Backward_footprint_monthly_final = Backward_footprint_monthly_final + Backward_footprint_monthly

    ## write monthly footprints into dataset
    Forward_footprint_monthly_sum = xr.DataArray(Forward_footprint_monthly_final, coords=[evap.month.values, evap.lat.values, evap.lon.values],
                 dims=['month', 'lat', 'lon'], name = "fw_evap_footprint_monthly", attrs=dict(description="Forward Monthly Evaporation Footprint", units="mm/month"))
    
    Backward_footprint_monthly_sum = xr.DataArray(Backward_footprint_monthly_final, coords=[evap.month.values, evap.lat.values, evap.lon.values],
                 dims=['month', 'lat', 'lon'], name = "bw_evap_footprint_monthly", attrs=dict(description="Backward Monthly Evaporation Footprint", units="mm/month"))
    
    return Forward_footprint_monthly_sum, Backward_footprint_monthly_sum

## TESTS
if __name__ == '__main__':
    ## TEST :
    x = 1234 # PLACEHOLDER
