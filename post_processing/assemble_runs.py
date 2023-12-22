import os
import glob
import pandas as pd
import numpy as np 
import xarray as xr 
import sys
import settings

OUT_FORWARD = os.path.join(settings.PATH_MOISTURE_FOOTPRINTS,'forward_complete') # output directory for forward footprints
OUT_BACKWARD = os.path.join(settings.PATH_MOISTURE_FOOTPRINTS,'backward_complete') # output directory for backward footprints

params = pd.read_csv(settings.PARAMS, index_col=0)
zone_id = params.zone.drop_duplicates().reset_index(drop=True)

for count,id in enumerate(zone_id):
    
    print(f'{id} -> [{count+1}/{len(zone_id)}]')
    
    fw_files = glob.glob(os.path.join(settings.OUT_ESHEDS,f'forward_footprint_{id}_*.nc'))
    bw_files = glob.glob(os.path.join(settings.OUT_PSHEDS,f'backward_footprint_{id}_*.nc'))

    if not fw_files:
        continue

    fw_sum = np.zeros(shape=(12,360,720))
    fw_sum_final = np.zeros(shape=(12,360,720))

    bw_sum = np.zeros(shape=(12,360,720))
    bw_sum_final = np.zeros(shape=(12,360,720))

    for i,j in zip(fw_files,bw_files):
        fw = xr.open_dataset(i, engine='netcdf4')['fw_evap_footprint_monthly']
        bw = xr.open_dataset(j, engine='netcdf4')['bw_evap_footprint_monthly']
        
        fw = xr.where(fw==np.nan,0,fw)
        bw = xr.where(bw==np.nan,0,bw)

        fw_sum_final = fw_sum_final + fw.values
        bw_sum_final = bw_sum_final + bw.values

    fw_assembled = xr.DataArray(fw_sum_final, coords=[fw.month.values, fw.lat.values, fw.lon.values],
                dims=['month', 'lat', 'lon'], name = "fw_evap_footprint_monthly", attrs=dict(description="Forward Monthly Evaporation Footprint", units="mm/month"))
    bw_assembled = xr.DataArray(bw_sum_final, coords=[fw.month.values, fw.lat.values, fw.lon.values],
                dims=['month', 'lat', 'lon'], name = "bw_evap_footprint_monthly", attrs=dict(description="Backward Monthly Evaporation Footprint", units="mm/month"))

    fw_assembled.to_netcdf(os.path.join(OUT_FORWARD, f'forward_footprint_{id}.nc'), engine='netcdf4')
    bw_assembled.to_netcdf(os.path.join(OUT_BACKWARD, f'backward_footprint_{id}.nc'), engine='netcdf4')
