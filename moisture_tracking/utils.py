# -*- coding: utf-8 -*-
import os
import sys
import glob
import settings
import numpy as np
import pandas as pd

## Screening function to isolate tracking zone coordinates in mask
def screening(mask):
    lat = []
    lon = []
    #screened_data = []
    for y in range(mask.shape[0]):
        for x in range(mask.shape[1]):
            if np.isnan(mask[y,x].values) == True:
                continue
            else:
                lat.append(round(float(mask[y,x].lat.values),2))
                lon.append(round(float(mask[y,x].lon.values),2))
                #screened_data.append(mask[y,x].values)
    
    # Extraxt the data as dataframe
    screened_data = pd.DataFrame({'Lat': lat,'Lon': lon}) #,'Screened': Screened_data
    return screened_data

## function that creates job packages for utrack cluster runs --> saved in .csv
def params_creator(files, STEP_SIZE):
    def get_ID(zone_ID):
        return zone_ID.get('zone')

    params = []
    for file in files:
        zone_id = int(file.split('_')[-1].split('.')[0])
        df = pd.read_csv(file, index_col=0)
        current_row = 0
        while current_row + STEP_SIZE < len(df):
            params.append({'zone': zone_id, 'start': current_row, 'stop': current_row + STEP_SIZE})
            current_row += STEP_SIZE
        params.append({'zone': zone_id, 'start': current_row, 'stop': len(df)})

    params.sort(key=get_ID)
    pd.DataFrame(params).to_csv(os.path.join(settings.WRKDIR,'params.csv'))
    
    print('############')
    print(f'Number of needed jobs: {len(params)}')
    print('Dont forget to adjust the slurm script before submission!')
    print('############')

if __name__ == '__main__':
    
    PATH_mask = os.path.join() # netCDF mask of tracking zones [x(zone_IDs),nan]
    PATH_IDs = os.path.join() # csv file with zone IDs that correspond to the mask 
    PATH_OUT = os.path.join() # output path for screened data --> target cells path needs to be added to settings.py

    df_wetlands = pd.read_csv(PATH_IDs, encoding='ISO-8859-1')
    zone_ID = df_wetlands.ramsarid

    import xarray as xr
    from multiprocessing import Pool
    def screen_func(ID):
        print(f'processing zone: {ID}')
        mask = xr.open_dataset(os.path.join(PATH_mask,f'b_remap{ID}.nc')).Band1
        
        screened_data = screening(mask)
        screened_data.to_csv(os.path.join(PATH_OUT,f'target_cells_{ID}.csv'))
    
    print('Start parallel processing')
    KERNEL = 7
    with Pool(KERNEL) as pool:
        pool.map(screen_func, zone_ID)    
