# =============================================================================
# ANTHROMES LAND-COVER EVAPORATION PARTITIONING
# Simon Felix Fahrländer
# March 2022; updated March 2023
# GE9011 thesis project; Wetland Vulnerability Paper 
# Stockholm University; Potsdam Insititue for Climate Impact Research
# =============================================================================
import xarray as xr
import numpy as np 
import pandas as pd 
import warnings
import time
import sys
import settings
import os 
from tqdm import tqdm
# =============================================================================
## Function that returns index of closest latitude in the tracking array to the one provided by the tracking zone
def get_closest_index(lats,lat):
        import operator
        lat_index, min_value = min(enumerate(abs(lats-lat)), key=operator.itemgetter(1))
        return lat_index

# =============================================================================
# MASK SCREENING: --> lists masked locations
# Screening function: takes mask array with range (nan-1) and gives out locs[lon,lat] of cells containing values.
# input: array-mask 
# output: DataFrame with lon,lat of masked cells
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

# =============================================================================
# function to sum evaporation over land-cover categories for each precipitationshes. 40 psheds with 6 categories each! 
def cat_precip_sum(zone_ID, save=True):
        # start = time.time()
        idx = zone_ID
# =============================================================================
        cat1 = [11,12] # Dense settlement
        cat2 = [21,22,31] # Rice and irrigated cropland
        cat3 = [23,32,33,34] # Rainfed cropland
        cat4 = [24,41,42,43] # Rangeland 
        cat5 = [51,52,53,61] # Woodland
        cat6 = [54,62] # Barrenland
        #cat7 = [np.nan] # 'NaNs'
        #cat_list = [cat1,cat2,cat3,cat4,cat5,cat6,cat7] # create list for reformatting loop
        # id_list = ['Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean'] # list of value names
# =============================================================================
        # cat1_ls = [] # Dense settlement
        # cat2_ls = [] # Rice and irrigated cropland
        # cat3_ls = [] # Rainfed cropland
        # cat4_ls = [] # Rangeland 
        # cat5_ls = [] # Woodland
        # cat6_ls = [] # Barrenland
        # cat7_ls = [] #'NaNs' --> can be approximated as oceanic source cells
        # cat8_ls = [] # Ocean
        # mtrack_ls = [] # 'tracked moisture'
# =============================================================================
        anthromes = xr.open_dataset(settings.ANTHROMES_PATH, engine='netcdf4').Band1
        
        # # Use of land-sea-mask to differentiate between ocean and nans. BUT coastal areas don't have an overlap'
        lsm = xr.open_dataset(settings.LSM_PATH, engine='netcdf4').lsm.squeeze('time')
        lsm = xr.where(lsm>0,1,0)
        # anthromes = original_anthromes.where(lsm==0,1)

        cell_area = xr.open_dataset(settings.GRID_PATH, engine='netcdf4').cell_area

        target_zone = pd.read_csv(os.path.join(settings.PATH_TARGET_ZONES,f'target_cells_{idx}.csv'))
        target_cells = len(target_zone)
# =============================================================================
    #for idx in range(len(zone_ID)):
        # Load mask and screen for lat,lon 
        track_fp = os.path.join(settings.OUT_PSHED,f'pshed70_{idx}.nc') # 70% psheds
        Evap_footprint_annual_sum = xr.open_dataset(track_fp, engine='netcdf4')['70%_pshed_values'] 
       
        #mask = Evap_footprint_annual_sum
        #mask = mask.where((mask>0),np.nan)
        # print(f'start screening zone {idx}...')
        screened_data = screening(Evap_footprint_annual_sum)
# =============================================================================
        # tracker for evaporation from land-cover categories:
        tracker1 = 0 # Dense settlement
        tracker2 = 0 # Rice and irrigated cropland
        tracker3 = 0 # Rainfed cropland
        tracker4 = 0 # Rangeland 
        tracker5 = 0 # Woodland
        tracker6 = 0 # Barrenland
        tracker7 = 0 # 'NaNs'
        tracker8 = 0 # Ocean
        
        # cell counters: 
        cells1 = 0
        cells2 = 0
        cells3 = 0
        cells4 = 0
        cells5 = 0
        cells6 = 0
        cells7 = 0
        cells8 = 0
        total_cells = len(screened_data)
        warnings.filterwarnings("ignore")
# =============================================================================
        # print(f'start processing zone {idx}: {len(screened_data)} cells...')
        for i in range(len(screened_data)):
            ## Create lat,lon arrays equal to dimensions of tracking dataset
            lats=np.arange(90,-90,-0.5)
            lons=np.arange(0,360,0.5)
    
            latidx=get_closest_index(lats,screened_data.Lat[i])
            lonidx=get_closest_index(lons,screened_data.Lon[i])

            if anthromes[latidx,lonidx].values in cat1:
                tracker1 = np.nansum(tracker1 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells1 += cell_area[latidx,lonidx].values
            elif anthromes[latidx,lonidx].values in cat2:
                tracker2 = np.nansum(tracker2 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells2 += cell_area[latidx,lonidx].values
            elif anthromes[latidx,lonidx].values in cat3:
                tracker3 = np.nansum(tracker3 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells3 += cell_area[latidx,lonidx].values
            elif anthromes[latidx,lonidx].values in cat4:
                tracker4 = np.nansum(tracker4 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells4 += cell_area[latidx,lonidx].values
            elif anthromes[latidx,lonidx].values in cat5:
                tracker5 = np.nansum(tracker5 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells5 += cell_area[latidx,lonidx].values
            elif anthromes[latidx,lonidx].values in cat6:
                tracker6 = np.nansum(tracker6 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells6 += cell_area[latidx,lonidx].values
            elif (np.isnan(anthromes[latidx,lonidx].values)) and (lsm[latidx,lonidx].values == 1):
                tracker7 = np.nansum(tracker7 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells7 += cell_area[latidx,lonidx].values
            else:
                tracker8 = np.nansum(tracker8 + Evap_footprint_annual_sum[latidx,lonidx].values)
                cells8 += cell_area[latidx,lonidx].values
# =============================================================================
        ## Sum up backtracking footprint
        mtrack = np.nansum(Evap_footprint_annual_sum)
        total_area = np.array([cells1,cells2,cells3,cells4,cells5,cells6,cells7,cells8]).sum()
        ## correct spatial component of mm unit
        lc_type_P = np.array([tracker1,tracker2,tracker3,tracker4,tracker5,tracker6,tracker7,tracker8,mtrack])/target_cells

        lc_type_P = np.append(idx,lc_type_P)
        ## m2 to km2
        lc_type_A = np.array([cells1,cells2,cells3,cells4,cells5,cells6,cells7,cells8,total_area])/1_000_000
        lc_type_A = np.append(idx,lc_type_A)
# =============================================================================
        df1 = pd.DataFrame(lc_type_P)
        df1 = df1.transpose()
        df1.columns = ['ramsarid','Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean','total']
        df2 = pd.DataFrame(lc_type_A)
        df2 = df2.transpose()
        df2.columns = ['ramsarid','Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean','total']
        
        frames = [df1,df2]
        df = pd.concat(frames)
        if save==True:
            df.to_csv(os.path.join(settings.OUT_LC,f'categorical_precipitation_source_{idx}.csv'))
            print(f'done with zone {idx}')
# =============================================================================
        # end = time.time()
        # elapse = round((end - start)/3600,2)
        # print(elapse)
        else:
            print(f'done with zone {idx}')
            return(df)

# =============================================================================
## function to unpack results and merge to single dfs for precipitation volume and land-cover type area
def unpack_MP_results(zone_IDs,labels):
    lc_vol = [] # volume of evaporation from land-cover categories
    lc_area = [] # number of cells of each land-cover category in psheds
    # =============================================================================
    # loop over dfs and draw out respective columns and save in lists
    for i,idx in tqdm(enumerate(zone_IDs)):
        filepath = os.path.join(settings.OUT_LC,f'categorical_precipitation_source_{idx}.csv')
        df = pd.read_csv(filepath, encoding='ISO-8859-1', index_col=0,dtype=np.float64)
        f1 = pd.DataFrame(df.iloc[0,:]).transpose()
        f2 = pd.DataFrame(df.iloc[1,:]).transpose()
        lc_vol.append(f1)
        lc_area.append(f2)
    
    df_vol = pd.concat(lc_vol,ignore_index=True)
    df_vol['name'] = labels

    df_area = pd.concat(lc_area,ignore_index=True)
    df_area['name'] = labels

    # return lists 
    return df_vol, df_area
# =============================================================================
# Calculate percentage of evaporation from each LC in each pshed
def calc_percentage(data):
    var_list = ['Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean'] # variable list for loop #,'NaNs'
    
    tracker = {} # tracker dictionnary to organize output DataFrame
    tracker['ID'] = data['ramsarid'] # Ramsar IDs
    tracker['name'] = data['name'] # wetland abbreviated names
    
    for var in var_list:
        x = data[var].astype('float')
        y = data['total'].astype('float')
        z = x/y # count of variable i divided by sum of all cells --> percentage
        tracker[f'{var}'] = z # save in dictionnary
# =============================================================================
    frame = pd.DataFrame(tracker) # create DateFrame from dictionnary 
    percentages = frame[['ID','name','Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean']] # reorder colums #,'NaNs'
    
    check = percentages[['Dense settlements','Rice and irrigated cropland','Rainfed cropland', 'Rangeland','Woodland','Barrenland','NaNs','Ocean']].sum(axis=1) #,'NaNs'
    
    return percentages, check

# =============================================================================
## TEST:
if __name__ == '__main__':
    zone = 32
    a = cat_precip_sum(zone,save=False)
