# =============================================================================
## UTrack tracking post-processing script 
## prequisite --> assembled runs (assemble_runs.py)
## calculates recycling indices and atmospheric watersheds
# =============================================================================
import numpy as np
import pandas as pd
import xarray as xr
import rasterio as rio
import geopandas as gpd
import os
from tqdm import tqdm
# =============================================================================
## function to define the cells that fall into the X-pshed
## input: target percent of pshed, 100% precipitationshed as DataArray
## output: last counted cell falling into X-pshed
def pshed_limit_counter(percent, pshed):
    pshed100 = np.nansum(pshed) # total precipitation amount
    percentile = percent/100
    pshedX = pshed100*percentile # X% precipitation
    
    # flatten and sort pshed-array 
    pshed100_copy = np.copy(pshed)
    pshed100_flat = pshed100_copy.flatten()
    pshed100_sort = np.sort(pshed100_flat)
    pshed100_sort = np.flip(pshed100_sort) # flattend and sorted from high to low values
    
    i = 0 
    psum = 0
    while psum < pshedX: 
        psum = psum + pshed100_sort[i]
        i = i + 1
        #print(i) 
    return i
# =============================================================================
## function that creates pshed-mask as matrix of 0-1 
## input: precipitationshed as DataArray, output from pshed_limit_counter()
## output: mask of X-pshed [0-1]
def pshed_matrix(pshed, last_cell):
    
    pshed_copy = np.copy(pshed)
    pshed_flat = pshed_copy.flatten()
    pshed_sort = np.sort(pshed_flat)
    pshed_sort = np.flip(pshed_sort) 
    
    pshedX_matrix = np.copy(pshed_copy)
    
    lat = pshed.lat
    row = np.arange(0,len(lat),1)
    lon = pshed.lon
    col = np.arange(0,len(lon),1)
    
    for j in row: 
        for k in col:
            if pshedX_matrix[j,k] >= pshed_sort[last_cell]:
                pshedX_matrix[j,k] = 1
            else:
                pshedX_matrix[j,k] = float('nan')
    
    return pshedX_matrix
# =============================================================================
## function that uses the above two funtions to create and extract X-pshed mask and values 
## input: file_path to 100% pshed, percent X
## output: mask and value array of X-pshed --> (mask,values)
## direction: backward --> Pshed, forward --> Eshed
## SET TRANSFORM ARGUMENT HERE! 
def percentile_pshed(input_fp, percent, direction,transform=False):
    ## load evaporation footprint
    if direction == 'backward':
        shed = xr.open_dataset(input_fp)['bw_evap_footprint_monthly']
    elif direction == 'forward':
        shed = xr.open_dataset(input_fp)['fw_evap_footprint_monthly']
    ## aggregate to annual average
    shed = shed.sum(dim='month', skipna=True)
    ## define last counted cell for shed
    last_cell = pshed_limit_counter(percent, shed)
    ## create shed matrix
    shedX_matrix = pshed_matrix(shed, last_cell)
    ## create mask (0-1)
    if direction == 'backward':
        shedX_mask = xr.DataArray(shedX_matrix, coords=[shed.lat.values, shed.lon.values],
                    dims=['lat', 'lon'], name = f'{percent}%_pshed_mask',
                    attrs=dict(description="Mask of Percentile Precipitationshed", units="[0-1]"))

    if direction == 'forward':
        shedX_mask = xr.DataArray(shedX_matrix, coords=[shed.lat.values, shed.lon.values],
                    dims=['lat', 'lon'], name = f'{percent}%_eshed_mask',
                    attrs=dict(description="Mask of Percentile Evaporationtionshed", units="[0-1]"))
    ## create value array
    shedX_values = shed*shedX_mask

    if direction == 'backward':
        shedX_values = shedX_values.rename(f'{percent}%_pshed_values')
    
    if direction == 'forward':
        shedX_values = shedX_values.rename(f'{percent}%_eshed_values')

    ## transform coordinates from UTrack to (lon:(-180,180), lat(-89.5,90)) 
    if transform == True:
        # reset coordinates
        shedX_mask.coords['lon'] = (shedX_mask.coords['lon'] + 180) % 360 - 180
        shedX_mask = shedX_mask.sortby(shedX_mask.lon)
        shedX_mask = shedX_mask.reindex(lat=list(reversed(shedX_mask.lat)))
        # shedX_mask.rio.write_crs("epsg:4326", inplace=True)

        shedX_values.coords['lon'] = (shedX_values.coords['lon'] + 180) % 360 - 180
        shedX_values = shedX_values.sortby(shedX_values.lon)
        shedX_values = shedX_values.reindex(lat=list(reversed(shedX_values.lat)))
        # shedX_values.rio.write_crs("epsg:4326", inplace=True)
    ## combine mask and value arrays to dataset
    combined_ds = xr.merge([shedX_values, shedX_mask])

    return combined_ds #shedX_mask, shedX_values
# =============================================================================
## CALCULATE ATMOSPHERIC WATERSHEDS FUNCTION
## INPUT: list or df of ids, percent of total moisture accounted for in watersheds (e.g. 70 or 90) 
def atmos_watersheds(ID,percent):
    ## Precipitation- and Evaporationsheds
    print(f'calculating pshed & eshed of zone {ID} --> num of zones: {len(ID)}')
    
    for count, idx in enumerate(tqdm(ID)):
        try:
            bw_path = os.path.join(settings.BACKWARD,f'backward_footprint_{idx}.nc')
        except FileNotFoundError:
            print(f'missing bw_footprint: {idx} --> moving on...')    
        try:  
            fw_path = os.path.join(settings.FORWARD,f'forward_footprint_{idx}.nc')
        except FileNotFoundError:
            print(f'missing fw_footprint: {idx} --> moving on...')

        pshed = percentile_pshed(bw_path, percent, direction = 'backward', transform=False)
        eshed = percentile_pshed(fw_path, percent, direction = 'forward', transform=False)

        pshed.to_netcdf(os.path.join(settings.OUT_PSHED, f'pshed{percent}_{idx}.nc'), engine='netcdf4')
        eshed.to_netcdf(os.path.join(settings.OUT_ESHED, f'eshed{percent}_{idx}.nc'), engine='netcdf4')

        # print(f'zone: {idx} --> [{count+1}/{len(ID)}]')
    print('done!')
# =============================================================================
## function to calculate internal and terrestrial precipitation recycling
## zone_identifier array with [ID,name]!! e.g. [FID,country]
def moisture_recycling(zone_identifier):
    
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    ID = zone_identifier.ID
    name = zone_identifier.name
    
    ## Create lists for saving in loop
    cell_list = [] # count of masked cells
    
    ## PRECIPITATION QUANTITIES IN KM3 PER YEAR
    ptotal_list = [] # general moisture import in km3/yr
    int_Precy_list = [] # amount of internal precipitation recycling in km3//yr
    terr_Precy_list = [] # total moisture from terrestrial source in km3//yr
    ext_Precy_list = [] # moisture import from terrestiral source in km3//yr
    pocean_list = [] # moisture from oceanic sources in km3/yr 

    ## PRECIPITATION QUANTITIES IN MM PER YEAR
    ptotal_mm_list = [] # general moisture import in mm/yr
    int_Precy_mm_list = [] # amount of internal precipitation recycling in mm/yr
    terr_Precy_mm_list = [] # moisture import from terrestrial source in mm/yr
    ext_Precy_mm_list = [] # moisture import from terrestiral source in mm/yr
    pocean_mm_list = [] # moisture from oceanic sources in mm/yr 

    ## PRECIPITATION RECYCLING RATIOS
    int_Pratio_list = [] # internal precipitation recycling ratio
    terr_Pratio_list = [] # terrestrial precipitation recycling ratio
    ext_Pratio_list = [] # external terrestrial precipitation recycling

    ## EVAPORATION QUANTITIES IN KM3 PER YEAR
    etotal_list = [] # general moisture export in km3/yr
    int_Erecy_list = [] # amount of internal evaporation recycling in km3/yr
    terr_Erecy_list = [] # total moisture to terrestrial sink in km3/yr
    ext_Erecy_list = [] # moisture export to terrestrial sink in km3/yr
    eocean_list = [] # moisture export to oceanic sinks in km3/yr 

    ## EVAPORATION QUANTITIES IN MM PER YEAR
    etotal_mm_list = [] # general moisture export in mm/yr
    int_Erecy_mm_list = [] # amount of internal evaporation recycling in mm/yr
    terr_Erecy_mm_list = [] # moisture export from terrestrial source in mm/yr
    ext_Erecy_mm_list = [] # moisture export to terrestrial sink in mm/yr
    eocean_mm_list = [] # moisture export to oceanic sinks in mm/yr

    ## EVAPORATION RECYCLING RATIOS
    int_Eratio_list = [] # internal evaporation recycling ratio
    terr_Eratio_list = [] # terrestrial evaporation recycling ratio
    ext_Eratio_list = [] # external terrestrial evaporation recycling

    ## EAR5 PRECIPITATION 
    era5_total_list = []
    era5_total_mm_list = []
    p_error_list = []

    ## Load land-sea-mask for terrestrical recycling and source/sink zones
    grid_area = xr.open_dataset(settings.GRID_PATH, engine='netcdf4').cell_area
    lsm = xr.open_dataset(settings.LSM_PATH, engine='netcdf4').lsm.squeeze('time')
    lsm = xr.where(lsm>0,1,np.nan)

    ## ERA5 precipitation
    p = xr.open_dataset(settings.PRECIP_PATH).tp.sel(time=slice('2008-01-01','2017-12-01'))
    p = p.resample(time='Y').sum(skipna=True)
    p = p.mean(dim='time', skipna=True)

    # ## STEAM precipitation
    # p = xr.open_dataset(settings.PRECIP_PATH).p_cur
    # p = p.rename({'month':'time'})
    # p = p.resample(time='Y').sum(skipna=True)
    # p = p.mean(dim='time', skipna=True)

    ## loop over zones, calculates imported moisture for each zone
    for count, (idx, n) in enumerate(tqdm(zip(ID,name))):
        ## COUNTRY MASK
        mask = xr.open_dataset(os.path.join(settings.PATH_mask,f'b_remap{idx}.nc'), engine='netcdf4').Band1
        cell_num = np.nansum(mask)
        ## COUNTRY AREA IN M^2
        country_area = np.nansum(mask*grid_area)

        ## LOAD BACKWARD FOOTPRINT FOR ZONE IDX
        p_track_fp = os.path.join(settings.BACKWARD,f'backward_footprint_{idx}.nc')
        try:
            BW_footprint_monthly = xr.open_dataset(p_track_fp, engine='netcdf4').bw_evap_footprint_monthly
        except FileNotFoundError:
            print(f'missing bw_footprint: {idx} --> {n}')
            ID = ID.drop(count)
            name = name.drop(count)
            continue
        BW_footprint_annual_sum = BW_footprint_monthly.sum(dim='month',skipna=True)

        ## LOAD FORWARD FOOTPRINT FOR ZONE IDX
        e_track_fp = os.path.join(settings.FORWARD,f'forward_footprint_{idx}.nc')
        try:
            FW_footprint_monthly = xr.open_dataset(e_track_fp, engine='netcdf4').fw_evap_footprint_monthly
        except FileNotFoundError:
            print(f'missing fw_footprint: {idx} --> {n}')
            ID = ID.drop(count)
            name = name.drop(count)
            continue
        FW_footprint_annual_sum = FW_footprint_monthly.sum(dim='month',skipna=True)

        ## TOTAL IMPORTED MOISTURE  => Sum up bacKward footprint
        ptotal_agg = np.nansum(BW_footprint_annual_sum) # total imported moisture in aggregated mm/yr 
        try:
            ptotal_mm = (ptotal_agg/cell_num) # total imported moisture in mm/yr
        except ZeroDivisionError:
            ptotal_mm = np.nan
        ptotal = (ptotal_mm*country_area* 10**(-12)) # total imported moisture in km3/yr 

        ## TOTAL EXPORTED MOISTURE  => Sum up backward footprint
        etotal_agg = np.nansum(FW_footprint_annual_sum) # total exported moisture in aggregated mm/yr
        try:
            etotal_mm = (etotal_agg/cell_num) # total exported moisture in mm/yr
        except ZeroDivisionError:
            etotal_mm = np.nan
        etotal = (etotal_mm*country_area* 10**(-12)) # total exported moisture in km3/yr

        ## PRECIPITATION FROM TERRESTRIAL SOURCES
        terr_BW_footprint_annual_sum = BW_footprint_annual_sum * lsm
        terr_Precycling_agg = np.nansum(terr_BW_footprint_annual_sum) # footprint for masked area originating from land surface in aggregated mm/yr
        try:
            terr_Precycling_mm = (terr_Precycling_agg/cell_num) # footprint for masked area originating from land surface in mm/yr
        except ZeroDivisionError:
            terr_Precycling_mm = np.nan
        terr_Precycling = (terr_Precycling_mm*country_area* 10**(-12)) # footprint for masked area originating from land surface in km3/yr

        ## EVAPORATION TO TERRESTRIAL SINKS
        terr_FW_footprint_annual_sum = FW_footprint_annual_sum * lsm
        terr_Erecycling_agg = np.nansum(terr_FW_footprint_annual_sum) # footprint for masked area going to land surface in aggregated mm/yr
        try:
            terr_Erecycling_mm = (terr_Erecycling_agg/cell_num) # footprint for masked area going to land surface in mm/yr
        except ZeroDivisionError:
            terr_Erecycling_mm = np.nan
        terr_Erecycling = (terr_Erecycling_mm*country_area* 10**(-12)) # footprint for masked area going to land surface in km3/yr

        ## PRECIPITATION FROM INTERNAL SOURCES: --> Sums of precipitation and recycled moisture in masked area and sum of terrestrial moisture
        intP = BW_footprint_annual_sum*mask
        intE = FW_footprint_annual_sum*mask

        int_Precycling_agg = np.nansum(intP) # recycled precip within mask area in aggregated mm/yr
        try:
            int_Precycling_mm = int_Precycling_agg/cell_num # recycled precip within mask area in mm/yr
        except: 
            int_Precycling_mm = np.nan
        int_Precycling = (int_Precycling_mm*country_area* 10**(-12)) # recycled precip within mask area in km3/yr

        int_Erecycling_agg = np.nansum(intE) # recycled evap within mask area in aggregated mm/yr
        try:
            int_Erecycling_mm = int_Erecycling_agg/cell_num # recycled evap within mask area in mm/yr
        except ZeroDivisionError:
            int_Erecycling_mm = np.nan
        int_Erecycling = (int_Erecycling_mm*country_area* 10**(-12)) # recycled evap within mask area in km3/yr

        ## INTERNAL PRECIPITATION RECYCLING RATIO
        try:
            int_p_ratio = int_Precycling/ptotal
        except ZeroDivisionError:
            int_p_ratio = np.nan

        ## INTERNAL EVAPORATION RECYCLING RATIO
        try:
            int_e_ratio = int_Erecycling/etotal
        except ZeroDivisionError:
            int_e_ratio = np.nan

        ## TERRESTRIAL PRECIPITATION RECYCLING RATIO
        try:
            terr_p_ratio = terr_Precycling/ptotal
        except ZeroDivisionError:
            terr_p_ratio = np.nan

        ## TERRESTRIAL EVAPORATION RECYCLING RATIO
        try:
            terr_e_ratio = terr_Erecycling/etotal
        except ZeroDivisionError:
            terr_e_ratio = np.nan

        ## EXTERNAL TERRESTRIAL PRECIPITATION RECYCLING
        ext_terr_p_ratio = (terr_Precycling - int_Precycling)/ptotal # ratio of P that originates from terrestrial sources outside of sink area
        ext_terr_P_agg = (terr_Precycling_agg - int_Precycling_agg) # volume of P from external terrestrial sources in aggregated mm/yr
        ext_terr_P = (terr_Precycling - int_Precycling) # volume of P from external terrestrial sources in km3/yr
        ext_terr_P_mm = (ext_terr_P_agg/cell_num) # volume of P from external terrestrial sources in mm/yr

        ## EXTERNAL TERRESTRIAL PRECIPITATION RECYCLING
        ext_terr_e_ratio = (terr_Erecycling - int_Erecycling)/etotal # ratio of E that precipitates over terrestrial land outside of source area
        ext_terr_E_agg = (terr_Erecycling_agg - int_Erecycling_agg) # volume of E going to external terrestrial sources in aggregated mm/yr
        ext_terr_E = (terr_Erecycling - int_Erecycling) # volume of E going to external terrestrial sinks in km3/yr
        ext_terr_E_mm = (ext_terr_E_agg/cell_num) # volume of E going to external terrestrial sinks in mm/yr

        ## PRECIPITATION FROM OCEANIC SOURCES
        P_ocean_km3 = ptotal - terr_Precycling # precipitation coming from ocean evaporation in km3/yr
        P_ocean_mm = ptotal_mm - terr_Precycling_mm # precipitation coming from ocean evaporation in mm/yr

        ## EVAPORATION GOING TO THE OCEANS
        E_ocean_km3 = etotal - terr_Erecycling # evaporation going to the ocean in km3/yr
        E_ocean_mm = etotal_mm - terr_Erecycling_mm # evaporation going to the ocean in mm/yr

        ## ERA5 PRECIPITATION --> validation and control instance! 
        era5_p_agg =  np.nansum((p*mask))
        era5_p_mm = (era5_p_agg/cell_num)
        era5_p_km3 = (era5_p_mm*country_area*10**(-12))
        era5_p = np.nansum(p*mask)
        # p_error = ptotal/era5_p_km3 # error ratio between model based and data based precipitation (1:optimal, <1:underestimation, >1: overestimation)
        p_error = np.nansum(BW_footprint_annual_sum)/era5_p

        ## SAVE CALCULATED VALUES IN LISTS
        cell_list.append(cell_num)

        ## SAVING PRECIPITATION INDICES
        ptotal_list.append(ptotal) # P_total in km3/yr
        int_Precy_list.append(int_Precycling) # P_E_internal in km3/yr
        terr_Precy_list.append(terr_Precycling) # P_E_terrestrial in km3/yr
        ext_Precy_list.append(ext_terr_P) # P from external terrestiral sources in km3/yr

        ptotal_mm_list.append(ptotal_mm) # P_total in mm/yr
        int_Precy_mm_list.append(int_Precycling_mm) # P_E_internal in mm/yr
        terr_Precy_mm_list.append(terr_Precycling_mm) # P_E_terrestrial in mm/yr
        ext_Precy_mm_list.append(ext_terr_P_mm) # P from external terrestiral sources in mm/yr

        int_Pratio_list.append(int_p_ratio) # internal precipitation recycling ratio
        terr_Pratio_list.append(terr_p_ratio) # terrestrial precipitation recycling ratio
        ext_Pratio_list.append(ext_terr_p_ratio) # external terrestrial p ratio

        pocean_list.append(P_ocean_km3) # P from oceanic sources in km3
        pocean_mm_list.append(P_ocean_mm) # P from oceanic source in mm

        ## SAVING EVAPORATION INDICES
        etotal_list.append(etotal) # E_total in km3/yr
        int_Erecy_list.append(int_Erecycling) # E_P_internal in km3/yr
        terr_Erecy_list.append(terr_Erecycling) # E_P_terrestrial in km3/yr
        ext_Erecy_list.append(ext_terr_E) # E to external terrestrial sinks in km3/yr

        etotal_mm_list.append(etotal_mm) # E_total in mm/yr
        int_Erecy_mm_list.append(int_Erecycling_mm) # E_P_internal in mm/yr
        terr_Erecy_mm_list.append(terr_Erecycling_mm) # E_P_terrestrial in mm/yr
        ext_Erecy_mm_list.append(ext_terr_E_mm) # E to external terrestrial sinks in mm/yr

        int_Eratio_list.append(int_e_ratio) # internal evaporation recycling ratio
        terr_Eratio_list.append(terr_e_ratio) # terrestrial evaporation recycling ratio
        ext_Eratio_list.append(ext_terr_e_ratio) # external terrestrial e ratio

        eocean_list.append(E_ocean_km3) # E export to oceanic sinks in km3
        eocean_mm_list.append(E_ocean_mm) # E export to oceanic sinks in mm

        ## SAVING CONTROL INSTANCES
        era5_total_list.append(era5_p_km3) # ERA5 precipitation in km3
        era5_total_mm_list.append(era5_p_mm) # ERA5 precipitation in mm
        p_error_list.append(p_error) # error of tracked moisture estimation
    
    ID_list = ID.to_list()
    name_list = name.to_list()
    ## Create DataFrame with lists
    df = pd.DataFrame({'zone_ID':ID_list,
                        'name':name_list,
                        'cell_num':cell_list,
                        'ERA5_P_km3':era5_total_list,
                        'ERA5_P_mm':era5_total_mm_list,
                        'P_error': p_error_list,

                        'P_total_km3':ptotal_list,
                        'P_terr_km3':terr_Precy_list,
                        'P_int_km3':int_Precy_list,
                        'P_ext_km3':ext_Precy_list,
                        'P_ocean_km3':pocean_list,
                        'P_total_mm':ptotal_mm_list,
                        'P_terr_mm':terr_Precy_mm_list,
                        'P_int_mm':int_Precy_mm_list,
                        'P_ext_mm':ext_Precy_mm_list,
                        'P_ocean_mm':pocean_mm_list,
                        'rho_int':int_Pratio_list,
                        'rho_terr':terr_Pratio_list,
                        'rho_ext':ext_Pratio_list,
                        
                        'E_total_km3':etotal_list,
                        'E_terr_km3':terr_Erecy_list,
                        'E_int_km3':int_Erecy_list,
                        'E_ext_km3':ext_Erecy_list,
                        'E_ocean_km3':eocean_list,
                        'E_total_mm':etotal_mm_list,
                        'E_terr_mm':terr_Erecy_mm_list,
                        'E_int_mm':int_Erecy_mm_list,
                        'E_ext_mm':ext_Erecy_mm_list,
                        'E_ocean_mm':eocean_mm_list,
                        'epsilon_int':int_Eratio_list,
                        'epsilon_terr':terr_Eratio_list,
                        'epsilon_ext':ext_Eratio_list})

    return df

if __name__ == '__main__':
    import sys
    import settings

    ## LOAD BASIN IDENTIFIERS FROM WORKSHEET
    df_wetlands = pd.read_csv(settings.PATH_IDs, encoding='ISO-8859-1')
    zone_ID = df_wetlands.ramsarid
    ID = pd.DataFrame(df_wetlands.ramsarid)
    name = pd.DataFrame(df_wetlands.name_short)
    df = pd.concat([ID,name], axis=1).rename(columns={'ramsarid':'ID','name_short':'name'})

    ## RUN CALCULATION FUNCTION AND SAVE OUTPUT
    x = moisture_recycling(df)
    x.to_csv(os.path.join(settings.OUT_PATH,'moisture_recycling_wetland_basins.csv'))

    ## DELINEATE ATMOSPHERIC WATERSHEDS
    params = pd.read_csv(settings.PARAMS, index_col=0)
    zone_id = params.zone.drop_duplicates().reset_index(drop=True)
    atmos_watersheds(zone_id,percent=70)
