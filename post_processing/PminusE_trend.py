# =============================================================================
# Hydroclimatic Trend Analysis on wetland specific zones w Mann-Kendall-Test 
# Simon Felix Fahrländer
# =============================================================================
import xarray as xr
import pandas as pd
import numpy as np
import pymannkendall as mk
import warnings
from tqdm import tqdm
import os 
import sys
sys.path.insert(1,'/p/projects/open/simon/Wetland_vulnerability_paper/post_processing')
import settings
from importlib import reload # reload
reload(settings)
# =============================================================================
def zonal_mean(data,zone_ID):
    
    inf_values = np.isinf(data)
    inf_values = xr.where(inf_values==True,np.nan,1)
    data = data*inf_values
    
    ds_list1 = []
    ds_list2 = []
    ds_list3 = []
    
    time = data.time.values
    years = data.time.dt.year.values
    df = pd.DataFrame({'time':years})
    ds_list1.append(df)
    ds_list2.append(df)
    ds_list3.append(df)
    
    for i, idx in enumerate(tqdm(zone_ID)): 
        
        basin_mask = xr.open_dataset(os.path.join(settings.PATH_mask, f'b_remap{idx}.nc')).Band1
        
        new_lat = basin_mask.lat.values 
        new_lon = basin_mask.lon.values
        pshed_mask = xr.open_dataset(os.path.join(settings.OUT_PSHED, f'pshed70_{idx}.nc'))['70%_pshed_mask']
        pshed_mask = pshed_mask.reindex(lat=pshed_mask.lat[::-1])
        pshed_mask = pshed_mask.interp(lat=new_lat,lon=new_lon,method='nearest')
        
        import_mask = xr.where(basin_mask==1,np.nan,pshed_mask)
        
        runner1 = data*basin_mask
        runner2 = data*pshed_mask
        runner3 = data*import_mask
        
        tracker1 = []
        tracker2 = []
        tracker3 = []
        
        for t in time:
            x1 = runner1.sel(time=t).mean(skipna=True).values
            x2 = runner2.sel(time=t).mean(skipna=True).values
            x3 = runner3.sel(time=t).mean(skipna=True).values
            
            tracker1.append(x1)
            tracker2.append(x2)
            tracker3.append(x3)
            
            f1 = pd.DataFrame({f'{idx}':tracker1})
            f2 = pd.DataFrame({f'{idx}':tracker2})
            f3 = pd.DataFrame({f'{idx}':tracker3})
        
        ds_list1.append(f1)
        ds_list2.append(f2)
        ds_list3.append(f3)
        
        output_basin = pd.concat(ds_list1,axis=1)
        output_pshed = pd.concat(ds_list2,axis=1)
        output_import = pd.concat(ds_list3,axis=1)
        
    return(output_basin,output_pshed,output_import)
# =============================================================================
def mk_test(data, ignore_warnings=True):
    
    results_nested = []
    result_keys=['trend','h','p','z','Tau','s','var_s','slope','intercept']#
    f1 = pd.DataFrame({'keys':result_keys})
    results_nested.append(f1)
    
    for i in range(len(zone_ID)):
        
        if ignore_warnings:
            warnings.filterwarnings("ignore")
            
        result = mk.original_test(data[f'{zone_ID[i]}'])
        
        result_list = list(result)
        f2 = pd.DataFrame({f'{zone_ID[i]}':result_list})
        results_nested.append(f2)
        
        frame = pd.concat(results_nested, axis=1).transpose()
        
    return(frame)
# =============================================================================
def PminusE_runner(pre_cur,et_cur, pre_pv, et_pv): 
    # calculate P-E 
    PminusE_cur = pre_cur - et_cur
    PminusE_pv = pre_pv - et_pv
    # calculate zonal means
    basin_cur, pshed_cur, imp_zon_cur = zonal_mean(PminusE_cur, zone_ID)
    basin_pv, pshed_pv, imp_zon_pv = zonal_mean(PminusE_pv, zone_ID)
    
    basin = abs(basin_pv) - abs(basin_cur)
    pshed = abs(pshed_pv) - abs(pshed_cur)
    imp_zon = abs(imp_zon_pv) - abs(imp_zon_cur)
    
    # calculate Mann-Kendall Test statistics
    basin_test = mk_test(basin)
    pshed_test = mk_test(pshed)
    import_test = mk_test(imp_zon)

    # summary table 
    b_result = basin_test[[0,2,7]]
    b_result = b_result.tail(b_result.shape[0] -1)
    b_result.rename(columns={0:'trend_basin',2:'p_basin',7:'slope_basin'}, inplace=True)

    p_result = pshed_test[[0,2,7]]
    p_result = p_result.tail(p_result.shape[0] -1)
    p_result.rename(columns={0:'trend_pshed',2:'p_pshed',7:'slope_pshed'}, inplace=True)

    i_result = import_test[[0,2,7]]
    i_result = i_result.tail(i_result.shape[0] -1)
    i_result.rename(columns={0:'trend_import',2:'p_import',7:'slope_import'}, inplace=True)

    glued = pd.concat([b_result,p_result,i_result],axis=1)

    return glued
# =============================================================================
wetland_csv = pd.read_csv(settings.PATH_IDs, encoding='cp1252')
zone_ID = wetland_csv.ramsarid
name = wetland_csv.name_short
# =============================================================================
## STEAM PV P and E data
P_STEAM_PV_PATH = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','STEAM_evap_prec','STEAM_p_pv_utrack_fixedgrid_remap.nc') 
E_STEAM_PV_PATH = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','STEAM_evap_prec','STEAM_e_pv_utrack_fixedgrid_remap.nc')

pre_pv = xr.open_dataset(P_STEAM_PV_PATH).p_pv
pre_pv = pre_pv.resample(month='Y').sum(dim='month',skipna=True)
# pre_pv = pre_pv.mean(dim='month',skipna=True)
pre_pv = pre_pv.rename({'month':'time'})

et_pv = xr.open_dataset(E_STEAM_PV_PATH).e_pv
et_pv = et_pv.resample(month='Y').sum(dim='month',skipna=True)
# et_pv = et_pv.mean(dim='month',skipna=True)
et_pv = et_pv.rename({'month':'time'})
# =============================================================================
## STEAM CUR P and E data
P_STEAM_CUR_PATH = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','STEAM_evap_prec','STEAM_p_cur_utrack_fixedgrid_remap.nc') 
E_STEAM_CUR_PATH = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','STEAM_evap_prec','STEAM_e_cur_utrack_fixedgrid_remap.nc')
pre_cur = xr.open_dataset(P_STEAM_CUR_PATH).p_cur
pre_cur = pre_cur.resample(month='Y').sum(dim='month',skipna=True)
# pre_cur = pre_cur.mean(dim='month',skipna=True)
pre_cur = pre_cur.rename({'month':'time'})

et_cur = xr.open_dataset(E_STEAM_CUR_PATH).e_cur
et_cur = et_cur.resample(month='Y').sum(dim='month',skipna=True)
# et_cur = et_cur.mean(dim='month',skipna=True)
et_cur = et_cur.rename({'month':'time'})
# =============================================================================
## RUN FUNCTIONS
df = PminusE_runner(pre_cur, et_cur, pre_pv, et_pv)

## OUTPUT FILEPATH
DATA_PATH = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','post_processing','process_output','PminusE_cur_pv.csv')
df.to_csv(DATA_PATH)

# =============================================================================
# =============================================================================
# =============================================================================
## ERA5 AVERAGE P IN BASINS
PRECIP_PATH = os.path.join('/p','projects','open','simon','BGWater','CI_paper','data','ERA5','ERA5_precipitation','era5_monthly_tp_1979_2021_mm_month_remap.nc') 
pre_cru = xr.open_dataset(PRECIP_PATH,engine='netcdf4').tp
pre_cru = pre_cru.resample(time='Y').sum(skipna=True)
pre_basin, pre_pshed, pre_import = zonal_mean(pre_cru, zone_ID)
pre_basin_avg = pre_basin.mean()
pre_basin_avg = pre_basin_avg.iloc[1:].reset_index(drop=True)
df = pd.concat([zone_ID,name,pre_basin_avg],axis=1)
df.columns = ['ramsarid','name','P_mean_basin']
df.to_csv(os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','basin_mean_precip.csv'))
# =============================================================================
