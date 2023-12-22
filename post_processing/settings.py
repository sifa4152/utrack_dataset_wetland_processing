import os

LSM_PATH = os.path.join() # path to land-sea mask
GRID_PATH = os.path.join() # path to grid containing grid cell area (m2)
ANTHROMES_PATH = os.path.join() # path to anthromes dataset (land cover dataset)

## HYDROCLIM DATA --> adjust according to your data needs to be in same format as UTrack dataset
PRECIP_PATH = os.path.join() 
EVAP_PATH = os.path.join() 

## PARAMETERS
PARAMS = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','moisture_tracking','basins_fpxP_era5','params.csv') # stays the same!
PATH_mask = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','masks','basin_masks','b_mask_netcdf') # 
PATH_IDs = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','wetlands_selection_csv.csv')
PATH_TARGET_ZONES = os.path.join('/p','projects','open','simon','Wetland_vulnerability_paper','data','target_cells')

## FOOTPRINTS
PATH_MOISTURE_FOOTPRINTS = os.path.join() # directory with moisture footprints from tracking runs
FORWARD = os.path.join(PATH_MOISTURE_FOOTPRINTS, 'forward_complete')
BACKWARD = os.path.join(PATH_MOISTURE_FOOTPRINTS, 'backward_complete')

## WRITE OUT RESULTS
OUT = os.path.join() # folder for post-processing output
OUT_PATH = os.path.join(OUT,'recycling')
OUT_ESHED = os.path.join(OUT,'esheds')
OUT_PSHED = os.path.join(OUT,'psheds')
OUT_LC = os.path.join(OUT,'land_cover')


