# Post-processing of moisture footprints of 40 Ramsar wetlands
## Simon F. Fahrländer, Feb 2023
## Potsdam Institute for Climate Impact Research - (PIK)

#### =================================================================
### Scripts: 
#### settings.py --> WKDIRs!
#### assemble_runs.py --> function that adds up moisture footprints of the same regions based on the different processing batches (params.csv)
#### process_functions_updated.py --> calculates moisture recycling indices 
#### PminusE_trend.py --> calculates Mann-Kendall trends of P-E data for selected regions 
#### LC_Psource_functions.py --> connects moisture sources to land cover classification
#### =================================================================


### Example start of routine: 

Assemble process batches 
`python3 assemble_runs.py´

Calculate moisture recycling indices 
`python3 process_functions.py' 

Calculate Mann-Kendall trends
`python3 PminusE_trend.py´

Connect land cover to moisture sources --> needs to be looped for multiple regions or run individually for each region 
`python3 LC_Psource_functions.py´

#### =================================================================
