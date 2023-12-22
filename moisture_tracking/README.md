# UTrack Database Moisture Tracking Routine 
## Simon F. Fahrländer, Feb 2023
## Potsdam Institute for Climate Impact Research - (PIK)


#### =================================================================
### Scripts: 
#### utils.py --> functions to screen masks for coordinates and writing out job packages for core splits 
#### create_params.py --> script to create parameter file with job packages 
#### settings.py --> WKDIRs!
#### utrack_functions.py --> Tracking Functions 
#### main.py --> job runner
#### utrack_slurm.sh --> send jobs to HPC
#### =================================================================


### Example start of routine: 

Activate conda environment containing all necessary modules

`conda activate geo_ext`

Create parameters

`python3 create_params.py`

Test single job runner locally

`python3 main.py 1`

Run all on one CPU or split up accordingly (e.g., on an HPC or server via the parameters) 

`python3 main.py´

#### =================================================================
