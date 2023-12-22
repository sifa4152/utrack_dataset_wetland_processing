from utils import params_creator
import glob
import os
import sys
import settings


if __name__ == '__main__':
    files = glob.glob(os.path.join(settings.PATH_TARGET_ZONES, 'target_cells_*.csv')) # target cell coordinates (lat,lon)
    params_creator(files, STEP_SIZE=50) # adjust chunk size to your needs