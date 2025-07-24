# wrapper script 
source("R/packages.R")
source("R/functions.R")

## Run on a single set of samples
data_dir = 'example_data/08202024/data'
sal_mass_file = 'example_data/Mass_8_20_2024.csv'
date_cal = '8/20/24'
cal_file = 'example_data/pHCalibration.csv'
output = 'example_data/08202024/08202024_results.csv'

run_titrator(data_dir, sal_mass_file, date_cal, cal_file, output)

## Run on a project directory that contains only 
project_dir <- 'Path'

