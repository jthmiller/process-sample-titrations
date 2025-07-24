# wrapper script 
library(seacarb) #used to calculate TA
library(tidyverse)
library(fs)

source("R/functions.R")




run_titrator(data_dir, sal_mass_file, date_cal, cal_file)