# Alberto Rovellini
# 12/28/2023
# This code is to bias correct Atlantis physical forcing files
# To start we will focus on temperature and salinity
# Ideally this should use CDO / NCO, based on the concepts we have used for GOACLIM

# we need:
# Historical run
# Hindcast
# Projection

# The idea is to calculate a delta (by month) between historical run and hindcast, in each cell
# then it adds each delta to each time step of the projections

# As of 12/28/2023 we have not yet translated the historical runs, so let's use the same years from the two hindcast (revised and unrevised)
# The resulting delta will be very small but not 0

# to get a minimal working example we need 2 files each.
# Use 2080-2081 for projection
# Use 2004-2005 from old hindcast as historical run
# Use 2004-2005 from new hindcast as hindcast

# NB:
# Writing this has highlighted some outstanding issues with the forcing files that should have been solved by now:
# 1. We still have leap years in the forcign files out of the translation process
# 2. Annual nc forcing file for  2006 will include midnight of 01/01/2007

# NB:
# Some of the commands below get a little hard to read, but paste() functions seem to struggle due to the formatting (e.g., nested quotes, escaped characters, etc) so let's just bear with it

# it seems like the ymonmean (or sub) call works on one file at a time, so we need to merge the annual files first
# I am unsure how this will scale with ~25 annual files for the overlapping period but should be fine 
# merge the hindcast
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\hindcast\\*.nc data\\hindcast\\hindcast_merged.nc")
# at the moment, the hindcast files include the first time step of the following year
# this causes ymonmean to struggle later on, so a temporary solution is to drop the last time step
# get the number of time steps
library(ncdf4)
ncf <- nc_open("data/hindcast/hindcast_merged.nc")
ntime <- length(ncf$dim$t$vals)
nc_close(ncf)

# drop the lat time step
shell(paste0("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" seltimestep,1/", 
             ntime-1, 
             " data\\hindcast\\hindcast_merged.nc data\\hindcast\\hindcast_merged_cut.nc"))

# merge the historical
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\historical\\*.nc data\\historical\\historical_merged.nc")

# now get monthly means per unit space (is this the right unit? Should it be daily instead?)
# hindcast
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" ymonmean data\\hindcast\\hindcast_merged_cut.nc data\\hindcast\\hindcast_ymonmean.nc")

# historical
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" ymonmean data\\historical\\historical_merged.nc data\\historical\\historical_ymonmean.nc")

# now subtract historical means from hindcast means
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" sub data\\hindcast\\hindcast_ymonmean.nc data\\historical\\historical_ymonmean.nc data\\delta\\monthly_delta.nc")

# now add deltas to projections
# The projection files are in 12-hour time steps, but we have monthly deltas
# So, now we move to an R approach for some finer handling
# Step 1: go from monthly deltas to daily deltas. Have one for regular years and one for leap years
# Step 2: loop over each projection year. Add the deltas to the projection years

library(ncdf4)
library(tidync)
library(tidyverse)
library(lubridate)

# what variable are we handling?
this_variable <- "temperature"

monthly_delta_nc <- nc_open("data/delta/monthly_delta.nc")
monthly_delta_values <- ncvar_get(monthly_delta_nc, this_variable)

# Initialize the output array
daily_delta_values <- array(dim = c(7, 109, 730))
daily_delta_values_leap <- array(dim = c(7, 109, 732))

# days for each month
days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
days_leap <- c(31,29,31,30,31,30,31,31,30,31,30,31)

# Function to fill the output array
fill_array <- function(monthly_deltas, daily_deltas, leap = F) {
  
  if(leap){dd = days_leap} else {dd = days} # determine if this array will be for a leap year
  
  start_index <- 1
  
  for (i in 1:12) {
    print(i)
    # Calculate the number of 12-hour intervals in the month
    num_intervals <- dd[i] * 2
    end_index <- start_index + num_intervals - 1
    
    # Fill the array for this month
    daily_deltas[,,start_index:end_index] <- array(monthly_deltas[,,i], dim = c(7, 109, num_intervals))
    
    # Update the start index for the next month
    start_index <- end_index + 1
  }
  return(daily_deltas)
}

# Fill the output array with your data
daily_delta_values <- fill_array(monthly_deltas = monthly_delta_values, daily_deltas = daily_delta_values, leap = F)
daily_delta_values_leap <- fill_array(monthly_deltas = monthly_delta_values, daily_deltas = daily_delta_values_leap, leap = T)

# now loop over the projection files, and add the deltas with R's array operations
# This approach will be time-step agnostic
# be careful with NA's
# this works on annual forcing files
proj_files <- list.files("data/projection/", full.names = T)

for(i in 1:length(proj_files)){
  
  # point at pre-correction file
  original_proj_file <- proj_files[i]
  # just file name to re-write later
  fn <- gsub("data/projection/","",original_proj_file)
  
  # open pre-correction file
  original_proj_nc <- nc_open(original_proj_file)
  original_proj_values <- ncvar_get(original_proj_nc, this_variable)
  
  # define time steps
  
  
  # calculate values for the delta-corrected file
  if(dim(original_proj_values)[3]==730){
    corrected_proj_values <- original_proj_values + daily_delta_values
  } else if (dim(original_proj_values)[3]==732) {
    corrected_proj_values <- original_proj_values + daily_delta_values_leap
  } else {
    stop("The projection nc array is not 730 or 732") # this would indicate a bug upstream in the translation
  }
  
  # now re-pack this to nc forcing files in a different folder (bias-corrected)
  # create a new nc file with the dimensions, variables, and attributes of the original one
  # the easiest way to do this is to copy the original file to a new location, open it, modify it, and close it
  corrected_proj_file <- paste0("data/projection_corrected/",fn)
  
  # copy the original necdf file that we are correcting
  file.copy(original_proj_file, corrected_proj_file)
  
  # open the new file
  corrected_proj_nc <- nc_open(corrected_proj_file, write = T)
  
  # write the delta-corrected variable to it
  ncvar_put(corrected_proj_nc, varid = this_variable, vals = corrected_proj_values)
  
  # close both nc files
  nc_close(original_proj_nc)
  nc_close(corrected_proj_nc)
  
}
