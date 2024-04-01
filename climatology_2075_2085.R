# Alberto Rovellini
# 01/23/2024
# This code takes temp and salt nc files for 2075 and 2085 and creates a climatology
# then it writes the files out to one forcing file each 
# This should have n=been done with CDO, however the 12-hr time step and the presence of leap years complicates it

library(tidyverese)
library(ncdf4)
library(tidync)

# apply the function we have created for the delta correction as the averaging is the same process
# merge corrected projection files
shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\ssp245\\temp\\OY_2075_2084\\*.nc data\\ssp245\\temp\\OY_2075_2084\\proj_merged.nc")

# what variable are we handling?
this_variable <- "temperature"

# files
projection_file <- "data/ssp245/temp/OY_2075_2084/proj_merged.nc"

# apply function
mean_proj <- make_delta_array(variable = this_variable,
                              sim_period = "projection",
                              mean = T,
                              leap = F)

# pack to netcdf
# for the purpose of not having to rebuild the netcdf from scratch, put the new variable to an exisiting file
mean_proj_file <- paste0("data/ssp245/temp/OY_2075_2084/mean_projection_", this_variable, ".nc")

# copy the original necdf file that we are correcting
file.copy("data/ssp245/temp/OY_2075_2084/goa_roms_temp_2075.nc", mean_proj_file)

# open the new file
mean_proj_nc <- nc_open(mean_proj_file, write = T)

# write the delta-corrected variable to it
ncvar_put(mean_proj_nc, varid = this_variable, vals = mean_proj)

# close the file
nc_close(mean_proj_nc)
