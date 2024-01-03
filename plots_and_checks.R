# Alberto Rovellini
# 1/2/2024
# This code plots NetCDF forcings and contains sanity checks for the delta- correction method

# List of packages for session
.packages = c("dplyr", "ggplot2", "purrr", "tidyr", "ncdf4", "tidync", "sf", "ggh4x", "rbgm", "cowplot", "maps", "mapdata")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

select <- dplyr::select

# Atlantis model spatial domain
bgm.file <- "data/GOA_WGS84_V4_final.bgm" 
atlantis_bgm <- bgm.file %>% read_bgm()
atlantis_box <- atlantis_bgm %>% 
  box_sf() %>%
  st_transform(crs = 4326) %>%
  mutate(botz = -1*botz)

# define depths
cum.depth <- c(1,30,100,200,500,1000,3969)

# for hindcast and historical, make sure you read in the same year(s). Read in the merged files
# hindcast
hindcast_file <- "data/hindcast/hindcast_merged.nc"
hindcast_nc <- nc_open(hindcast_file)

# HINDCAST and HISTORICAL -------------------------------------------------
# historical
historical_file <- "data/historical/historical_merged.nc"
historical_nc <- nc_open(historical_file)

# pull variables with tidync
temp_hind <- tidync(hindcast_file) 
temp_hist <- tidync(historical_file)

these_vars <- hyper_grids(temp_hind) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    temp_hind %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

grids <- these_vars %>% filter(name=="temperature") %>% pluck('grd')

dat_temp_hind <- temp_hind %>% activate(grids) %>% hyper_tibble()
dat_temp_hist <- temp_hist %>% activate(grids) %>% hyper_tibble()

# numbering of b and z starts from 1 - change
dat_temp_hind <- dat_temp_hind %>% 
  mutate(b = y - 1, z = x - 1, run = "hind") %>%
  select(-x,-y)

dat_temp_hist <- dat_temp_hist %>% 
  mutate(b = y - 1, z = x - 1, run = "hist") %>%
  select(-x,-y)

dat <- rbind(dat_temp_hind, dat_temp_hist)
# view
# summarize in space
# should weight these averages by cell volume for a better representation, but this is only visual
dat %>%
  group_by(run,t) %>%
  summarize(meantemp = mean(temperature)) %>%
  ggplot()+
  geom_line(aes(x = t, y = meantemp, color = run), linewidth = 1)+
  theme_bw()

# let's plot some random boxes: 
b_toplot <- sample(unique(dat$b), size = 6, replace = T)
# remember: 6 is the sediment, 0 is the layer above the sediment, and then up to the surface

dat %>%
  filter(b %in% b_toplot) %>%
  ggplot()+
  geom_line(aes(x = t, y = temperature, color = factor(z), linetype = run), linewidth = 1)+
  theme_bw()+
  facet_wrap(~b)

# YMONMEAN ----------------------------------------------------------------
# plot monthly means and compare to plot above
hind_ymonmean_file <- "data/hindcast/hindcast_ymonmean.nc"
hind_ymonmean_nc <- nc_open(hind_ymonmean_file)

# historical
hist_ymonmean_file <- "data/historical/historical_ymonmean.nc"
hist_ymonmean_nc <- nc_open(hist_ymonmean_file)

# pull variables with tidync
temp_hind_ym <- tidync(hind_ymonmean_file) 
temp_hist_ym <- tidync(hist_ymonmean_file)

these_vars <- hyper_grids(temp_hind_ym) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    temp_hind_ym %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

grids <- these_vars %>% filter(name=="temperature") %>% pluck('grd')

dat_temp_hind_ym <- temp_hind_ym %>% activate(grids) %>% hyper_tibble()
dat_temp_hist_ym <- temp_hist_ym %>% activate(grids) %>% hyper_tibble()

# numbering of b and z starts from 1 - change
dat_temp_hind_ym <- dat_temp_hind_ym %>% 
  mutate(b = y - 1, z = x - 1, run = "hind") %>%
  select(-x,-y)

dat_temp_hist_ym <- dat_temp_hist_ym %>% 
  mutate(b = y - 1, z = x - 1, run = "hist") %>%
  select(-x,-y)

dat_ym <- rbind(dat_temp_hind_ym, dat_temp_hist_ym)
# view
# summarize in space
dat_ym %>%
  group_by(run,t) %>%
  summarize(meantemp = mean(temperature)) %>%
  ggplot()+
  geom_line(aes(x = t, y = meantemp, color = run), linewidth = 1)+
  theme_bw()

# let's plot some random boxes: 
# remember: 6 is the sediment, 0 is the layer above the sediment, and then up to the surface

dat_ym %>%
  filter(b %in% b_toplot) %>%
  ggplot()+
  geom_line(aes(x = t, y = temperature, color = factor(z), linetype = run), linewidth = 1)+
  theme_bw()+
  facet_wrap(~b)

# MONTHLY DELTAS ----------------------------------------------------------
monthly_delta_file <- "data/delta/monthly_delta.nc"
delta_nc <- nc_open(monthly_delta_file)

# pull variables with tidync
temp_delta <- tidync(monthly_delta_file) 

these_vars <- hyper_grids(temp_delta) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    temp_delta %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

grids <- these_vars %>% filter(name=="temperature") %>% pluck('grd')

dat_temp_delta <- temp_delta %>% activate(grids) %>% hyper_tibble()

# numbering of b and z starts from 1 - change
dat_temp_delta <- dat_temp_delta %>% 
  mutate(b = y - 1, z = x - 1) %>%
  select(-x,-y)

# view
# summarize in space
dat_temp_delta %>%
  group_by(t) %>%
  summarize(meandelta = mean(temperature), linewidth = 1) %>%
  ggplot()+
  geom_bar(aes(x = t, y = meandelta), stat = "identity")+
  theme_bw()

# let's plot some random boxes: 
# remember: 6 is the sediment, 0 is the layer above the sediment, and then up to the surface

dat_temp_delta %>%
  filter(b %in% b_toplot) %>%
  ggplot()+
  geom_line(aes(x = t, y = temperature, color = factor(z)), linewidth = 1)+
  theme_bw()+
  facet_wrap(~b)

# PROJECTIONS -------------------------------------------------------------

# projection files do not get merged
# pick one year and plot it before and after the delta correction
# changes should be consistent with the delta shown above
# pre-correction
original_proj_file <- "data/projection/goa_roms_temp_2080.nc"
original_proj_nc <- nc_open(original_proj_file)

# corrected
corrected_proj_file <- "data/projection_corrected/goa_roms_temp_2080.nc"
corrected_proj_nc <- nc_open(corrected_proj_file)

# pull variables with tidync
temp_o <- tidync(original_proj_file) 
temp_c <- tidync(corrected_proj_file)

these_vars <- hyper_grids(temp_o) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    temp_o %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

grids <- these_vars %>% filter(name=="temperature") %>% pluck('grd')

dat_temp_o <- temp_o %>% activate(grids) %>% hyper_tibble()
dat_temp_c <- temp_c %>% activate(grids) %>% hyper_tibble()

# numbering of b and z starts from 1 - change
dat_temp_o <- dat_temp_o %>% 
  mutate(b = b - 1, z = z - 1, run = "original")

dat_temp_c <- dat_temp_c %>% 
  mutate(b = b - 1, z = z - 1, run = "corrected")

dat_proj <- rbind(dat_temp_o, dat_temp_c)
# view
# summarize in space
# should weight these averages by cell volume for a better representation, but this is only visual
dat_proj %>%
  group_by(run,t) %>%
  summarize(meantemp = mean(temperature), linewidth = 1) %>%
  ggplot()+
  geom_line(aes(x = t, y = meantemp, color = run))+
  theme_bw()

# remember: 6 is the sediment, 0 is the layer above the sediment, and then up to the surface
dat_proj %>%
  filter(b %in% b_toplot) %>%
  ggplot()+
  geom_line(aes(x = t, y = temperature, color = factor(z), linetype = run), linewidth = 1)+
  theme_bw()+
  facet_wrap(~b)

# Everything seems to be working from a coding / algebra standpoint
# however, it is apparent that this causes some jagged jumps
# This is because we are calculating monthly deltas that are applied equally to each 12-hourly time step
# So, if the time series in the projection is smooth, jumping from a delta to the next is very visible.
# This is not ideal for Atlantis, we can have unintended jumps in temperature that may trickle down to the biology
# Let's explore some options:
# Calculate daily deltas over the overlapping period. Does this even make sense conceptually?
# 1. Smooth deltas by assigning the value of the monthly delta to the middle of each month, and then interpolate for the days
# 2. This option makes sense to me: if in month A the delta was +1 (i.e., the hindcast was 1 C warmer than the historical),
# and in month A+1 the delta was -1, one may assume that the delta has changed smoothly over time. This may or may not be true
