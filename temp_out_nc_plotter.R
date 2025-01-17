# # Alberto Rovellini
library(tidync)
library(tidyverse)
library(ncdf4)
library(rbgm)
library(sf)
# 
# # GOAL: plot temperature time series from out.nc
# # Get model-wide means weighted by volume
# # Do surface and bottom and all
# # Do spatial plots
# # organize as a function to then compare different runs to one another
# 
# this_no <- '1771'
# run_end <- 30
# layer = "bottom"
ssp_key <- data.frame("run" = c('1770','1771','1772','1773'),
                      ssp = c("base_1999","ssp126","ssp245","ssp585"))


# Function ----------------------------------------------------------------
#' Calculate Surface or Bottom Temperature from Atlantis Output
#' 
#' @description 
#' Extracts and processes temperature data from Atlantis model output files, calculating
#' either surface or bottom temperatures for specified model runs. The function handles
#' BGM (box geometry) files and NetCDF output files, accounting for boundary boxes and
#' volume-weighted averaging.
#' 
#' @param this_no numeric. The run number to process
#' @param run_end numeric. The end time of the run in years (default: 5)
#' @param layer character. Which vertical layer to analyze: "surface" or "bottom" (default: "surface")
#' @param goa_wide logical. If TRUE, averages across all boxes; if FALSE, keeps box-specific values (default: TRUE)
#' 
#' @return A tibble containing the following columns:
#'   \itemize{
#'     \item run: The run number
#'     \item time: Time in years from start of simulation
#'     \item b: Box ID (only present if goa_wide = FALSE)
#'     \item mean_temp: Volume-weighted mean temperature
#'     \item yr: Year of simulation (floor of time)
#'     \item doy: Day of year (0-365, in 73-day intervals)
#'     \item toy: Time of year as formatted date (e.g., "Jan-01")
#'   }
#'
#' @details
#' The function performs the following steps:
#'   1. Reads and processes the BGM file for box geometry
#'   2. Loads NetCDF output file containing temperature data
#'   3. Excludes boundary boxes from calculations
#'   4. Processes temperature data using volume-weighted means
#'   5. Aggregates results by time and optionally by box
#'   6. Adds temporal information (year, day of year, formatted date)
#'
#' @note
#' - Assumes NetCDF output is provided every 73 days
#' - First timestep (t=0) is dropped as it contains initialization values
#' - Requires box geometry (BGM) file in specified directory
#' - Surface layer is assumed to be the highest vertical layer (excluding sediment)
#' - Bottom layer is the lowest vertical layer with positive volume
#'
#' @importFrom dplyr filter mutate select group_by summarise ungroup left_join
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dfr
#' @importFrom sf st_crs
#' @importFrom ncdf4 nc_open ncvar_get
#' @importFrom tidync tidync hyper_filter hyper_tibble
#'
#' @examples
#' # Calculate surface temperatures for run 1, averaged across all boxes
#' surf_temps <- calculate_surface_temp(1)
#'
#' # Calculate bottom temperatures for run 2, keeping box-specific values
#' bottom_temps <- calculate_surface_temp(2, layer = "bottom", goa_wide = FALSE)
#' 
calculate_surface_temp <- function(this_no, 
                                   run_end = 5, 
                                   layer = "surface",
                                   goa_wide = T) {
  # Set directory path
  this.dir <- paste0('C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_', this_no, '_SSP')
  
  # Load and process BGM file
  fl <- paste(this.dir,'GOA_WGS84_V4_final.bgm',sep='/')
  bgm <- read_bgm(fl)
  goa_sf <- box_sf(bgm)
  goa_sf <- goa_sf %>% mutate(box_id = box_id+1)
  st_crs(goa_sf) <- st_crs(attr(goa_sf$geometry, "crs")$proj)
  
  # Load and process NC file
  out_fl <- paste(this.dir, paste0('outputGOA', '0', this_no, '_test.nc'),sep='/')
  out <- tidync(out_fl)
  this.nc <- ncdf4::nc_open(out_fl)
  
  # Get volumes and time dimensions
  volumes <- out %>% hyper_filter(t=t==0) %>% hyper_tibble(select_var="volume") %>% dplyr::select(-t)
  ts <- ncdf4::ncvar_get(this.nc,varid = "t") %>% as.numeric
  tyrs <- ts/(60*60*24*365)
  
  # Get boundary boxes
  boundary_boxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id)
  
  # Process temperature data
  temp <- ncvar_get(this.nc, "Temp")
  temp[, boundary_boxes, ] <- NA
  run_length <- dim(temp)[3]
  
  # Create temperature dataframe
  temp_df <- map_dfr(1:run_length, function(i) {
    matrix_i <- temp[, , i]
    as.data.frame(matrix_i) %>%
      mutate(z = 1:7,
             time = tyrs[i])
  })
  
  # Clean and reshape data
  temp_df <- temp_df %>%
    pivot_longer(-c(time,z), values_to = "temp", names_to = "b") %>%
    mutate(b = as.numeric(gsub("V","",b))) %>%
    dplyr::select(time,b,z,temp) %>%
    arrange(time,b,z) %>%
    filter(time <= run_end) %>%
    left_join(volumes, by = c("b","z")) %>%
    filter(z < 7) %>% # 7 is sediment in the output
    filter(volume > 0) %>%
    mutate(run = this_no)
  
  # get temp for this slice
  temp_slice <- temp_df %>% 
    group_by(run, time, b) %>%
    filter(z == ifelse(layer == "surface", max(z), min(z))) %>% #surface is always 6 in the output, bottom is lowest number
    group_by(!!!syms(c("run", "time", if(!isTRUE(goa_wide)) "b"))) %>%
    summarise(mean_temp = weighted.mean(temp, volume, na.rm = T)) %>%
    ungroup()
  
  # add which day of the year we are at
  temp_slice <- temp_slice %>%
    mutate(yr = floor(time)) %>%
    group_by(!!!syms(c("run", "yr", if(!isTRUE(goa_wide)) "b"))) %>%
    mutate(doy = (row_number()-1)*73) %>% # nc output is every 73 days
    ungroup() %>%
    mutate(toy = as.Date(doy),
           toy = format(toy, "%b-%d"))
  
  # order toy chronologically
  temp_slice$toy <- factor(temp_slice$toy, levels = c("Jan-01","Mar-15","May-27","Aug-08","Oct-20"))
  
  # drop t0 as that is the temp fillvalue from init.nc
  temp_slice <- temp_slice %>% filter(time > 0)
  
  return(temp_slice)
}


# Model-wide averages -----------------------------------------------------
this_temp_df <- bind_rows(
  lapply(
    list('1770','1771','1772','1773'), 
    calculate_surface_temp, 
    run_end = 80, 
    layer = "surface",
    goa_wide = T
  )
)

# add ssp info
this_temp_df <- this_temp_df %>%
  left_join(ssp_key, by = "run")

# view
# plot
p_goa_wide <- this_temp_df %>%
  ggplot(aes(x = time+2015, y = mean_temp, color = ssp))+
  geom_line(linewidth = 1.2, alpha = 0.75)+
  theme_bw()+
  scale_color_manual(values = c("#4575B4", "#FFD700", "#FF8C00", "#D73027"))+
  labs(x = "Year", y = "Temperature (C)", color = "SSP")+
  facet_wrap(~toy)
p_goa_wide

# By box ------------------------------------------------------------------
this_temp_df <- bind_rows(
  lapply(
    list('1770','1771','1772','1773'), 
    calculate_surface_temp, 
    run_end = 80, 
    layer = "surface",
    goa_wide = F
  )
)

# add ssp info
this_temp_df <- this_temp_df %>%
  left_join(ssp_key, by = "run")

# view
# plot
p2 <- this_temp_df %>%
  filter(!is.nan(mean_temp)) %>%
  ggplot(aes(x = time+2015, y = mean_temp, color = b, group = b))+
  geom_line()+
  scale_color_viridis_c(option = "inferno") +
  theme_bw()+
  labs(x = "Year", y = "Temperature (C)", color = "SSP")+
  facet_grid(ssp~toy)
p2
