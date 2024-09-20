

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue, ggspatial, climateR, geodata, rnaturalearthdata, rnaturalearth, readxl)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50, path = './tmpr')

## Trial points
pnts <- read_excel('//catalogue/workspace-cluster9/2024/WCR/V1/r/tbl/points/IMLVT Site info.xlsx')
pnts <- dplyr::select(pnts, id, Longitude, Latitude)
pnts <- mutate(pnts, id = gsub('Bulegeni\t\t_UGA', 'Bulegeni_UGA', id))

## World points
# pnts.wrld <- read_csv('./tble/Arabica_cleaned_up_points.csv')[,-1]
pnts.wrld <- read_csv('../V6_worldpoints_bios/tbl/CoffeePoints_2.csv')
pnts.wrld <- mutate(pnts.wrld, source = 'Bunn', id = 1:nrow(pnts.wrld))
pnts <- mutate(pnts, source = 'WCR')

pnts <- rbind(pnts %>% dplyr::select(source, id, Longitude, Latitude), pnts.wrld %>% dplyr::select(source, id, Longitude, Latitude))

# Current worldclim -------------------------------------------------------
prec.wrld <- worldclim_global(res = 0.5, var = 'prec', path = './tmpr') 
tmin.wrld <- worldclim_global(res = 0.5, var = 'tmin', path = './tmpr') 
tmax.wrld <- worldclim_global(res = 0.5, var = 'tmax', path = './tmpr') 

## To change the names
names(prec.wrld) <- glue('prec_{1:12}')
names(tmin.wrld) <- glue('tmin_{1:12}')
names(tmax.wrld) <- glue('tmax_{1:12}')

# Future worldclim --------------------------------------------------------
root <- '//ALLIANCEDFS.ALLIANCE.CGIAR.ORG/data_cluster19/GLOBAL/climate/Agroclimas/data/ipcc_6ar_wcl_downscaled'
dirs <- dir_ls(root, type = 'directory') %>% as.character()

# Extract values from the points ------------------------------------------

## Current -----------------------------------------------------------------
pnts.crnt.prec <- as_tibble(cbind(pnts, as_tibble(terra::extract(prec.wrld, pnts[,c('Longitude', 'Latitude')]))[,-1]))
pnts.crnt.tmin <- as_tibble(cbind(pnts, as_tibble(terra::extract(tmin.wrld, pnts[,c('Longitude', 'Latitude')]))[,-1]))
pnts.crnt.tmax <- as_tibble(cbind(pnts, as_tibble(terra::extract(tmax.wrld, pnts[,c('Longitude', 'Latitude')]))[,-1]))

pnts.crnt <- list(pnts.crnt.prec, pnts.crnt.tmin, pnts.crnt.tmax) %>% reduce(., inner_join)
write.csv(pnts.crnt, './tble/trials-bunn_current.csv', row.names = FALSE)

## Future ------------------------------------------------------------------

## Function to use
extr.vles.ftre <- function(sspe, year){
  
  # sspe <- dirs[1]
  # year <- '2050s'
  
  ## To list the files 
  cat('To process: ', sspe, '\n')
  dire <- dir_ls(sspe) %>% grep(year, ., value = T) %>% dir_ls()
  
  ## To extract the values 
  vles <- map(.x = 1:length(dire), .f = function(i){
    
    try(
      expr = {
        
        ## To list the files
        cat('To start the process!\n')
        dir <- dire[i]
        fls <- dir_ls(dir) %>% dir_ls(.) %>% as.character()
        
        ## To read as a raster files
        ppt <- grep('prec', fls, value = T) %>% rast()
        tmn <- grep('tmin', fls, value = T) %>% grep('.tif$', ., value = T) %>% rast()
        tmx <- grep('tmax', fls, value = T) %>% grep('.tif$', ., value = T) %>% rast()
        
        ## To extract the values
        ppt.vls <- as_tibble(cbind(pnts, terra::extract(ppt, pnts[,c('Longitude', 'Latitude')])[,-1]))
        tmn.vls <- as_tibble(cbind(pnts, terra::extract(tmn, pnts[,c('Longitude', 'Latitude')])[,-1]))
        tmx.vls <- as_tibble(cbind(pnts, terra::extract(tmx, pnts[,c('Longitude', 'Latitude')])[,-1]))
        
        ## To join the three tables into only one 
        vls <- list(ppt.vls, tmn.vls, tmx.vls) %>% reduce(., inner_join)
        vls <- mutate(vls, ssp = basename(sspe), gcm = basename(dir), .before = 'source')
        
        ## Finish
        cat('Done!\n')
        return(vls)
        
      }
    )
    
  })

  ## Join all the tables into only one
  vles <- vles[map_lgl(vles, is.data.frame)]
  vles <- bind_rows(vles)
  vles <- mutate(vles, period = year, .before = 'ssp')
  
  ## Finish
  return(vles)
  
}

## To apply the function
pnts.370.50 <- extr.vles.ftre(sspe = dirs[3], year = '2050s')# pnts.370 <- map(.x = c('2030s', '2050s', '2070s', '2090s'), .f = function(i){extr.vles.ftre(sspe = dirs[3], year = i)})
pnts.370.50 <- pnts.370.50[,1:43]
colnames(pnts.370.50) <- gsub('wc2.1_30s_', '', colnames(pnts.370.50))
write.csv(pnts.370.50, './tble/trials-bunn_ssp370-2050s-gcms.csv', row.names = FALSE)

## To make the summarise 
smmr.370.50 <- pnts.370.50 %>% 
  group_by(period, ssp, source, id, Longitude, Latitude) %>% 
  summarise_if(is.numeric, mean) %>% 
  ungroup()

# Join both tables into only one ------------------------------------------
pnts.crnt
pnts.crnt <- pnts.crnt %>% mutate(period = 'current', ssp = 'current')
pnts.crnt <- pnts.crnt %>% dplyr::select(period, ssp, everything())
colnames(smmr.370.50) <- gsub('_0', '_', colnames(smmr.370.50))

colnames(pnts.crnt)
colnames(smmr.370.50)

unique(pnts.370.50$gcm)
# "access_cm2"       "access_esm1_5"    "awi_cm_1_1_mr"    "bcc_csm2_mr"      "canesm5"          "cmcc_esm2"       
# "cnrm_esm2_1"      "ec_earth3_veg_lr" "giss_e2_1_g"      "giss_e2_1_h"      "inm_cm4_8"        "inm_cm5_0"       
# "ipsl_cm6a_lr"     "miroc6"           "miroc_es2l" 

vles <- rbind(pnts.crnt, smmr.370.50)
write.csv(vles, './tble/trials-bunn_current-ssp370-2050s-avrg.csv', row.names = FALSE)

