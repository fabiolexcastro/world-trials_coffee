
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue, hablar, dismo, gtools, ggspatial, climateR, geodata, rnaturalearthdata, rnaturalearth, readxl)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
wrld <- ne_countries(returnclass = 'sf', scale = 50)
vles <- read_csv('./tble/trials-bunn_current-ssp370-2050s-avrg.csv')

## Solar radiation
srad <- '//catalogue/workspace-cluster9/DATA/ET_SolRad'
srad <- dir_ls(srad, regexp = 'et_solrad_') %>% as.character() %>% mixedsort()
srad <- rast(srad)

# To extract srad variables from the points -------------------------------
srad.vles <- as_tibble(terra::extract(srad, vles[,c('Longitude', 'Latitude')])[,-1])
colnames(srad.vles) <- glue('srad_{1:12}')
allv <- cbind(vles, srad.vles)
allv <- as_tibble(allv)

# To estimate bioclimatic variables ---------------------------------------
bios <- dismo::biovars(prec = as.matrix(allv[,7:18]), tmin = as.matrix(allv[,19:30]), tmax = as.matrix(allv[,31:42]))
bios <- as.data.frame(bios)
bios <- as_tibble(bios)
colnames(bios) <- glue('bioc_{1:19}')
bios <- as_tibble(cbind(allv[,1:6], bios))

# To estimate PET variables -----------------------------------------------
dfrm <- allv %>% gather(var, value, -c(period, ssp, source, id, Longitude, Latitude)) %>% separate(., col = 'var', into = c('variable', 'month'), sep = '_')
dfrm <- dfrm %>% spread(variable, value)
dfrm <- dfrm %>% mutate(tavg = (tmax + tmin) / 2)

dfrm <- dfrm %>% mutate(etp = 0.0013 * 0.408 * srad * (tavg + 17) * (tmax - tmin - 0.0123 * prec) ^ 0.76)
dfrm <- inner_join(dfrm, tibble(month = as.character(1:12), mult = c(31,29,31,30,31,30,31,31,30,31,30,31)), by = 'month')
dfrm <- mutate(dfrm, pet = etp * mult)
dfrm <- dplyr::select(dfrm, -etp, -mult)
dfrm <- mutate(dfrm, pet = ifelse(is.na(pet), 0, pet))

# To estimate def ---------------------------------------------------------
dfrm <- mutate(dfrm, def = prec - pet)

## Driest def
dfrm.def <- dplyr::select(dfrm, period, ssp, source, id, Longitude, Latitude, month, def)
dfrm.ppt <- dplyr::select(dfrm, period, ssp, source, id, Longitude, Latitude, month, prec) 

## Seasons
ssns <- tibble(seasons = LETTERS[1:12], month_1 = c(1:12), month_2 = c(2:12, 1), month_3 = c(3:12, 1:2))
gids <- unique(dfrm.ppt$id)

# Function to calculate driest quarter by each coffee presence
extr.drst <- function(gid){
  
  ## Filtering
  cat('To process: ', gid, '\n')
  ppt <- dfrm.ppt %>% filter(id == gid & period == 'current')
  def <- dfrm.def %>% filter(id == gid)
  
  ## Accumulate by each month 
  drs <- map_dfr(.x = 1:12, .f = function(s){
    ssn <- ssns[s,] %>% dplyr::select(2:4) %>% as.numeric()
    rsl <- ppt %>% filter(month %in% as.character(ssn))
    rsl <- tibble(id = gid, month_1 = ssn[1], month_2 = ssn[2], month_3 = ssn[3], prec_ssn = sum(pull(rsl, prec)))
    return(rsl)
  }) %>% top_n(x = ., n = -1, wt = prec_ssn) %>% .[1,]
  
  dfr.mnt <- drs %>% dplyr::select(2:4) %>% as.numeric()
  
  rsl <- def %>% 
    filter(month %in% as.character(dfr.mnt)) %>% 
    group_by(period, ssp, source, id, Longitude, Latitude) %>%
    dplyr::summarise(def_driest = sum(def, na.rm = T)) %>% ungroup()

  return(rsl)
  
}

## To apply the function 
drst <- map(gids, extr.drst)
drst <- bind_rows(drst)

# Summarise main dfrm -----------------------------------------------------
smmr <- dfrm %>% 
  group_by(period, ssp, source, id, Longitude, Latitude) %>% 
  dplyr::summarise(
    prec = sum(prec), 
    tmax = mean(tmax),
    tmin = mean(tmin),
    pet = sum(pet)
  ) %>% 
  ungroup()
smmr <- inner_join(smmr, drst, by = c('period', 'ssp', 'source', 'id', 'Longitude', 'Latitude'))
smmr <- inner_join(smmr, bios, by = c('period', 'ssp', 'source', 'id', 'Longitude', 'Latitude'))
colnames(smmr)[7] <- 'ppt'

write.csv(smmr, './tble/send/bunn-trials_wc-ssp370.csv', row.names = FALSE)

