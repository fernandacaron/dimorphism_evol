rm(list = ls())

## Code to compile NPP data 

library(raster)

setwd("/Volumes/Apps")

files1 <- list.files(path = "NPP2009", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)
files2 <- list.files(path = "NPP2010", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files3 <- list.files(path = "NPP2011", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files4 <- list.files(path = "NPP2012", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files5 <- list.files(path = "NPP2013", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files6 <- list.files(path = "NPP2014", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files7 <- list.files(path = "NPP2015", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files8 <- list.files(path = "NPP2016", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files9 <- list.files(path = "NPP2017", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files10 <- list.files(path = "NPP2018", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

all_files <- c(files1, files2, files3, files4, files5, files6, files7, 
               files8, files9, files10)
for (i in 1:length(all_files)) {
  ras <- raster(all_files[i])
  if (i == 1) {
    raster_stack_all <- ras
  } else {
    raster_stack_all <- addLayer(raster_stack_all, ras)
  }
}

npp <- stackApply(raster_stack_all, indices = rep(1, length(all_files)), 
                  fun = mean, na.rm = T)

r_low <- raster(ncols = 2160, nrows = 900, ymn = -60, vals = 1:(2160*900))

npp_low <- resample(npp, r_low)

save(npp_low, file = "~/Documents/lab/dimorph_evol/data/spatial/npp_2009-2018_2160_900.RData")