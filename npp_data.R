rm(list = ls())

library(raster)

# Código rodado no cluster da UFG

setwd("dimorphism_evol")

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

all_files <- c(files1, files2, files3, files4, files5, files6, files7, files8,
               files9, files10)
for (i in 1:length(all_files)) {
  ras <- raster(all_files[i])
  if (i == 1) {
    raster_stack_all <- ras
  } else {
    raster_stack_all <- addLayer(raster_stack_all, ras)
  }
}
full <- stackApply(raster_stack_all, indices = rep(1, length(all_files)), fun = mean, na.rm = T)

save(full, file = "raster_stack_npp_Wang2021.RData")
