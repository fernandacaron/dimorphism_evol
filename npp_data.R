rm(list = ls())

library(raster)

setwd("Documents/lab/dimorph_evol")

files <- list.files(path = "data/npp_Wang&al2021/NPP1981", pattern = '.tif$', 
                    all.files = TRUE, full.names = TRUE)

for (i in 1:length(files)) {
  ras <- raster(files[i])
  if (i == 1) {
    raster_stack <- ras
  } else {
    raster_stack <- addLayer(raster_stack, ras)
  }
}

files2 <- list.files(path = "data/npp_Wang&al2021/NPP1982", pattern = '.tif$', 
                    all.files = TRUE, full.names = TRUE)

for (i in 1:length(files2)) {
  ras <- raster(files2[i])
  raster_stack <- addLayer(raster_stack, ras)
}

files3 <- list.files(path = "data/npp_Wang&al2021/WMJ/for_Qi/V2_NPPprog_CI_EF0.5_1/wmjtestdata/1983", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files3)) {
  ras <- raster(files3[i])
  raster_stack <- addLayer(raster_stack, ras)
}

files4 <- list.files(path = "data/npp_Wang&al2021/NPP1984", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files4)) {
  ras <- raster(files4[i])
  raster_stack <- addLayer(raster_stack, ras)
}

files5 <- list.files(path = "data/npp_Wang&al2021/NPP1985", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files5)) {
  ras <- raster(files5[i])
  raster_stack <- addLayer(raster_stack, ras)
}

files6 <- list.files(path = "data/npp_Wang&al2021/NPP1986", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files6)) {
  ras <- raster(files6[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files7 <- list.files(path = "data/npp_Wang&al2021/NPP1987", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files7)) {
  ras <- raster(files7[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files8 <- list.files(path = "data/npp_Wang&al2021/NPP1988", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files8)) {
  ras <- raster(files8[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files9 <- list.files(path = "data/npp_Wang&al2021/NPP1989", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files9)) {
  ras <- raster(files9[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files10 <- list.files(path = "data/npp_Wang&al2021/NPP1990", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)

for (i in 1:length(files10)) {
  ras <- raster(files10[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files11 <- list.files(path = "data/npp_Wang&al2021/NPP1991", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files11)) {
  ras <- raster(files11[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files12 <- list.files(path = "data/npp_Wang&al2021/NPP1992", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files12)) {
  ras <- raster(files12[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files13 <- list.files(path = "data/npp_Wang&al2021/NPP1993", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files13)) {
  ras <- raster(files13[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files14 <- list.files(path = "data/npp_Wang&al2021/NPP1994", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files14)) {
  ras <- raster(files14[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files15 <- list.files(path = "data/npp_Wang&al2021/NPP1995", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files15)) {
  ras <- raster(files15[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files16 <- list.files(path = "data/npp_Wang&al2021/NPP1996", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files16)) {
  ras <- raster(files16[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files17 <- list.files(path = "data/npp_Wang&al2021/NPP1997", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files17)) {
  ras <- raster(files17[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files18 <- list.files(path = "data/npp_Wang&al2021/NPP1998", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files18)) {
  ras <- raster(files18[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files19 <- list.files(path = "data/npp_Wang&al2021/NPP1999", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files19)) {
  ras <- raster(files19[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files20 <- list.files(path = "data/npp_Wang&al2021/NPP2000", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files20)) {
  ras <- raster(files20[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files21 <- list.files(path = "data/npp_Wang&al2021/NPP2001", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files21)) {
  ras <- raster(files21[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files22 <- list.files(path = "data/npp_Wang&al2021/NPP2002", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files22)) {
  ras <- raster(files22[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files23 <- list.files(path = "data/npp_Wang&al2021/NPP2003", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files23)) {
  ras <- raster(files23[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files24 <- list.files(path = "data/npp_Wang&al2021/NPP2004", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files24)) {
  ras <- raster(files24[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files25 <- list.files(path = "data/npp_Wang&al2021/NPP2005", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files25)) {
  ras <- raster(files25[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files26 <- list.files(path = "data/npp_Wang&al2021/NPP2006", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files26)) {
  ras <- raster(files26[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files27 <- list.files(path = "data/npp_Wang&al2021/NPP2007", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files27)) {
  ras <- raster(files27[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files28 <- list.files(path = "data/npp_Wang&al2021/NPP2008", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files28)) {
  ras <- raster(files28[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files29 <- list.files(path = "data/npp_Wang&al2021/NPP2009", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files29)) {
  ras <- raster(files29[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files30 <- list.files(path = "data/npp_Wang&al2021/NPP2010", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files30)) {
  ras <- raster(files30[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files31 <- list.files(path = "data/npp_Wang&al2021/NPP2011", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files31)) {
  ras <- raster(files31[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files32 <- list.files(path = "data/npp_Wang&al2021/NPP2012", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files32)) {
  ras <- raster(files32[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files33 <- list.files(path = "data/npp_Wang&al2021/NPP2013", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files33)) {
  ras <- raster(files33[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files34 <- list.files(path = "data/npp_Wang&al2021/NPP2014", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files34)) {
  ras <- raster(files34[i])
  raster_stack <- addLayer(raster_stack, ras)
}



files35 <- list.files(path = "data/npp_Wang&al2021/NPP2015", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files35)) {
  ras <- raster(files35[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files36 <- list.files(path = "data/npp_Wang&al2021/NPP2016", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files36)) {
  ras <- raster(files36[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files37 <- list.files(path = "data/npp_Wang&al2021/NPP2017", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

for (i in 1:length(files37)) {
  ras <- raster(files37[i])
  raster_stack <- addLayer(raster_stack, ras)
}


files38 <- list.files(path = "data/npp_Wang&al2021/NPP2018", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)

files_2 <- c(files37, files38)
for (i in 1:length(files_2)) {
  ras <- raster(files_2[i])
  if (i == 1) {
    raster_stack2 <- ras
  } else {
    raster_stack2 <- addLayer(raster_stack2, ras)
  }
}
full <- stackApply(raster_stack2, indices = rep(1, length(files_2)), fun = mean, na.rm = T)

save(raster_stack, file = "data/npp_Wang&al2021/raster_stack.RData")
full <- stackApply(raster_stack, indices = rep(1, 46*38), fun = mean, na.rm = T)


####################

setwd("~")

files1 <- list.files(path = "/Volumes/Apps/NPP2009", pattern = '.tif$', 
                      all.files = TRUE, full.names = TRUE)
files2 <- list.files(path = "/Volumes/Apps/NPP2010", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files3 <- list.files(path = "/Volumes/Apps/NPP2011", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files4 <- list.files(path = "/Volumes/Apps/NPP2012", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files5 <- list.files(path = "/Volumes/Apps/NPP2013", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files6 <- list.files(path = "/Volumes/Apps/NPP2014", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files7 <- list.files(path = "/Volumes/Apps/NPP2015", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files8 <- list.files(path = "/Volumes/Apps/NPP2016", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files9 <- list.files(path = "/Volumes/Apps/NPP2017", pattern = '.tif$', 
                     all.files = TRUE, full.names = TRUE)
files10 <- list.files(path = "/Volumes/Apps/NPP2018", pattern = '.tif$', 
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

save(raster_stack, file = "data/npp_Wang&al2021/raster_stack.RData")
