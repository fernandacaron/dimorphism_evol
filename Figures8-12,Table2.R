rm(list = ls())

setwd("Documents/lab/dimorph_evol")

library(phytools)
library(caper)
library(sf)
library(raster)
library(dplyr)
library(fasterize)
library(viridis)
library(rgeos)

dat <- read.csv("data/BodySizeAves_30may22_edit.csv", row.names = 1)

dat_red <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
                 complete.cases(dat$Body_mass_g_F_mean), ]

SDI <- function (male = male, female = female, cutoff = FALSE, 
                 cut.value = 0.10) {
  if (cutoff == TRUE) {
    if (male <= female) {
      ifelse((female - male) >= (male * cut.value),
             SDI <- ((female/male) - 1),
             SDI <- 0)
    } else {
      if (female < male) {
        ifelse((male - female) >= (female * cut.value),
               SDI <- -((male/female) - 1),
               SDI <- 0)
      }
    }
  }
  
  if (cutoff == FALSE) {
    if (male <= female) {
      SDI <- (female/male) - 1
    } else {
      if (female < male) {
        SDI <- -((male/female) - 1)
      } 
    }
    
  }
  
  return(SDI)
}

sdi <- numeric()
for (i in 1:nrow(dat_red)) {
  sdi[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
                female = dat_red$Body_mass_g_F_mean[i],
                cutoff = FALSE)
}
names(sdi) <- dat_red$Scientific_name
sdi <- sdi[complete.cases(sdi)]

load(file = "data/spatial/acc.RData")
load(file = "data/spatial/ans.RData")
load(file = "data/spatial/apo.RData")
load(file = "data/spatial/cha.RData")
load(file = "data/spatial/col.RData")
load(file = "data/spatial/gal.RData")
#load(file = "data/spatial/pas.RData")
load(file = "data/spatial/pic.RData")
load(file = "data/spatial/psi.RData")

sdi_acc <- sdi[names(sdi) %in% acc$sci_name]
sdi_ans <- sdi[names(sdi) %in% ans$sci_name]
sdi_apo <- sdi[names(sdi) %in% apo$sci_name]
sdi_cha <- sdi[names(sdi) %in% cha$sci_name]
sdi_col <- sdi[names(sdi) %in% col$sci_name]
sdi_gal <- sdi[names(sdi) %in% gal$sci_name]
#sdi_pas <- sdi[names(sdi) %in% pas$sci_name]
sdi_pic <- sdi[names(sdi) %in% pic$sci_name]
sdi_psi <- sdi[names(sdi) %in% psi$sci_name]

maps_acc <- subset(acc, acc$presence %in% c(1) & acc$origin %in% c(1,2) & 
                     acc$seasonal %in% c(1,2) & acc$sci_name %in% names(sdi_acc))
maps_ans <- subset(ans, ans$presence %in% c(1) & ans$origin %in% c(1,2) & 
                     ans$seasonal %in% c(1,2) & ans$sci_name %in% names(sdi_ans))
maps_apo <- subset(apo, apo$presence %in% c(1) & apo$origin %in% c(1,2) & 
                     apo$seasonal %in% c(1,2) & apo$sci_name %in% names(sdi_apo))
maps_cha <- subset(cha, cha$presence %in% c(1) & cha$origin %in% c(1,2) & 
                     cha$seasonal %in% c(1,2) & cha$sci_name %in% names(sdi_cha))
maps_col <- subset(col, col$presence %in% c(1) & col$origin %in% c(1,2) & 
                     col$seasonal %in% c(1,2) & col$sci_name %in% names(sdi_col))
maps_gal <- subset(gal, gal$presence %in% c(1) & gal$origin %in% c(1,2) & 
                     gal$seasonal %in% c(1,2) & gal$sci_name %in% names(sdi_gal))
#maps_pas <- subset(pas, pas$presence %in% c(1) & pas$origin %in% c(1,2) & 
#                     pas$seasonal %in% c(1,2) & pas$sci_name %in% names(sdi_pas))
maps_pic <- subset(pic, pic$presence %in% c(1) & pic$origin %in% c(1,2) & 
                     pic$seasonal %in% c(1,2) & pic$sci_name %in% names(sdi_pic))
maps_psi <- subset(psi, psi$presence %in% c(1) & psi$origin %in% c(1,2) & 
                     psi$seasonal %in% c(1,2) & psi$sci_name %in% names(sdi_psi))

sdi_acc <- sdi_acc[names(sdi_acc) %in% maps_acc$sci_name]
sdi_ans <- sdi_ans[names(sdi_ans) %in% maps_ans$sci_name]
sdi_apo <- sdi_apo[names(sdi_apo) %in% maps_apo$sci_name]
sdi_cha <- sdi_cha[names(sdi_cha) %in% maps_cha$sci_name]
sdi_col <- sdi_col[names(sdi_col) %in% maps_col$sci_name]
sdi_gal <- sdi_gal[names(sdi_gal) %in% maps_gal$sci_name]
#sdi_pas <- sdi_pas[names(sdi_pas) %in% maps_pas$sci_name]
sdi_pic <- sdi_pic[names(sdi_pic) %in% maps_pic$sci_name]
sdi_psi <- sdi_psi[names(sdi_psi) %in% maps_psi$sci_name]

getVarCells <- function(maps, sdi, res = 2.5) {
  
  r <- raster(ncols = 2160, nrows = 900, ymn = -60)
  raster_stack <- r
  raster_stack_rich <- r
  
  for (i in 1:length(sdi)) {

    species_i <- names(sdi)[i]
    
    maps_i <- subset(maps, maps$sci_name == species_i) 
    
    for (j in 1:length(maps_i$Shape)) { 
      try(maps_i$Shape[[j]] <- st_cast(maps_i$Shape[[j]], 'MULTIPOLYGON'))
    }
    try(maps_i$Shape <- st_cast(maps_i$Shape, 'MULTIPOLYGON'))
    
    raster_i <- fasterize(st_as_sf(maps_i$Shape), r)
    raster_i2 <- fasterize(st_as_sf(maps_i$Shape), r)
    
    rastercells <- which(getValues(!is.na(raster_i))) 
    raster_i[rastercells] <- sdi[i]
    
    if (i == 1) {
      raster_stack <- raster_i
      raster_stack_rich <- raster_i2
    } else {
      raster_stack <- addLayer(raster_stack, raster_i)
      raster_stack_rich <- addLayer(raster_stack_rich, raster_i2)
    }
  }
  
  av_full <- stackApply(raster_stack, indices = rep(1, length(sdi)), 
                            fun = median, na.rm = T) 
  av_full_rich <- stackApply(raster_stack_rich, indices = rep(1, length(sdi)), 
                                 fun = sum, na.rm = T) 
  
  envar <- getData("worldclim", var = "bio", res = res) 
  envars <- stack(envar[[4]], envar[[15]])
  
  npp <- raster("data/npp-geotiff/npp_geotiff.tif")
  
  r <- raster(ncols = 2160, nrows = 900, ymn = -60)
  r[] <- NA 
  r[] <- 1:length(values(r))
  
  env_vals <- raster::extract(envars, r[])
  npp_vals <- raster::extract(npp, r[])
  sdi_vals <- raster::extract(av_full, r[])
  ric_vals <- raster::extract(av_full_rich, r[])
  
  cells <- which(getValues(av_full_rich) > 0)
  
  lon <- raster::xFromCell(av_full_rich, cell = cells)
  lat <- raster::yFromCell(av_full_rich, cell = cells)
  
  cells_stats <- as.data.frame(cbind(sdi_vals[cells],
                                     ric_vals[cells],
                                     lat, 
                                     lon,
                                     env_vals[cells, ],
                                     npp_vals[cells]))
  colnames(cells_stats) <- c("SDI", "Richness", "Latitude", "Longitude", "bio4", 
                             "bio15", "NPP")
  
  res <- list()
  res$stats <- cells_stats
  res$sdi <- av_full
  res$rich <- av_full_rich
  
  return(res)
  
}

getVarSpp <- function(maps, sdi, taxon, env, npp) {
  
  r <- raster(ncols = 8640, nrows = 3600, ymn = -60, vals = 1:(8640*3600))

  env_vals <- env
  npp_vals <- npp
  
  for(i in 1:length(unique(maps$sci_name))) try({
    
    print(paste0(i, "/", length(unique(maps$sci_name))))

    species_i <- unique(maps$sci_name)[i]
  
    maps_i <- maps %>% filter(sci_name == species_i)
  
    for (j in 1:length(maps_i$Shape)) { 
      try(maps_i$Shape[j] <- st_make_valid(maps_i$Shape[j]))
      try(maps_i$Shape[[j]] <- st_cast(maps_i$Shape[[j]], 'MULTIPOLYGON'))
    }
    try(maps_i$Shape <- st_cast(maps_i$Shape, 'MULTIPOLYGON'))
  
    raster_i <- fasterize(st_as_sf(maps_i$Shape), r) 

    lat <- gCentroid(as(maps_i$Shape, "Spatial"))@coords[2]
    
    rastercells <- which(getValues(raster_i) > 0)
  
    values_env <- env_vals[rastercells, ]
    values_npp <- npp_vals[rastercells]
    
    if(nrow(values_env) == 1) {
      ifelse(is.na(values_env[1]), 
             bio4 <- NA, 
             bio4 <- values_env[1])
      ifelse(is.na(values_env[2]), 
             bio15 <- NA, 
             bio15 <- values_env[2])
    } else { 
      ifelse(all(is.na(values_env[, 1])),
             bio4 <- NA,
             bio4 <- mean(values_env[, 1], na.rm = TRUE))
      ifelse(all(is.na(values_env[, 2])), 
             bio15 <- NA,
             bio15 <- mean(values_env[, 2], na.rm = TRUE))
    }
  
    if(length(values_npp) == 1) {
      ifelse(is.na(values_npp[1]), 
             npp <- NA, 
             npp <- values_npp[1])
    } else { 
      ifelse(all(is.na(values_npp)),
             npp <- NA,
             npp <- mean(values_npp, na.rm = TRUE))
    }
  
    study_results <- data.frame(species = maps_i$sci_name[1],
                                sdi[names(sdi) == species_i][[1]],
                                lat,
                                bio4,
                                bio15,
                                npp)
  
    write.table(study_results, paste0("data/data_env_per_spp_", taxon, ".csv"),
                sep = ",", col.names = FALSE, append = TRUE, row.names = F)
  })
  
  res <- read.csv(paste0("data/data_env_per_spp_", taxon, ".csv"),
                  row.names = 1, header = F)
  colnames(res) <- c("SDI", "Latitude", "bio4", "bio15", "NPP")
  write.csv(res, paste0("data/data_env_per_spp_", taxon, ".csv"))
  return(res)
}

envar <- getData("worldclim", var = "bio", res = 2.5)
envars <- stack(envar[[4]], envar[[15]])

npp <- raster("data/npp-geotiff/npp_geotiff.tif")

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

av_full <- stackApply(raster_stack, indices = rep(1, length(sdi)), 
                      fun = median, na.rm = T) 

npp <- stack(list.files(path = "data/npp_Wang&al2021/NPP1981", 
                        pattern = '.tif$', all.files = TRUE, full.names = TRUE))
save(npp, file = "data/npp_Wang&al2021/test.RData")

r <- raster(ncols = 8640, nrows = 3600, ymn = -60, vals = 1:(8640*3600))

env_vals <- raster::extract(envars, r[])
npp_vals <- raster::extract(npp, r[])

env_acc <- getVarSpp(maps_acc, sdi_acc, "Accipitriformes", env = env_vals, 
                     npp = npp_vals)
env_ans <- getVarSpp(maps_ans, sdi_ans, "Anseriformes", env = env_vals, 
                     npp = npp_vals)
env_apo <- getVarSpp(maps_apo, sdi_apo, "Apodiformes", env = env_vals, 
                     npp = npp_vals)
env_cha <- getVarSpp(maps_cha, sdi_cha, "Charadriiformes", env = env_vals, 
                     npp = npp_vals)
env_col <- getVarSpp(maps_col, sdi_col, "Columbiformes", env = env_vals, 
                     npp = npp_vals)
env_gal <- getVarSpp(maps_gal, sdi_gal, "Galliformes", env = env_vals, 
                     npp = npp_vals)
#env_pas <- getVarSpp(maps_pas, sdi_pas, "Passeriformes", env = env_vals, 
#                     npp = npp_vals)
env_pic <- getVarSpp(maps_pic, sdi_pic, "Piciformes", env = env_vals, 
                     npp = npp_vals)
env_psi <- getVarSpp(maps_psi, sdi_psi, "Psittaciformes", env = env_vals, 
                     npp = npp_vals)

var_acc <- getVarCells(maps_acc, sdi_acc)
var_ans <- getVarCells(maps_ans, sdi_ans)
var_apo <- getVarCells(maps_apo, sdi_apo)
var_cha <- getVarCells(maps_cha, sdi_cha)
var_col <- getVarCells(maps_col, sdi_col)
var_gal <- getVarCells(maps_gal, sdi_gal)
#var_pas <- getVarCells(maps_pas, sdi_pas)
var_pic <- getVarCells(maps_pic, sdi_pic)
var_psi <- getVarCells(maps_psi, sdi_psi)

#save(var_acc, file = "data/var_acc.RData")
#save(var_ans, file = "data/var_ans.RData")
#save(var_apo, file = "data/var_apo.RData")
#save(var_cha, file = "data/var_cha.RData")
#save(var_col, file = "data/var_col.RData")
#save(var_gal, file = "data/var_gal.RData")
#save(var_pas, file = "data/var_pas.RData")
#save(var_pic, file = "data/var_pic.RData")
#save(var_psi, file = "data/var_psi.RData")

#load(file = "data/var_acc.RData")
#load(file = "data/var_ans.RData")
#load(file = "data/var_apo.RData")
#load(file = "data/var_cha.RData")
#load(file = "data/var_col.RData")
#load(file = "data/var_gal.RData")
#load(file = "data/var_pas.RData")
#load(file = "data/var_pic.RData")
#load(file = "data/var_psi.RData")

pdf("figures/Figure8.pdf", width = 8)

layout(matrix(1:9, ncol = 3, byrow = TRUE))
par(mar = c(4, 4, 2, 1))
cols1 <- rgb(68/255, 1/255, 84/255, 0.2)
plot(abs(SDI) ~ Richness, data = var_acc$stats, pch = 16, col = cols1,
     ylab = "Absolute SDI", xlab = "")
title("Accipitriformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_ans$stats, pch = 16, col = cols1,
     ylab = "", xlab = "")
title("Anseriformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_apo$stats, pch = 16, col = cols1,
     ylab = "", xlab = "")
title("Apodiformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_cha$stats, pch = 16, col = cols1,
     ylab = "Absolute SDI", xlab = "")
title("Charadriiformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_col$stats, pch = 16, col = cols1,
     ylab = "", xlab = "")
title("Columbiformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_gal$stats, pch = 16, col = cols1,
     ylab = "", xlab = "")
title("Galliformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_pas$stats, pch = 16, col = cols1,
     ylab = "Absolute SDI", xlab = "Richness")
title("Passeriformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_pic$stats, pch = 16, col = cols1,
     ylab = "", xlab = "Richness")
title("Piciformes", adj = 0)
plot(abs(SDI) ~ Richness, data = var_psi$stats, pch = 16, col = cols1,
     ylab = "", xlab = "Richness")
title("Psittaciformes", adj = 0)

dev.off()

pdf("figures/Figure9.pdf", width = 8)

layout(matrix(1:9, ncol = 3, byrow = TRUE))
par(mar = c(4, 4, 2, 1))
cols2 <- rgb(49/255, 104/255, 142/255, 0.2)
plot(abs(SDI) ~ NPP, data = var_acc$stats, pch = 16, col = cols2,
     ylab = "Absolute SDI", xlab = "")
title("Accipitriformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_ans$stats, pch = 16, col = cols2,
     ylab = "", xlab = "")
title("Anseriformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_apo$stats, pch = 16, col = cols2,
     ylab = "", xlab = "")
title("Apodiformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_cha$stats, pch = 16, col = cols2,
     ylab = "Absolute SDI", xlab = "")
title("Charadriiformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_col$stats, pch = 16, col = cols2,
     ylab = "", xlab = "")
title("Columbiformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_gal$stats, pch = 16, col = cols2,
     ylab = "", xlab = "")
title("Galliformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_pas$stats, pch = 16, col = cols2,
     ylab = "Absolute SDI", xlab = "NPP")
title("Passeriformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_pic$stats, pch = 16, col = cols2,
     ylab = "", xlab = "NPP")
title("Piciformes", adj = 0)
plot(abs(SDI) ~ NPP, data = var_psi$stats, pch = 16, col = cols2,
     ylab = "", xlab = "NPP")
title("Psittaciformes", adj = 0)

dev.off()

pdf("figures/Figure10.pdf", width = 8)

layout(matrix(1:9, ncol = 3, byrow = TRUE))
par(mar = c(4, 4, 2, 1))
cols3 <- rgb(31/255, 158/255, 137/255, 0.2)
plot(SDI ~ bio4, data = var_acc$stats, pch = 16, col = cols3,
     ylab = "SDI", xlab = "")
title("Accipitriformes", adj = 0)
plot(SDI ~ bio4, data = var_ans$stats, pch = 16, col = cols3,
     ylab = "", xlab = "")
title("Anseriformes", adj = 0)
plot(SDI ~ bio4, data = var_apo$stats, pch = 16, col = cols3,
     ylab = "", xlab = "")
title("Apodiformes", adj = 0)
plot(SDI ~ bio4, data = var_cha$stats, pch = 16, col = cols3,
     ylab = "SDI", xlab = "")
title("Charadriiformes", adj = 0)
plot(SDI ~ bio4, data = var_col$stats, pch = 16, col = cols3,
     ylab = "", xlab = "")
title("Columbiformes", adj = 0)
plot(SDI ~ bio4, data = var_gal$stats, pch = 16, col = cols3,
     ylab = "", xlab = "")
title("Galliformes", adj = 0)
plot(SDI ~ bio4, data = var_pas$stats, pch = 16, col = cols3,
     ylab = "SDI", xlab = "Temperature Seasonality")
title("Passeriformes", adj = 0)
plot(SDI ~ bio4, data = var_pic$stats, pch = 16, col = cols3,
     ylab = "", xlab = "Temperature Seasonality")
title("Piciformes", adj = 0)
plot(SDI ~ bio4, data = var_psi$stats, pch = 16, col = cols3,
     ylab = "", xlab = "Temperature Seasonality")
title("Psittaciformes", adj = 0)

dev.off()

pdf("figures/Figure11.pdf", width = 8)

layout(matrix(1:9, ncol = 3, byrow = TRUE))
par(mar = c(4, 4, 2, 1))
cols4 <- rgb(180/255, 222/225, 44/255, 0.2)
plot(SDI ~ bio15, data = var_acc$stats, pch = 16, col = cols4,
     ylab = "SDI", xlab = "")
title("Accipitriformes", adj = 0)
plot(SDI ~ bio15, data = var_ans$stats, pch = 16, col = cols4,
     ylab = "", xlab = "")
title("Anseriformes", adj = 0)
plot(SDI ~ bio15, data = var_apo$stats, pch = 16, col = cols4,
     ylab = "", xlab = "")
title("Apodiformes", adj = 0)
plot(SDI ~ bio15, data = var_cha$stats, pch = 16, col = cols4,
     ylab = "SDI", xlab = "")
title("Charadriiformes", adj = 0)
plot(SDI ~ bio15, data = var_col$stats, pch = 16, col = cols4,
     ylab = "", xlab = "")
title("Columbiformes", adj = 0)
plot(SDI ~ bio15, data = var_gal$stats, pch = 16, col = cols4,
     ylab = "", xlab = "")
title("Galliformes", adj = 0)
plot(SDI ~ bio15, data = var_pas$stats, pch = 16, col = cols4,
     ylab = "SDI", xlab = "Precipitation Seasonality")
title("Passeriformes", adj = 0)
plot(SDI ~ bio15, data = var_pic$stats, pch = 16, col = cols4,
     ylab = "", xlab = "Precipitation Seasonality")
title("Piciformes", adj = 0)
plot(SDI ~ bio15, data = var_psi$stats, pch = 16, col = cols4,
     ylab = "", xlab = "Precipitation Seasonality")
title("Psittaciformes", adj = 0)

dev.off()

mod1_acc <- lm(abs(SDI) ~ Richness + NPP, data = var_acc$stats)
mod1_ans <- lm(abs(SDI) ~ Richness + NPP, data = var_ans$stats)
mod1_apo <- lm(abs(SDI) ~ Richness + NPP, data = var_apo$stats)
mod1_cha <- lm(abs(SDI) ~ Richness + NPP, data = var_cha$stats)
mod1_col <- lm(abs(SDI) ~ Richness + NPP, data = var_col$stats)
mod1_gal <- lm(abs(SDI) ~ Richness + NPP, data = var_gal$stats)
#mod1_pas <- lm(abs(SDI) ~ Richness + NPP, data = var_pas$stats)
mod1_pic <- lm(abs(SDI) ~ Richness + NPP, data = var_pic$stats)
mod1_psi <- lm(abs(SDI) ~ Richness + NPP, data = var_psi$stats)

mod2_acc <- lm(SDI ~ bio4 + bio15, data = var_acc$stats)
mod2_ans <- lm(SDI ~ bio4 + bio15, data = var_ans$stats)
mod2_apo <- lm(SDI ~ bio4 + bio15, data = var_apo$stats)
mod2_cha <- lm(SDI ~ bio4 + bio15, data = var_cha$stats)
mod2_col <- lm(SDI ~ bio4 + bio15, data = var_col$stats)
mod2_gal <- lm(SDI ~ bio4 + bio15, data = var_gal$stats)
#mod2_pas <- lm(SDI ~ bio4 + bio15, data = var_pas$stats)
mod2_pic <- lm(SDI ~ bio4 + bio15, data = var_pic$stats)
mod2_psi <- lm(SDI ~ bio4 + bio15, data = var_psi$stats)

mod3_acc <- lm(SDI ~ Latitude, data = var_acc$stats)
mod3_ans <- lm(SDI ~ Latitude, data = var_ans$stats)
mod3_apo <- lm(SDI ~ Latitude, data = var_apo$stats)
mod3_cha <- lm(SDI ~ Latitude, data = var_cha$stats)
mod3_col <- lm(SDI ~ Latitude, data = var_col$stats)
mod3_gal <- lm(SDI ~ Latitude, data = var_gal$stats)
#mod3_pas <- lm(SDI ~ Latitude, data = var_pas$stats)
mod3_pic <- lm(SDI ~ Latitude, data = var_pic$stats)
mod3_psi <- lm(SDI ~ Latitude, data = var_psi$stats)

subset_mal_acc <- var_acc$stats[var_acc$stats$SDI < 0, ]
subset_mal_ans <- var_ans$stats[var_ans$stats$SDI < 0, ]
subset_mal_apo <- var_apo$stats[var_apo$stats$SDI < 0, ]
subset_mal_cha <- var_cha$stats[var_cha$stats$SDI < 0, ]
subset_mal_col <- var_col$stats[var_col$stats$SDI < 0, ]
subset_mal_gal <- var_gal$stats[var_gal$stats$SDI < 0, ]
#subset_mal_pas <- var_pas$stats[var_pas$stats$SDI < 0, ]
subset_mal_pic <- var_pic$stats[var_pic$stats$SDI < 0, ]
subset_mal_psi <- var_psi$stats[var_psi$stats$SDI < 0, ]

subset_fem_acc <- var_acc$stats[var_acc$stats$SDI > 0, ]
subset_fem_ans <- var_ans$stats[var_ans$stats$SDI > 0, ]
subset_fem_apo <- var_apo$stats[var_apo$stats$SDI > 0, ]
subset_fem_cha <- var_cha$stats[var_cha$stats$SDI > 0, ]
subset_fem_col <- var_col$stats[var_col$stats$SDI > 0, ]
subset_fem_gal <- var_gal$stats[var_gal$stats$SDI > 0, ]
#subset_fem_pas <- var_pas$stats[var_pas$stats$SDI > 0, ]
subset_fem_pic <- var_pic$stats[var_pic$stats$SDI > 0, ]
subset_fem_psi <- var_psi$stats[var_psi$stats$SDI > 0, ]

mod4_acc <- lm(SDI ~ Latitude, data = subset_mal_acc)
mod4_ans <- lm(SDI ~ Latitude, data = subset_mal_ans)
mod4_apo <- lm(SDI ~ Latitude, data = subset_mal_apo)
mod4_cha <- lm(SDI ~ Latitude, data = subset_mal_cha)
mod4_col <- lm(SDI ~ Latitude, data = subset_mal_col)
mod4_gal <- lm(SDI ~ Latitude, data = subset_mal_gal)
#mod4_pas <- lm(SDI ~ Latitude, data = subset_mal_pas)
mod4_pic <- lm(SDI ~ Latitude, data = subset_mal_pic)
mod4_psi <- lm(SDI ~ Latitude, data = subset_mal_psi)

mod5_acc <- lm(SDI ~ Latitude, data = subset_fem_acc)
mod5_ans <- lm(SDI ~ Latitude, data = subset_fem_ans)
mod5_apo <- lm(SDI ~ Latitude, data = subset_fem_apo)
mod5_cha <- lm(SDI ~ Latitude, data = subset_fem_cha)
mod5_col <- lm(SDI ~ Latitude, data = subset_fem_col)
mod5_gal <- lm(SDI ~ Latitude, data = subset_fem_gal)
#mod5_pas <- lm(SDI ~ Latitude, data = subset_fem_pas)
mod5_pic <- lm(SDI ~ Latitude, data = subset_fem_pic)
mod5_psi <- lm(SDI ~ Latitude, data = subset_fem_psi)

res_lm <- matrix (nrow = 9, ncol = 14)
colnames(res_lm) <- c("Slope1_mod1", "p1_mod1", "Slope2_mod1", "p2_mod1",
                      "Slope1_mod2", "p1_mod2", "Slope2_mod2", "p2_mod2",
                      "Slope_mod3", "p_mod3",
                      "Slope_mod4", "p_mod4",
                      "Slope_mod5", "p_mod5")
rownames(res_lm) <- c("Accipitriformes", "Anseriformes", "Apodiformes",
                      "Charadriiformes", "Columbiformes", "Galliformes",
                      "Passeriformes", "Piciformes", "Psittaciformes")

res_lm[1, 1:2] <- round(summary(mod1_acc)$coefficients[2, c(1, 4)], 3)
res_lm[1, 3:4] <- round(summary(mod1_acc)$coefficients[3, c(1, 4)], 3)
res_lm[1, 5:6] <- round(summary(mod2_acc)$coefficients[2, c(1, 4)], 3)
res_lm[1, 7:8] <- round(summary(mod2_acc)$coefficients[3, c(1, 4)], 3)
res_lm[1, 9:10] <- round(summary(mod3_acc)$coefficients[2, c(1, 4)], 3)
res_lm[1, 11:12] <- round(summary(mod4_acc)$coefficients[2, c(1, 4)], 3)
res_lm[1, 13:14] <- round(summary(mod5_acc)$coefficients[2, c(1, 4)], 3)

res_lm[2, 1:2] <- round(summary(mod1_ans)$coefficients[2, c(1, 4)], 3)
res_lm[2, 3:4] <- round(summary(mod1_ans)$coefficients[3, c(1, 4)], 3)
res_lm[2, 5:6] <- round(summary(mod2_ans)$coefficients[2, c(1, 4)], 3)
res_lm[2, 7:8] <- round(summary(mod2_ans)$coefficients[3, c(1, 4)], 3)
res_lm[2, 9:10] <- round(summary(mod3_ans)$coefficients[2, c(1, 4)], 3)
res_lm[2, 11:12] <- round(summary(mod4_ans)$coefficients[2, c(1, 4)], 3)
res_lm[2, 13:14] <- round(summary(mod5_ans)$coefficients[2, c(1, 4)], 3)

res_lm[3, 1:2] <- round(summary(mod1_apo)$coefficients[2, c(1, 4)], 3)
res_lm[3, 3:4] <- round(summary(mod1_apo)$coefficients[3, c(1, 4)], 3)
res_lm[3, 5:6] <- round(summary(mod2_apo)$coefficients[2, c(1, 4)], 3)
res_lm[3, 7:8] <- round(summary(mod2_apo)$coefficients[3, c(1, 4)], 3)
res_lm[3, 9:10] <- round(summary(mod3_apo)$coefficients[2, c(1, 4)], 3)
res_lm[3, 11:12] <- round(summary(mod4_apo)$coefficients[2, c(1, 4)], 3)
res_lm[3, 13:14] <- round(summary(mod5_apo)$coefficients[2, c(1, 4)], 3)

res_lm[4, 1:2] <- round(summary(mod1_cha)$coefficients[2, c(1, 4)], 3)
res_lm[4, 3:4] <- round(summary(mod1_cha)$coefficients[3, c(1, 4)], 3)
res_lm[4, 5:6] <- round(summary(mod2_cha)$coefficients[2, c(1, 4)], 3)
res_lm[4, 7:8] <- round(summary(mod2_cha)$coefficients[3, c(1, 4)], 3)
res_lm[4, 9:10] <- round(summary(mod3_cha)$coefficients[2, c(1, 4)], 3)
res_lm[4, 11:12] <- round(summary(mod4_cha)$coefficients[2, c(1, 4)], 3)
res_lm[4, 13:14] <- round(summary(mod5_cha)$coefficients[2, c(1, 4)], 3)

res_lm[5, 1:2] <- round(summary(mod1_col)$coefficients[2, c(1, 4)], 3)
res_lm[5, 3:4] <- round(summary(mod1_col)$coefficients[3, c(1, 4)], 3)
res_lm[5, 5:6] <- round(summary(mod2_col)$coefficients[2, c(1, 4)], 3)
res_lm[5, 7:8] <- round(summary(mod2_col)$coefficients[3, c(1, 4)], 3)
res_lm[5, 9:10] <- round(summary(mod3_col)$coefficients[2, c(1, 4)], 3)
res_lm[5, 11:12] <- round(summary(mod4_col)$coefficients[2, c(1, 4)], 3)
res_lm[5, 13:14] <- round(summary(mod5_col)$coefficients[2, c(1, 4)], 3)

res_lm[6, 1:2] <- round(summary(mod1_gal)$coefficients[2, c(1, 4)], 3)
res_lm[6, 3:4] <- round(summary(mod1_gal)$coefficients[3, c(1, 4)], 3)
res_lm[6, 5:6] <- round(summary(mod2_gal)$coefficients[2, c(1, 4)], 3)
res_lm[6, 7:8] <- round(summary(mod2_gal)$coefficients[3, c(1, 4)], 3)
res_lm[6, 9:10] <- round(summary(mod3_gal)$coefficients[2, c(1, 4)], 3)
res_lm[6, 11:12] <- round(summary(mod4_gal)$coefficients[2, c(1, 4)], 3)
res_lm[6, 13:14] <- round(summary(mod5_gal)$coefficients[2, c(1, 4)], 3)

#res_lm[7, 1:2] <- round(summary(mod1_pas)$coefficients[2, c(1, 4)], 3)
#res_lm[7, 3:4] <- round(summary(mod1_pas)$coefficients[3, c(1, 4)], 3)
#res_lm[7, 5:6] <- round(summary(mod2_pas)$coefficients[2, c(1, 4)], 3)
#res_lm[7, 7:8] <- round(summary(mod2_pas)$coefficients[3, c(1, 4)], 3)
#res_lm[7, 9:10] <- round(summary(mod3_pas)$coefficients[2, c(1, 4)], 3)
#res_lm[7, 11:12] <- round(summary(mod4_pas)$coefficients[2, c(1, 4)], 3)
#res_lm[7, 13:14] <- round(summary(mod5_pas)$coefficients[2, c(1, 4)], 3)

res_lm[8, 1:2] <- round(summary(mod1_pic)$coefficients[2, c(1, 4)], 3)
res_lm[8, 3:4] <- round(summary(mod1_pic)$coefficients[3, c(1, 4)], 3)
res_lm[8, 5:6] <- round(summary(mod2_pic)$coefficients[2, c(1, 4)], 3)
res_lm[8, 7:8] <- round(summary(mod2_pic)$coefficients[3, c(1, 4)], 3)
res_lm[8, 9:10] <- round(summary(mod3_pic)$coefficients[2, c(1, 4)], 3)
res_lm[8, 11:12] <- round(summary(mod4_pic)$coefficients[2, c(1, 4)], 3)
res_lm[8, 13:14] <- round(summary(mod5_pic)$coefficients[2, c(1, 4)], 3)

res_lm[9, 1:2] <- round(summary(mod1_psi)$coefficients[2, c(1, 4)], 3)
res_lm[9, 3:4] <- round(summary(mod1_psi)$coefficients[3, c(1, 4)], 3)
res_lm[9, 5:6] <- round(summary(mod2_psi)$coefficients[2, c(1, 4)], 3)
res_lm[9, 7:8] <- round(summary(mod2_psi)$coefficients[3, c(1, 4)], 3)
res_lm[9, 9:10] <- round(summary(mod3_psi)$coefficients[2, c(1, 4)], 3)
res_lm[9, 11:12] <- round(summary(mod4_psi)$coefficients[2, c(1, 4)], 3)
res_lm[9, 13:14] <- round(summary(mod5_psi)$coefficients[2, c(1, 4)], 3)

write.csv(res_lm, "tables/Table2_unformatted.csv")

