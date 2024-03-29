rm(list = ls())

setwd("~/Documents/lab/dimorph_evol")

library(sf)
library(raster)
library(dplyr)
library(fasterize)
library(rgeos)
library(phytools)
library(geiger)
library(caper)
library(wesanderson)

tr <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
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
load(file = "data/spatial/pas.RData")
load(file = "data/spatial/pic.RData")
load(file = "data/spatial/psi.RData")

sdi_acc <- sdi[names(sdi) %in% acc$sci_name]
sdi_ans <- sdi[names(sdi) %in% ans$sci_name]
sdi_apo <- sdi[names(sdi) %in% apo$sci_name]
sdi_cha <- sdi[names(sdi) %in% cha$sci_name]
sdi_col <- sdi[names(sdi) %in% col$sci_name]
sdi_gal <- sdi[names(sdi) %in% gal$sci_name]
sdi_pas <- sdi[names(sdi) %in% pas$sci_name]
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
maps_pas <- subset(pas, pas$presence %in% c(1) & pas$origin %in% c(1,2) & 
                     pas$seasonal %in% c(1,2) & pas$sci_name %in% names(sdi_pas))
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
sdi_pas <- sdi_pas[names(sdi_pas) %in% maps_pas$sci_name]
sdi_pic <- sdi_pic[names(sdi_pic) %in% maps_pic$sci_name]
sdi_psi <- sdi_psi[names(sdi_psi) %in% maps_psi$sci_name]

envar <- raster::getData("worldclim", var = "bio", res = 10)
envars <- stack(envar[[4]], envar[[15]])

r <- raster(ncols = 2160, nrows = 900, ymn = -60, vals = 1:(2160*900))

env_vals <- raster::extract(envars, 1:length(r[]))

load("data/spatial/npp_vals_2009-2018_2160x900.RData")

npp_vals <- raster::extract(npp_low, 1:length(r[]))

getVarSpp <- function(maps, sdi, taxon, env, npp) {

     r <- raster(ncols = 2160, nrows = 900, ymn = -60, vals = rep(1, 2160*900))

     raster_rich <- r

     for (i in 1:length(sdi)) {

          print(paste0("1. ", i, "/", length(unique(maps$sci_name))))

          species_i <- names(sdi)[i]
    
          maps_i <- subset(maps, maps$sci_name == species_i) 
    
          for (j in 1:length(maps_i$Shape)) { 
               try(maps_i$Shape[[j]] <- st_cast(maps_i$Shape[[j]], 'MULTIPOLYGON'))
          }
          try(maps_i$Shape <- st_cast(maps_i$Shape, 'MULTIPOLYGON'))
    
          raster_i <- fasterize(st_as_sf(maps_i$Shape), r)
    
          rastercells <- which(getValues(!is.na(raster_i)))

          raster_i[rastercells] <- 1

          if (i == 1) {
               raster_rich <- raster_i
          } else {
               raster_rich <- addLayer(raster_rich, raster_i)
          }
          removeTmpFiles(h=0)
     }
  
     full_rich <- stackApply(raster_rich, indices = rep(1, length(sdi)),
                             fun = sum, na.rm = T)
  
     ric_vals <- raster::extract(full_rich, 1:length(r[]))
  
     env_vals <- env
     npp_vals <- npp
  
     for (i in 1:length(unique(maps$sci_name))) try({
    
          print(paste0("2. ", i, "/", length(unique(maps$sci_name))))

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
          values_ric <- ric_vals[rastercells]

          if (nrow(values_env) == 1) {
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
  
          if (length(values_npp) == 1) {
               ifelse(is.na(values_npp[1]),
                      npp <- NA,
                      npp <- values_npp[1])
          } else { 
               ifelse(all(is.na(values_npp)),
                      npp <- NA,
                      npp <- mean(values_npp, na.rm = TRUE))
          }
  
          if (length(values_ric) == 1) {
               ifelse(is.na(values_ric[1]),
                      ric <- NA,
                      ric <- values_ric[1])
          } else { 
               ifelse(all(is.na(values_ric)),
                      ric <- NA,
                      ric <- max(values_ric, na.rm = TRUE))
          }
          
          study_results <- data.frame(species = maps_i$sci_name[1],
                                      sdi[names(sdi) == species_i][[1]],
                                      lat,
                                      ric,
                                      bio4,
                                      bio15,
                                      npp)
  
          write.table(study_results, 
                      paste0("data/EnvSpp_2160x900_", taxon, ".csv"),
                      sep = ",", col.names = FALSE, append = TRUE, 
                      row.names = F)
     })
  
     res <- read.csv(paste0("data/EnvSpp_2160x900_", taxon, ".csv"),
                     row.names = 1, header = F)
     colnames(res) <- c("SDI", "Latitude", "Richness", "bio4", "bio15", "NPP")
     write.csv(res, paste0("data/EnvSpp_2160x900_", taxon, ".csv"))
     return(res)

}

env_acc <- getVarSpp(maps_acc, sdi_acc, "Accipitriformes", env = env_vals, 
                     npp = npp_vals)
#save(env_acc, file = "data/results/EnvSpp_2160x900_acc.RData")

env_ans <- getVarSpp(maps_ans, sdi_ans, "Anseriformes", env = env_vals, 
                     npp = npp_vals)
#save(env_ans, file = "data/results/EnvSpp_2160x900_ans.RData")

env_apo <- getVarSpp(maps_apo, sdi_apo, "Apodiformes", env = env_vals, 
                     npp = npp_vals)
#save(env_apo, file = "data/results/EnvSpp_2160x900_apo.RData")

env_cha <- getVarSpp(maps_cha, sdi_cha, "Charadriiformes", env = env_vals, 
                     npp = npp_vals)
#save(env_cha, file = "data/results/EnvSpp_2160x900_cha.RData")

env_col <- getVarSpp(maps_col, sdi_col, "Columbiformes", env = env_vals, 
                     npp = npp_vals)
#save(env_col, file = "data/results/EnvSpp_2160x900_col.RData")

env_gal <- getVarSpp(maps_gal, sdi_gal, "Galliformes", env = env_vals, 
                     npp = npp_vals)
#save(env_gal, file = "data/results/EnvSpp_2160x900_gal.RData")

env_pas <- getVarSpp(maps_pas, sdi_pas, "Passeriformes", env = env_vals, 
                     npp = npp_vals)
#save(env_pas, file = "data/results/EnvSpp_2160x900_pas.RData")

env_pic <- getVarSpp(maps_pic, sdi_pic, "Piciformes", env = env_vals, 
                     npp = npp_vals)
#save(env_pic, file = "data/results/EnvSpp_2160x900_pic.RData")

env_psi <- getVarSpp(maps_psi, sdi_psi, "Psittaciformes", env = env_vals, 
                     npp = npp_vals)
#save(env_psi, file = "data/results/EnvSpp_2160x900_psi.RData")

#load(file = "data/results/EnvSpp_2160x900_acc.RData")
#load(file = "data/results/EnvSpp_2160x900_ans.RData")
#load(file = "data/results/EnvSpp_2160x900_apo.RData")
#load(file = "data/results/EnvSpp_2160x900_cha.RData")
#load(file = "data/results/EnvSpp_2160x900_col.RData")
#load(file = "data/results/EnvSpp_2160x900_gal.RData")
#load(file = "data/results/EnvSpp_2160x900_pas.RData")
#load(file = "data/results/EnvSpp_2160x900_pic.RData")
#load(file = "data/results/EnvSpp_2160x900_psi.RData")

env_acc <- env_acc[complete.cases(env_acc), ]
env_ans <- env_ans[complete.cases(env_ans), ]
env_apo <- env_apo[complete.cases(env_apo), ]
env_cha <- env_cha[complete.cases(env_cha), ]
env_col <- env_col[complete.cases(env_col), ]
env_gal <- env_gal[complete.cases(env_gal), ]
env_pas <- env_pas[complete.cases(env_pas), ]
env_pic <- env_pic[complete.cases(env_pic), ]
env_psi <- env_psi[complete.cases(env_psi), ]

env_acc <- cbind(rownames(env_acc), env_acc)
env_ans <- cbind(rownames(env_ans), env_ans)
env_apo <- cbind(rownames(env_apo), env_apo)
env_cha <- cbind(rownames(env_cha), env_cha)
env_col <- cbind(rownames(env_col), env_col)
env_gal <- cbind(rownames(env_gal), env_gal)
env_pas <- cbind(rownames(env_pas), env_pas)
env_pic <- cbind(rownames(env_pic), env_pic)
env_psi <- cbind(rownames(env_psi), env_psi)

colnames(env_acc)[1] <- colnames(env_ans)[1] <- colnames(env_apo)[1] <-
    colnames(env_cha)[1] <- colnames(env_col)[1] <- colnames(env_gal)[1] <-
        colnames(env_pas)[1] <- colnames(env_pic)[1] <- colnames(env_psi)[1] <-
            "Species"

res_acc <- res_ans <- res_apo <- res_cha <- res_col <- res_gal <- 
    res_pas <- res_pic <- res_psi <- matrix(nrow = 1000, ncol = 12)

colnames(res_acc) <- colnames(res_ans) <- colnames(res_apo) <- 
    colnames(res_cha) <- colnames(res_col) <- colnames(res_gal) <- 
        colnames(res_pas) <- colnames(res_pic) <- colnames(res_psi) <- 
            c("Intercept", "Richness", "p_Richness", "Latitude", "p_Latitude", 
              "bio4", "p_bio4", "bio15", "p_bio15", "NPP", "p_NPP", 
              "R_squared")

pgls_acc <- pgls_ans <- pgls_apo <- pgls_cha <- pgls_col <- pgls_gal <-
  pgls_pas <- pgls_pic <- pgls_psi <- list()

for (i in 1:1000) {
    print(i)
    pruned_acc <- treedata(tr[[i]], env_acc, warnings = F)$phy
    comp_acc <- comparative.data(pruned_acc, env_acc, Species)
    pgls_acc[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_acc)
    
    summ_acc <- summary(pgls_acc[[i]])
    res_acc[i, 1] <- round(summ_acc$coefficients[1, 1], 3)
    res_acc[i, 2] <- round(summ_acc$coefficients[2, 1], 3)
    res_acc[i, 3] <- round(summ_acc$coefficients[2, 4], 3)
    res_acc[i, 4] <- round(summ_acc$coefficients[3, 1], 3)
    res_acc[i, 5] <- round(summ_acc$coefficients[3, 4], 3)
    res_acc[i, 6] <- round(summ_acc$coefficients[4, 1], 3)
    res_acc[i, 7] <- round(summ_acc$coefficients[4, 4], 3)
    res_acc[i, 8] <- round(summ_acc$coefficients[5, 1], 3)
    res_acc[i, 9] <- round(summ_acc$coefficients[5, 4], 3)
    res_acc[i, 10] <- round(summ_acc$coefficients[6, 1], 3)
    res_acc[i, 11] <- round(summ_acc$coefficients[6, 4], 3)
    res_acc[i, 12] <- round(summ_acc$r.squared, 3)
}
#write.csv(res_acc, "tables/Table4_unformatted_acc.csv")
#save(pgls_acc, file = "data/results/pgls_acc.RData")

for (i in 1:1000) {
    print(i)
    pruned_ans <- treedata(tr[[i]], env_ans, warnings = F)$phy
    comp_ans <- comparative.data(pruned_ans, env_ans, Species)
    pgls_ans[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_ans)
    
    summ_ans <- summary(pgls_ans[[i]])
    res_ans[i, 1] <- round(summ_ans$coefficients[1, 1], 3)
    res_ans[i, 2] <- round(summ_ans$coefficients[2, 1], 3)
    res_ans[i, 3] <- round(summ_ans$coefficients[2, 4], 3)
    res_ans[i, 4] <- round(summ_ans$coefficients[3, 1], 3)
    res_ans[i, 5] <- round(summ_ans$coefficients[3, 4], 3)
    res_ans[i, 6] <- round(summ_ans$coefficients[4, 1], 3)
    res_ans[i, 7] <- round(summ_ans$coefficients[4, 4], 3)
    res_ans[i, 8] <- round(summ_ans$coefficients[5, 1], 3)
    res_ans[i, 9] <- round(summ_ans$coefficients[5, 4], 3)
    res_ans[i, 10] <- round(summ_ans$coefficients[6, 1], 3)
    res_ans[i, 11] <- round(summ_ans$coefficients[6, 4], 3)
    res_ans[i, 12] <- round(summ_ans$r.squared, 3)
}
#write.csv(res_ans, "tables/Table4_unformatted_ans.csv")
#save(pgls_ans, file = "data/results/pgls_ans.RData")

for (i in 1:1000) {
    print(i)
    pruned_apo <- treedata(tr[[i]], env_apo, warnings = F)$phy
    comp_apo <- comparative.data(pruned_apo, env_apo, Species)
    pgls_apo[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_apo)
    
    summ_apo <- summary(pgls_apo[[i]])
    res_apo[i, 1] <- round(summ_apo$coefficients[1, 1], 3)
    res_apo[i, 2] <- round(summ_apo$coefficients[2, 1], 3)
    res_apo[i, 3] <- round(summ_apo$coefficients[2, 4], 3)
    res_apo[i, 4] <- round(summ_apo$coefficients[3, 1], 3)
    res_apo[i, 5] <- round(summ_apo$coefficients[3, 4], 3)
    res_apo[i, 6] <- round(summ_apo$coefficients[4, 1], 3)
    res_apo[i, 7] <- round(summ_apo$coefficients[4, 4], 3)
    res_apo[i, 8] <- round(summ_apo$coefficients[5, 1], 3)
    res_apo[i, 9] <- round(summ_apo$coefficients[5, 4], 3)
    res_apo[i, 10] <- round(summ_apo$coefficients[6, 1], 3)
    res_apo[i, 11] <- round(summ_apo$coefficients[6, 4], 3)
    res_apo[i, 12] <- round(summ_apo$r.squared, 3)
}
#write.csv(res_apo, "tables/Table4_unformatted_apo.csv")
#save(pgls_apo, file = "data/results/pgls_apo.RData")

for (i in 1:1000) {
    print(i)
    pruned_cha <- treedata(tr[[i]], env_cha, warnings = F)$phy
    comp_cha <- comparative.data(pruned_cha, env_cha, Species)
    pgls_cha[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_cha)
    
    summ_cha <- summary(pgls_cha[[i]])
    res_cha[i, 1] <- round(summ_cha$coefficients[1, 1], 3)
    res_cha[i, 2] <- round(summ_cha$coefficients[2, 1], 3)
    res_cha[i, 3] <- round(summ_cha$coefficients[2, 4], 3)
    res_cha[i, 4] <- round(summ_cha$coefficients[3, 1], 3)
    res_cha[i, 5] <- round(summ_cha$coefficients[3, 4], 3)
    res_cha[i, 6] <- round(summ_cha$coefficients[4, 1], 3)
    res_cha[i, 7] <- round(summ_cha$coefficients[4, 4], 3)
    res_cha[i, 8] <- round(summ_cha$coefficients[5, 1], 3)
    res_cha[i, 9] <- round(summ_cha$coefficients[5, 4], 3)
    res_cha[i, 10] <- round(summ_cha$coefficients[6, 1], 3)
    res_cha[i, 11] <- round(summ_cha$coefficients[6, 4], 3)
    res_cha[i, 12] <- round(summ_cha$r.squared, 3)
}
#write.csv(res_cha, "tables/Table4_unformatted_cha.csv")
#save(pgls_cha, file = "data/results/pgls_cha.RData")

for (i in 1:1000) {
    print(i)
    pruned_col <- treedata(tr[[i]], env_col, warnings = F)$phy
    comp_col <- comparative.data(pruned_col, env_col, Species)
    pgls_col[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_col)
    
    summ_col <- summary(pgls_col[[i]])
    res_col[i, 1] <- round(summ_col$coefficients[1, 1], 3)
    res_col[i, 2] <- round(summ_col$coefficients[2, 1], 3)
    res_col[i, 3] <- round(summ_col$coefficients[2, 4], 3)
    res_col[i, 4] <- round(summ_col$coefficients[3, 1], 3)
    res_col[i, 5] <- round(summ_col$coefficients[3, 4], 3)
    res_col[i, 6] <- round(summ_col$coefficients[4, 1], 3)
    res_col[i, 7] <- round(summ_col$coefficients[4, 4], 3)
    res_col[i, 8] <- round(summ_col$coefficients[5, 1], 3)
    res_col[i, 9] <- round(summ_col$coefficients[5, 4], 3)
    res_col[i, 10] <- round(summ_col$coefficients[6, 1], 3)
    res_col[i, 11] <- round(summ_col$coefficients[6, 4], 3)
    res_col[i, 12] <- round(summ_col$r.squared, 3)
}
#write.csv(res_col, "tables/Table4_unformatted_col.csv")
#save(pgls_col, file = "data/results/pgls_col.RData")

for (i in 1:1000) {
    print(i)
    pruned_gal <- treedata(tr[[i]], env_gal, warnings = F)$phy
    comp_gal <- comparative.data(pruned_gal, env_gal, Species)
    pgls_gal[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_gal)
    
    summ_gal <- summary(pgls_gal[[i]])
    res_gal[i, 1] <- round(summ_gal$coefficients[1, 1], 3)
    res_gal[i, 2] <- round(summ_gal$coefficients[2, 1], 3)
    res_gal[i, 3] <- round(summ_gal$coefficients[2, 4], 3)
    res_gal[i, 4] <- round(summ_gal$coefficients[3, 1], 3)
    res_gal[i, 5] <- round(summ_gal$coefficients[3, 4], 3)
    res_gal[i, 6] <- round(summ_gal$coefficients[4, 1], 3)
    res_gal[i, 7] <- round(summ_gal$coefficients[4, 4], 3)
    res_gal[i, 8] <- round(summ_gal$coefficients[5, 1], 3)
    res_gal[i, 9] <- round(summ_gal$coefficients[5, 4], 3)
    res_gal[i, 10] <- round(summ_gal$coefficients[6, 1], 3)
    res_gal[i, 11] <- round(summ_gal$coefficients[6, 4], 3)
    res_gal[i, 12] <- round(summ_gal$r.squared, 3)
}
#write.csv(res_gal, "tables/Table4_unformatted_gal.csv")
#save(pgls_gal, file = "data/results/pgls_gal.RData")

for (i in 1:1000) {
    print(i)
    pruned_pas <- treedata(tr[[i]], env_pas, warnings = F)$phy
    comp_pas <- comparative.data(pruned_pas, env_pas, Species)
    pgls_pas[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_pas)
    
    summ_pas <- summary(pgls_pas[[i]])
    res_pas[i, 1] <- round(summ_pas$coefficients[1, 1], 3)
    res_pas[i, 2] <- round(summ_pas$coefficients[2, 1], 3)
    res_pas[i, 3] <- round(summ_pas$coefficients[2, 4], 3)
    res_pas[i, 4] <- round(summ_pas$coefficients[3, 1], 3)
    res_pas[i, 5] <- round(summ_pas$coefficients[3, 4], 3)
    res_pas[i, 6] <- round(summ_pas$coefficients[4, 1], 3)
    res_pas[i, 7] <- round(summ_pas$coefficients[4, 4], 3)
    res_pas[i, 8] <- round(summ_pas$coefficients[5, 1], 3)
    res_pas[i, 9] <- round(summ_pas$coefficients[5, 4], 3)
    res_pas[i, 10] <- round(summ_pas$coefficients[6, 1], 3)
    res_pas[i, 11] <- round(summ_pas$coefficients[6, 4], 3)
    res_pas[i, 12] <- round(summ_pas$r.squared, 3)
}
#write.csv(res_pas, "tables/Table4_unformatted_pas.csv")
save(pgls_pas, file = "data/results/pgls_pas.RData")

for (i in 1:1000) {
    print(i)
    pruned_pic <- treedata(tr[[i]], env_pic, warnings = F)$phy
    comp_pic <- comparative.data(pruned_pic, env_pic, Species)
    pgls_pic[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_pic)
    
    summ_pic <- summary(pgls_pic[[i]])
    res_pic[i, 1] <- round(summ_pic$coefficients[1, 1], 3)
    res_pic[i, 2] <- round(summ_pic$coefficients[2, 1], 3)
    res_pic[i, 3] <- round(summ_pic$coefficients[2, 4], 3)
    res_pic[i, 4] <- round(summ_pic$coefficients[3, 1], 3)
    res_pic[i, 5] <- round(summ_pic$coefficients[3, 4], 3)
    res_pic[i, 6] <- round(summ_pic$coefficients[4, 1], 3)
    res_pic[i, 7] <- round(summ_pic$coefficients[4, 4], 3)
    res_pic[i, 8] <- round(summ_pic$coefficients[5, 1], 3)
    res_pic[i, 9] <- round(summ_pic$coefficients[5, 4], 3)
    res_pic[i, 10] <- round(summ_pic$coefficients[6, 1], 3)
    res_pic[i, 11] <- round(summ_pic$coefficients[6, 4], 3)
    res_pic[i, 12] <- round(summ_pic$r.squared, 3)
}
#write.csv(res_pic, "tables/Table4_unformatted_pic.csv")
#save(pgls_pic, file = "data/results/pgls_pic.RData")

for (i in 1:1000) {
    print(i)
    pruned_psi <- treedata(tr[[i]], env_psi, warnings = F)$phy
    comp_psi <- comparative.data(pruned_psi, env_psi, Species)
    pgls_psi[[i]] <- pgls(SDI ~ Richness + abs(Latitude) + bio4 + bio15 + NPP,
                          data = comp_psi)
    
    summ_psi <- summary(pgls_psi[[i]])
    res_psi[i, 1] <- round(summ_psi$coefficients[1, 1], 3)
    res_psi[i, 2] <- round(summ_psi$coefficients[2, 1], 3)
    res_psi[i, 3] <- round(summ_psi$coefficients[2, 4], 3)
    res_psi[i, 4] <- round(summ_psi$coefficients[3, 1], 3)
    res_psi[i, 5] <- round(summ_psi$coefficients[3, 4], 3)
    res_psi[i, 6] <- round(summ_psi$coefficients[4, 1], 3)
    res_psi[i, 7] <- round(summ_psi$coefficients[4, 4], 3)
    res_psi[i, 8] <- round(summ_psi$coefficients[5, 1], 3)
    res_psi[i, 9] <- round(summ_psi$coefficients[5, 4], 3)
    res_psi[i, 10] <- round(summ_psi$coefficients[6, 1], 3)
    res_psi[i, 11] <- round(summ_psi$coefficients[6, 4], 3)
    res_psi[i, 12] <- round(summ_psi$r.squared, 3)
}
#write.csv(res_psi, "tables/Table4_unformatted_psi.csv")
#save(pgls_psi, file = "data/results/pgls_psi.RData")

#load(file = "data/results/pgls_acc.RData")
#load(file = "data/results/pgls_ans.RData")
#load(file = "data/results/pgls_apo.RData")
#load(file = "data/results/pgls_cha.RData")
#load(file = "data/results/pgls_col.RData")
#load(file = "data/results/pgls_gal.RData")
#load(file = "data/results/pgls_pas.RData")
#load(file = "data/results/pgls_pic.RData")
#load(file = "data/results/pgls_psi.RData")

stats <- function(x) {
    mean <- round(mean(x), 3)
    min <- round(range(x)[1], 3)
    max <- round(range(x)[2], 3)

    res <- paste0(mean, " (", min, "-", max, ")")
    res
}

res_pgls <- as.data.frame(matrix(ncol = 12, nrow = 9))
colnames(res_pgls) <- colnames(res_acc)
rownames(res_pgls) <- c("Accipitriformes", "Anseriformes", "Apodiformes",
                        "Charadriiformes", "Columbiformes", "Galliformes",
                        "Passeriformes", "Piciformes", "Psittaciformes")

res_pgls[1, 1:12] <- apply(res_acc, 2, stats)
res_pgls[2, 1:12] <- apply(res_ans, 2, stats)
res_pgls[3, 1:12] <- apply(res_apo, 2, stats)
res_pgls[4, 1:12] <- apply(res_cha, 2, stats)
res_pgls[5, 1:12] <- apply(res_col, 2, stats)
res_pgls[6, 1:12] <- apply(res_gal, 2, stats)
res_pgls[7, 1:12] <- apply(res_pas, 2, stats)
res_pgls[8, 1:12] <- apply(res_pic, 2, stats)
res_pgls[9, 1:12] <- apply(res_psi, 2, stats)

write.csv(res_pgls, "tables/Table4_unformatted.csv")

col_lin <- rgb(0.7, 0.7, 0.7, 0.1)
col_ric <- paste0(wes_palette("Zissou1")[1], 80)
col_lat <- paste0(wes_palette("Zissou1")[2], 80)
col_b04 <- paste0(wes_palette("Zissou1")[3], 80)
col_b15 <- paste0(wes_palette("Zissou1")[4], 80)
col_npp <- paste0(wes_palette("Zissou1")[5], 60)

# Figure S3

pdf("figures/FigureS3.pdf")

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(4, 4, 2, 2))

plot(SDI ~ Richness, data = env_acc, pch = 16, col = col_ric, xlab = "")
#abline(a = res_acc[1, 1], b = res_acc[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_acc[i, 1], b = res_acc[i, 2], col = col_lin)
title("Accipitriformes", adj = 0)

plot(SDI ~ Richness, data = env_ans, pch = 16, col = col_ric, xlab = "", 
     ylab = "")
#abline(a = res_ans[1, 1], b = res_ans[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_ans[i, 1], b = res_ans[i, 2], col = col_lin)
title("Anseriformes", adj = 0)

plot(SDI ~ Richness, data = env_apo, pch = 16, col = col_ric, xlab = "", 
     ylab = "")
#abline(a = res_apo[1, 1], b = res_apo[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_apo[i, 1], b = res_apo[i, 2], col = col_lin)
title("Apodiformes", adj = 0)

plot(SDI ~ Richness, data = env_cha, pch = 16, col = col_ric, xlab = "")
#abline(a = res_cha[1, 1], b = res_cha[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_cha[i, 1], b = res_cha[i, 2], col = col_lin)
title("Charadriiformes", adj = 0)

plot(SDI ~ Richness, data = env_col, pch = 16, col = col_ric, xlab = "", 
     ylab = "")
#abline(a = res_col[1, 1], b = res_col[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_col[i, 1], b = res_col[i, 2], col = col_lin)
title("Columbiformes", adj = 0)

plot(SDI ~ Richness, data = env_gal, pch = 16, col = col_ric, xlab = "", 
     ylab = "")
#abline(a = res_gal[1, 1], b = res_gal[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_gal[i, 1], b = res_gal[i, 2], col = col_lin)
title("Galliformes", adj = 0)

plot(SDI ~ Richness, data = env_pas, pch = 16, col = col_ric, 
     xlab = "Richness")
#abline(a = res_pas[1, 1], b = res_pas[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_pas[i, 1], b = res_pas[i, 2], col = col_lin)
title("Passeriformes", adj = 0)

plot(SDI ~ Richness, data = env_pic, pch = 16, col = col_ric, ylab = "",
     xlab = "Richness")
#abline(a = res_pic[1, 1], b = res_pic[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_pic[i, 1], b = res_pic[i, 2], col = col_lin)
title("Piciformes", adj = 0)

plot(SDI ~ Richness, data = env_psi, pch = 16, col = col_ric, ylab = "",
     xlab = "Richness")
#abline(a = res_psi[1, 1], b = res_psi[1, 2], col = col_lin)
#for (i in 2:1000) abline(a = res_psi[i, 1], b = res_psi[i, 2], col = col_lin)
title("Psittaciformes", adj = 0)

dev.off()

# Figure S4

pdf("figures/FigureS4.pdf")

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(4, 4, 2, 2))

plot(SDI ~ abs(Latitude), data = env_acc, pch = 16, col = col_lat, xlab = "")
title("Accipitriformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_ans, pch = 16, col = col_lat, xlab = "", 
     ylab = "")
title("Anseriformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_apo, pch = 16, col = col_lat, xlab = "", 
     ylab = "")
title("Apodiformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_cha, pch = 16, col = col_lat, xlab = "")
title("Charadriiformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_col, pch = 16, col = col_lat, xlab = "", 
     ylab = "")
title("Columbiformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_gal, pch = 16, col = col_lat, xlab = "", 
     ylab = "")
title("Galliformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_pas, pch = 16, col = col_lat,
     xlab = "Absolute latitude")
title("Passeriformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_pic, pch = 16, col = col_lat, ylab = "",
     xlab = "Absolute latitude")
title("Piciformes", adj = 0)

plot(SDI ~ abs(Latitude), data = env_psi, pch = 16, col = col_lat, ylab = "",
     xlab = "Absolute latitude")
title("Psittaciformes", adj = 0)

dev.off()

# Figure S5

pdf("figures/FigureS5.pdf")

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(4, 4, 2, 2))

plot(SDI ~ bio4, data = env_acc, pch = 16, col = col_b04, xlab = "")
title("Accipitriformes", adj = 0)

plot(SDI ~ bio4, data = env_ans, pch = 16, col = col_b04, xlab = "", ylab = "")
title("Anseriformes", adj = 0)

plot(SDI ~ bio4, data = env_apo, pch = 16, col = col_b04, xlab = "", ylab = "")
title("Apodiformes", adj = 0)

plot(SDI ~ bio4, data = env_cha, pch = 16, col = col_b04, xlab = "")
title("Charadriiformes", adj = 0)

plot(SDI ~ bio4, data = env_col, pch = 16, col = col_b04, xlab = "", ylab = "")
title("Columbiformes", adj = 0)

plot(SDI ~ bio4, data = env_gal, pch = 16, col = col_b04, xlab = "", ylab = "")
title("Galliformes", adj = 0)

plot(SDI ~ bio4, data = env_pas, pch = 16, col = col_b04,
     xlab = "Temperature Seasonality")
title("Passeriformes", adj = 0)

plot(SDI ~ bio4, data = env_pic, pch = 16, col = col_b04, ylab = "",
     xlab = "Temperature Seasonality")
title("Piciformes", adj = 0)

plot(SDI ~ bio4, data = env_psi, pch = 16, col = col_b04, ylab = "",
     xlab = "Temperature Seasonality")
title("Psittaciformes", adj = 0)

dev.off()

# Figure S6

pdf("figures/FigureS6.pdf")

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(4, 4, 2, 2))

plot(SDI ~ bio15, data = env_acc, pch = 16, col = col_b15, xlab = "")
title("Accipitriformes", adj = 0)

plot(SDI ~ bio15, data = env_ans, pch = 16, col = col_b15, xlab = "", 
     ylab = "")
title("Anseriformes", adj = 0)

plot(SDI ~ bio15, data = env_apo, pch = 16, col = col_b15, xlab = "", 
     ylab = "")
title("Apodiformes", adj = 0)

plot(SDI ~ bio15, data = env_cha, pch = 16, col = col_b15, xlab = "")
title("Charadriiformes", adj = 0)

plot(SDI ~ bio15, data = env_col, pch = 16, col = col_b15, xlab = "", 
     ylab = "")
title("Columbiformes", adj = 0)

plot(SDI ~ bio15, data = env_gal, pch = 16, col = col_b15, xlab = "", 
     ylab = "")
title("Galliformes", adj = 0)

plot(SDI ~ bio15, data = env_pas, pch = 16, col = col_b15,
     xlab = "Precipitation Seasonality")
title("Passeriformes", adj = 0)

plot(SDI ~ bio15, data = env_pic, pch = 16, col = col_b15, ylab = "",
     xlab = "Precipitation Seasonality")
title("Piciformes", adj = 0)

plot(SDI ~ bio15, data = env_psi, pch = 16, col = col_b15, ylab = "",
     xlab = "Precipitation Seasonality")
title("Psittaciformes", adj = 0)

dev.off()

# Figure S7

pdf("figures/FigureS7.pdf")

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(4, 4, 2, 2))

plot(SDI ~ NPP, data = env_acc, pch = 16, col = col_npp, xlab = "")
title("Accipitriformes", adj = 0)

plot(SDI ~ NPP, data = env_ans, pch = 16, col = col_npp, xlab = "", ylab = "")
title("Anseriformes", adj = 0)

plot(SDI ~ NPP, data = env_apo, pch = 16, col = col_npp, xlab = "", ylab = "")
title("Apodiformes", adj = 0)

plot(SDI ~ NPP, data = env_cha, pch = 16, col = col_npp, xlab = "")
title("Charadriiformes", adj = 0)

plot(SDI ~ NPP, data = env_col, pch = 16, col = col_npp, xlab = "", ylab = "")
title("Columbiformes", adj = 0)

plot(SDI ~ NPP, data = env_gal, pch = 16, col = col_npp, xlab = "", ylab = "")
title("Galliformes", adj = 0)

plot(SDI ~ NPP, data = env_pas, pch = 16, col = col_npp, 
     xlab = "Net Primary Productivity")
title("Passeriformes", adj = 0)

plot(SDI ~ NPP, data = env_pic, pch = 16, col = col_npp, ylab = "", 
     xlab = "Net Primary Productivity")
title("Piciformes", adj = 0)

plot(SDI ~ NPP, data = env_psi, pch = 16, col = col_npp, ylab = "", 
     xlab = "Net Primary Productivity")
title("Psittaciformes", adj = 0)

dev.off()
