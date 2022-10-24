rm(list = ls())

setwd("~/dimorphism_evol")

library(sf)
#library(epm)
library(stringi)
#library(rgdal)
library(raster)
library(dplyr)
#library(gdalUtilities)
library(fasterize)
library(maptools)
library(MetBrewer)
library(ggplot2)

file <- "BOTW/BOTW.gdb"

birds <- st_read(dsn = file, layer = "All_Species")

birds$sci_name <- stri_replace_all_fixed(birds$sci_name, " ", "_")

dat <- read.csv("BodySizeAves_30may22_edit.csv", row.names = 1)

SDI <- function (male = male, female = female) {
  if (male <= female) {
    SDI <- (female/male) - 1
  } else {
    if (female < male) {
      SDI <- -((male/female) - 1)
    } 
  }
  
  return(SDI)
  
}

dat_red <- dat[, c(5, 9, 13)]
dat_red <- dat_red[complete.cases(dat_red$Body_mass_g_M_mean) & 
                     complete.cases(dat_red$Body_mass_g_F_mean), ]

sdi <- numeric()
for (i in 1:nrow(dat_red)) {
  sdi[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
                female = dat_red$Body_mass_g_F_mean[i])
}
names(sdi) <- dat_red$Scientific_name

r <- raster(ncols = 2160, nrows = 900, ymn = -60)
raster_stack_fem <- r
raster_stack_mal <- r

birds2 <- birds %>% filter(st_geometry_type(Shape) != "MULTISURFACE")

sdi <- sdi[names(sdi) %in% birds2$sci_name]

sdi_fem <- sdi[sdi > 0]
sdi_mal <- sdi[sdi < 0]

for (i in 1:length(sdi_fem)) {
  print(i) 
  
  s <- as.character(names(sdi_fem)[i]) 
  map_i <- subset(birds2, birds2$sci_name == s)
  
  for (j in 1:length(map_i$Shape)) { 
    try(map_i$Shape[[j]] <- st_cast(map_i$Shape[[j]], 'MULTIPOLYGON'))
  }
  try(map_i$Shape <- st_cast(map_i$Shape, 'MULTIPOLYGON'))
  
  raster_i <- fasterize(st_as_sf(map_i$Shape), r)
  
  rastercells <- which(getValues(!is.na(raster_i))) 
  raster_i[rastercells] <- sdi_fem[i]
  
  if(i == 1) {
    raster_stack_fem <- raster_i
  } else {
    raster_stack_fem <- addLayer(raster_stack_fem, raster_i)
  }
}

av_full_fem <- stackApply(raster_stack_fem, indices = rep(1, length(sdi_fem)), 
                          fun = median, na.rm = T) 

for (i in 1:length(sdi_mal)) {
  print(i) 
  
  s <- as.character(names(sdi_mal)[i]) 
  map_i <- subset(birds2, birds2$sci_name == s)
  
  for (j in 1:length(map_i$Shape)) { 
    try(map_i$Shape[[j]] <- st_cast(map_i$Shape[[j]], 'MULTIPOLYGON'))
  }
  try(map_i$Shape <- st_cast(map_i$Shape, 'MULTIPOLYGON'))
  
  raster_i <- fasterize(st_as_sf(map_i$Shape), r)
  
  rastercells <- which(getValues(!is.na(raster_i))) 
  raster_i[rastercells] <- sdi_mal[i]
  
  if(i == 1) {
    raster_stack_mal <- raster_i
  } else {
    raster_stack_mal <- addLayer(raster_stack_mal, raster_i)
  }
}

av_full_mal <- stackApply(raster_stack_mal, indices = rep(1, length(sdi_mal)), 
                          fun = median, na.rm = T) 

data(wrld_simpl)
rem_fem1 <- extract(av_full_fem, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_fem2 <- do.call(rbind.data.frame, rem_fem1)[, 1]
values(av_full_fem)[-rem_fem2] <- NA

rem_mal1 <- extract(av_full_mal, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_mal2 <- do.call(rbind.data.frame, rem_mal1)[, 1]
values(av_full_mal)[-rem_mal2] <- NA

brks <- quantile(values(av_full_fem)[order(values(av_full_fem))], probs = seq(0, 1, 0.02), na.rm = T)
av_full_fem_p <- rasterToPoints(av_full_fem, spatial = TRUE)
av_full_fem_df  <- data.frame(av_full_fem_p)
av_full_fem_df <- av_full_fem_df %>% mutate(index_2 = cut(index_1, breaks = brks))

brks <- quantile(values(av_full_mal)[order(values(av_full_mal))], probs = seq(0, 1, 0.02), na.rm = T)
av_full_mal_p <- rasterToPoints(av_full_mal, spatial = TRUE)
av_full_mal_df  <- data.frame(av_full_mal_p)
av_full_mal_df <- av_full_mal_df %>% mutate(index_2 = cut(index_1, breaks = brks))

colors_fem <- c("#1E466E", "#224C75", "#27527C", "#2B5883", "#305E8A", "#346491",
                "#396B97", "#3E729B", "#4379A0", "#4881A4", "#4D88A9", "#528FAD", 
                "#5898B5", "#5EA0BC", "#64A8C3", "#6AB0CB", "#70B9D2", "#78BFD6", 
                "#83C5D8", "#8DCBDA", "#97D1DC", "#A2D7DE", "#ADDCDE", "#BDDED6", 
                "#CCE0CF", "#DCE1C7", "#EBE3C0", "#FBE5B8", "#FFE2AC", "#FFDE9F", 
                "#FFDA92", "#FFD685", "#FFD277", "#FECD6D", "#FDC669", "#FBBF65", 
                "#FAB860", "#F8B15C", "#F7AA58", "#F5A455", "#F49E52", "#F2994E", 
                "#F1934B", "#EF8D48", "#EE8648", "#EC7F4A", "#EB784C", "#E9704F", 
                "#E86951", "#E76254")
colors_mal <- c("#E76254", "#E86951", "#E9704F", "#EB784C", "#EC7F4A", "#EE8648", 
                "#EF8D48", "#F1934B", "#F2994E", "#F49E52", "#F5A455", "#F7AA58", 
                "#F8B15C", "#FAB860", "#FBBF65", "#FDC669", "#FECD6D", "#FFD277", 
                "#FFD685", "#FFDA92", "#FFDE9F", "#FFE2AC", "#FBE5B8", "#EBE3C0", 
                "#DCE1C7", "#CCE0CF", "#BDDED6", "#ADDCDE", "#A2D7DE", "#97D1DC", 
                "#8DCBDA", "#83C5D8", "#78BFD6", "#70B9D2", "#6AB0CB", "#64A8C3", 
                "#5EA0BC", "#5898B5", "#528FAD", "#4D88A9", "#4881A4", "#4379A0", 
                "#3E729B", "#396B97", "#346491", "#305E8A", "#2B5883", "#27527C", 
                "#224C75", "#1E466E")

ggplot_fem <- ggplot() +
  geom_raster(data = av_full_fem_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_fem) +
  theme_void() 

inches <- 4.5
res <- 600
plot_fem <- "Figure3_female.tiff"
tiff(plot_fem, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_fem)  
dev.off()  

ggplot_mal <- ggplot() +
  geom_raster(data = av_full_mal_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mal) +
  theme_void() 

inches <- 4.5
res <- 600
plot_mal <- "Figure3_male.tiff"
tiff(plot_mal, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_mal)  
dev.off() 

plot_sca_fem <- "Figure3_ScaleBar.tiff"
tiff(plot_sca_fem, width = 1*res, h = 0.1*res, units = "px")
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))
brks <- seq(0, 1, 0.02)
breaks <- seq(0, 100, length.out = length(brks))

ix <- 1:2
iy <- breaks
nBreaks <- length(breaks)
midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",ylab = "", 
      col = colors_fem, breaks = breaks)
axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                  cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                  labels = rep("", 5), tck = 0.2)
do.call("axis", axis.args)
box()
dev.off()

as.numeric(quantile(values(av_full_fem)[order(values(av_full_fem))], probs = seq(0, 1, 0.01),
                    na.rm = T)[seq(1, 101, 25)])

as.numeric(quantile(values(av_full_mal)[order(values(av_full_mal))], probs = seq(0, 1, 0.01),
                    na.rm = T)[seq(1, 101, 25)])

# Environmental data

envar <- getData("worldclim", var = "bio", res = 2.5) 

bio_1 <- envar[[1]]
bio_2 <- envar[[2]]
bio_4 <- envar[[4]]
bio_5 <- envar[[5]]
bio_6 <- envar[[6]]
bio_7 <- envar[[7]]
bio_12 <- envar[[12]]
bio_13 <- envar[[13]]
bio_14 <- envar[[14]]
bio_15 <- envar[[15]]

envars <- stack(bio_1, bio_2, bio_4, bio_5, bio_6, bio_7, bio_12, bio_13, bio_14, bio_15)

r <- raster(ncol = 8640, nrow = 3600, ymn = -60)
r[] <- NA 
template_raster <- r
r[] <- 1:length(values(r))

var_vals <- raster::extract(envars, r[])

species <- names(sdi)
  
maps_cut <- subset(birds, birds$presence %in% c(1) & birds$origin %in% c(1,2) & birds$seasonal %in% c(1,2) & birds$sci_name %in% species)

for(i in 1:length(unique(maps_cut$sci_name))) try({
  print(i)
  
  species_i <- unique(maps_cut$sci_name)[i]
  
  maps_cut_i <- maps_cut %>% filter(sci_name == species_i)
  
  ##Convert the maps into multipolygon format
  for (j in 1:length(maps_cut_i$Shape)) { 
    try(maps_cut_i$Shape[j] <- st_make_valid(maps_cut_i$Shape[j]))
    try(maps_cut_i$Shape[[j]] <- st_cast(maps_cut_i$Shape[[j]], 'MULTIPOLYGON'))  #changes the shapes in the feature
  }
  try(maps_cut_i$Shape <- st_cast(maps_cut_i$Shape, 'MULTIPOLYGON'))
  
  #Rasterize the polygons
  raster_i <- fasterize(st_as_sf(maps_cut_i$Shape), r) #fasterize polygon onto raster
  
  #Select the projected cells of the raster and take the values at these cells by comparison to the var_vals array
  rastercells <- which(getValues(raster_i) > 0) #select the cell
  
  ## you could extract it here but it was faster to preextract
  ## now we just take the values for the cells where our raster is present
  values <- var_vals[rastercells, ]
  
  #some only have one rastercell and therefore only one value - take this value
  #cannot divide by 10 if value is NA therefore ifelse 
  # Divide by 10 as temperature rasters are stored as *10 in worldclim.org
  if(nrow(values) == 1) {
    ifelse(is.na(values[1]), 
           bio1 <- NA, 
           bio1 <- values[1]/10)
    ifelse(is.na(values[2]), 
           bio2 <- NA, 
           bio2 <-values[2]/10)
    ifelse(is.na(values[3]), 
           bio4 <- NA, 
           bio4 <-values[3])
    ifelse(is.na(values[4]), 
           bio5 <- NA, 
           bio5 <-values[4]/10)
    ifelse(is.na(values[5]), 
           bio6 <- NA, 
           bio6 <-values[5]/10)
    ifelse(is.na(values[6]), 
           bio7 <- NA, 
           bio7 <-values[6]/10)
    ifelse(is.na(values[7]), 
           bio12 <- NA, 
           bio12 <-values[7])
    ifelse(is.na(values[8]), 
           bio13 <- NA, 
           bio13 <-values[8])
    ifelse(is.na(values[9]), 
           bio14 <- NA, 
           bio14 <-values[9])
    ifelse(is.na(values[10]), 
           bio15 <- NA, 
           bio15 <-values[10])
    
    
  } else { 
    #most have multiple so take the mean of all the values
    #However cannot take the mean or divide by 10 if all vals are NA hense ifelse()
    
    ifelse(all(is.na(values[,1])), 
           bio1 <- NA,
           bio1 <- mean(values[,1], na.rm = TRUE)/10) #write them in
    ifelse(all(is.na(values[,2])),
           bio2 <- NA,
           bio2 <- mean(values[,2], na.rm = TRUE)/10) #write them in
    ifelse(all(is.na(values[,3])), 
           bio4 <- NA,
           bio4 <- mean(values[,3], na.rm = TRUE)) #write them in
    ifelse(all(is.na(values[,4])),
           bio5 <- NA,
           bio5 <- mean(values[,4], na.rm = TRUE)/10) #write them in
    ifelse(all(is.na(values[,5])), 
           bio6 <- NA,
           bio6 <- mean(values[,5], na.rm = TRUE)/10) #write them in
    ifelse(all(is.na(values[,6])), 
           bio7 <- NA,
           bio7 <-mean(values[,6], na.rm = TRUE)/10) #write them in
    ifelse(all(is.na(values[,7])), 
           bio12 <- NA,
           bio12 <- mean(values[,7], na.rm = TRUE)) #write them in
    ifelse(all(is.na(values[,8])), 
           bio13 <- NA,
           bio13 <- mean(values[,8], na.rm = TRUE)) #write them in
    ifelse(all(is.na(values[,9])), 
           bio14 <- NA,
           bio14 <- mean(values[,9], na.rm = TRUE)) #write them in
    ifelse(all(is.na(values[,10])), 
           bio15 <- NA,
           bio15 <- mean(values[,10], na.rm = TRUE)) #write them in
  }
  
  #get the study_results
  study_results <- data.frame(species = maps_cut_i$SCINAME[1],
                              bio1,
                              bio2,
                              bio4,
                              bio5,
                              bio6,
                              bio7,
                              bio12,
                              bio13,
                              bio14,
                              bio15)
  
  ## Notice I am appending the csv as the loop continues. This is just incase it crashes and we need to restart but keep same amount of species
  write.table(study_results, "BL_environmental.csv", sep = ",", col.names = FALSE, append = TRUE)
  
})

