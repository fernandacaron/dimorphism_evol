rm(list = ls())

setwd("~/Documents/lab/dimorph_evol")

library(phytools)
library(geiger)
library(MetBrewer)
library(colorspace)
library(stringi)
library(sf)
library(raster)
library(dplyr)
library(fasterize)
library(maptools)
library(ggplot2)

dat <- read.csv("data/BodySizeAves_30may22_edit.csv", row.names = 1)
tr <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")

male <- "#9966FF"
monom <- "gray"
female <- "#E69F00"
male_al <- rgb(153/255, 102/255, 255/255, alpha = 0.5)
female_al <- rgb(230/255, 159/255, 0/255, alpha = 0.5)

dat_red <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
                 complete.cases(dat$Body_mass_g_F_mean), ]

Body_mass_g_mean <- numeric()
for (i in 1:nrow(dat_red)) {
  Body_mass_g_mean[i] <- mean(c(dat_red$Body_mass_g_M_mean[i], 
                                dat_red$Body_mass_g_F_mean[i]))
}
names(Body_mass_g_mean) <- dat_red$Scientific_name 

dat_red <- cbind(dat_red, Body_mass_g_mean)

## Figure 1

pdf("figures/Figure1_org.pdf", height = 10, width = 13)
layout(matrix(1:9, ncol = 3, byrow = TRUE))

par(mar = c(4, 4, 2, 2))

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Accipitriformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 20, col = male_al, 
     border = male)
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Accipitriformes"]),
     breaks = 20, col = female_al, border = female, add = T)
title("Accipitriformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Anseriformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 20, col = male_al, 
     border = male, ylab = "")
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Anseriformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Anseriformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Apodiformes"]),
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 20, col = male_al,
     border = male, ylab = "", ylim = c(0, 60))
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Apodiformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Apodiformes", adj = 0)
legend("topright", pch = 15, bty = 'n', col = c(male, female), 
       legend  = c("Male", "Female"))

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Charadriiformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 20, col = male_al,
     border = male)
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Charadriiformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Charadriiformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Columbiformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 10, col = male_al,
     border = male, ylab = "")
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Columbiformes"]), 
     breaks = 10, col = female_al, border = female, add = T)
title("Columbiformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Galliformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "", breaks = 20, col = male_al,
     border = male, ylim = c(0, 25), ylab = "")
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Galliformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Galliformes", adj = 0)

par(mar = c(5, 4, 2, 2))

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Passeriformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "log body size (g)", breaks = 20, 
     col = male_al, border = male)
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Passeriformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Passeriformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Piciformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "log body size (g)", breaks = 20, 
     col = male_al, border = male, ylab = "")
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Piciformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Piciformes", adj = 0)

hist(log10(dat_red$Body_mass_g_M_mean[dat_red$Order == "Psittaciformes"]), 
     xlim = c(0, 4.5), main = "", xlab = "log body size (g)", breaks = 20, 
     col = male_al, border = male, ylab = "")
hist(log10(dat_red$Body_mass_g_F_mean[dat_red$Order == "Psittaciformes"]), 
     breaks = 20, col = female_al, border = female, add = T)
title("Psittaciformes", adj = 0)

dev.off()

## Figure 2

## SDI = (female size/male size - 1) when female larger
## SDI = -(male size/female size - 1) when male larger

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

tr_map <- treedata(tr[[1]], sdi)$phy

sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- 1
sdi_disc[sdi_disc == 0] <- 0
sdi_disc[sdi_disc < 0] <- -1

map <- make.simmap(tr_map, sdi_disc)

cols_pal <- c(male, monom, female)
set_cols <- setNames(c(male, monom, female), c(-1, 0, 1))

col1_lig <- lighten(cols_pal[1], amount = 0.4)
col2_lig <- lighten(cols_pal[2], amount = 0.4)
col3_lig <- lighten(cols_pal[3], amount = 0.4)

col_body <- character()
for (i in 1:length(Body_mass_g_mean)) {
	col_body[i] <- 
		ifelse(sdi_disc[names(sdi_disc) == names(Body_mass_g_mean)[i]] == 1, 
		       col3_lig, 
		       ifelse(sdi_disc[names(sdi_disc) == names(Body_mass_g_mean)[i]] == -1,
		              col1_lig, col2_lig))
}
names(col_body) <- names(Body_mass_g_mean)
col_body <- col_body[map$tip.label]

rbPal <- colorRampPalette(cols_pal)
cols <- rbPal(101)[as.numeric(cut(1:101, breaks = 101))]

pdf("figures/Figure2_org.pdf", width = 9)

par(mar = c(0, 0, 1, 0))

plotTree.wBars(map, log10(Body_mass_g_mean), scale = 4, tip.labels = F,
    type = "fan", method = "plotSimmap", colors = set_cols, lwd = 1, 
    border = NA, col = col_body, mar = c(0, 0, 1, 0), part = 0.5)

legend("topright", col = cols_pal, bty = "n", pch = 15,
       legend = c("Male-biased SSD", "Monomorphism", "Female-biased SSD"))

par(xpd = TRUE)

arc.cladelabels(text = "Accipitriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Accipitriformes"][dat$Scientific_name[dat$Order == 
                                  "Accipitriformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Anseriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Anseriformes"][dat$Scientific_name[dat$Order == 
                                  "Anseriformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Apodiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Apodiformes"][dat$Scientific_name[dat$Order == 
                                  "Apodiformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Charadriiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"][dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Columbiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Columbiformes"][dat$Scientific_name[dat$Order == 
                                  "Columbiformes"] %in% names(sdi_disc)])))

arc.cladelabels(text = "Galliformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Galliformes"][dat$Scientific_name[dat$Order == 
                                  "Galliformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Passeriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Passeriformes"][dat$Scientific_name[dat$Order == 
                                  "Passeriformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Piciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Piciformes"][dat$Scientific_name[dat$Order == 
                                  "Piciformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Psittaciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"][dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"] %in% names(sdi_disc)])))

dev.off()

### Scale bar

pdf("figures/Figure2_scalebar.pdf", width = 9)

par(mar = c(0, 0, 1, 0))

plotTree.barplot(map, log10(Body_mass_g_mean), scale = 4, tip.labels = F,
                 lwd = 1, border = NA, mar = c(0, 0, 1, 0), part = 0.5, 
                 args.axis = list(cex.axis = 0.8, at = c(0, 2.5, 5)))

dev.off()

# Figure 3

## This part of the code was provided by Thomas Weeks

birds <- st_read(dsn = "~/Documents/lab/data/spatial/BOTW/BOTW.gdb", layer = "All_Species")

birds$sci_name <- stri_replace_all_fixed(birds$sci_name, " ", "_")

dat <- read.csv("data/BodySizeAves_30may22_edit.csv", row.names = 1)

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
raster_stack_all <- r
raster_stack_fem <- r
raster_stack_mal <- r

birds2 <- birds %>% filter(st_geometry_type(Shape) != "MULTISURFACE")

sdi <- sdi[names(sdi) %in% birds2$sci_name]

sdi_all <- abs(sdi)
sdi_fem <- sdi[sdi > 0]
sdi_mal <- sdi[sdi < 0]

for (i in 1:length(sdi_all)) {
  print(i) 
  
  s <- as.character(names(sdi_all)[i]) 
  map_i <- subset(birds2, birds2$sci_name == s)
  
  for (j in 1:length(map_i$Shape)) { 
    try(map_i$Shape[[j]] <- st_cast(map_i$Shape[[j]], 'MULTIPOLYGON'))
  }
  try(map_i$Shape <- st_cast(map_i$Shape, 'MULTIPOLYGON'))
  
  raster_i <- fasterize(st_as_sf(map_i$Shape), r)
  
  rastercells <- which(getValues(!is.na(raster_i))) 
  raster_i[rastercells] <- sdi_all[i]
  
  if(i == 1) {
    raster_stack_all <- raster_i
  } else {
    raster_stack_all <- addLayer(raster_stack_all, raster_i)
  }
  removeTmpFiles(h = 0)
}

av_full_all <- stackApply(raster_stack_all, indices = rep(1, length(sdi_all)), 
                          fun = median, na.rm = T)

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

rem_all1 <- extract(av_full_all, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_all2 <- do.call(rbind.data.frame, rem_all1)[, 1]
values(av_full_all)[-rem_all2] <- NA

rem_fem1 <- extract(av_full_fem, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_fem2 <- do.call(rbind.data.frame, rem_fem1)[, 1]
values(av_full_fem)[-rem_fem2] <- NA

rem_mal1 <- extract(av_full_mal, wrld_simpl, cellnumbers = T, weights = T, 
                    small = T)
rem_mal2 <- do.call(rbind.data.frame, rem_mal1)[, 1]
values(av_full_mal)[-rem_mal2] <- NA

brks_all <- quantile(values(av_full_all)[order(values(av_full_all))], 
                     probs = seq(0, 1, 0.02), na.rm = T)
av_full_all_p <- rasterToPoints(av_full_all, spatial = TRUE)
av_full_all_df  <- data.frame(av_full_all_p)
av_full_all_df <- av_full_all_df %>% mutate(index_2 = cut(index_1, 
                                                          breaks = brks_all))

brks_fem <- quantile(values(av_full_fem)[order(values(av_full_fem))], 
                 probs = seq(0, 1, 0.02), na.rm = T)
av_full_fem_p <- rasterToPoints(av_full_fem, spatial = TRUE)
av_full_fem_df  <- data.frame(av_full_fem_p)
av_full_fem_df <- av_full_fem_df %>% mutate(index_2 = cut(index_1, 
                                                          breaks = brks_fem))

brks_mal <- quantile(values(av_full_mal)[order(values(av_full_mal))], 
                 probs = seq(0, 1, 0.02), na.rm = T)
av_full_mal_p <- rasterToPoints(av_full_mal, spatial = TRUE)
av_full_mal_df  <- data.frame(av_full_mal_p)
av_full_mal_df <- av_full_mal_df %>% mutate(index_2 = cut(index_1, 
                                                          breaks = brks_mal))

colors_all <- met.brewer(name = "Hiroshige", n = 50, direction = -1)[1:50]
colors_fem <- met.brewer(name = "Hiroshige", n = 50, direction = -1)[1:50]
colors_mal <- met.brewer(name = "Hiroshige", n = 50, direction = 1)[1:50]

inches <- 4.5
res <- 600

ggplot_all <- ggplot() +
  geom_raster(data = av_full_all_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_all) +
  theme_void() 

plot_all <- "Figure3_all.tiff"
tiff(plot_all, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_all)  
dev.off()  

ggplot_fem <- ggplot() +
  geom_raster(data = av_full_fem_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_fem) +
  theme_void() 

plot_fem <- "Figure3_female.tiff"
tiff(plot_fem, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_fem)  
dev.off()  

ggplot_mal <- ggplot() +
  geom_raster(data = av_full_mal_df, aes(x = x, y = y, fill = index_2), show.legend = FALSE) +
  scale_fill_manual(values = colors_mal) +
  theme_void() 

plot_mal <- "Figure3_male.tiff"
tiff(plot_mal, width = inches*res, height = inches*res/2, units = "px")
print(ggplot_mal)  
dev.off() 

plot_sca<- "Figure3_ScaleBar.tiff"
tiff(plot_sca, width = 1*res, h = 0.1*res, units = "px")
par(mfrow = c(1, 1))
par(mar = c(1, 1, 1, 1))
brks <- seq(0, 1, 0.02)
breaks <- seq(0, 100, length.out = length(brks))

ix <- 1:2
iy <- breaks
nBreaks <- length(breaks)
midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
      col = colors_fem, breaks = breaks)
axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                  cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                  labels = rep("", 5), tck = 0.2)
do.call("axis", axis.args)
box()
dev.off()

round(as.numeric(quantile(values(av_full_all)[order(values(av_full_all))], 
                    probs = seq(0, 1, 0.01), na.rm = T)[seq(1, 101, 25)]), 3)

round(as.numeric(quantile(values(av_full_fem)[order(values(av_full_fem))], 
                    probs = seq(0, 1, 0.01), na.rm = T)[seq(1, 101, 25)]), 3)

round(as.numeric(quantile(values(av_full_mal)[order(values(av_full_mal))], 
                    probs = seq(0, 1, 0.01), na.rm = T)[seq(1, 101, 25)]), 3)
