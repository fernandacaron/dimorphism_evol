rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(MetBrewer)
library(colorspace)
library(stringi)
library(rgeos)
library(sf)
library(ggplot2)

dat <- read.csv("data/aves/BodySizeAves_30may22_edit.csv", row.names = 1)
tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")

# definir cores
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

#########################

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

#########################

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

#########################

# Figure S1 - cutoff

sdi_cut <- numeric()
for (i in 1:nrow(dat_red)) {
  sdi_cut[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
                    female = dat_red$Body_mass_g_F_mean[i],
                    cutoff = TRUE, cut.value = 0.10)
}
names(sdi_cut) <- dat_red$Scientific_name
sdi_cut <- sdi_cut[complete.cases(sdi_cut)]

tr_map_cut <- treedata(tr[[1]], sdi_cut)$phy

sdi_cut_disc <- sdi_cut
sdi_cut_disc[sdi_cut_disc > 0] <- 1
sdi_cut_disc[sdi_cut_disc == 0] <- 0
sdi_cut_disc[sdi_cut_disc < 0] <- -1

map_cut <- make.simmap(tr_map_cut, sdi_cut_disc)

col_body_cut <- character()
for (i in 1:length(Body_mass_g_mean)) {
  col_body_cut[i] <- 
    ifelse(sdi_cut_disc[names(sdi_cut_disc) == names(Body_mass_g_mean)[i]] == 1,
           col3_lig, 
           ifelse(sdi_cut_disc[names(sdi_cut_disc) == names(Body_mass_g_mean)[i]] == -1,
                  col1_lig, col2_lig))
}
names(col_body_cut) <- names(Body_mass_g_mean)
col_body_cut <- col_body_cut[map_cut$tip.label]

rbPal <- colorRampPalette(cols_pal)
cols <- rbPal(101)[as.numeric(cut(1:101, breaks = 101))]

pdf("figures/FigureS1_org.pdf", width = 9)

par(mar = c(0, 0, 1, 0))

plotTree.wBars(map_cut, log10(Body_mass_g_mean), scale = 4, tip.labels = F,
    type = "fan", method = "plotSimmap", colors = set_cols, lwd = 1, 
    border = NA, col = col_body_cut, mar = c(0, 0, 1, 0), part = 0.5)

legend("topright", col = cols_pal, bty = "n", pch = 15,
       legend = c("Male-biased SSD", "Monomorphism", "Female-biased SSD"))

par(xpd = TRUE)

arc.cladelabels(text = "Accipitriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Accipitriformes"][dat$Scientific_name[dat$Order == 
                                  "Accipitriformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Anseriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Anseriformes"][dat$Scientific_name[dat$Order == 
                                  "Anseriformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Apodiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Apodiformes"][dat$Scientific_name[dat$Order == 
                                  "Apodiformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Charadriiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"][dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Columbiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Columbiformes"][dat$Scientific_name[dat$Order == 
                                  "Columbiformes"] %in% names(sdi_cut_disc)])))

arc.cladelabels(text = "Galliformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Galliformes"][dat$Scientific_name[dat$Order == 
                                  "Galliformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Passeriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Passeriformes"][dat$Scientific_name[dat$Order == 
                                  "Passeriformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Piciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Piciformes"][dat$Scientific_name[dat$Order == 
                                  "Piciformes"] %in% names(sdi_cut_disc)])))
arc.cladelabels(text = "Psittaciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(map_cut, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"][dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"] %in% names(sdi_cut_disc)])))

dev.off()

#########################

# Figure 3

map.SSD <- function(data, func = mean, cols = NULL, figFolder, fileName) { 
  # de Tobias et al. (2022)
  
  # Behrmann equal area (96 x 96km) grid shapefile
  grid <- rgdal::readOGR("data/aves/spatial/BehrmannMeterGrid_WGS84_land.shp")
  
  # Country borders shapefile
  countriesGeo <- rgdal::readOGR("data/aves/spatial/all_countries.shp")
  
  # gridded species geographic range data - Birdlife taxonomy 
  rangeData <- read.csv("data/aves/spatial/AllSpeciesBirdLifeMaps2019.csv")

  rangeData$Species <- stri_replace_all_fixed(rangeData$Species, " ", "_")

  rangeData <- rangeData[rangeData$Species %in% names(data), ]
  
  # data processing and cleaning

  # set grid and country shapefile to Behrmann projection
  P4S.Behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
  gridB <- spTransform(grid, P4S.Behr)
  countries <- spTransform(countriesGeo, P4S.Behr)
  # simplify country shapefile - needed for plotting   
  countriesS <- gSimplify(countries, tol = 10000, topologyPreserve = TRUE)
  # convert to simple feature for plotting with ggplot
  countriesS2 <- st_as_sf(countriesS)

  # Maps

  # assign trait data to each species in range database
  rangeData$SDI <- data[match(rangeData$Species, names(data))]

  # remove rows (i.e. cells x species) with no trait data
  rangeData <- na.omit(rangeData)
  length(unique(na.omit(rangeData$Species)))

  # calculate median trait value per cell
  SDIperCell <- split(rangeData$SDI, rangeData$WorldID)
  SDIperCell <- lapply(SDIperCell, function(x) x[!is.na(x)])
  SDI_median_perCell <- sapply(SDIperCell, func)

  # assign values to grid shapefile
  gridB@data$SDI_median_perCell <- NA
  gridB@data$SDI_median_perCell[match(names(SDI_median_perCell), 
                                      gridB@data$WorldID)] <- 
                                                as.numeric(SDI_median_perCell)

  # set scale and colors
  brks <- quantile(gridB@data$SDI_median_perCell, probs = seq(0, 1, 0.02),
                   na.rm = T)
  gridB@data$col_SDI_median <- NA
  gridB@data$col_SDI_median <- findInterval(gridB@data$SDI_median_perCell, brks,
                                            all.inside = TRUE)

  if (is.null(cols)) {
    colors <- c(brewer.pal(9, "Blues")[2:4], brewer.pal(9, "YlGnBu")[5:9])
  } else {
    colors <- cols
  }
  colors <- colorRampPalette(colors)(50)

  # plot maps
  gridB2 <- st_as_sf(gridB)

  inches <- 4.5
  res <- 600

  ggplot.ssd <- ggplot(gridB2) +
                 geom_sf(aes(fill = col_SDI_median, color = col_SDI_median)) +
                 scale_colour_gradientn(colors = colors) +
                 scale_fill_gradientn(colors = colors) +
                 theme_void() 
  
  plot.ssd1 <- paste0(figFolder, fileName, ".tiff")
  tiff(plot.ssd1, width = inches*res, height = inches*res/2, units = "px")
  print(ggplot.ssd)  
  dev.off()  

  plot.ssd2 <- paste0(figFolder, fileName, "_ScaleBar", ".tiff")
  tiff(plot.ssd2, width = 1*res, h = 0.1*res, units = "px")
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
        col = colors, breaks = breaks)
  axis.args <- list(side = 1, padj = -1, mgp = c(3, 1, 0), las = 0, 
                    cex.axis = 0.5, mgp = c(1, 0, 0), at = seq(0, 100, 25),
                    labels = rep("", 5), tck = 0.2)
  do.call("axis", axis.args)
  box()
  dev.off()
  
  # tick values
  as.numeric(quantile(gridB@data$SDI_median_perCell, probs = seq(0, 1, 0.01),
                      na.rm = T)[seq(1, 101, 25)])

}

sdi_fem <- sdi[sdi > 0]
sdi_mal <- sdi[sdi < 0]

map.SSD(data = sdi_mal, func = median, figFolder = "figures/", 
        fileName = "Figure3_male", 
        cols = met.brewer(name = "Hiroshige", n = 100, direction = 1))

map.SSD(data = sdi_fem, func = median, figFolder = "figures/", 
        fileName = "Figure3_female", 
        cols = met.brewer(name = "Hiroshige", n = 100, direction = -1))


#########################

# Figure S2

sdi_cut_fem <- sdi_cut[sdi_cut > 0]
sdi_cut_mal <- sdi_cut[sdi_cut < 0]

map.SSD(data = sdi_cut_mal, func = median, figFolder = "figures/", 
        fileName = "FigureS2_male", 
        cols = met.brewer(name = "Hiroshige", n = 100, direction = 1))

map.SSD(data = sdi_cut_fem, func = median, figFolder = "figures/", 
        fileName = "FigureS2_female", 
        cols = met.brewer(name = "Hiroshige", n = 100, direction = -1))


