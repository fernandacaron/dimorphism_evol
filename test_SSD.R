rm(list = ls())

setwd("Documents/lab/dimorph_evol")

library(phytools)
library(stringi)
library(rgeos)
library(sf)
library(RColorBrewer)
library(ggplot2)
library(viridis)

####################################
# Histograms difference across sources
####################################
# Ocampo et al. 2021
ref1 <- read.csv("data/Ocampo&al2020_DataSπ1_BodyMass.csv")
males1 <- ref1$Average.weight..gr...Male.[
  complete.cases(ref1$Average.weight..gr...Male.)]
names(males1) <- ref1$Species[complete.cases(ref1$Average.weight..gr...Male.)]
females1 <- ref1$Average.weight..gr...Female.[
  complete.cases(ref1$Average.weight..gr...Female.)]
names(females1) <- ref1$Species[complete.cases(ref1$Average.weight..gr...Female.)]

# Lislevand et al. 2007
ref2 <- read.csv("data/Lislevand2007_30may22.csv") 
ref2$Species_name <- stri_replace_all_fixed(ref2$Species_name, " ", "_")
ref2$M_mass[ref2$M_mass == -999] <- NA
ref2$F_mass[ref2$F_mass == -999] <- NA
males2 <- ref2$M_mass[complete.cases(ref2$M_mass)]
names(males2) <- ref2$Species_name[complete.cases(ref2$M_mass)]
females2 <- ref2$F_mass[complete.cases(ref2$F_mass)]
names(females2) <- ref2$Species_name[complete.cases(ref2$F_mass)]

# Myhrvold et al. 2015 (Amniote Database)
ref3 <- read.csv("data/AmnioteDatabase2015_02set2021.csv") 
ref3["Species"] <- NA
for (i in 1:nrow(ref3)) {
  ref3$Species[i] <- paste0(ref3$genus[i], "_", ref3$species[i])
}
ref3$female_body_mass_g[ref3$female_body_mass_g == -999] <- NA
ref3$male_body_mass_g[ref3$male_body_mass_g == -999] <- NA
males3 <- ref3$male_body_mass_g[complete.cases(ref3$male_body_mass_g)]
names(males3) <- ref3$Species[complete.cases(ref3$male_body_mass_g)]
females3 <- ref3$female_body_mass_g[complete.cases(ref3$female_body_mass_g)]
names(females3) <- ref3$Species[complete.cases(ref3$female_body_mass_g)]

###########
# Male

# Difference (male) source 1 (Ocampo et al. 2021) and 2 (Lislevand et al. 2007)
names_m12 <- names(males1)[names(males1) %in% names(males2)]
prop_m12 <- numeric()
for (i in 1:length(names_m12)) {
  m1 <- males1[names(males1) == names_m12[i]]
  m2 <- males2[names(males2) == names_m12[i]]
  prop_m12[i] <- abs(m1-m2)/min(m1, m2)
}
names(prop_m12) <- names_m12

# Difference (male) source 1 (Ocampo et al. 2021) and 3 (Myhrvold et al. 2015)
names_m13 <- names(males1)[names(males1) %in% names(males3)]
prop_m13 <- numeric()
for (i in 1:length(names_m13)) {
  m1 <- males1[names(males1) == names_m13[i]]
  m3 <- males3[names(males3) == names_m13[i]]
  prop_m13[i] <- abs(m1-m3)/min(m1, m3)
}
names(prop_m13) <- names_m13

# Difference (male) source 2 (Lislevand et al. 2007) and 3 (Myhrvold et al. 2015)
names_m23 <- names(males2)[names(males2) %in% names(males3)]
prop_m23 <- numeric()
for (i in 1:length(names_m23)) {
  m2 <- males2[names(males2) == names_m23[i]]
  m3 <- males3[names(males3) == names_m23[i]]
  prop_m23[i] <- abs(m2-m3)/min(m2, m3)
}
names(prop_m23) <- names_m23

##########
# Female

# Difference (female) source 1 (Ocampo et al. 2021) and 2 (Lislevand et al. 2007)
names_f12 <- names(females1)[names(females1) %in% names(females2)]
prop_f12 <- numeric()
for (i in 1:length(names_f12)) {
  f1 <- females1[names(females1) == names_f12[i]]
  f2 <- females2[names(females2) == names_f12[i]]
  prop_f12[i] <- abs(f1-f2)/min(f1, f2)
}
names(prop_f12) <- names_f12

# Difference (female) source 1 (Ocampo et al. 2021) and 3 (Myhrvold et al. 2015)
names_f13 <- names(females1)[names(females1) %in% names(females3)]
prop_f13 <- numeric()
for (i in 1:length(names_f13)) {
  f1 <- females1[names(females1) == names_f13[i]]
  f3 <- females3[names(females3) == names_f13[i]]
  prop_f13[i] <- abs(f1-f3)/min(f1, f3)
}
names(prop_f13) <- names_f13

# Difference (female) source 2 (Lislevand et al. 2007) and 3 (Myhrvold et al. 2015)
names_f23 <- names(females2)[names(females2) %in% names(females3)]
prop_f23 <- numeric()
for (i in 1:length(names_f23)) {
  f2 <- females2[names(females2) == names_f23[i]]
  f3 <- females3[names(females3) == names_f23[i]]
  prop_f23[i] <- abs(f2-f3)/min(f2, f3)
}
names(prop_f23) <- names_f23

layout(matrix(1:6, ncol = 3, byrow = T))

hist(prop_m12, main = "Male: Ocampo et al. 2021 vs Lislevand et al. 2007", 
     xlab = "", freq = F)
hist(prop_m13, main = "Male: Ocampo et al. 2021 vs Myhrvold et al. 2015", 
     xlab = "", ylab = "", freq = F)
hist(prop_m23, main = "Male: Lislevand et al. 2007 vs Myhrvold et al. 2015", 
     xlab = "", ylab = "", freq = F)
hist(prop_f12, main = "Female: Ocampo et al. 2021 vs Lislevand et al. 2007", 
     xlab = "", freq = F)
hist(prop_f13, main = "Female: Ocampo et al. 2021 vs Myhrvold et al. 2015", 
     xlab = "Proporção que a maior medida é maior em relação a pequena",
     freq = F, ylab = "")
hist(prop_f23, main = "Female: Lislevand et al. 2007 vs Myhrvold et al. 2015", 
     xlab = "", ylab = "", freq = F)

names12 <- c(names(prop_m12[prop_m12 >= 1]), names(prop_f12[prop_f12 >= 1]))
names13 <- c(names(prop_m13[prop_m13 >= 1]), names(prop_f13[prop_f13 >= 1]))
names23 <- c(names(prop_m23[prop_m23 >= 1]), names(prop_f23[prop_f23 >= 1]))
names12 <- unique(names12)
names13 <- unique(names13)
names23 <- unique(names23)

all_names <- c(names12, names13, names23)
all_names <- unique(all_names)

dat_nosso <- read.csv("data/BodySizeAves_30may22.csv", row.names = 1)

dat_nosso <- dat_nosso[complete.cases(dat_nosso$Body_mass_g_M_mean) & 
                       complete.cases(dat_nosso$Body_mass_g_F_mean), ]

dat_nosso <- dat_nosso[!(dat_nosso$Scientific_name %in% all_names), ]

write.csv(dat_nosso, "data/BodySizeAves_30may22_edit.csv")

sdi_nosso1 <- log10(dat_nosso$Body_mass_g_M_mean/dat_nosso$Body_mass_g_F_mean)
names(sdi_nosso1) <- dat_nosso$Scientific_name
sdi_nosso1 <- sdi_nosso1[complete.cases(sdi_nosso1)]

####################################
# Comparing SDI
####################################

# Gonzalez-Voyer et al. (2022)

dat1 <- read.csv("data/aves/SexroleEcologyFinal.csv")
dat1$species[dat1$species == "Brachypteracias_squamigera"] <- 
                                                    "Brachypteracias_squamiger"
dat1$species[dat1$species=="Nectarinia_neergaardi"] <- "Nectarinia_neergardi"

# Eles usaram uma métrica diferente (log10(Male/Female)) e fizeram mapa de SSD
# com média das medidas do corpo e asa

# Código do trabalho deles com algumas modificações (https://github.com/AlejandroG-V/AvianSexRoles/blob/main/Sex%20Roles%20Analyses%20Code.R) :

mass.ssd <- log10(dat1$m.mass/dat1$f.mass)
wing.ssd <- log10(dat1$m.wing/dat1$f.wing)
tarsus.ssd <- log10(dat1$m.tarsus/dat1$f.tarsus)
ssd <- as.data.frame(cbind(mass.ssd, wing.ssd, tarsus.ssd))
ssd$species <- dat1$species
ssd <- ssd[, c(4, 1, 2, 3)]
ssdMean <- rowMeans(ssd[, 2:4], na.rm = TRUE)
ssdMean <- as.data.frame(cbind(ssd$species, ssdMean), stringsAsFactors = FALSE)
colnames(ssdMean) <- c("species", "ssdM")
ssdMean$ssdM <- as.numeric(as.character(ssdMean$ssdM))
SSDMean <- ssdMean$ssdM/sd(na.omit(ssdMean$ssdM))

sdi1 <- ssdMean$ssdM
names(sdi1) <- ssdMean$species
sdi1 <- sdi1[complete.cases(sdi1)]

sdi1_stand <- SSDMean
names(sdi1_stand) <- ssdMean$species
sdi1_stand <- sdi1_stand[complete.cases(sdi1_stand)]

sdi_nosso_stand <- sdi_nosso1/sd(sdi_nosso1)
names(sdi_nosso_stand) <- names(sdi_nosso1)

range(sdi1) # -0.3685255  0.6005903
range(sdi_nosso1) #  -0.3979400  0.7979596
range(sdi1_stand) # -6.918551 11.275243
range(sdi_nosso_stand) # -5.660937 11.351459

#########

# Friedman & Remes (2015)

# Eles não tem código pq usaram o SAM

# Dados apenas de Lislevand et al. 2007
dat2 <- read.csv("data/Lislevand2007_30may22.csv") 
dat2$Species_name <- stri_replace_all_fixed(dat2$Species_name, " ", "_")
dat2$M_mass[dat2$M_mass == -999] <- NA
dat2$F_mass[dat2$F_mass == -999] <- NA

dat2 <- dat2[complete.cases(dat2$M_mass) & complete.cases(dat2$F_mass), ]

# Eles falam que conseguiram 2581 spp, mas só tem 1245 nessa base de dados 
# com valor de tamanho de corpo tanto de macho quanto de fêmea

# Função para calcular SDI como eles calcularam
SDI <- function (male = male, female = female) {
  if (male <= female) {
      SDI <- -((female/male) - 1)
    } else {
      if (female < male) {
        SDI <- (male/female) - 1
      } 
    }
  
  return(SDI)
}

sdi2 <- numeric()
for (i in 1:nrow(dat2)) {
  sdi2[i] <- SDI(male = dat2$M_mass[i], female = dat2$F_mass[i])
}
names(sdi2) <- dat2$Species_name
sdi2 <- sdi2[complete.cases(sdi2)]

sdi_nosso2 <- numeric()
for (i in 1:nrow(dat_nosso)) {
  sdi_nosso2[i] <- SDI(male = dat_nosso$Body_mass_g_M_mean[i], 
                       female = dat_nosso$Body_mass_g_F_mean[i])
}
names(sdi_nosso2) <- dat_nosso$Scientific_name
sdi_nosso2 <- sdi_nosso2[complete.cases(sdi_nosso2)]

range(sdi2) # -1.022472  2.138103
range(sdi_nosso2) # -1.50  5.28

####################################
# Testing the maps
####################################

# Mapa com letsR não tá funcionando, peguei o código que usaram no Tobias et al.
# (2022), do trabalho do AVONET. Não peguei dos trabalhos que mapearam SSD, 
# porque eles não tem os códigos dessa parte.

# Tobias et al. (2022) - Usaram a mesma base de dados de ranges das aves, mas 
# plotaram de uma maneira diferente. Eles primeiro convertaram os ranges em uma 
# matriz de presença-ausência, mas em relação ao shapefile do mapa
# mundi, o que eu não sei fazer. Então eu peguei essa matriz deles. Não tem 
# código explicando isso:

# "For the BirdLife dataset, species geographical range maps were obtained from 
# Birdlife International. We selected breeding and resident ranges in areas 
# where the species is coded as extant and either native or reintroduced.
# Mapping species under the BirdTree taxonomy dataset often required the 
# combination of maps for multiple BirdLife species to form an expanded range
# map for a single BirdTree species (Supplementary Dataset 1). We did not 
# generate new maps for the eBird dataset which can be integrated directly with 
# citizen-science locality data and a suite of spatial research tools (Sullivan 
# et al., 2014). [...] To map morphological traits, we extracted species ranges 
# onto an equal area grid (Behrmann projection) with a resolution of 110 km
# (≈1° at the equator)."

map.SSD <- function(data, func = mean, cols = NULL, figFolder, fileName) { 
  # de Tobias et al. (2022)
  
  # Behrmann equal area (96 x 96km) grid shapefile
  grid <- rgdal::readOGR("data/spatial/BehrmannMeterGrid_WGS84_land.shp")
  
  # Country borders shapefile
  countriesGeo <- rgdal::readOGR("data/spatial/all_countries.shp")
  
  # gridded species geographic range data - Birdlife taxonomy 
  rangeData <- read.csv("data/spatial/AllSpeciesBirdLifeMaps2019.csv")

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
  
  plot.ssd1 <- paste0(figFolder, "SSD_Map_", fileName, ".tiff")
  tiff(plot.ssd1, width = inches*res, height = inches*res/2, units = "px")
  print(ggplot.ssd)  
  dev.off()  

  plot.ssd2 <- paste0(figFolder, "SSD_Map_ScaleBar_", fileName, ".tiff")
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

########################

# Repetindo mapas de outros trabalhos para ver se código está funcionando

########################

# Gonzalez-Voyer et al. (2022)

# Eles dizem isso, então acredito que devem ter feito mais algo com os dados do 
# que está nesse script:
# "Spatial distributions of avian sex role components. (a–e) show the mean 
# values of each sex role component among species per 100 km grid cell. 
# In (a–d) the values are standardised and centred on zero (no bias in sex 
# role) with diverging colour palette to identify regions with male-biased
# (green-blue) or female-biased (yellow-red) sex roles. Panel e shows the log 
# of the extent of sex role bias. The colour ramps are scaled from the 1st to 
# 99th percentiles of the data to minimise the effects of outliers on 
# visualisation of variation."

# Legenda fica com valores errados, tem que fazer separada
# Figura salva em arquivo "figFolder/SSD_Map_fileName.tiff" e legenda em
# "figFolder/SSD_Map_ScaleBar_fileName.tiff" mas precisa colocar valores 
# nos ticks da legenda retornados da função abaixo

# Mas o mapa fica parecido, mesmo com os valores máximos e mínimos não sendo
# iguais
map.SSD(sdi1, cols = c(viridis(100), "orange")[101:1], figFolder = "figures/",
        fileName = "GV")

# Comparando com nosso mapa mesmo com aqueles valores errados:
map.SSD(sdi_nosso1, cols = c(viridis(100), "orange")[101:1], 
        figFolder = "figures/", fileName = "GV_nosso")


########################

# Friedman & Remes (2015)

# Eles não tem código pq usaram o SAM

# Dados apenas de Lislevand et al. 2007
dat2 <- read.csv("data/Lislevand2007_19aug21.csv") 
dat2$Species_name <- stri_replace_all_fixed(dat2$Species_name, " ", "_")
dat2$M_mass[dat2$M_mass == -999] <- NA
dat2$F_mass[dat2$F_mass == -999] <- NA

dat2 <- dat2[complete.cases(dat2$M_mass) & complete.cases(dat2$F_mass), ]

# Eles falam que conseguiram 2581 spp, mas só tem 1245 nessa base de dados 
# com valor de tamanho de corpo tanto de macho quanto de fêmea

# Função para calcular SDI como eles calcularam
SDI <- function (male = male, female = female) {
  if (male <= female) {
      SDI <- -((female/male) - 1)
    } else {
      if (female < male) {
        SDI <- (male/female) - 1
      } 
    }
  
  return(SDI)
}

sdi2 <- numeric()
for (i in 1:nrow(dat2)) {
	sdi2[i] <- SDI(male = dat2$M_mass[i], female = dat2$F_mass[i])
}
names(sdi2) <- dat2$Species_name
sdi2 <- sdi2[complete.cases(sdi2)]

range(sdi2)

map.SSD(sdi2, cols = c("white", "orange", "red"), figFolder = "figures/",
        fileName = "FR")

map.SSD(sdi_nosso, cols = c("white", "orange", "red"), figFolder = "figures/",
        fileName = "FR_nosso")

