rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(sp)
library(raster)
library(letsR)
library(stringi)
library(caper)

# Presence-absence matrix
pam <- readRDS("data/aves/pam.rds")
colnames(pam$Presence_and_Absence_Matrix) <-
  stri_replace_all_fixed(colnames(pam$Presence_and_Absence_Matrix), " ", "_")
pam$Species_name <- stri_replace_all_fixed(pam$Species_name, " ", "_")

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

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

sdi <- sdi[names(sdi) %in% pam$Species_name]
sdi <- sdi[pam$Species_name]

# Latitude midpoint

# Midpoint da longitude (x) e latitude (y) - Não sei se o método é esse mesmo, 
# mas visualizando parece ok
pam1 <- lets.subsetPAM(pam, pam$Species_name[1:500], remove.cells = T)
pam2 <- lets.subsetPAM(pam, pam$Species_name[501:1000], remove.cells = T)
pam3.1 <- lets.subsetPAM(pam, pam$Species_name[1001:1125], remove.cells = T)
pam3.2 <- lets.subsetPAM(pam, pam$Species_name[1126:1250], remove.cells = T)
pam3.3 <- lets.subsetPAM(pam, pam$Species_name[1250:1350], remove.cells = T)
pam3.4 <- lets.subsetPAM(pam, pam$Species_name[1351:1500], remove.cells = T)
pam4.1 <- lets.subsetPAM(pam, pam$Species_name[1501:1625], remove.cells = T)
pam4.2 <- lets.subsetPAM(pam, pam$Species_name[1626:1750], remove.cells = T)
pam4.3 <- lets.subsetPAM(pam, pam$Species_name[1751:2000], remove.cells = T)
pam5 <- lets.subsetPAM(pam, pam$Species_name[2001:2500], remove.cells = T)
pam6.1 <- lets.subsetPAM(pam, pam$Species_name[2501:2750], remove.cells = T)
pam6.2 <- lets.subsetPAM(pam, pam$Species_name[2751:3000], remove.cells = T)
pam7 <- lets.subsetPAM(pam, pam$Species_name[3001:3500], remove.cells = T)
pam8 <- lets.subsetPAM(pam, pam$Species_name[3501:4000], remove.cells = T)
pam9 <- lets.subsetPAM(pam, pam$Species_name[4001:4080], remove.cells = T)

mid1 <- lets.midpoint(pam1, method = "CMD")
mid2 <- lets.midpoint(pam2, method = "CMD")
mid3.1 <- lets.midpoint(pam3.1, method = "CMD")
mid3.2 <- lets.midpoint(pam3.2, method = "CMD")
mid3.3 <- lets.midpoint(pam3.3, method = "CMD")
mid3.3 <- lets.midpoint(pam3.4, method = "CMD")
mid4.1 <- lets.midpoint(pam4.1, method = "CMD")
mid4.2 <- lets.midpoint(pam4.2, method = "CMD")
mid4.3 <- lets.midpoint(pam4.3, method = "CMD")
mid5 <- lets.midpoint(pam5, method = "CMD")
mid6.1 <- lets.midpoint(pam6.1, method = "CMD")
mid6.2 <- lets.midpoint(pam6.2, method = "CMD")
mid7 <- lets.midpoint(pam7, method = "CMD")
mid8 <- lets.midpoint(pam8, method = "CMD")
mid9 <- lets.midpoint(pam9, method = "CMD")

pam1 <- lets.subsetPAM(pam, pam$Species_name[1:2040], remove.cells = T)

mid <- lets.midpoint(pam, method = "CMD")

# Climatic data
bio <- getData('worldclim', var = 'bio', res = 10)

# https://pubs.usgs.gov/ds/691/ds691.pdf

# bio1: Annual Mean Temp
bio1 <- bio@layers[[1]]
bio1 <- bio1/10
pam_bio1 <- as.data.frame(lets.addvar(pam, bio1))

# bio12: Annual Prec
bio12 <- bio@layers[[12]]
pam_bio12 <- as.data.frame(lets.addvar(pam, bio12))

# bio4: Season Temp
bio4 <- bio@layers[[4]]
pam_bio4 <- as.data.frame(lets.addvar(pam, bio4))

# bio15: Season Prec
bio15 <- bio@layers[[15]]
pam_bio15 <- as.data.frame(lets.addvar(pam, bio15))

npp <- raster("data/npp-geotiff/npp_geotiff.tif")
pam_npp <- as.data.frame(lets.addvar(pam, npp))

birds <- data.frame(Species = pam$Species_name, bio1 = NA, bio12 = NA, 
                    bio4 = NA, bio15 = NA, npp = NA, SDI = NA)
for (i in 1:length(pam$Species_name)) {
  spp <- pam$Species_name[[i]]
  
  spp_bio1 <- pam_bio1$bio1_mean[pam_bio1[, spp] == 1]
  spp_bio12 <- pam_bio12$bio12_mean[pam_bio12[, spp] == 1]
  spp_bio4 <- pam_bio4$bio4_mean[pam_bio4[, spp] == 1]
  spp_bio15 <- pam_bio15$bio15_mean[pam_bio15[, spp] == 1]
  spp_npp <- pam_npp$npp_geotiff_mean[pam_npp[, spp] == 1]
  
  birds$bio1[i] <- mean(spp_bio1, na.rm = TRUE)
  birds$bio12[i] <- mean(spp_bio12, na.rm = TRUE)
  birds$bio4[i] <- mean(spp_bio4, na.rm = TRUE)
  birds$bio15[i] <- mean(spp_bio15, na.rm = TRUE)
  birds$npp[i] <- mean(spp_npp, na.rm = TRUE)
  
  birds$SDI[i] <- sdi[names(sdi) == spp]
}

# We excluded highly correlated variables (assessed by variance inflation 
# factor values >10) to minimize multicollinearity.
for (i in 1:10) {
  comp_data <- comparative.data(tr[[i]], birds, names.col = Species)
  pgls <- pgls(SDI ~ bio1 + bio12 + bio4 + bio15 + npp, comp_data)
}
