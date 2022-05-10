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

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")[[1]]

# Climatic data
bio <- getData('worldclim', var = 'bio', res = 10)

# bio1: Annual Mean Temp
bio1 <- bio@layers[[1]]
bio1 <- bio1/10
pam_bio1 <- as.data.frame(lets.addvar(pam, bio1))

# bio12: Annual Prec
bio12 <- bio@layers[[12]]
bio12 <- bio12/10
pam_bio12 <- as.data.frame(lets.addvar(pam, bio12))

# bio4: Season Temp
bio4 <- bio@layers[[4]]
bio4 <- bio4/10
pam_bio4 <- as.data.frame(lets.addvar(pam, bio4))

# bio15: Season Prec
bio15 <- bio@layers[[15]]
bio15 <- bio15/10
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

#for (i in 1:3) {
  comp_data <- comparative.data(tr, birds, names.col = Species)
  pgls <- pgls(SDI ~ bio1 + bio12 + bio4 + bio15 + npp, comp_data)
#}
