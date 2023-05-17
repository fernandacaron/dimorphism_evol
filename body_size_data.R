rm(list = ls())

setwd("Documents/lab/dimorph_evol")

# Code to compile datasets

library(phytools)
library(stringi)

## This function will automatically place the data of the respective traits in 
## the data of the indicated taxa

getSppTraits <- function(data, ref, traitCol, trait) {

	for (i in 1:nrow(data)) {
		data[colnames(data) == traitCol][i, ] <- 
			ifelse(length(ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ]) == 0,
			       data[colnames(data) == traitCol][i, ],
             ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ])
		}
  return(data)
}


dat_av <- read.csv("data/taxonomyAve.csv")

#Ocampo et al. 2021 - Body_mass_g_M_1 e Body_mass_g_F_1
ref_av_1 <- read.csv("data/Ocampo&al2020_DataS1_BodyMass.csv")

#Lislevand et al. 2007 - Body_mass_g_M_2 e Body_mass_g_F_2
ref_av_2 <- read.csv("data/Lislevand2007_30may22.csv") 
colnames(ref_av_2)[colnames(ref_av_2) == "Species_name"] <- "Species"
ref_av_2$Species <- stri_replace_all_fixed(ref_av_2$Species, " ", "_")
ref_av_2$M_mass[ref_av_2$M_mass == -999] <- NA
ref_av_2$F_mass[ref_av_2$F_mass == -999] <- NA

#Myhrvold et al. 2015 (Amniote Database) - Body_mass_g_M_3 e Body_mass_g_F_3
ref_av_3 <- read.csv("data/AmnioteDatabase2015_02set2021.csv") 
ref_av_3["Species"] <- NA
for (i in 1:nrow(ref_av_3)) {
  ref_av_3$Species[i] <- paste0(ref_av_3$genus[i], "_", ref_av_3$species[i])
}
ref_av_3$female_body_mass_g[ref_av_3$female_body_mass_g == -999] <- NA
ref_av_3$male_body_mass_g[ref_av_3$male_body_mass_g == -999] <- NA

dat_av["Body_mass_g_F_mean"] <- dat_av["Body_mass_g_F_3"] <- 
  dat_av["Body_mass_g_F_2"] <- dat_av["Body_mass_g_F_1"] <- 
  dat_av["Body_mass_g_M_mean"] <- dat_av["Body_mass_g_M_3"] <- 
  dat_av["Body_mass_g_M_2"] <- dat_av["Body_mass_g_M_1"] <-  NA

dat_av <- getSppTraits(dat_av, ref_av_1, "Body_mass_g_M_1", 
                    "Average.weight..gr...Male.")
dat_av <- getSppTraits(dat_av, ref_av_1, "Body_mass_g_F_1", 
                    "Average.weight..gr...Female.")
dat_av <- getSppTraits(dat_av, ref_av_2, "Body_mass_g_M_2", "M_mass")
dat_av <- getSppTraits(dat_av, ref_av_2, "Body_mass_g_F_2", "F_mass")
dat_av <- getSppTraits(dat_av, ref_av_3, "Body_mass_g_M_3", "male_body_mass_g")
dat_av <- getSppTraits(dat_av, ref_av_3, "Body_mass_g_F_3", 
                       "female_body_mass_g")

for (i in 1:nrow(dat_av)) {
  y <- c(dat_av$Body_mass_g_M_1[i], dat_av$Body_mass_g_M_2[i], 
       dat_av$Body_mass_g_M_3[i])
  y <- as.numeric(y)
  if (length(is.na(y)[is.na(y) == "TRUE"]) != 3) {
    dat_av$Body_mass_g_M_mean[i] <- mean(y, na.rm = T)
  } else {
    dat_av$Body_mass_g_M_mean[i] <- NA
  }
  
  z <- c(dat_av$Body_mass_g_F_1[i], dat_av$Body_mass_g_F_2[i], 
       dat_av$Body_mass_g_F_3[i])
  z <- as.numeric(z)
  if (length(is.na(z)[is.na(z) == "TRUE"]) != 3) {
    dat_av$Body_mass_g_F_mean[i] <- mean(z, na.rm = T)
  } else {
    dat_av$Body_mass_g_F_mean[i] <- NA
  }
}

write.csv(dat_av, "data/BodySizeAve_30may22.csv")

