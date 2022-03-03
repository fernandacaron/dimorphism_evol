rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(stringi)

## Função para colocar automaticamente os dados dos traits respectivos nos 
## dados (data) do grupo indicado. "data" - dataframe onde serão colocados os 
## dados com coluna "Scientific_name" e coluna do trait como "traitCol", "ref" 
## - de onde serão retirados os dados com coluna trait de nome "trait" a "data" 
## e com coluna "Species"

getSppTraits <- function(data, ref, traitCol, trait) {
	for (i in 1:nrow(data)) {
		data[colnames(data) == traitCol][i, ] <- 
			ifelse(length(ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ]) == 0,
			       data[colnames(data) == traitCol][i,],
			       ifelse(is.na(ref[colnames(ref) == 
                          trait][ref$Species == data$Scientific_name[i], ]),
			              data[colnames(data) == traitCol][i, ], 
			              ref[colnames(ref) == 
                    trait][ref$Species == data$Scientific_name[i], ]))
		}
  return(data)
}

#Reptilia

dat_sq <- read.csv("data/squamata/taxonomySqu.csv")

# Myhrvold et al. 2015 (Amniote Database) - Body_mass_g_F_1, Body_mass_g_M_1,
## SVL_mm_F_1 e SVL_mm_M_1
ref_sq_1 <- read.csv("data/squamata/AmnioteDatabase2015_02set2021.csv") 
ref_sq_1["Species"] <- NA
for (i in 1:nrow(ref_sq_1)) {
  ref_sq_1$Species[i] <- paste0(ref_sq_1$genus[i], "_", ref_sq_1$species[i])
}
ref_sq_1$female_body_mass_g[ref_sq_1$female_body_mass_g == -999] <- NA
ref_sq_1$male_body_mass_g[ref_sq_1$male_body_mass_g == -999] <- NA
ref_sq_1$male_svl_cm[ref_sq_1$male_svl_cm == -999] <- NA
ref_sq_1$male_svl_cm <- (ref_sq_1$male_svl_cm)*10 #SVL em cm, transformar em mm
ref_sq_1$female_svl_cm[ref_sq_1$female_svl_cm == -999] <- NA
ref_sq_1$female_svl_cm <- (ref_sq_1$female_svl_cm)*10 #SVL em cm, transformar em mm

dat_sq["SVL_mm_F_1"] <- dat_sq["SVL_mm_M_1"] <- dat_sq["Body_mass_g_F_1"] <- 
  dat_sq["Body_mass_g_M_1"] <- NA

dat_sq <- getSppTraits(dat_sq, ref_sq_1, "Body_mass_g_F_1", 
                       "female_body_mass_g")
dat_sq <- getSppTraits(dat_sq, ref_sq_1, "Body_mass_g_M_1", "male_body_mass_g")
dat_sq <- getSppTraits(dat_sq, ref_sq_1, "SVL_mm_F_1", "female_svl_cm")
dat_sq <- getSppTraits(dat_sq, ref_sq_1, "SVL_mm_M_1", "male_svl_cm")

write.csv(dat_sq, "data/squamata/BodySizeSqu_18jan22.csv")

##################

#Aves

dat_av <- read.csv("data/aves/taxonomyAve.csv")

#Ocampo et al. 2021 - Body_mass_g_M_1 e Body_mass_g_F_1
ref_av_1 <- read.csv("data/aves/Ocampo&al2020_DataS1_BodyMass.csv")

#Lislevand et al. 2007 - Body_mass_g_M_2 e Body_mass_g_F_2
ref_av_2 <- read.csv("data/aves/Lislevand2007_19aug21.csv") 
colnames(ref_av_2)[colnames(ref_av_2) == "Species_name"] <- "Species"
ref_av_2$Species <- stri_replace_all_fixed(ref_av_2$Species, " ", "_")
ref_av_2$M_mass[ref_av_2$M_mass == -999] <- NA
ref_av_2$F_mass[ref_av_2$F_mass == -999] <- NA

#Myhrvold et al. 2015 (Amniote Database) - Body_mass_g_M_3 e Body_mass_g_F_3
ref_av_3 <- read.csv("data/aves/AmnioteDatabase2015_02set2021.csv") 
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

write.csv(dat_av, "data/aves/BodySizeAve_18jan22.csv")

## Mammals

dat_ma <- read.csv("data/mammalia/taxonomyMam.csv")

#Ocampo et al. 2021 - Body_mass_g_M_1 e Body_mass_g_F_1
ref_ma_1 <- read.csv("data/mammalia/Ocampo&al2020_DataS1_BodyMass.csv") 

#Sibly et al. 2012 - Body_mass_g_M_1 e Body_mass_g_F_1
ref_ma_2 <- read.csv("data/mammalia/Sibly&al2012AmNat_BodySizeHerb.csv")
ref_ma_2["Species"] <- NA
for (i in 1:nrow(ref_ma_2)) {
  ref_ma_2$Species[i] <- paste0(ref_ma_2$Genus[i], "_", ref_ma_2$species[i])
}
ref_ma_2$Male.mass..kg. <- ref_ma_2$Male.mass..kg. * 1000
ref_ma_2$Female.mass..kg. <- ref_ma_2$Female.mass..kg. * 1000

dat_ma["Body_mass_g_M_mean"] <- dat_ma["Body_mass_g_M_2"] <- 
  dat_ma["Body_mass_g_M_1"] <- dat_ma["Body_mass_g_F_mean"] <- 
  dat_ma["Body_mass_g_F_2"] <- dat_ma["Body_mass_g_F_1"] <- NA

dat_ma <- getSppTraits(dat_ma, ref_ma_1, "Body_mass_g_M_1", 
                       "Average.weight..gr...Male.")
dat_ma <- getSppTraits(dat_ma, ref_ma_1, "Body_mass_g_F_1", 
                       "Average.weight..gr...Female.")
dat_ma <- getSppTraits(dat_ma, ref_ma_2, "Body_mass_g_M_2", 
                       "Male.mass..kg.")
dat_ma <- getSppTraits(dat_ma, ref_ma_2, "Body_mass_g_F_2", 
                       "Female.mass..kg.")

for (i in 1:nrow(dat_ma)) {
  y <- c(dat_ma$Body_mass_g_M_1[i], dat_ma$Body_mass_g_M_2[i])
  y <- as.numeric(y)
  if (length(is.na(y)[is.na(y) == "TRUE"]) != 2) {
    dat_ma$Body_mass_g_M_mean[i] <- mean(y, na.rm = T)
  } else {
    dat_ma$Body_mass_g_M_mean[i] <- NA
  }
  
  z <- c(dat_ma$Body_mass_g_F_1[i], dat_ma$Body_mass_g_F_2[i])
  z <- as.numeric(z)
  if (length(is.na(z)[is.na(z) == "TRUE"]) != 2) {
    dat_ma$Body_mass_g_F_mean[i] <- mean(z, na.rm = T)
  } else {
    dat_ma$Body_mass_g_F_mean[i] <- NA
  }
}

write.csv(dat_ma, "data/mammalia/BodySizeMam_26jan22.csv")
