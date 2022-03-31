rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(letsR)
library(stringi)
library(viridis)
library(maptools)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

pam <- readRDS("data/aves/pam.rds")
colnames(pam$Presence_and_Absence_Matrix) <-
    stri_replace_all_fixed(colnames(pam$Presence_and_Absence_Matrix), " ", "_")
pam$Species_name <- stri_replace_all_fixed(pam$Species_name, " ", "_")

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

dat_red <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
               complete.cases(dat$Body_mass_g_F_mean), c(5, 9, 13)]

sdi <- numeric()
for (i in 1:nrow(dat_red)) {
	sdi[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
	              female = dat_red$Body_mass_g_F_mean[i],
	              cutoff = FALSE)
}
names(sdi) <- dat_red$Scientific_name

sdi <- sdi[names(sdi) %in% pam$Species_name]

sdi <- sdi[pam$Species_name]
sdi <- sdi[sdi <= quantile(sdi, probs = 0.975) & 
             sdi >= quantile(sdi, probs = 0.025)]

abs_sdi <- abs(sdi)

res <- lets.maplizer(pam, abs_sdi, pam[[3]], ras = TRUE)

data(wrld_simpl)

# Tirar mapa da água, deixar só terra
res_crop <- lets.pamcrop(res, wrld_simpl, remove.sp = TRUE)
pam_crop <- lets.pamcrop(pam, wrld_simpl, remove.sp = TRUE)

# https://sedac.ciesin.columbia.edu/data/set/hanpp-net-primary-productivity/data-download
npp <- raster("data/npp-geotiff/npp_geotiff.tif")

pdf("figures/map.pdf")

layout(matrix(1:3, ncol = 1))

par(mar = c(4, 0, 1, 0))

pal <- colorRampPalette(viridis(20))

map("world", fill = TRUE, col = "white", bg = "white", border = NA)
plot(pam_crop$Richness_Raster, add = TRUE, col = pal(20))
title("Richness")

map("world", fill = TRUE, col = "white", bg = "white", border = NA)
plot(res_crop$Raster, add = TRUE, col = pal(200))
title("Absolute SSD")

map("world", fill = TRUE, col = "white", bg = "white", border = NA)
plot(npp, add = TRUE, col = pal(200))
title("NPP")

dev.off()

