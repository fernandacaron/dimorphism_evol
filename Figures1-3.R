rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(viridis)
library(colorspace)
library(letsR)
library(stringi)
library(maptools)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

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
pdf("figures/Figure1.pdf", height = 9)
layout(matrix(1:2, ncol = 1))

cols <- mako(7)
cols_al <- rgb(t(col2rgb(cols)), alpha = 150, maxColorValue = 255)

hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Apodiformes"]), 
     ylim = c(0, 1), xlim = c(0, 10), main = "", 
     xlab = "log body mass (g) males", breaks = 30, 
     col = cols_al[1], freq = F, border = cols[1])
hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Charadriiformes"]), 
     add = T, breaks = 40, col = cols_al[2], freq = F, border = cols[2])
hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Columbiformes"]), add = T, 
     breaks = 40, col = cols_al[3], freq = F, border = cols[3])
hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Passeriformes"]), add = T, 
     breaks = 40, col = cols_al[4], freq = F, border = cols[4])
hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Piciformes"]), add = T, 
     breaks = 40, col = cols_al[5], freq = F, border = cols[5])
hist(log(dat_red$Body_mass_g_M_mean[dat_red$Order == "Psittaciformes"]), 
     add = T, breaks = 40, col = cols_al[6], freq = F, border = cols[6])
title("A", adj = 0)

hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Apodiformes"]),
     ylim = c(0, 1), xlim = c(0, 10), main = "", 
     xlab = "log body mass (g) females", breaks = 30, 
     col = cols_al[1], freq = F, border = cols[1])
hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Charadriiformes"]), 
     add = T, breaks = 40, col = cols_al[2], freq = F, border = cols[2])
hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Columbiformes"]), add = T, 
     breaks = 40, col = cols_al[3], freq = F, border = cols[3])
hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Passeriformes"]), add = T, 
     breaks = 40, col = cols_al[4], freq = F, border = cols[4])
hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Piciformes"]), add = T, 
     breaks = 40, col = cols_al[5], freq = F, border = cols[5])
hist(log(dat_red$Body_mass_g_F_mean[dat_red$Order == "Psittaciformes"]), 
     add = T, breaks = 40, col = cols_al[6], freq = F, border = cols[6])

title("B", adj = 0)
legend("topright", pch = 15, bty = 'n', col = cols_al[1:6], 
       legend  = c("Apodiformes", "Charadriiformes", "Columbiformes",
                   "Passeriformes", "Piciformes", "Psittaciformes"))

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
	              cutoff = TRUE)
}
names(sdi) <- dat_red$Scientific_name
sdi <- sdi[complete.cases(sdi)]

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")

tr_map <- treedata(tr[[1]], sdi)$phy

sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- 1
sdi_disc[sdi_disc == 0] <- 0
sdi_disc[sdi_disc < 0] <- -1

cols_pal <- c("#9966FF", "gray", "#E69F00")

map <- make.simmap(tr_map, sdi_disc)
set_cols <- setNames(c("#9966FF", "gray", "#E69F00"), c(-1, 0, 1))

col1_lig <- lighten(cols_pal[1], amount = 0.4)
col2_lig <- lighten(cols_pal[2], amount = 0.4)
col3_lig <- lighten(cols_pal[3], amount = 0.4)

col_body <- character()
for (i in 1:length(body_size)) {
	col_body[i] <- 
		ifelse(sdi_disc[names(sdi_disc) == names(body_size)[i]] == 1, 
		       col3_lig, 
		       ifelse(sdi_disc[names(sdi_disc) == names(body_size)[i]] == -1,
		              col1_lig, col2_lig))
}
names(col_body) <- names(body_size)
col_body <- col_body[map$tip.label]

rbPal <- colorRampPalette(cols_pal)
cols <- rbPal(101)[as.numeric(cut(1:101, breaks = 101))]

pdf("figures/Figure2_cutoff.pdf", width = 9)

par(mar = c(0, 0, 1, 0))

plotTree.wBars(map, log(body_size), scale = 2, tip.labels = F,
    type = "fan", method = "plotSimmap", colors = set_cols, lwd = 1, 
    border = NA, col = col_body, mar = c(0, 0, 1, 0), part = 0.5)

legend("topright", col = cols_pal, bty = "n", pch = 15,
       legend = c("Male-biased SSD", "Monomorphism", "Female-biased SSD"))

par(xpd = TRUE)

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

# Figure 3

pam <- readRDS("data/aves/pam.rds")
colnames(pam$Presence_and_Absence_Matrix) <-
  stri_replace_all_fixed(colnames(pam$Presence_and_Absence_Matrix), " ", "_")
pam$Species_name <- stri_replace_all_fixed(pam$Species_name, " ", "_")

subsdi <- sdi_disc[names(sdi_disc) %in% pam$Species_name]

subsdi <- subsdi[pam$Species_name]

pam_mal <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% 
                                          names(subsdi[subsdi == -1])])
pam_fem <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% 
                                          names(subsdi[subsdi == 1])])

data(wrld_simpl)

pam_mal <- lets.pamcrop(pam_mal, wrld_simpl, remove.sp = TRUE)
pam_fem <- lets.pamcrop(pam_fem, wrld_simpl, remove.sp = TRUE)

pdf("figures/Figure3.pdf", height = 18, width = 15)

layout(matrix(1:2, ncol = 1))

par(mar = c(0, 0, 0, 4))

pal <- colorRampPalette(viridis(20))

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(pam_mal, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
title("Richness male-biased SSD", adj = 0, line = -3, cex.main = 2)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(pam_fem, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
title("Richness female-biased SSD", adj = 0, line = -3, cex.main = 2)

dev.off()


