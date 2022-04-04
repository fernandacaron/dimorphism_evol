rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(letsR)
library(stringi)
library(viridis)
library(maptools)
library(colorspace)

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

## Figure 2

pdf("figures/Figure2.pdf", width = 14)

layout(matrix(1:2, ncol = 2))

col1 <- "#9966FF"
col3 <- "#E69F00"

col1_lig <- lighten(col1, amount = 0.5)
col3_lig <- lighten(col3, amount = 0.5)

q <- c(-1, -0.1, 0, 0.1, 1)

h1 <- hist(sdi, breaks = 30, plot = FALSE)
cuts <- cut(h1$breaks[-1], c(-Inf, 0, Inf))
plot(h1, col = c(col1_lig, col3_lig)[cuts], border = c(col1, col3)[cuts],
     xlab = "SDI", main = "")
abline(v = q[1], col = "black")
abline(v = q[2], col = "black")
abline(v = q[3], col = "black")
abline(v = q[4], col = "black")
abline(v = q[5], col = "black")

sdi_95 <- sdi[sdi <= quantile(sdi, probs = 0.975) &
                sdi >= quantile(sdi, probs = 0.025)]
h2 <- hist(sdi_95, breaks = 30, plot = FALSE)
cuts <- cut(h2$breaks[-1], c(-Inf, 0, Inf))
plot(h2, col = c(col1_lig, col3_lig)[cuts], border = c(col1, col3)[cuts],
     xlab = "95% SDI", main = "")
abline(v = q[1], col = "black")
abline(v = q[2], col = "black")
abline(v = q[3], col = "black")
abline(v = q[4], col = "black")
abline(v = q[5], col = "black")

legend("topright", legend = c("Male-biased SSD", "Female-biased SSD"),
       fill = c(col1_lig, col3_lig), border = c(col1, col3), bty = "n")

dev.off()

# Figure 3

sdi <- sdi[names(sdi) %in% pam$Species_name]

sdi <- sdi[pam$Species_name]

q1 <- names(subset(sdi, sdi<q[1]))
q2 <- names(subset(sdi, sdi>=q[1] & sdi<=q[2]))
q3 <- names(subset(sdi, sdi>q[2] & sdi<q[3]))
q0 <- names(subset(sdi, sdi==q[3]))
q4 <- names(subset(sdi, sdi>q[3] & sdi<q[4]))
q5 <- names(subset(sdi, sdi>=q[4] & sdi<=q[5]))
q6 <- names(subset(sdi, sdi>q[5]))

res1 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q1])
res2 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q2])
res3 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q3])
res0 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q0])
res4 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q4])
res5 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q5])
res6 <- lets.subsetPAM(pam, pam[[3]][pam[[3]] %in% q6])

data(wrld_simpl)

res <- lets.pamcrop(pam, wrld_simpl, remove.sp = TRUE)
res1 <- lets.pamcrop(res1, wrld_simpl, remove.sp = TRUE)
res2 <- lets.pamcrop(res2, wrld_simpl, remove.sp = TRUE)
res3 <- lets.pamcrop(res3, wrld_simpl, remove.sp = TRUE)
res0 <- lets.pamcrop(res0, wrld_simpl, remove.sp = TRUE)
res4 <- lets.pamcrop(res4, wrld_simpl, remove.sp = TRUE)
res5 <- lets.pamcrop(res5, wrld_simpl, remove.sp = TRUE)
res6 <- lets.pamcrop(res6, wrld_simpl, remove.sp = TRUE)

pdf("figures/Figure3.pdf", width = 20, height = 8)

layout(matrix(1:9, ncol = 3))

par(mar = c(0, 0, 0, 4))

pal <- colorRampPalette(viridis(20))

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext("Richness", side = 3, line = -1.5, adj = 0, cex = 0.9)

npp <- raster("data/npp-geotiff/npp_geotiff.tif")
map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(npp, add = TRUE, col = pal(200))
mtext("NPP", side = 3, line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 5))
plot(res1, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(min(sdi), 2), " - ", round(q[1], 2)), side = 3,
      line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res2, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(q[1], 2), " - ", round(q[2], 2)), side = 3,
      line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res3, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(q[2], 2), " - ", round(q[3], 2)), side = 3,
      line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 5))
plot(res0, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext("SDI = 0", side = 3, line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res4, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(q[3], 2), " - ", round(q[4], 2)), side = 3, 
      line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res5, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(q[4], 2), " - ", round(q[5],2)), side = 3,
      line = -1.5, adj = 0, cex = 0.9)

map("world", fill = TRUE, col = "gray", bg = "white", border = NA, 
    mar = c(0, 0, 0, 4))
plot(res6, axes = FALSE, box = FALSE, col_rich = pal, world = FALSE, 
     plot = FALSE, add = TRUE)
mtext(paste0("SDI = ", round(q[5], 2), " - ", round(max(sdi), 2)), side = 3,
      line = -1.5, adj = 0, cex = 0.9)

dev.off()

####################################








title("Absolute SSD")

map("world", fill = TRUE, col = "gray", bg = "white", border = NA)
plot(npp, add = TRUE, col = pal(200))




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

