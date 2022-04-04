rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(viridis)
library(colorspace)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

## Ref Blanckenhorn et al. 2007 (Livro Fairbairn et al. 2007, Chapter 6, p. 65)
## SDI = (female size/male size - 1) when female larger
## SDI = -(male size/female size - 1) when male larger

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

## 
sdi <- sdi[complete.cases(sdi)]
tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")

tr_map <- treedata(tr[[1]], sdi)$phy

sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- 1
sdi_disc[sdi_disc == 0] <- 0
sdi_disc[sdi_disc < 0] <- -1

col1 <- "#9966FF"
col2 <- "gray"
col3 <- "#E69F00"
cols_pal <- c(col1, col2, col3)

map <- contMap(tr_map, sdi_disc, plot = F)
plot <- setMap(map, cols_pal)

body_size <- numeric()
for (i in 1:nrow(dat_red)) {
	body_size[i] <- mean(c(dat_red$Body_mass_g_M_mean[i], 
	                       dat_red$Body_mass_g_F_mean[i]))
}
names(body_size) <- dat_red$Scientific_name 

col1_lig <- lighten(col1, amount = 0.4)
#rgb(153/255, 102/255, 255/255, 0.5)
col2_lig <- lighten(col2, amount = 0.4)
#rgb(190/255, 190/255, 190/255, 0.5)
col3_lig <- lighten(col3, amount = 0.4)
#rgb(230/255, 159/255, 0/255, 0.5)
col_body <- character()
for (i in 1:length(body_size)) {
	col_body[i] <- 
		ifelse(sdi_disc[names(sdi_disc) == names(body_size)[i]] == 1, 
		       col3_lig, 
		       ifelse(sdi_disc[names(sdi_disc) == names(body_size)[i]] == -1,
		              col1_lig, col2_lig))
}
names(col_body) <- names(body_size)
col_body <- col_body[plot$tree$tip.label]

rbPal <- colorRampPalette(cols_pal)
cols <- rbPal(101)[as.numeric(cut(1:101, breaks = 101))]


## Figure 1

pdf("figures/Figure1.pdf", width = 9)

par(mar = c(0, 0, 1, 0))

plotTree.wBars(plot$tree, log(body_size), scale = 2, tip.labels = F,
    type = "fan", method = "plotSimmap", colors = plot$cols, lwd = 1, 
    border = NA, col = col_body, mar = c(0, 0, 1, 0), part = 0.5)

add.color.bar(40, cols = cols, title = "SSD", c(-1, 1), prompt = F, 
             x = 0.9*par()$usr[1], y = -10, subtitle = "")

par(xpd = TRUE)

arc.cladelabels(text = "Apodiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Apodiformes"][dat$Scientific_name[dat$Order == 
                                  "Apodiformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Charadriiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"][dat$Scientific_name[dat$Order == 
                                  "Charadriiformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Columbiformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Columbiformes"][dat$Scientific_name[dat$Order == 
                                  "Columbiformes"] %in% names(sdi_disc)])))

arc.cladelabels(text = "Passeriformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Passeriformes"][dat$Scientific_name[dat$Order == 
                                  "Passeriformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Piciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.25, lab.offset = 1.3, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Piciformes"][dat$Scientific_name[dat$Order == 
                                  "Piciformes"] %in% names(sdi_disc)])))
arc.cladelabels(text = "Psittaciformes", mark.node = F, cex = 0.6, 
                col = "black", lwd = 1, ln.offset = 1.35, lab.offset = 1.4, 
                node = findMRCA(plot$tree, 
                                c(dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"][dat$Scientific_name[dat$Order == 
                                  "Psittaciformes"] %in% names(sdi_disc)])))

dev.off()
