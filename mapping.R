rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(letsR)

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

#https://lpdaac.usgs.gov/product_search/?keyword=Net+Primary+Productivity+%28NPP%29&view=cards&sort=title
#https://lpdaac.usgs.gov/product_search/?keyword=Gross+Primary+Productivity+%28GPP%29&view=cards&sort=title
squamataq5 <- readOGR("data/REPTILES_IUCN", layer="REPTILES", dropNULLGeometries = TRUE)[readOGR("data/REPTILES_IUCN", layer="REPTILES", dropNULLGeometries = TRUE)@data$family %in% qS5,]
load("data/AMPHIBIANS_IUCN/amphibians2.Rdata")
pamA <- lets.presab(amphibians, crs=CRS(amphibians@proj4string@projargs), resol = 1)


map("world",fill=TRUE, col="grey", bg="white", border=NA)
title("Amphibia", line = 2)
plot(pamA, axes = FALSE, box = FALSE, col_rich = inferno, world = FALSE, plot = FALSE, add = TRUE)
