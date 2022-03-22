rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(letsR)
library(rgdal)
library(raster)
library(ecostructure)

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

#https://lpdaac.usgs.gov/product_search/?keyword=Net+Primary+Productivity+%28NPP%29&view=cards&sort=title
#https://lpdaac.usgs.gov/product_search/?keyword=Gross+Primary+Productivity+%28GPP%29&view=cards&sort=title


bird <- sf::st_read(dsn = "data/aves/BOTW/BOTW.gdb", layer = "All_Species",
                    type = 6)

## Presence-absence matrix
presab <- dsp_create_from_gdb(gdb_object = bird, raster_resolution = 5,
                              thresh = 3, raster_latlim = c(-90,90),
                              raster_longlim = c(-180,180),
                              species_feature = "sci_name")

save(presab, file = "data/aves/presab.RData")





#save(shapes, file="data/aves/shapes.Rdata")
#load("data/aves/shapes.Rdata")

head(shapes@data)

maps <- lets.presab(shapes, resol = 3, cover = 0.01, 
                    remove.cells = F) #evita a remocao de celulas sem info
#save(maps, file="data/aves/maps.Rdata")
#load("data/aves/maps.Rdata")

plot(maps, axes = F, main = "Birds Richness")

summary(maps)
str(maps)

pa <- maps$Presence_and_Absence_Matrix[, -c(1:2)]

# mudando nome de acordo com a filogenia funcao
muda_nome <- function(x) {
	unlist(lapply(strsplit(x = as.character(x), " "), function(x) { 
		paste(x[1], x[2], sep = "_")
	}
	)
	)
}

# mudando nome da matrix local vs espécie
colnames(pa) = muda_nome(colnames(pa))

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

sdi <- sdi[complete.cases(sdi)]

sp <- intersect(names(sdi), colnames(pa))

shapes <- shapes[shapes$binomialMod %in% sp, ]
dim(shapes)

pa <- pa[, colnames(pa) %in% sp] 
dim(pa)

# basta mudar os valores do raster original pela metrica desejada PDiv
maps2 <- maps

values(maps2$XXX) <- SDI
plot(maps2, main = "SDI")

library(gdalUtils)
setwd("~Documents/NPP")
files <- dir(pattern = ".hdf")

filename <- substr(files, 12, 44)
filename <- paste0("NPP", filename, ".tif")

for (i in 1:length(filename)) {
  sds <- get_subdatasets(files[i])
  gdal_translate(sds[1], dst_dataset = filename[i])
}

shapes_npp <- readOGR("~Documents/NPP")

npp <- lets.presab(shapes_npp, resol = 3, cover = 0.01, remove.cells = F)
#save(npp, file="data/npp.Rdata")
#load("data/npp.Rdata")
plot(npp, axes = F, main = "NPP")


pamA <- lets.presab(amphibians, crs=CRS(amphibians@proj4string@projargs), resol = 1)


map("world",fill=TRUE, col="grey", bg="white", border=NA)
title("Amphibia", line = 2)
plot(pamA, axes = FALSE, box = FALSE, col_rich = inferno, world = FALSE, plot = FALSE, add = TRUE)
