rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

# load packages
library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(letsR)
library(stringi)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)
birds <- st_read(dsn = "data/aves/BOTW/BOTW.gdb", layer = "All_Species")

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

sel_spp <- birds %>% filter(sci_name %in% 
                            stri_replace_all_fixed(names(sdi), "_", " ")) 

bird_ranges <- sel_spp %>% dplyr::select(id_no, sci_name, presence, Shape)

bird_ranges <- bird_ranges %>%
               filter(st_geometry_type(Shape) == "MULTIPOLYGON")

#saveRDS(bird_ranges, "data/aves/bird_ranges.rds")
#bird_ranges <- readRDS("data/aves/bird_ranges.rds")

bird_ranges_spatial <- as_Spatial(bird_ranges)

#saveRDS(bird_ranges_spatial, "data/aves/bird_ranges_spatial.rds")
#bird_ranges_spatial <- readRDS("data/aves/bird_ranges_spatial.rds")

# Código abaixo modificado com base na função lets.presab do pacote letsR 

xmn <- -180
xmx <- 180
ymn <- -90
ymx <- 90
resol <- 1
remove.cells <- TRUE
remove.sp <- TRUE
show.matrix <- FALSE
crs <- CRS(bird_ranges_spatial@proj4string@projargs)
crs.grid <- crs
cover <- 0

ras <- raster(xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs.grid)
res(ras) <- resol
values(ras) <- 1
ncellras <- ncell(ras)
coord <- xyFromCell(object = ras, cell = 1:ncellras)
colnames(coord) <- c("Longitude(x)", "Latitude(y)")
areashape <- NULL
areagrid <- NULL

nomes <- sort(unique(bird_ranges_spatial$sci_name))
n <- length(bird_ranges_spatial$sci_name)
nomes2 <- bird_ranges_spatial$sci_name
nomes <- nomes[nomes %in% nomes2]
matriz <- matrix(0, ncol = length(nomes), nrow = ncellras)
colnames(matriz) <- nomes


extractpos <- function (ras, shapepol, nomes, nomes2, cover, areashape, 
                        areagrid, i) { 
    shapepol1 <- shapepol
    celulas <- try(celulas <- extract(ras, shapepol1, cellnumbers = TRUE, 
        weights = TRUE, small = TRUE, normalizeWeights = FALSE), 
        silent = T)
    if (class(celulas) != "try-error") {
        celulas2 <- extract(ras, shapepol1, cellnumbers = TRUE, 
            weights = FALSE, small = TRUE, normalizeWeights = FALSE)
        keep_this <- celulas[[1]][, 1] %in% celulas2[[1]][, 1]
        celulas[[1]] <- celulas[[1]][keep_this, , drop = FALSE]
        if (length(celulas[[1]][, 1]) == 0) {
          celulas <- try(celulas <- extract(ras, shapepol1, cellnumbers = TRUE,
                                            weights = TRUE, small = TRUE,
                                            normalizeWeights = FALSE),
                        silent = T)
          celulas <- extract(ras, shapepol1, cellnumbers = TRUE)
          nen <- sapply(celulas, nrow)
          for (ky in 1:length(nen)) {
              weigth <- rep(0, nen[ky])
              celulas[[ky]] <- cbind(celulas[[ky]], weigth)
          }
      }
    }
    if (class(celulas) == "try-error") {
        celulas <- extract(ras, shapepol1, cellnumbers = TRUE)
        nen <- sapply(celulas, nrow)
        for (ky in 1:length(nen)) {
            weigth <- rep(0, nen[ky])
            celulas[[ky]] <- cbind(celulas[[ky]], weigth)
        }
    }
    celulas <- celulas[!sapply(celulas, is.null)]
    if (length(celulas) > 0) {
        .rename <- function(x) {
            colnames(x) <- 1:3
            return(x)
        }
        celulas <- lapply(celulas, .rename)
    }
    pos <- which(nomes2[i] == nomes)
    pos2 <- do.call(rbind.data.frame, celulas)
    if (cover > 0) {
        prop <- round((pos2[, 3] * areashape[i])/areagrid[pos2[, 
            1]], 2)
        if (any(prop > 1)) {
            prop[prop > 1] <- 1
        }
        pos2 <- pos2[which(prop >= cover), , drop = FALSE]
    }
    return(list(pos = pos, pos2 = pos2))
}



for (i in 1:n) {
  print(i)
  pospos2 <- extractpos(ras, bird_ranges_spatial[i, ], nomes, nomes2, cover,
                        areashape, areagrid, i)
  matriz[pospos2$pos2[, 1], pospos2$pos] <- 1
}

Resultado <- cbind(coord, matriz)

if (remove.cells) Resultado <- letsR:::.removeCells(Resultado)

if (remove.sp) Resultado <- letsR:::.removeSp(Resultado)

values(ras) <- rowSums(matriz)
pam <- list(Presence_and_Absence_Matrix = Resultado,
              Richness_Raster = ras, 
              Species_name = colnames(Resultado)[-(1:2)])
class(pam) <- "PresenceAbsence"

saveRDS(pam, "data/aves/pam.rds")
