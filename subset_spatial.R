rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(stringi)
library(sf)
library(raster)
library(dplyr)
library(fasterize)
library(rgeos)

dat <- read.csv("data/aves/BodySizeAves_30may22_edit.csv", row.names = 1)

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

birds <- st_read(dsn = "data/aves/spatial/BOTW/BOTW.gdb", layer = "All_Species")
birds$sci_name <- stri_replace_all_fixed(birds$sci_name, " ", "_")
birds2 <- birds %>% filter(st_geometry_type(Shape) != "MULTISURFACE")

acc <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Accipitriformes"])
ans <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Anseriformes"])
apo <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Apodiformes"])
cha <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Charadriiformes"])
col <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Columbiformes"])
gal <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Galliformes"])
pas <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Passeriformes"])
pic <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Piciformes"])
psi <- birds2 %>% filter(sci_name %in% dat_red$Scientific_name[dat_red$Order ==
                            "Psittaciformes"])

save(acc, file = "data/aves/spatial/acc.RData")
save(ans, file = "data/aves/spatial/ans.RData")
save(apo, file = "data/aves/spatial/apo.RData")
save(cha, file = "data/aves/spatial/cha.RData")
save(col, file = "data/aves/spatial/col.RData")
save(gal, file = "data/aves/spatial/gal.RData")
save(pas, file = "data/aves/spatial/pas.RData")
save(pic, file = "data/aves/spatial/pic.RData")
save(psi, file = "data/aves/spatial/psi.RData")

