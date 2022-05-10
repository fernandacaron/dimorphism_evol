rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(RRphylo)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)
tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")[[1]]

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

mass <- numeric()
for (i in 1:nrow(dat_red)) {
  mass[i] <- mean(c(dat_red$Body_mass_g_M_mean[i],
                    dat_red$Body_mass_g_F_mean[i]))
}
names(mass) <- dat_red$Scientific_name 
mass <- log(mass)

mass_male <- log(dat_red$Body_mass_g_M_mean)
names(mass_male) <- dat_red$Scientific_name

mass_female <- log(dat_red$Body_mass_g_F_mean)
names(mass_female) <- dat_red$Scientific_name

## 1. Exploração geral do SSD, auto-detectando shifts 
a <- Sys.time()

tr_sdi <- treedata(tr, sdi, warnings = F)$phy
subsdi <- sdi[names(sdi) %in% tr_sdi$tip.label]
subsdi <- subsdi[match(tr_sdi$tip.label, names(subsdi))]
subsdi <- as.matrix(subsdi)
RR_sdi <- RRphylo(tree = tr_sdi, y = subsdi)

b <- Sys.time()

tr_mass <- treedata(tr, mass, warnings = F)$phy
submass <- mass[names(mass) %in% tr_mass$tip.label]
submass <- submass[match(tr_mass$tip.label, names(submass))]
RR_mass <- RRphylo(tree = tr_mass, y = submass)

cov <- c(RR_mass$ace, submass)

search.shift(RR_sdi, status.type = "clade", cov = cov, 
             filename = paste(tempdir(), "RRphylo1", sep = "/"))



## 2. Usando SSD como state
sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- 1
sdi_disc[sdi_disc == 0] <- 0
sdi_disc[sdi_disc < 0] <- -1
sdi_disc <- sdi_disc[match(pruned_tr1$tip.label, names(sdi_disc))]

tr_male <- treedata(tr, mass_male, warnings = F)$phy
submale <- mass_male[names(mass_male) %in% tr_male$tip.label]
submale <- submale[match(tr_male$tip.label, names(submale))]
RR_male <- RRphylo(tree = tr_male, y = submale)

tr_female <- treedata(tr, mass_female, warnings = F)$phy
subfemale <- mass_female[names(mass_female) %in% tr_female$tip.label]
subfemale <- subfemale[match(tr_female$tip.label, names(subfemale))]
RR_female <- RRphylo(tree = tr_female, y = subfemale)

search.shift(RR_male, status.type = "sparse", state = sdi_disc, cov = submale,
             filename = paste(tempdir(), "rrphylo_male", sep = "/"))

search.shift(RR_female, status.type = "sparse", state = sdi_disc, 
             cov = subfemale, filename = paste(tempdir(), "rrphylo_female", 
                                               sep = "/"))

