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

## 1. Exploração geral do SSD, auto-detectando shifts 
pruned_tr1 <- treedata(tr, sdi, warnings = F)$phy
subsdi <- treedata(tr, sdi, warnings = F)$data
subsdi <- subsdi[match(pruned_tr1$tip.label, names(subsdi))]
RR_sdi <- RRphylo(tree = pruned_tr1, y = subsdi)

pruned_tr2 <- treedata(tr, mass, warnings = F)$phy
submass <- treedata(tr, mass, warnings = F)$data
submass <- submass[match(pruned_tr2$tip.label, names(subdat))]
RR_mass <- RRphylo(tree = pruned_tr2, y = submass)

cov <- c(RR_mass$ace, submass)

search.shift(RR, auto.recognize = "yes", covariate = "TRUE", cov = cov)

## 2. Usando SSD como state
sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- 1
sdi_disc[sdi_disc == 0] <- 0
sdi_disc[sdi_disc < 0] <- -1
sdi_disc <- sdi_disc[match(pruned_tr1$tip.label, names(sdi_disc))]

search.shift(RR, status.type = "sparse", state = sdi_disc, 
             auto.recognize = "yes", covariate = "TRUE", cov = cov)


