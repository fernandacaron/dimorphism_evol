rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(RRphylo)

dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)
tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")

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
search.shift.birds1 <- function(phy, taxon) {
  
  if (!is.null(taxon)) {
    subsdi <- sdi[names(sdi) %in% dat$Scientific_name[dat$Order == taxon]]
  } else subsdi <- sdi

  pruned_tr <- treedata(phy, subsdi, warnings = F)$phy
  subsdi <- subsdi[names(subsdi) %in% pruned_tr$tip.label]
  submass <- mass[names(mass) %in% names(subsdi)]

  subsdi <- subsdi[match(pruned_tr$tip.label, names(subsdi))]
  submass <- submass[match(pruned_tr$tip.label, names(submass))]
  
  RR_sdi <- RRphylo(tree = pruned_tr, y = subsdi)
  
  shift <- search.shift(RR = RR_sdi, status.type = "clade", cov = submass,
                        filename = paste(tempdir(), "RRphylo1", taxon,
                                         sep = "/"))
  
  res <- list()
  res$shift <- shift
  res$phy <- pruned_tr
  
  res
}

plot.shift <- function(RR_shift, ncol = 1) {
  layout(matrix(1:length(RR_shift), ncol = 5, byrow = TRUE))

  par(mar = c(0, 0, 0, 0))

  shift_dec <- shift_inc <- list()
  for (i in 1:length(RR_shift)){
    shift_dec[[i]] <- as.numeric(rownames(RR_shift[[i]]$shift$all.clades[
                                          RR_shift[[i]]$shift$all.clades[, 2] <=
                                          0.025, ]))
    names(shift_dec[[i]]) <- round(RR_shift[[i]]$shift$all.clades[
                                   RR_shift[[i]]$shift$all.clades[, 2] <=
                                   0.025, 1], 3)

    shift_inc[[i]] <- as.numeric(rownames(RR_shift[[i]]$shift$all.clades[
                                          RR_shift[[i]]$shift$all.clades[, 2] >=
                                          0.975, ]))
    names(shift_inc[[i]]) <- round(RR_shift[[i]]$shift$all.clades[
                                   RR_shift[[i]]$shift$all.clades[, 2] >=
                                   0.975, 1], 3)
  }

  for (i in 1:length(RR_shift)) {
    plot(RR_shift[[i]]$phy)

    if (length(shift_dec[[i]]) > 0) {
      for (j in 1:length(shift_dec[[i]])) nodelabels(node = shift_dec[[i]][j],
                                                     pch = 16, cex = 2,
                                                     col = "blue")
      for (j in 1:length(shift_dec[[i]])) nodelabels(text = 
                                                     names(shift_dec[[i]][j]),
                                                     node = shift_dec[[i]][j],
                                                     col = "black",
                                                     bg = "white",
                                                     cex = 0.7,
                                                     adj = c(0.2, -1.5))
    }

    if (length(shift_inc[[i]]) > 0) {
        for (j in 1:length(shift_inc[[i]])) nodelabels(node = shift_inc[[i]][j],
                                                       pch = 16, cex = 2,
                                                       col = "red")
        for (j in 1:length(shift_inc[[i]])) nodelabels(text = 
                                                        names(shift_inc[[i]][j]),
                                                        node = shift_inc[[i]][j],
                                                        col = "black",
                                                        bg = "white",
                                                        cex = 0.7, 
                                                        adj = c(0.2, -1.5))
      
    }
  }
}
# $all.clades: for each detected node, the data-frame includes the average 
# rate difference (computed as the mean rate over all branches subtended by the 
# node minus the average rate for the rest of the tree) and the probability
# that it do represent a real shift. Probabilities are contrasted to simulations
# shuffling the rates across the tree branches for a number of replicates
# specified by the argument nrep. Note that the p-values refer to the number of 
# times the real average rates are larger (or smaller) than the rates averaged
# over the rest of the tree, divided by the number of simulations. Hence, large
# rates are significantly larger than the rest of the tree (at alpha = 0.05), 
# when the probability is > 0.975; and small rates are significantly small 
# for p < 0.025.
# $single.clades the same as with ’all.clades’ but restricted to the largest/
# smallest rate values along a single clade (i.e. nested clades with smaller 
# rate shifts are excluded). Large rates are significantly larger than the rest 
# of the tree (at alpha = 0.05), when the probability is > 0.975; and small 
# rates are significantly small for p < 0.025.

RR1_col <- lapply(tr[1:10], search.shift.birds1, "Columbiformes")
RR1_psi <- lapply(tr[1:10], search.shift.birds1, "Psittaciformes")
RR1_ans <- lapply(tr[1:10], search.shift.birds1, "Anseriformes")
RR1_acc <- lapply(tr[1:10], search.shift.birds1, "Accipitriformes")
RR1_gal <- lapply(tr[1:10], search.shift.birds1, "Galliformes")
RR1_pic <- lapply(tr[1:10], search.shift.birds1, "Piciformes")
RR1_apo <- lapply(tr[1:10], search.shift.birds1, "Apodiformes")
RR1_cha <- lapply(tr[1:10], search.shift.birds1, "Charadriiformes")
RR1_pas <- lapply(tr[1:10], search.shift.birds1, "Passeriformes")
RR1_ave <- lapply(tr[1:10], search.shift.birds1)

pdf("figures/RR1.pdf", width = 14, height = 10)
plot.shift(RR1_col, ncol = 5)
plot.shift(RR1_psi, ncol = 5)
plot.shift(RR1_ans, ncol = 5)
plot.shift(RR1_acc, ncol = 5)
plot.shift(RR1_gal, ncol = 5)
plot.shift(RR1_pic, ncol = 5)
plot.shift(RR1_apo, ncol = 5)
plot.shift(RR1_cha, ncol = 5)
plot.shift(RR1_pas, ncol = 5)
dev.off()


## 2. Usando SSD como state
sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- "Female-biased SSD"
sdi_disc[sdi_disc == 0] <- "Monomorphism"
sdi_disc[sdi_disc < 0] <- "Male-biased SSD"
sdi_disc <- as.factor(sdi_disc)

# Acho que tem que mudar o 0 do subsdi (ou não)
search.shift.birds2 <- function(phy, taxon = NULL) {
  
  if (!is.null(taxon)) {
    subsdi_disc <- sdi_disc[names(sdi_disc) %in% dat$Scientific_name[dat$Order 
                                                                      == taxon]]
    sub_mal <- mass_male[names(mass_male) %in% dat$Scientific_name[dat$Order 
                                                                      == taxon]]
    sub_fem <- mass_female[names(mass_female) %in% dat$Scientific_name[dat$Order
                                                                      == taxon]]

  } else {
    subsdi_disc <- sdi_disc
    sub_mal <- mass_male
    sub_fem <- mass_female
  }

  pruned_tr <- treedata(phy, subsdi_disc, warnings = F)$phy
  subsdi_disc <- subsdi_disc[names(subsdi_disc) %in% pruned_tr$tip.label]
  sub_mal <- sub_mal[names(sub_mal) %in% names(subsdi_disc)]
  sub_fem <- sub_fem[names(sub_fem) %in% names(subsdi_disc)]

  subsdi_disc <- subsdi_disc[match(pruned_tr$tip.label, names(subsdi_disc))]  
  sub_mal <- sub_mal[match(pruned_tr$tip.label, names(sub_mal))]
  sub_fem <- sub_fem[match(pruned_tr$tip.label, names(sub_fem))]
  
  RR_male <- RRphylo(tree = pruned_tr, y = sub_mal)

  RR_female <- RRphylo(tree = pruned_tr, y = sub_fem)

  shift_male <- search.shift(RR_male, status.type = "sparse",
                             state = subsdi_disc, cov = sub_mal,
                             filename = paste(tempdir(), 
                                              paste0("RRphylo2", "_", taxon, 
                                                     "_", "M"),
                                              sep = "/"))

  shift_female <- search.shift(RR_female, status.type = "sparse",
                               state = subsdi_disc, cov = sub_fem,
                               filename = paste(tempdir(), 
                                              paste0("RRphylo2", "_", taxon, 
                                                     "_", "F"),
                                              sep = "/"))

  res <- list()
  res$shift_male <- shift_male
  res$shift_female <- shift_female
  res$phy <- pruned_tr
  res
}

RR2_col <- lapply(tr[1:100], search.shift.birds2, "Columbiformes")
save(RR2_col, file = "data/aves/RR2_col.RData")

RR2_psi <- lapply(tr[1:100], search.shift.birds2, "Psittaciformes")
save(RR2_psi, file = "data/aves/RR2_psi.RData")

RR2_ans <- lapply(tr[1:100], search.shift.birds2, "Anseriformes")
save(RR2_ans, file = "data/aves/RR2_ans.RData")

RR2_acc <- lapply(tr[1:100], search.shift.birds2, "Accipitriformes")
save(RR2_acc, file = "data/aves/RR2_acc.RData")

RR2_gal <- lapply(tr[1:100], search.shift.birds2, "Galliformes")
save(RR2_gal, file = "data/aves/RR2_gal.RData")

RR2_pic <- lapply(tr[1:100], search.shift.birds2, "Piciformes")
save(RR2_pic, file = "data/aves/RR2_pic.RData")

RR2_apo <- lapply(tr[1:100], search.shift.birds2, "Apodiformes")
save(RR2_apo, file = "data/aves/RR2_apo.RData")

RR2_cha <- lapply(tr[1:100], search.shift.birds2, "Charadriiformes")
save(RR2_cha, file = "data/aves/RR2_cha.RData")

RR2_pas <- lapply(tr[1:100], search.shift.birds2, "Passeriformes")
save(RR2_pas, file = "data/aves/RR2_pas.RData")

RR2_ave <- lapply(tr[1:100], search.shift.birds2)
save(RR2_ave, file = "data/aves/RR2_ave.RData")
