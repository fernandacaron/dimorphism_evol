rm(list = ls())

setwd("Documents/lab/dimorph_evol")

library(phytools)
library(geiger)
library(plotrix)
library(AICcmodavg)

########## MK + DR ##########

dat <- read.csv("data/BodySizeAves_30may22_edit.csv", row.names = 1)
tr <- read.nexus("data/aves_Ericson_VertLife_27JUL20.nex")
tr <- tr[1:100]

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

sdi[sdi > 0] <- "1"
sdi[sdi == 0] <- NA
sdi[sdi < 0] <- "-1"
sdi <- sdi[complete.cases(sdi)]
sdi <- as.factor(sdi)

fitMK <- function(phy, sdi, data, taxon) {
  subsdi <- sdi[names(sdi) %in% dat$Scientific_name[dat$Order == taxon]]
  N_M <- length(subsdi[subsdi == -1])
  N_F <- length(subsdi[subsdi == 1])
  
  aicc <- aiccw <- est1 <- est2 <- list()
  for (i in 1:length(phy)) {
    pruned_tr <- treedata(phy[[i]], subsdi, warnings = F)$phy
    subsdi <- subsdi[names(subsdi) %in% pruned_tr$tip.label]
    subsdi <- subsdi[match(pruned_tr$tip.label, names(subsdi))]  
    
    fitSYM <- fitMk(pruned_tr, subsdi, model = "SYM", pi = "fitzjohn")
    fitARD <- fitMk(pruned_tr, subsdi, model = "ARD", pi = "fitzjohn")
    modEst <- c(fitSYM$rates, fitARD$rates)
    names(modEst) <- c("SYM", "ARD1", "ARD2")
    modL <- c(fitSYM$logLik, fitARD$logLik)
    names(modL) <- c("SYM", "ARD")
    modK <- c(1, 2)
    outTab <- aictabCustom(logL = modL, K = modK, modnames = c("SYM", "ARD"),
                           nobs = length(subsdi))

    AICc <- outTab$AICc
    names(AICc) <- outTab$Modnames
    aicc[[i]] <- c(AICc[names(AICc) == "SYM"], AICc[names(AICc) == "ARD"])
    names(aicc[[i]]) <- c("AICcSYM", "AICcARD")
    AICcWt <- outTab$AICcWt
    names(AICcWt) <- outTab$Modnames
    aiccw[[i]] <- c(AICcWt[names(AICcWt) == "SYM"], AICcWt[names(AICcWt) == "ARD"])
    names(aiccw[[i]]) <- c("AICcWtSYM", "AICcWtARD")

    modsEst1 <- c(modEst[1], modEst[2])
    est1[[i]] <- modavgCustom(logL = modL, K = modK, 
                              modnames = c("SYM", "ARD"),
                              nobs = length(subsdi),
                              estimate = modsEst1,
                              se = c(NA, NA))
    modsEst2 <- c(modEst[1], modEst[3])
    est2[[i]] <- modavgCustom(logL = modL, K = modK, 
                              modnames = c("SYM", "ARD"), 
                              nobs = length(subsdi), 
                              estimate = modsEst2,
                              se = c(NA, NA))
  }
  
  aiccw_res <- matrix(nrow = 1, ncol = 4)
  colnames(aiccw_res) <- c("AICc_SYM", "AICc_ARD", "AICcW_SYM", "AICcW_ARD")
  rownames(aiccw_res) <- taxon

  stats1 <- function(x) {
    stat_sym <- stat_ard <- numeric()
    for (i in 1:length(phy)) {
      stat_sym[i] <- x[[i]][[1]]
      stat_ard[i] <- x[[i]][[2]]
    }
    me_sym <- round(mean(stat_sym), 3)
    me_ard <- round(mean(stat_ard), 3)
    range_sym <- range(stat_sym)
    range_ard <- range(stat_ard)
    return(c(paste0(me_sym, " (", round(range_sym[1], 3), "-", 
                    round(range_sym[2], 3), ")"),
             paste0(me_ard, " (", round(range_ard[1], 3), "-", 
                    round(range_ard[2], 3), ")")))
  }
  
  aiccw_res[1, 1:2] <- stats1(aicc)
  aiccw_res[1, 3:4] <- stats1(aiccw)
    
  stats2 <- function(x) {
    stat_est <- numeric()
    for (i in 1:length(phy)) {
      stat_est[i] <- x[[i]]$Mod.avg.est
    }
    me_est <- round(mean(stat_est), 3)
    range_est <- range(stat_est)
    return(c(paste0(me_est, " (", round(range_est[1], 3), "-", 
                    round(range_est[2], 3), ")")))
  }

  avgest_res <- matrix(nrow = 1, ncol = 2)
  colnames(avgest_res) <- c("rate_1to-1", "rate_-1to1")
  rownames(avgest_res) <- taxon

  avgest_res[1, 1] <- stats2(est1)
  avgest_res[1, 2] <- stats2(est2)

  res <- list()
  res$fitSYM1 <- fitSYM
  res$fitARD1 <- fitARD
  res$NM <- N_M
  res$NF <- N_F
  res$aiccw <- aiccw_res
  res$rates <- avgest_res

  return(res)
}

# Análises feitas no cluster da UFG

fitMK_col <- fitMK(tr, sdi = sdi, data = dat, taxon = "Columbiformes")
#save(fitMK_col, file = "data/aves/results/fitMK_col_AICcmmodavg.RData")

fitMK_psi <- fitMK(tr, sdi = sdi, data = dat, taxon = "Psittaciformes")
#save(fitMK_psi, file = "data/aves/results/fitMK_psi_AICcmmodavg.RData")

fitMK_ans <- fitMK(tr, sdi = sdi, data = dat, taxon = "Anseriformes")
#save(fitMK_ans, file = "data/aves/results/fitMK_ans_AICcmmodavg.RData")

fitMK_acc <- fitMK(tr, sdi = sdi, data = dat, taxon = "Accipitriformes")
#save(fitMK_acc, file = "data/aves/results/fitMK_acc_AICcmmodavg.RData")

fitMK_gal <- fitMK(tr, sdi = sdi, data = dat, taxon = "Galliformes")
#save(fitMK_gal, file = "data/aves/results/fitMK_gal_AICcmmodavg.RData")

fitMK_pic <- fitMK(tr, sdi = sdi, data = dat, taxon = "Piciformes")
#save(fitMK_pic, file = "data/aves/results/fitMK_pic_AICcmmodavg.RData")

fitMK_apo <- fitMK(tr, sdi = sdi, data = dat, taxon = "Apodiformes")
#save(fitMK_apo, file = "data/aves/results/fitMK_apo_AICcmmodavg.RData")

fitMK_cha <- fitMK(tr, sdi = sdi, data = dat, taxon = "Charadriiformes")
#save(fitMK_cha, file = "data/aves/results/fitMK_cha_AICcmmodavg.RData")

fitMK_pas <- fitMK(tr, sdi = sdi, data = dat, taxon = "Passeriformes")
#save(fitMK_pas, file = "data/aves/results/fitMK_pas_AICcmmodavg.RData")

#load(file = "data/aves/results/fitMK_col_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_psi_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_ans_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_acc_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_gal_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_pic_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_apo_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_cha_AICcmmodavg.RData")
#load(file = "data/aves/results/fitMK_pas_AICcmmodavg.RData")

aiccw <- rbind(fitMK_acc$aiccw,
               fitMK_ans$aiccw,
               fitMK_apo$aiccw,
               fitMK_cha$aiccw,
               fitMK_col$aiccw,
               fitMK_gal$aiccw,
               fitMK_pas$aiccw,
               fitMK_pic$aiccw,
               fitMK_psi$aiccw)

write.csv(aiccw, "tables/Table1_unformatted_AICcmmodavg.csv")

male <- "#9966FF"
monom <- "gray"
female <- "#E69F00"

pdf("figures/Figure4_AICcmmodavg.pdf", height = 9, width = 9)

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(0.5, 0.5, 0.5, 0.5))
cols2 <- setNames(c(male, female), c("Male-biased SSD", "Female-biased SSD"))

fitMK_acc$fitARD1$rates <- fitMK_acc$rates
plot_acc <- plot(fitMK_acc$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
mapply(draw.circle, plot_acc$x, plot_acc$y, col = cols2, border = NA, 
       MoreArgs = list(radius = 0.30))
text(plot_acc$x, plot_acc$y, c(paste("M\nN = ", fitMK_acc$NM, sep = ""),
                               paste("F\nN = ", fitMK_acc$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Accipitriformes", adj = 0, line = -1)

fitMK_ans$fitARD1$rates <- fitMK_ans$rates
plot_ans <- plot(fitMK_ans$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_ans$x, plot_ans$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_ans$x, plot_ans$y, c(paste("M\nN = ", fitMK_ans$NM, sep = ""),
                               paste("F\nN = ", fitMK_ans$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Anseriformes", adj = 0, line = -1)

fitMK_apo$fitARD1$rates <- fitMK_apo$rates
plot_apo <- plot(fitMK_apo$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_apo$x, plot_apo$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_apo$x, plot_apo$y, c(paste("M\nN = ", fitMK_apo$NM, sep = ""),
                               paste("F\nN = ", fitMK_apo$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Apodiformes", adj = 0, line = -1)

fitMK_cha$fitARD1$rates <- fitMK_cha$rates
plot_cha <- plot(fitMK_cha$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_cha$x, plot_cha$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_cha$x, plot_cha$y, c(paste("M\nN = ", fitMK_cha$NM, sep = ""),
                               paste("F\nN = ", fitMK_cha$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Charadriiformes", adj = 0, line = -1)

fitMK_col$fitARD1$rates <- fitMK_col$rates
plot_col <- plot(fitMK_col$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_col$x, plot_col$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_col$x, plot_col$y, c(paste("M\nN = ", fitMK_col$NM, sep = ""),
                               paste("F\nN = ", fitMK_col$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Columbiformes", adj = 0, line = -1)

fitMK_gal$fitARD1$rates <- fitMK_gal$rates
plot_gal <- plot(fitMK_gal$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_gal$x, plot_gal$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_gal$x, plot_gal$y, c(paste("M\nN = ", fitMK_gal$NM, sep = ""),
                               paste("F\nN = ", fitMK_gal$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Galliformes", adj = 0, line = -1)

fitMK_pas$fitARD1$rates <- fitMK_pas$rates
plot_pas <- plot(fitMK_pas$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_pas$x, plot_pas$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_pas$x, plot_pas$y, c(paste("M\nN = ", fitMK_pas$NM, sep = ""),
                               paste("F\nN = ", fitMK_pas$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Passeriformes", adj = 0, line = -1)

fitMK_pic$fitARD1$rates <- fitMK_pic$rates
plot_pic <- plot(fitMK_pic$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_pic$x, plot_pic$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_pic$x, plot_pic$y, c(paste("M\nN = ", fitMK_pic$NM, sep = ""),
                               paste("F\nN = ", fitMK_pic$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Piciformes", adj = 0, line = -1)

fitMK_psi$fitARD1$rates <- fitMK_psi$rates
plot_psi <- plot(fitMK_psi$fitARD1, show.zeros = FALSE, spacer = 0.4, 
                 mar = rep(0.5, 4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_psi$x, plot_psi$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.3)))
text(plot_psi$x, plot_psi$y, c(paste("M\nN = ", fitMK_psi$NM, sep = ""),
                               paste("F\nN = ", fitMK_psi$NF, sep = "")), 
     cex = 1.3, col = "white", font = 2)
title("Psittaciformes", adj = 0, line = -1)

dev.off()

## Figure 5 - Regressão DR contra grau de SSD

source("https://raw.githubusercontent.com/mgharvey/ES-sim/master/R/essim.R")

sdi <- numeric()
for (i in 1:nrow(dat_red)) {
  sdi[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
                female = dat_red$Body_mass_g_F_mean[i],
                cutoff = FALSE)
}
names(sdi) <- dat_red$Scientific_name
sdi <- sdi[complete.cases(sdi)]

# complete tree: calculting es statistic 
es <- subset_es <- list()
for(i in 1:length(tr)) {
  es[[i]] <- compute_es(tr[[i]])
  subset_es[[i]] <- es[[i]][names(es[[i]]) %in% names(sdi)]
  subset_es[[i]] <- subset_es[[i]][match(names(sdi), names(subset_es[[i]]))]
}

essim_ave <- lapply(tr[1:100], essim, trait = sdi, nsim = 100)
save(essim_ave, file = "data/aves/results/essim_ave.RData")

male_al <- rgb(153/255, 102/255, 255/255, 0.3)
monom_al <- rgb(190/255, 190/255, 190/255, 0.3)
female_al <- rgb(230/255, 159/255, 0/255, 0.3)

cols <- rep(NA, length(sdi))
names(cols) <- sdi
for (i in 1:length(sdi)) {
  cols[i] <- ifelse(sdi[i] > 0, female_al, ifelse(sdi[i] < 0, male_al,
                                                  monom_al))
}

pdf("figures/Figure5.pdf")

plot(sdi ~ subset_es[[1]], xlab = expression(lambda["DR"]), ylab = "SDI", main = "",
     pch = 16, col = cols)
for (i in 2:100) {
  points(sdi ~ subset_es[[i]], pch = 16, col = cols)
}

legend("bottomright", pch = 16, bty = 'n', 
       col = c(male, monom, female), 
       legend  = c("Male-biased SSD", "Monomorphism", "Female-biased SSD"))

dev.off()
