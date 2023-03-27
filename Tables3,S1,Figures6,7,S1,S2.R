rm(list = ls())

setwd("~/Documents/lab/dimorph_evol")

library(mvMORPH)
library(phytools)
library(geiger)
library(RRphylo)
library(viridis)

## mvBM + RRphylo

tr <- read.nexus("~/Documents/lab/data/trees/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/BodySizeAves_30may22_edit.csv", row.names = 1)

subdat <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
               complete.cases(dat$Body_mass_g_F_mean), ]
rownames(subdat) <- dat$Scientific_name[complete.cases(dat$Body_mass_g_M_mean) &
										complete.cases(dat$Body_mass_g_F_mean)]

subdat$log_Body_mass_g_M_mean <- log(subdat$Body_mass_g_M_mean)
subdat$log_Body_mass_g_F_mean <- log(subdat$Body_mass_g_F_mean)

## mvBM

multiModelComp <- function(phy, data, taxon = NULL, ntree = 1) {

	if (!is.null(taxon)) {
		subdat <- subset(data, Order == taxon, c(log_Body_mass_g_M_mean,
		                                         log_Body_mass_g_F_mean))

		rownames(subdat) <- subset(data, Order == taxon, Scientific_name)[, 1]

	} else {
		subdat <- data[, c("log_Body_mass_g_M_mean", "log_Body_mass_g_F_mean")]

		rownames(subdat) <- data$Scientific_name

	}

	fit <- list()
	pval <- numeric()
	for (i in 1:length(phy)) {

		try <- tryCatch({

			xx <- treedata(phy[[i]], subdat, warnings = FALSE)
			tree <- xx$phy
			data <- xx$data
			data <- data[tree$tip.label, ]

			# BMM multi-rate/multi-selective regimes, BM1 unique rate of evolution/trait
			fit1 <- mvBM(tree, data, model = "BM1", 
			             optimization = "L-BFGS-B")

			if (fit1$convergence > 0) {
				fit1 <- mvBM(tree, data, model = "BM1", 
				             optimization = "Nelder-Mead")
			}

			if (fit1$convergence > 0) {
				fit1 <- mvBM(tree, data, model = "BM1", 
				             optimization = "subplex")
			}

			if (fit1$convergence > 0) {
				stop("convergence of the optimizer has not been reached, try simpler model")
			}

			# Constraint faz que com sigma fique igual para os dois caracteres
			fit2 <- mvBM(tree, data, model = "BM1", 
			             param = list(constraint = TRUE),
			             optimization = c("L-BFGS-B"))

			if (fit2$convergence > 0) {
				fit2 <- mvBM(tree, data, model = "BM1", 
				             param = list(constraint = TRUE), 
				             optimization = "Nelder-Mead")
			}

			if (fit2$convergence > 0) {
				fit2 <- mvBM(tree, data, model = "BM1", 
				             param = list(constraint = TRUE), 
				             optimization = "subplex")
			}

			if (fit2$convergence > 0) {
				stop("convergence of the optimizer has not been reached, try simpler model")
			}

			res <- list()
			res$res <- matrix(nrow = 2, ncol = 6)
			colnames(res$res) <- c("LogLik", "AICc", "sigma_M", "sigma_F", 
			                       "nparam", "convergence")
			rownames(res$res) <- c("fit1", "fit2_constraint")
			res$res[1, ] <- c(fit1$LogLik, fit1$AICc, fit1$sigma[1, 1], 
			                  fit1$sigma[2, 2], fit1$param$nparam, 
			                  fit1$convergence)
			res$res[2, ] <- c(fit2$LogLik, fit2$AICc, fit2$sigma[1, 1], 
			                  fit2$sigma[2, 2], fit2$param$nparam, 
			                  fit2$convergence)
			res$fit1 <- fit1
			res$fit2 <- fit2
			res

		}, error = function(e) e)

		if(inherits(try, "error")) next

		else {
			
			fit[[length(fit) + 1]] <- try$res

			pval[length(pval) + 1] <- LRT(try$fit1, try$fit2, 
			                              echo = FALSE)$pval

			if (length(fit) == ntree) break

		}
	}

	res <- list()
	res$fit <- fit
	res$pval <- pval
	res
}

# optimization = "L-BFGS-B"    # 3.852359 hours
# optimization = "Nelder-Mead" # 3.257686 hours
# optimization = "subplex"     # 12.4687 hours

# Análises feitas no cluster da UFG

## Columbiformes = 105 spp
fit_col <- multiModelComp(phy = tr, data = subdat, taxon = "Columbiformes",
                          ntree = 100)
#save(fit_col, file = "data/results/fitBM_col.RData")

## Psittaciformes = 131 spp
fit_psi <- multiModelComp(phy = tr, data = subdat, taxon = "Psittaciformes",
                          ntree = 100)
#save(fit_psi, file = "data/results/fitBM_psi.RData")

## Anseriformes = 150 spp
fit_ans <- multiModelComp(phy = tr, data = subdat, taxon = "Anseriformes",
                          ntree = 100)
#save(fit_ans, file = "data/results/fitBM_ans.RData")

## Accipitriformes = 169 spp
fit_acc <- multiModelComp(phy = tr, data = subdat, taxon = "Accipitriformes",
                          ntree = 100)
#save(fit_acc, file = "data/results/fitBM_acc.RData")

## Galliformes = 179 spp
fit_gal <- multiModelComp(phy = tr, data = subdat, taxon = "Galliformes",
                          ntree = 100)
#save(fit_gal, file = "data/results/fitBM_gal.RData")

## Piciformes = 217 spp
fit_pic <- multiModelComp(phy = tr, data = subdat, taxon = "Piciformes",
                          ntree = 100)
#save(fit_pic, file = "data/results/fitBM_pic.RData")

## Apodiformes = 247 spp
fit_apo <- multiModelComp(phy = tr, data = subdat, taxon = "Apodiformes",
                          ntree = 100)
#save(fit_apo, file = "data/results/fitBM_apo.RData")

## charadriiformes = 252 spp
fit_cha <- multiModelComp(phy = tr, data = subdat, taxon = "charadriiformes",
                          ntree = 100)
#save(fit_cha, file = "data/results/fitBM_cha.RData")

## Passeriformes = 2290 spp
fit_pas <- multiModelComp(phy = tr, data = subdat, taxon = "Passeriformes",
                          ntree = 100)
#save(fit_pas, file = "data/results/fitBM_pas.RData")

#load(file = "data/results/fitBM_col.RData")
#load(file = "data/results/fitBM_psi.RData")
#load(file = "data/results/fitBM_ans.RData")
#load(file = "data/results/fitBM_acc.RData")
#load(file = "data/results/fitBM_gal.RData")
#load(file = "data/results/fitBM_pic.RData")
#load(file = "data/results/fitBM_apo.RData")
#load(file = "data/results/fitBM_cha.RData")
#load(file = "data/results/fitBM_pas.RData")

## Table 2

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
for (i in 1:nrow(subdat)) {
	sdi[i] <- SDI(male = subdat$Body_mass_g_M_mean[i],
	              female = subdat$Body_mass_g_F_mean[i],
	              cutoff = FALSE)
}
subdat_sdi <- cbind(subdat, sdi)

## Sigmas do fit 1 (sigmas diferentes)
sigmaM_ave <- sigmaM_acc <- sigmaM_ans <- sigmaM_apo <- sigmaM_cha <- 
  sigmaM_col <- sigmaM_gal <- sigmaM_pas <- sigmaM_pic <- sigmaM_psi <- 
  numeric()
sigmaF_ave <- sigmaF_acc <- sigmaF_ans <- sigmaF_apo <- sigmaF_cha <- 
  sigmaF_col <- sigmaF_gal <- sigmaF_pas <- sigmaF_pic <- sigmaF_psi <- 
  numeric()
for (i in 1:length(fit_acc)) {
  sigmaM_acc[i] <- fit_acc$fit[[i]][1, 3]
	sigmaM_ans[i] <- fit_ans$fit[[i]][1, 3]
	sigmaM_apo[i] <- fit_apo$fit[[i]][1, 3]
	sigmaM_cha[i] <- fit_cha$fit[[i]][1, 3]
	sigmaM_col[i] <- fit_col$fit[[i]][1, 3]
	sigmaM_gal[i] <- fit_gal$fit[[i]][1, 3]
	sigmaM_pas[i] <- fit_pas$fit[[i]][1, 3]
	sigmaM_pic[i] <- fit_pic$fit[[i]][1, 3]
	sigmaM_psi[i] <- fit_psi$fit[[i]][1, 3]

	sigmaF_acc[i] <- fit_acc$fit[[i]][1, 4]
	sigmaF_ans[i] <- fit_ans$fit[[i]][1, 4]
	sigmaF_apo[i] <- fit_apo$fit[[i]][1, 4]
	sigmaF_cha[i] <- fit_cha$fit[[i]][1, 4]
	sigmaF_col[i] <- fit_col$fit[[i]][1, 4]
	sigmaF_gal[i] <- fit_gal$fit[[i]][1, 4]
	sigmaF_pas[i] <- fit_pas$fit[[i]][1, 4]
	sigmaF_pic[i] <- fit_pic$fit[[i]][1, 4]
	sigmaF_psi[i] <- fit_psi$fit[[i]][1, 4]
}

fit_birds <- data.frame(nspp = c(167, 150, 247, 255, 103, 196, 2510, 212, 130),
                        SDI = rep(NA, 9),
                        sigmaM = rep(NA, 9), 
                        sigmaF = rep(NA, 9), 
                        pval = rep(NA, 9),
                        row.names = c("Accipitriformes", "Anseriformes", 
                                      "Apodiformes", "Charadriiformes", 
                                      "Columbiformes", "Galliformes", 
                                      "Passeriformes", "Piciformes", 
                                      "Psittaciformes"))

summ <- function(x, sdi = FALSE) {
  if (sdi) {
    return(paste0(
           "Median = ",
           round(median(x), 3),
           " (Min = ",
           round(min(x), 3),
           "; Max = ",
           round(max(x), 3), 
           ")"))
  } else {
    return(paste0(
      round(median(x), 3),
      " (",
      round(min(x), 3),
      "-",
      round(max(x), 3), 
      ")"))
  }
}

fit_birds[1, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order ==
                                             "Accipitriformes"], TRUE),
                       summ(sigmaM_acc, FALSE),
                       summ(sigmaF_acc, FALSE),
                       summ(fit_acc$pval, FALSE)
)

fit_birds[2, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Anseriformes"], 
                            TRUE),
                       summ(sigmaM_ans, FALSE),
                       summ(sigmaF_ans, FALSE),
                       summ(fit_ans$pval, FALSE)
)

fit_birds[3, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Apodiformes"], 
                            TRUE),
                       summ(sigmaM_apo, FALSE),
                       summ(sigmaF_apo, FALSE),
                       summ(fit_apo$pval, FALSE)
)

fit_birds[4, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == 
                                             "Charadriiformes"], TRUE),
                       summ(sigmaM_cha, FALSE),
                       summ(sigmaF_cha, FALSE),
                       summ(fit_cha$pval, FALSE)
)

fit_birds[5, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Columbiformes"], 
                            TRUE),
                       summ(sigmaM_col, FALSE),
                       summ(sigmaF_col, FALSE),
                       summ(fit_col$pval, FALSE)
)

fit_birds[6, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Galliformes"], 
                            TRUE),
                       summ(sigmaM_gal, FALSE),
                       summ(sigmaF_gal, FALSE),
                       summ(fit_gal$pval, FALSE)
)

fit_birds[7, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Passeriformes"], 
                            TRUE),
                       summ(sigmaM_pas, FALSE),
                       summ(sigmaF_pas, FALSE),
                       summ(fit_pas$pval, FALSE)
)

fit_birds[8, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Piciformes"], 
                            TRUE),
                       summ(sigmaM_pic, FALSE),
                       summ(sigmaF_pic, FALSE),
                       summ(fit_pic$pval, FALSE)
)

fit_birds[9, 2:5] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == 
                                              "Psittaciformes"], TRUE),
                       summ(sigmaM_psi, FALSE),
                       summ(sigmaF_psi, FALSE),
                       summ(fit_psi$pval, FALSE)
)

write.csv(fit_birds, "tables/Table3_unformatted.csv")

## RRPhylo

dat_red <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
                 complete.cases(dat$Body_mass_g_F_mean), ]

sdi <- numeric()
for (i in 1:nrow(dat_red)) {
  sdi[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
                female = dat_red$Body_mass_g_F_mean[i],
                cutoff = FALSE)
}
names(sdi) <- dat_red$Scientific_name 

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

## 1. General exploration of SSD evolution (auto-detecting shifts)

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

# Análises feitas no cluster da UFG
RR1_col <- lapply(tr[1:100], search.shift.birds1, "Columbiformes")
#save(RR1_col, file = "data/results/RR1_col.RData")

RR1_psi <- lapply(tr[1:100], search.shift.birds1, "Psittaciformes")
#save(RR1_psi, file = "data/results/RR1_psi.RData")

RR1_ans <- lapply(tr[1:100], search.shift.birds1, "Anseriformes")
#save(RR1_ans, file = "data/results/RR1_ans.RData")

RR1_acc <- lapply(tr[1:100], search.shift.birds1, "Accipitriformes")
#save(RR1_acc, file = "data/results/RR1_acc.RData")

RR1_gal <- lapply(tr[1:100], search.shift.birds1, "Galliformes")
#save(RR1_gal, file = "data/results/RR1_gal.RData")

RR1_pic <- lapply(tr[1:100], search.shift.birds1, "Piciformes")
#save(RR1_pic, file = "data/results/RR1_pic.RData")

RR1_apo <- lapply(tr[1:100], search.shift.birds1, "Apodiformes")
#save(RR1_apo, file = "data/results/RR1_apo.RData")

RR1_cha <- lapply(tr[1:100], search.shift.birds1, "Charadriiformes")
#save(RR1_cha, file = "data/results/RR1_cha.RData")

RR1_pas <- lapply(tr[1:100], search.shift.birds1, "Passeriformes")
#save(RR1_pas, file = "data/results/RR1_pas.RData")

#load(file = "data/results/RR1_col.RData")
#load(file = "data/results/RR1_psi.RData")
#load(file = "data/results/RR1_ans.RData")
#load(file = "data/results/RR1_acc.RData")
#load(file = "data/results/RR1_gal.RData")
#load(file = "data/results/RR1_pic.RData")
#load(file = "data/results/RR1_apo.RData")
#load(file = "data/results/RR1_cha.RData")
#load(file = "data/results/RR1_pas.RData")

sig_col <- lapply(RR1_col, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                                       0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                                       0.975, ])
sig_psi <- lapply(RR1_psi, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_ans <- lapply(RR1_ans, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_acc <- lapply(RR1_acc, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_gal <- lapply(RR1_gal, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_pic <- lapply(RR1_pic, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_apo <- lapply(RR1_apo, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_cha <- lapply(RR1_cha, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

sig_pas <- lapply(RR1_pas, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
                                                            0.025 | 
                                                            x$shift$all.clades$p.value >= 
                                                            0.975, ])

data.boxplot <- function(data) {
  pos <- lapply(data, function(x) x$rate.difference[x$rate.difference > 0])
  neg <- lapply(data, function(x) x$rate.difference[x$rate.difference < 0])
  mat_pos <- lapply(pos, function(x) c("Positive", length(x)))
  mat_neg <- lapply(neg, function(x) c("Negative", length(x)))
  mat_pos <- do.call(rbind, mat_pos)
  mat_neg <- do.call(rbind, mat_neg)
  mat <- as.data.frame(rbind(mat_pos, mat_neg), row.names = NA)
  mat[, 2] <- as.numeric(mat[, 2])
  colnames(mat) <- c("State", "Quantity")
  mat
}

dat_box_col <- data.boxplot(sig_col)
dat_box_psi <- data.boxplot(sig_psi)
dat_box_ans <- data.boxplot(sig_ans)
dat_box_acc <- data.boxplot(sig_acc)
dat_box_gal <- data.boxplot(sig_gal)
dat_box_pic <- data.boxplot(sig_pic)
dat_box_apo <- data.boxplot(sig_apo)
dat_box_cha <- data.boxplot(sig_cha)
dat_box_pas <- data.boxplot(sig_pas)

dat_hist_col <- do.call(rbind, sig_col)
dat_hist_psi <- do.call(rbind, sig_psi)
dat_hist_ans <- do.call(rbind, sig_ans)
dat_hist_acc <- do.call(rbind, sig_acc)
dat_hist_gal <- do.call(rbind, sig_gal)
dat_hist_pic <- do.call(rbind, sig_pic)
dat_hist_apo <- do.call(rbind, sig_apo)
dat_hist_cha <- do.call(rbind, sig_cha)
dat_hist_pas <- do.call(rbind, sig_pas)

# Figure 6

pdf("figures/Figure6.pdf", width = 9, height = 9)

layout(matrix(1:9, ncol = 3, byrow = TRUE))

par(mar = c(4, 4, 4, 2))

col_alp <- c(rgb(63/255, 57/255, 148/255, 0.7),
             rgb(165/255, 19/255, 1/255, 0.7))

col <- c(turbo(20)[2], turbo(20)[19])

boxplot(Quantity ~ State, data = dat_box_acc, ylab = "Number of shifts/tree", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Accipitriformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_ans, ylab = "", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Anseriformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_apo, ylab = "", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Apodiformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_cha, ylab = "Number of shifts/tree", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Charadriiformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_col, ylab = "", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Columbiformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_gal, ylab = "", 
        xlab = "", names = NA, col = col_alp, border = col)
title("Galliformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_pas, ylab = "Number of shifts/tree", 
        xlab = "Rate difference", col = col_alp, border = col)
title("Passeriformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_pic, ylab = "", 
        xlab = "Rate difference", col = col_alp, border = col)
title("Piciformes", adj = 0)

boxplot(Quantity ~ State, data = dat_box_psi, ylab = "", 
        xlab = "Rate difference", col = col_alp, border = col)
title("Psittaciformes", adj = 0)

#boxplot(Quantity ~ State, data = dat_box_ave)

dev.off()

# Figure 7

pdf("figures/Figure7.pdf", width = 9, height = 9)

layout(matrix(1:9, nrow = 3, byrow = TRUE))

par(mar = c(4, 4, 4, 2))

hist(dat_hist_acc$rate.difference[dat_hist_acc$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", 
     xlab = "", breaks = 20)
hist(abs(dat_hist_acc$rate.difference[dat_hist_acc$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 15)
title("Accipitriformes", adj = 0)

hist(dat_hist_ans$rate.difference[dat_hist_ans$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "", breaks = 30, ylim = c(0, 40))
hist(abs(dat_hist_ans$rate.difference[dat_hist_ans$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 30)
title("Anseriformes", adj = 0)

hist(dat_hist_apo$rate.difference[dat_hist_apo$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "", breaks = 20)
hist(abs(dat_hist_apo$rate.difference[dat_hist_apo$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Apodiformes", adj = 0)
legend("topright", pch = 15, bty = 'n', col = col_alp, 
       legend  = c("Decrease", "Increase"))

hist(dat_hist_cha$rate.difference[dat_hist_cha$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", 
     xlab = "", breaks = 40)
hist(abs(dat_hist_cha$rate.difference[dat_hist_cha$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Charadriiformes", adj = 0)

hist(dat_hist_col$rate.difference[dat_hist_col$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "", breaks = 40)
hist(abs(dat_hist_col$rate.difference[dat_hist_col$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Columbiformes", adj = 0)

hist(dat_hist_gal$rate.difference[dat_hist_gal$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "", breaks = 30, ylim = c(0, 20))
hist(abs(dat_hist_gal$rate.difference[dat_hist_gal$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Galliformes", adj = 0)

hist(dat_hist_pas$rate.difference[dat_hist_pas$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", 
     xlab = "Absolute rate difference", breaks = 10)
hist(abs(dat_hist_pas$rate.difference[dat_hist_pas$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 10)
title("Passeriformes", adj = 0)

hist(dat_hist_pic$rate.difference[dat_hist_pic$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "Absolute rate difference", breaks = 20, ylim = c(0, 30))
hist(abs(dat_hist_pic$rate.difference[dat_hist_pic$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Piciformes", adj = 0)

hist(dat_hist_psi$rate.difference[dat_hist_psi$rate.difference > 0], 
     col = col_alp[2], border = col[2], xlim = c(0, 0.86), main = "", ylab = "",
     xlab = "Absolute rate difference", breaks = 20, ylim = c(0, 8))
hist(abs(dat_hist_psi$rate.difference[dat_hist_psi$rate.difference < 0]), 
     col = col_alp[1], border = col[1], add = T, breaks = 20)
title("Psittaciformes", adj = 0)

#hist(dat_hist_ave$rate.difference)

dev.off()

## 2. Using SSD as a state

sdi_disc <- sdi
sdi_disc[sdi_disc > 0] <- "Female-biased SSD"
sdi_disc[sdi_disc == 0] <- "Monomorphism"
sdi_disc[sdi_disc < 0] <- "Male-biased SSD"
sdi_disc <- as.factor(sdi_disc)

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
#save(RR2_col, file = "data/results/RR2_col.RData")

RR2_psi <- lapply(tr[1:100], search.shift.birds2, "Psittaciformes")
#save(RR2_psi, file = "data/results/RR2_psi.RData")

RR2_ans <- lapply(tr[1:100], search.shift.birds2, "Anseriformes")
#save(RR2_ans, file = "data/results/RR2_ans.RData")

RR2_acc <- lapply(tr[1:100], search.shift.birds2, "Accipitriformes")
#save(RR2_acc, file = "data/results/RR2_acc.RData")

RR2_gal <- lapply(tr[1:100], search.shift.birds2, "Galliformes")
#save(RR2_gal, file = "data/results/RR2_gal.RData")

RR2_pic <- lapply(tr[1:100], search.shift.birds2, "Piciformes")
#save(RR2_pic, file = "data/results/RR2_pic.RData")

RR2_apo <- lapply(tr[1:100], search.shift.birds2, "Apodiformes")
#save(RR2_apo, file = "data/results/RR2_apo.RData")

RR2_cha <- lapply(tr[1:100], search.shift.birds2, "Charadriiformes")
#save(RR2_cha, file = "data/results/RR2_cha.RData")

RR2_pas <- lapply(tr[1:100], search.shift.birds2, "Passeriformes")
#save(RR2_pas, file = "data/results/RR2_pas.RData")

#load(file = "data/results/RR2_col.RData")
#load(file = "data/results/RR2_psi.RData")
#load(file = "data/results/RR2_ans.RData")
#load(file = "data/results/RR2_acc.RData")
#load(file = "data/results/RR2_gal.RData")
#load(file = "data/results/RR2_pic.RData")
#load(file = "data/results/RR2_apo.RData")
#load(file = "data/results/RR2_cha.RData")
#load(file = "data/results/RR2_pas.RData")

RR2_acc_M_p <- RR2_acc_F_p <- RR2_ans_M_p <- RR2_ans_F_p <- RR2_col_M_p <-
  RR2_col_F_p <- RR2_gal_M_p <- RR2_gal_F_p <- numeric()

RR2_apo_M_p <- RR2_apo_F_p <- RR2_cha_M_p <- RR2_cha_F_p <- RR2_pas_M_p <- 
  RR2_pas_F_p <- RR2_pic_M_p <- RR2_pic_F_p <- RR2_psi_M_p <- RR2_psi_F_p <- 
  matrix(nrow = 100, ncol = 6)

colnames(RR2_apo_M_p) <- colnames(RR2_apo_F_p) <- colnames(RR2_cha_M_p) <- 
  colnames(RR2_cha_F_p) <- colnames(RR2_pas_M_p) <- colnames(RR2_pas_F_p) <- 
  colnames(RR2_pic_M_p) <- colnames(RR2_pic_F_p) <- colnames(RR2_psi_M_p) <- 
  colnames(RR2_psi_F_p) <- rownames(RR2_apo[[1]]$shift_male$state.results)
for (i in 1:100) {
  RR2_acc_M_p[i] <- RR2_acc[[i]]$shift_male$state.results[[2]]
  RR2_acc_F_p[i] <- RR2_acc[[i]]$shift_female$state.results[[2]]
  
  RR2_ans_M_p[i] <- RR2_ans[[i]]$shift_male$state.results[[2]]
  RR2_ans_F_p[i] <- RR2_ans[[i]]$shift_female$state.results[[2]]
  
  RR2_apo_M_p[i, ] <- RR2_apo[[i]]$shift_male$state.results[[2]]
  RR2_apo_F_p[i, ] <- RR2_apo[[i]]$shift_female$state.results[[2]]

  RR2_cha_M_p[i, ] <- RR2_cha[[i]]$shift_male$state.results[[2]]
  RR2_cha_F_p[i, ] <- RR2_cha[[i]]$shift_female$state.results[[2]]

  RR2_col_M_p[i] <- RR2_col[[i]]$shift_male$state.results[[2]]
  RR2_col_F_p[i] <- RR2_col[[i]]$shift_female$state.results[[2]]
  
  RR2_gal_M_p[i] <- RR2_gal[[i]]$shift_male$state.results[[2]]
  RR2_gal_F_p[i] <- RR2_gal[[i]]$shift_female$state.results[[2]]

  RR2_pas_M_p[i, ] <- RR2_pas[[i]]$shift_male$state.results[[2]]
  RR2_pas_F_p[i, ] <- RR2_pas[[i]]$shift_female$state.results[[2]]

  RR2_pic_M_p[i, ] <- RR2_pic[[i]]$shift_male$state.results[[2]]
  RR2_pic_F_p[i, ] <- RR2_pic[[i]]$shift_female$state.results[[2]]

  RR2_psi_M_p[i, ] <- RR2_psi[[i]]$shift_male$state.results[[2]]
  RR2_psi_F_p[i, ] <- RR2_psi[[i]]$shift_female$state.results[[2]]
}

res_RR2 <- matrix(nrow = 9, ncol = 12)
rownames(res_RR2) <- c("Accipitriformes", "Anseriformes", "Apodiformes", 
                       "Charadriiformes", "Columbiformes", "Galliformes", 
                       "Passeriformes", "Piciformes", "Psittaciformes")
colnames(res_RR2) <- c(paste0("M_", rownames(RR2_apo[[1]]$shift_male$state.results)),
                       paste0("F_", rownames(RR2_apo[[1]]$shift_female$state.results)))

res_RR2[1, 1] <- summ(RR2_acc_M_p)
res_RR2[1, 7] <- summ(RR2_acc_F_p)

res_RR2[2, 1] <- summ(RR2_ans_M_p)
res_RR2[2, 7] <- summ(RR2_ans_F_p)

res_RR2[3, 1:6] <- c(summ(RR2_apo_M_p[, 1]), summ(RR2_apo_M_p[, 2]), 
                     summ(RR2_apo_M_p[, 3]), summ(RR2_apo_M_p[, 4]), 
                     summ(RR2_apo_M_p[, 5]), summ(RR2_apo_M_p[, 6]))
res_RR2[3, 7:12] <- c(summ(RR2_apo_F_p[, 1]), summ(RR2_apo_F_p[, 2]), 
                     summ(RR2_apo_F_p[, 3]), summ(RR2_apo_F_p[, 4]), 
                     summ(RR2_apo_F_p[, 5]), summ(RR2_apo_F_p[, 6]))

res_RR2[4, 1:6] <- c(summ(RR2_cha_M_p[, 1]), summ(RR2_cha_M_p[, 2]), 
                     summ(RR2_cha_M_p[, 3]), summ(RR2_cha_M_p[, 4]), 
                     summ(RR2_cha_M_p[, 5]), summ(RR2_cha_M_p[, 6]))
res_RR2[4, 7:12] <- c(summ(RR2_cha_F_p[, 1]), summ(RR2_cha_F_p[, 2]), 
                      summ(RR2_cha_F_p[, 3]), summ(RR2_cha_F_p[, 4]), 
                      summ(RR2_cha_F_p[, 5]), summ(RR2_cha_F_p[, 6]))

res_RR2[5, 1] <- summ(RR2_col_M_p)
res_RR2[5, 7] <- summ(RR2_col_F_p)

res_RR2[6, 1] <- summ(RR2_gal_M_p)
res_RR2[6, 7] <- summ(RR2_gal_F_p)
                      
res_RR2[7, 1:6] <- c(summ(RR2_pas_M_p[, 1]), summ(RR2_pas_M_p[, 2]), 
                     summ(RR2_pas_M_p[, 3]), summ(RR2_pas_M_p[, 4]), 
                     summ(RR2_pas_M_p[, 5]), summ(RR2_pas_M_p[, 6]))
res_RR2[7, 7:12] <- c(summ(RR2_pas_F_p[, 1]), summ(RR2_pas_F_p[, 2]), 
                      summ(RR2_pas_F_p[, 3]), summ(RR2_pas_F_p[, 4]), 
                      summ(RR2_pas_F_p[, 5]), summ(RR2_pas_F_p[, 6]))           
  
res_RR2[8, 1:6] <- c(summ(RR2_pic_M_p[, 1]), summ(RR2_pic_M_p[, 2]), 
                     summ(RR2_pic_M_p[, 3]), summ(RR2_pic_M_p[, 4]), 
                     summ(RR2_pic_M_p[, 5]), summ(RR2_pic_M_p[, 6]))
res_RR2[8, 7:12] <- c(summ(RR2_pic_F_p[, 1]), summ(RR2_pic_F_p[, 2]), 
                      summ(RR2_pic_F_p[, 3]), summ(RR2_pic_F_p[, 4]), 
                      summ(RR2_pic_F_p[, 5]), summ(RR2_pic_F_p[, 6]))  

res_RR2[9, 1:6] <- c(summ(RR2_psi_M_p[, 1]), summ(RR2_psi_M_p[, 2]), 
                     summ(RR2_psi_M_p[, 3]), summ(RR2_psi_M_p[, 4]), 
                     summ(RR2_psi_M_p[, 5]), summ(RR2_psi_M_p[, 6]))
res_RR2[9, 7:12] <- c(summ(RR2_psi_F_p[, 1]), summ(RR2_psi_F_p[, 2]), 
                      summ(RR2_psi_F_p[, 3]), summ(RR2_psi_F_p[, 4]), 
                      summ(RR2_psi_F_p[, 5]), summ(RR2_psi_F_p[, 6]))

write.csv(res_RR2, "tables/TableS1_unformatted.csv")

acc_F_fem <- acc_F_mal <- acc_M_fem <- acc_M_mal <- list()
ans_F_fem <- ans_F_mal <- ans_M_fem <- ans_M_mal <- list()
apo_F_fem <- apo_F_mal <- apo_F_mon <- apo_M_fem <- apo_M_mal <- apo_M_mon <- list()
cha_F_fem <- cha_F_mal <- cha_F_mon <- cha_M_fem <- cha_M_mal <- cha_M_mon <- list()
col_F_fem <- col_F_mal <- col_M_fem <- col_M_mal <- list()
gal_F_fem <- gal_F_mal <- gal_M_fem <- gal_M_mal <- list()
pas_F_fem <- pas_F_mal <- pas_F_mon <- pas_M_fem <- pas_M_mal <- pas_M_mon <- list()
pic_F_fem <- pic_F_mal <- pic_F_mon <- pic_M_fem <- pic_M_mal <- pic_M_mon <- list()
psi_F_fem <- psi_F_mal <- psi_F_mon <- psi_M_fem <- psi_M_mal <- psi_M_mon <- list()
for (i in 1:100) {
  acc_M_mal[[i]] <- RR2_acc[[i]]$shift_male$rates[rownames(RR2_acc[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  acc_M_fem[[i]] <- RR2_acc[[i]]$shift_male$rates[rownames(RR2_acc[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  acc_F_mal[[i]] <- RR2_acc[[i]]$shift_female$rates[rownames(RR2_acc[[i]]$shift_female$rates) %in%
                                                  names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  acc_F_fem[[i]] <- RR2_acc[[i]]$shift_female$rates[rownames(RR2_acc[[i]]$shift_female$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]

  ans_M_mal[[i]] <- RR2_ans[[i]]$shift_male$rates[rownames(RR2_ans[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  ans_M_fem[[i]] <- RR2_ans[[i]]$shift_male$rates[rownames(RR2_ans[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  ans_F_mal[[i]] <- RR2_ans[[i]]$shift_female$rates[rownames(RR2_ans[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  ans_F_fem[[i]] <- RR2_ans[[i]]$shift_female$rates[rownames(RR2_ans[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
 
  apo_M_mal[[i]] <- RR2_apo[[i]]$shift_male$rates[rownames(RR2_apo[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  apo_M_fem[[i]] <- RR2_apo[[i]]$shift_male$rates[rownames(RR2_apo[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  apo_M_mon[[i]] <- RR2_apo[[i]]$shift_male$rates[rownames(RR2_apo[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]
  apo_F_mal[[i]] <- RR2_apo[[i]]$shift_female$rates[rownames(RR2_apo[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  apo_F_fem[[i]] <- RR2_apo[[i]]$shift_female$rates[rownames(RR2_apo[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  apo_F_mon[[i]] <- RR2_apo[[i]]$shift_female$rates[rownames(RR2_apo[[i]]$shift_female$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]

  cha_M_mal[[i]] <- RR2_cha[[i]]$shift_male$rates[rownames(RR2_cha[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  cha_M_fem[[i]] <- RR2_cha[[i]]$shift_male$rates[rownames(RR2_cha[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  cha_M_mon[[i]] <- RR2_cha[[i]]$shift_male$rates[rownames(RR2_cha[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]
  cha_F_mal[[i]] <- RR2_cha[[i]]$shift_female$rates[rownames(RR2_cha[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  cha_F_fem[[i]] <- RR2_cha[[i]]$shift_female$rates[rownames(RR2_cha[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  cha_F_mon[[i]] <- RR2_cha[[i]]$shift_female$rates[rownames(RR2_cha[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Monomorphism"]]

  col_M_mal[[i]] <- RR2_col[[i]]$shift_male$rates[rownames(RR2_col[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  col_M_fem[[i]] <- RR2_col[[i]]$shift_male$rates[rownames(RR2_col[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  col_F_mal[[i]] <- RR2_col[[i]]$shift_female$rates[rownames(RR2_col[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  col_F_fem[[i]] <- RR2_col[[i]]$shift_female$rates[rownames(RR2_col[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]

  gal_M_mal[[i]] <- RR2_gal[[i]]$shift_male$rates[rownames(RR2_gal[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  gal_M_fem[[i]] <- RR2_gal[[i]]$shift_male$rates[rownames(RR2_gal[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  gal_F_mal[[i]] <- RR2_gal[[i]]$shift_female$rates[rownames(RR2_gal[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  gal_F_fem[[i]] <- RR2_gal[[i]]$shift_female$rates[rownames(RR2_gal[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]

  pas_M_mal[[i]] <- RR2_pas[[i]]$shift_male$rates[rownames(RR2_pas[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  pas_M_fem[[i]] <- RR2_pas[[i]]$shift_male$rates[rownames(RR2_pas[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  pas_M_mon[[i]] <- RR2_pas[[i]]$shift_male$rates[rownames(RR2_pas[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]
  pas_F_mal[[i]] <- RR2_pas[[i]]$shift_female$rates[rownames(RR2_pas[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  pas_F_fem[[i]] <- RR2_pas[[i]]$shift_female$rates[rownames(RR2_pas[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  pas_F_mon[[i]] <- RR2_pas[[i]]$shift_female$rates[rownames(RR2_pas[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Monomorphism"]]

  pic_M_mal[[i]] <- RR2_pic[[i]]$shift_male$rates[rownames(RR2_pic[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  pic_M_fem[[i]] <- RR2_pic[[i]]$shift_male$rates[rownames(RR2_pic[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  pic_M_mon[[i]] <- RR2_pic[[i]]$shift_male$rates[rownames(RR2_pic[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]
  pic_F_mal[[i]] <- RR2_pic[[i]]$shift_female$rates[rownames(RR2_pic[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  pic_F_fem[[i]] <- RR2_pic[[i]]$shift_female$rates[rownames(RR2_pic[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  pic_F_mon[[i]] <- RR2_pic[[i]]$shift_female$rates[rownames(RR2_pic[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Monomorphism"]]

  psi_M_mal[[i]] <- RR2_psi[[i]]$shift_male$rates[rownames(RR2_psi[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  psi_M_fem[[i]] <- RR2_psi[[i]]$shift_male$rates[rownames(RR2_psi[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  psi_M_mon[[i]] <- RR2_psi[[i]]$shift_male$rates[rownames(RR2_psi[[i]]$shift_male$rates) %in%
                                                    names(sdi_disc)[sdi_disc == "Monomorphism"]]
  psi_F_mal[[i]] <- RR2_psi[[i]]$shift_female$rates[rownames(RR2_psi[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Male-biased SSD"]]
  psi_F_fem[[i]] <- RR2_psi[[i]]$shift_female$rates[rownames(RR2_psi[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Female-biased SSD"]]
  psi_F_mon[[i]] <- RR2_psi[[i]]$shift_female$rates[rownames(RR2_psi[[i]]$shift_female$rates) %in%
                                                      names(sdi_disc)[sdi_disc == "Monomorphism"]]
}
acc_M_mal <- do.call(c, acc_M_mal)
acc_M_fem <- do.call(c, acc_M_fem)
acc_F_mal <- do.call(c, acc_F_mal)
acc_F_fem <- do.call(c, acc_F_fem)

ans_M_mal <- do.call(c, ans_M_mal)
ans_M_fem <- do.call(c, ans_M_fem)
ans_F_mal <- do.call(c, ans_F_mal)
ans_F_fem <- do.call(c, ans_F_fem)

apo_M_mal <- do.call(c, apo_M_mal)
apo_M_fem <- do.call(c, apo_M_fem)
apo_M_mon <- do.call(c, apo_M_mon)
apo_F_mal <- do.call(c, apo_F_mal)
apo_F_fem <- do.call(c, apo_F_fem)
apo_F_mon <- do.call(c, apo_F_mon)

cha_M_mal <- do.call(c, cha_M_mal)
cha_M_fem <- do.call(c, cha_M_fem)
cha_M_mon <- do.call(c, cha_M_mon)
cha_F_mal <- do.call(c, cha_F_mal)
cha_F_fem <- do.call(c, cha_F_fem)
cha_F_mon <- do.call(c, cha_F_mon)

col_M_mal <- do.call(c, col_M_mal)
col_M_fem <- do.call(c, col_M_fem)
col_F_mal <- do.call(c, col_F_mal)
col_F_fem <- do.call(c, col_F_fem)

gal_M_mal <- do.call(c, gal_M_mal)
gal_M_fem <- do.call(c, gal_M_fem)
gal_F_mal <- do.call(c, gal_F_mal)
gal_F_fem <- do.call(c, gal_F_fem)

pas_M_mal <- do.call(c, pas_M_mal)
pas_M_fem <- do.call(c, pas_M_fem)
pas_M_mon <- do.call(c, pas_M_mon)
pas_F_mal <- do.call(c, pas_F_mal)
pas_F_fem <- do.call(c, pas_F_fem)
pas_F_mon <- do.call(c, pas_F_mon)

pic_M_mal <- do.call(c, pic_M_mal)
pic_M_fem <- do.call(c, pic_M_fem)
pic_M_mon <- do.call(c, pic_M_mon)
pic_F_mal <- do.call(c, pic_F_mal)
pic_F_fem <- do.call(c, pic_F_fem)
pic_F_mon <- do.call(c, pic_F_mon)

psi_M_mal <- do.call(c, psi_M_mal)
psi_M_fem <- do.call(c, psi_M_fem)
psi_M_mon <- do.call(c, psi_M_mon)
psi_F_mal <- do.call(c, psi_F_mal)
psi_F_fem <- do.call(c, psi_F_fem)
psi_F_mon <- do.call(c, psi_F_mon)

## Figure 8
pdf("figures/FigureS1.pdf", height = 9, width = 9)

layout(matrix(1:9, ncol = 3, byrow = TRUE))

male <- "#9966FF"
monom <- "gray"
female <- "#E69F00"

plot(density(abs(acc_M_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 1))
lines(density(abs(acc_M_mal)), col = male)
title("Accipitriformes", adj = 0)

plot(density(abs(ans_M_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.5))
lines(density(abs(ans_M_mal)), col = male)
title("Anseriformes", adj = 0)

plot(density(abs(apo_M_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.6))
lines(density(abs(apo_M_mal)), col = male)
lines(density(abs(apo_M_mon)), col = monom)
title("Apodiformes", adj = 0)

plot(density(abs(cha_M_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 0.6))
lines(density(abs(cha_M_mal)), col = male)
lines(density(abs(cha_M_mon)), col = monom)
title("Charadriiformes", adj = 0)

plot(density(abs(col_M_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 1))
lines(density(abs(col_M_mal)), col = male)
title("Columbiformes", adj = 0)

plot(density(abs(gal_M_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.8))
lines(density(abs(gal_M_mal)), col = male)
title("Galliformes", adj = 0)

plot(density(abs(pas_M_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 0.7))
lines(density(abs(pas_M_mal)), col = male)
lines(density(abs(pas_M_mon)), col = monom)
title("Passeriformes", adj = 0)

plot(density(abs(pic_M_fem)), type = "l", main = "", ylab = "", col = female,
     xlab = "Absolute male body mass rate residuals", ylim = c(0, 0.6))
lines(density(abs(pic_M_mal)), col = male)
lines(density(abs(pic_M_mon)), col = monom)
title("Piciformes", adj = 0)

plot(density(abs(psi_M_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.6))
lines(density(abs(psi_M_mal)), col = male)
lines(density(abs(psi_M_mon)), col = monom)
title("Psittaciformes", adj = 0)

dev.off()

## Figure 9 

pdf("figures/FigureS2.pdf", height = 9, width = 9)

layout(matrix(1:9, ncol = 3, byrow = TRUE))

plot(density(abs(acc_F_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 1.5))
lines(density(abs(acc_F_mal)), col = male)
title("Accipitriformes", adj = 0)
title("*", adj = 1, cex = 2)
plot(density(abs(ans_F_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.5))
lines(density(abs(ans_F_mal)), col = male)
title("Anseriformes", adj = 0)

plot(density(abs(apo_F_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.7))
lines(density(abs(apo_F_mal)), col = male)
lines(density(abs(apo_F_mon)), col = monom)
title("Apodiformes", adj = 0)

plot(density(abs(cha_F_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 0.6))
lines(density(abs(cha_F_mal)), col = male)
lines(density(abs(cha_F_mon)), col = monom)
title("Charadriiformes", adj = 0)

plot(density(abs(col_F_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 1))
lines(density(abs(col_F_mal)), col = male)
title("Columbiformes", adj = 0)

plot(density(abs(gal_F_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.8))
lines(density(abs(gal_F_mal)), col = male)
title("Galliformes", adj = 0)

plot(density(abs(pas_F_fem)), type = "l", main = "", xlab = "", 
     ylab = "Density", col = female, ylim = c(0, 0.7))
lines(density(abs(pas_F_mal)), col = male)
lines(density(abs(pas_F_mon)), col = monom)
title("Passeriformes", adj = 0)

plot(density(abs(pic_F_fem)), type = "l", main = "", ylab = "", col = female,
     xlab = "Absolute female body mass rate residuals", ylim = c(0, 0.6))
lines(density(abs(pic_F_mal)), col = male)
lines(density(abs(pic_F_mon)), col = monom)
title("Piciformes", adj = 0)

plot(density(abs(psi_F_fem)), type = "l", main = "", xlab = "", ylab = "", 
     col = female, ylim = c(0, 0.7))
lines(density(abs(psi_F_mal)), col = male)
lines(density(abs(psi_F_mon)), col = monom)
title("Psittaciformes", adj = 0)

dev.off()