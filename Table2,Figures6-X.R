rm(list = ls())

setwd("Documents/lab/dimorph_evol")

library(mvMORPH)
library(phytools)
library(geiger)
library(RRphylo)

########### mvBM + RRphylo ###########

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_30may22_edit.csv", row.names = 1)

subdat <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
               complete.cases(dat$Body_mass_g_F_mean), ]
rownames(subdat) <- dat$Scientific_name[complete.cases(dat$Body_mass_g_M_mean) &
										complete.cases(dat$Body_mass_g_F_mean)]

subdat$Body_mass_g_M_mean <- log(subdat$Body_mass_g_M_mean)
subdat$Body_mass_g_F_mean <- log(subdat$Body_mass_g_F_mean)

############

## mvBM

multiModelComp <- function(phy, data, taxon = NULL, ntree = 1) {

	if (!is.null(taxon)) {
		subdat <- subset(data, Order == taxon, c(Body_mass_g_M_mean,
		                                         Body_mass_g_F_mean))

		rownames(subdat) <- subset(data, Order == taxon, Scientific_name)[, 1]

	} else {
		subdat <- data[, c("Body_mass_g_M_mean", "Body_mass_g_F_mean")]

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
# Resultados iguais nos três

# Análises feitas no cluster da UFG

## Columbiformes = 105 spp
fit_col <- multiModelComp(phy = tr, data = subdat, taxon = "Columbiformes",
                          ntree = 100)
#save(fit_col, file = "data/aves/results/fitBM_col.RData")

## Psittaciformes = 131 spp
fit_psi <- multiModelComp(phy = tr, data = subdat, taxon = "Psittaciformes",
                          ntree = 100)
#save(fit_psi, file = "data/aves/results/fitBM_psi.RData")

## Anseriformes = 150 spp
fit_ans <- multiModelComp(phy = tr, data = subdat, taxon = "Anseriformes",
                          ntree = 100)
#save(fit_ans, file = "data/aves/results/fitBM_ans.RData")

## Accipitriformes = 169 spp
fit_acc <- multiModelComp(phy = tr, data = subdat, taxon = "Accipitriformes",
                          ntree = 100)
#save(fit_acc, file = "data/aves/results/fitBM_acc.RData")

## Galliformes = 179 spp
fit_gal <- multiModelComp(phy = tr, data = subdat, taxon = "Galliformes",
                          ntree = 100)
#save(fit_gal, file = "data/aves/results/fitBM_gal.RData")

## Piciformes = 217 spp
fit_pic <- multiModelComp(phy = tr, data = subdat, taxon = "Piciformes",
                          ntree = 100)
#save(fit_pic, file = "data/aves/results/fitBM_pic.RData")

## Apodiformes = 247 spp
fit_apo <- multiModelComp(phy = tr, data = subdat, taxon = "Apodiformes",
                          ntree = 100)
#save(fit_apo, file = "data/aves/results/fitBM_apo.RData")

## Charadriiformes = 252 spp
fit_cha <- multiModelComp(phy = tr, data = subdat, taxon = "Charadriiformes",
                          ntree = 100)
#save(fit_cha, file = "data/aves/results/fitBM_cha.RData")

## Passeriformes = 2290 spp
fit_pas <- multiModelComp(phy = tr, data = subdat, taxon = "Passeriformes",
                          ntree = 100)
#save(fit_pas, file = "data/aves/results/fitBM_pas.RData")

#load(file = "data/aves/results/fitBM_col.RData")
#load(file = "data/aves/results/fitBM_psi.RData")
#load(file = "data/aves/results/fitBM_ans.RData")
#load(file = "data/aves/results/fitBM_acc.RData")
#load(file = "data/aves/results/fitBM_gal.RData")
#load(file = "data/aves/results/fitBM_pic.RData")
#load(file = "data/aves/results/fitBM_apo.RData")
#load(file = "data/aves/results/fitBM_cha.RData")
#load(file = "data/aves/results/fitBM_pas.RData")

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

write.csv(fit_birds, "tables/Table2_unformatted.csv")

#################

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

# Análises feitas no cluster da UFPR
RR1_col <- lapply(tr[1:100], search.shift.birds1, "Columbiformes")
save(RR1_col, file = "data/aves/results/RR1_col.RData")

RR1_psi <- lapply(tr[1:100], search.shift.birds1, "Psittaciformes")
save(RR1_psi, file = "data/aves/results/RR1_psi.RData")

RR1_ans <- lapply(tr[1:100], search.shift.birds1, "Anseriformes")
save(RR1_ans, file = "data/aves/results/RR1_ans.RData")

RR1_acc <- lapply(tr[1:100], search.shift.birds1, "Accipitriformes")
save(RR1_acc, file = "data/aves/results/RR1_acc.RData")

RR1_gal <- lapply(tr[1:100], search.shift.birds1, "Galliformes")
save(RR1_gal, file = "data/aves/results/RR1_gal.RData")

RR1_pic <- lapply(tr[1:100], search.shift.birds1, "Piciformes")
save(RR1_pic, file = "data/aves/results/RR1_pic.RData")

RR1_apo <- lapply(tr[1:100], search.shift.birds1, "Apodiformes")
save(RR1_apo, file = "data/aves/results/RR1_apo.RData")

RR1_cha <- lapply(tr[1:100], search.shift.birds1, "Charadriiformes")
save(RR1_cha, file = "data/aves/results/RR1_cha.RData")

RR1_pas <- lapply(tr[1:100], search.shift.birds1, "Passeriformes")
save(RR1_pas, file = "data/aves/results/RR1_pas.RData")

#load(file = "data/aves/results/RR1_col.RData")
#load(file = "data/aves/results/RR1_psi.RData")
#load(file = "data/aves/results/RR1_ans.RData")
#load(file = "data/aves/results/RR1_acc.RData")
#load(file = "data/aves/results/RR1_gal.RData")
#load(file = "data/aves/results/RR1_pic.RData")
#load(file = "data/aves/results/RR1_apo.RData")
#load(file = "data/aves/results/RR1_cha.RData")
#load(file = "data/aves/results/RR1_pas.RData")

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

#sig_ave <- lapply(RR1_ave, function(x) x$shift$all.clades[x$shift$all.clades$p.value <= 
#                                                            0.025 | 
#                                                            x$shift$all.clades$p.value >= 
#                                                            0.975, ])

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
#dat_box_ave <- data.boxplot(sig_ave)

dat_hist_col <- do.call(rbind, sig_col)
dat_hist_psi <- do.call(rbind, sig_psi)
dat_hist_ans <- do.call(rbind, sig_ans)
dat_hist_acc <- do.call(rbind, sig_acc)
dat_hist_gal <- do.call(rbind, sig_gal)
dat_hist_pic <- do.call(rbind, sig_pic)
dat_hist_apo <- do.call(rbind, sig_apo)
dat_hist_cha <- do.call(rbind, sig_cha)
dat_hist_pas <- do.call(rbind, sig_pas)
#dat_hist_ave <- do.call(rbind, sig_ave)

# Figure 4

pdf("figures/Figure4.pdf", width = 9, height = 9)

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

# Figure 5 

pdf("figures/Figure5.pdf", width = 9, height = 9)

layout(matrix(1:9, nrow = 3, byrow = TRUE))

par(mar = c(4, 4, 4, 2))

h_acc <- hist(dat_hist_acc$rate.difference, breaks = 25, plot = F)
cuts_acc <- cut(h_acc$breaks, c(min(dat_hist_acc$rate.difference), 0,
                            max(dat_hist_acc$rate.difference)))
plot(h_acc, main = "", xlab = "", col = col_alp[cuts_acc], border = col[cuts_acc], 
     xlim = c(-0.79, 0.86))
title("Accipitriformes", adj = 0)

h_ans <- hist(dat_hist_ans$rate.difference, breaks = 25, plot = F)
cuts_ans <- cut(h_ans$breaks, c(min(dat_hist_ans$rate.difference), 0,
                                max(dat_hist_ans$rate.difference)))
plot(h_ans, main = "", xlab = "", ylab = "", col = col_alp[cuts_ans], 
     border = col[cuts_ans], xlim = c(-0.79, 0.86))
title("Anseriformes", adj = 0)

h_apo <- hist(dat_hist_apo$rate.difference, breaks = 20, plot = F)
cuts_apo <- cut(h_apo$breaks, c(min(dat_hist_apo$rate.difference), 0,
                                max(dat_hist_apo$rate.difference)))
plot(h_apo, main = "", xlab = "", ylab = "", col = col_alp[cuts_apo], 
     border = col[cuts_apo], xlim = c(-0.79, 0.86))
title("Apodiformes", adj = 0)

h_cha <- hist(dat_hist_cha$rate.difference, breaks = 20, plot = F)
cuts_cha <- cut(h_cha$breaks, c(min(dat_hist_cha$rate.difference), 0,
                                max(dat_hist_cha$rate.difference)))
plot(h_cha, main = "", xlab = "", col = col_alp[cuts_cha], 
     border = col[cuts_cha], xlim = c(-0.79, 0.86))
title("Charadriiformes", adj = 0)

h_col <- hist(dat_hist_col$rate.difference, breaks = 30, plot = F)
cuts_col <- cut(h_col$breaks, c(min(dat_hist_col$rate.difference), 0,
                                max(dat_hist_col$rate.difference)))
plot(h_col, main = "", xlab = "", ylab = "", col = col_alp[cuts_col], 
     border = col[cuts_col], xlim = c(-0.79, 0.86))
title("Columbiformes", adj = 0)

h_gal <- hist(dat_hist_gal$rate.difference, breaks = 20, plot = F)
cuts_gal <- cut(h_gal$breaks, c(min(dat_hist_gal$rate.difference), 0,
                                max(dat_hist_gal$rate.difference)))
plot(h_gal, main = "", xlab = "", ylab = "", col = col_alp[cuts_gal], 
     border = col[cuts_gal], xlim = c(-0.79, 0.86))
title("Galliformes", adj = 0)

h_pas <- hist(dat_hist_pas$rate.difference, breaks = 10, plot = F)
cuts_pas <- cut(h_pas$breaks, c(min(dat_hist_pas$rate.difference), 0,
                                max(dat_hist_pas$rate.difference)))
plot(h_pas, main = "", xlab = "Rate difference", col = col_alp[cuts_pas], 
     border = col[cuts_pas], xlim = c(-0.79, 0.86))
title("Passeriformes", adj = 0)

h_pic <- hist(dat_hist_pic$rate.difference, breaks = 20, plot = F)
cuts_pic <- cut(h_pic$breaks, c(min(dat_hist_pic$rate.difference), 0,
                                max(dat_hist_pic$rate.difference)))
plot(h_pic, main = "", xlab = "Rate difference", ylab = "", col = col_alp[cuts_pic], 
     border = col[cuts_pic], xlim = c(-0.79, 0.86))
title("Piciformes", adj = 0)

h_psi <- hist(dat_hist_psi$rate.difference, breaks = 30, plot = F)
cuts_psi <- cut(h_psi$breaks, c(min(dat_hist_psi$rate.difference), 0,
                                max(dat_hist_psi$rate.difference)))
plot(h_psi, main = "", xlab = "Rate difference", ylab = "", col = col_alp[cuts_psi], 
     border = col[cuts_psi], xlim = c(-0.79, 0.86))
title("Psittaciformes", adj = 0)

#hist(dat_hist_ave$rate.difference)

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
save(RR2_col, file = "data/aves/results/RR2_col.RData")

RR2_psi <- lapply(tr[1:100], search.shift.birds2, "Psittaciformes")
save(RR2_psi, file = "data/aves/results/RR2_psi.RData")

RR2_ans <- lapply(tr[1:100], search.shift.birds2, "Anseriformes")
save(RR2_ans, file = "data/aves/results/RR2_ans.RData")

RR2_acc <- lapply(tr[1:100], search.shift.birds2, "Accipitriformes")
save(RR2_acc, file = "data/aves/results/RR2_acc.RData")

RR2_gal <- lapply(tr[1:100], search.shift.birds2, "Galliformes")
save(RR2_gal, file = "data/aves/results/RR2_gal.RData")

RR2_pic <- lapply(tr[1:100], search.shift.birds2, "Piciformes")
save(RR2_pic, file = "data/aves/results/RR2_pic.RData")

RR2_apo <- lapply(tr[1:100], search.shift.birds2, "Apodiformes")
save(RR2_apo, file = "data/aves/results/RR2_apo.RData")

RR2_cha <- lapply(tr[1:100], search.shift.birds2, "Charadriiformes")
save(RR2_cha, file = "data/aves/results/RR2_cha.RData")

RR2_pas <- lapply(tr[1:100], search.shift.birds2, "Passeriformes")
save(RR2_pas, file = "data/aves/results/RR2_pas.RData")

