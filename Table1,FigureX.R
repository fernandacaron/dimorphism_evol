rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(mvMORPH)
library(phytools)
library(geiger)
library(colorspace)

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

subdat <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
               complete.cases(dat$Body_mass_g_F_mean), ]
rownames(subdat) <- dat$Scientific_name[complete.cases(dat$Body_mass_g_M_mean) &
										complete.cases(dat$Body_mass_g_F_mean)]

subdat$Body_mass_g_M_mean <- log(subdat$Body_mass_g_M_mean)
subdat$Body_mass_g_F_mean <- log(subdat$Body_mass_g_F_mean)

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


## Aves = 105 spp
fit_ave <- multiModelComp(phy = tr, data = subdat, ntree = 5)
save(fit_ave, file = "data/aves/fit_ave.RData")

## Columbiformes = 105 spp
fit_col <- multiModelComp(phy = tr, data = subdat, taxon = "Columbiformes",
                          ntree = 100)
save(fit_col, file = "data/aves/fit_col.RData")

## Psittaciformes = 131 spp
fit_psi <- multiModelComp(phy = tr, data = subdat, taxon = "Psittaciformes",
                          ntree = 100)
save(fit_psi, file = "data/aves/fit_psi.RData")

## Anseriformes = 150 spp
fit_ans <- multiModelComp(phy = tr, data = subdat, taxon = "Anseriformes",
                          ntree = 100)
save(fit_ans, file = "data/aves/fit_ans.RData")

## Accipitriformes = 169 spp
fit_acc <- multiModelComp(phy = tr, data = subdat, taxon = "Accipitriformes",
                          ntree = 100)
save(fit_acc, file = "data/aves/fit_acc.RData")

## Galliformes = 179 spp
fit_gal <- multiModelComp(phy = tr, data = subdat, taxon = "Galliformes",
                          ntree = 100)
save(fit_gal, file = "data/aves/fit_gal.RData")

## Piciformes = 217 spp
fit_pic <- multiModelComp(phy = tr, data = subdat, taxon = "Piciformes",
                          ntree = 100)
save(fit_pic, file = "data/aves/fit_pic.RData")

## Apodiformes = 247 spp
fit_apo <- multiModelComp(phy = tr, data = subdat, taxon = "Apodiformes",
                          ntree = 100)
save(fit_apo, file = "data/aves/fit_apo.RData")

## Charadriiformes = 252 spp
fit_cha <- multiModelComp(phy = tr, data = subdat, taxon = "Charadriiformes",
                          ntree = 100)
save(fit_cha, file = "data/aves/fit_cha.RData")

## Passeriformes = 2290 spp
fit_pas <- multiModelComp(phy = tr, data = subdat, taxon = "Passeriformes",
                          ntree = 100)
save(fit_pas, file = "data/aves/fit_pas.RData")


## Table 1

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
sigmaM_acc <- sigmaM_ans <- sigmaM_apo <- sigmaM_cha <- sigmaM_col <-
	sigmaM_gal <- sigmaM_pas <- sigmaM_pic <- sigmaM_psi <- numeric()
sigmaF_acc <- sigmaF_ans <- sigmaF_apo <- sigmaF_cha <- sigmaF_col <-
	sigmaF_gal <- sigmaF_pas <- sigmaF_pic <- sigmaF_psi <- numeric()
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

fit_birds <- data.frame(nspp = c(169, 150, 247, 252, 105, 179, 2290, 217, 131),
                        nsim = c(100, 100, 100, 100, 100, 100, 5, 100, 100),
                        SDI = rep(NA, 9),
                        sigmaM = rep(NA, 9), 
                        sigmaF = rep(NA, 9), 
                        pval = rep(NA, 9),
                        row.names = c("Accipitriformes", "Anseriformes",
                                      "Apodiformes", "Charadriiformes",
                                      "Columbiformes", "Galliformes",
                                      "Passeriformes", "Piciformes",
                                      "Psittaciformes"))

fit_birds[1, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Accipitriformes"]), 3),
                                        " (",
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Accipitriformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Accipitriformes"]), 3), 
                                        ")"),
                           paste0(round(median(sigmaM_acc), 3),
                                  " (",
                                  round(min(sigmaM_acc), 3),
                                  "-",
                                  round(max(sigmaM_acc), 3),
                                  ")"),
                           paste0(round(median(sigmaF_acc), 3),
                                  " (",
                                  round(min(sigmaF_acc), 3),
                                  "-",
                                  round(max(sigmaF_acc), 3),
                                  ")"),
                           paste0(round(median(fit_acc$pval), 3),
                                  " (",
                                  round(min(fit_acc$pval), 3),
                                  "-",
                                  round(max(fit_acc$pval), 3),
                                  ")")
)

fit_birds[2, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Anseriformes"]), 3),
                                        " (",
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Anseriformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Anseriformes"]), 3), 
                                        ")"),
                           paste0(round(median(sigmaM_ans), 3),
                                  " (",
                                  round(min(sigmaM_ans), 3),
                                  "-",
                                  round(max(sigmaM_ans), 3),
                                  ")"),
                           paste0(round(median(sigmaF_ans), 3),
                                  " (",
                                  round(min(sigmaF_ans), 3),
                                  "-",
                                  round(max(sigmaF_ans), 3),
                                  ")"),
                           paste0(round(median(fit_ans$pval), 3),
                                  " (",
                                  round(min(fit_ans$pval), 3),
                                  "-",
                                  round(max(fit_ans$pval), 3),
                                  ")")
)


fit_birds[3, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Apodiformes"]), 3),
                                        " (", 
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Apodiformes"]), 3), 
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Apodiformes"]), 3),
                                        ")"),
                           paste0(round(median(sigmaM_apo), 3),
                                  " (",
                                  round(min(sigmaM_apo), 3),
                                  "-",
                                  round(max(sigmaM_apo), 3),
                                  ")"),
                           paste0(round(median(sigmaF_apo), 3),
                                  " (",
                                  round(min(sigmaF_apo), 3),
                                  "-",
                                  round(max(sigmaF_apo), 3),
                                  ")"),
                           paste0(round(median(fit_apo$pval), 3),
                                  " (",
                                  round(min(fit_apo$pval), 3),
                                  "-",
                                  round(max(fit_apo$pval), 3),
                                  ")")
)


fit_birds[4, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Charadriiformes"]), 3),
                                        " (", 
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Charadriiformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Charadriiformes"]), 3),
                                        ")"),
                           paste0(round(median(sigmaM_cha), 3),
                                  " (",
                                  round(min(sigmaM_cha), 3),
                                  "-",
                                  round(max(sigmaM_cha), 3),
                                  ")"),
                           paste0(round(median(sigmaF_cha), 3),
                                  " (",
                                  round(min(sigmaF_cha), 3),
                                  "-",
                                  round(max(sigmaF_cha), 3),
                                  ")"),
                           paste0(round(median(fit_cha$pval), 3),
                                  " (",
                                  round(min(fit_cha$pval), 3),
                                  "-",
                                  round(max(fit_cha$pval), 3),
                                  ")")
)


fit_birds[5, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Columbiformes"]), 3),
                                        " (",
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Columbiformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Columbiformes"]), 3),
                                        ")"),
                           paste0(round(median(sigmaM_col), 3),
                                  " (",
                                  round(min(sigmaM_col), 3),
                                  "-",
                                  round(max(sigmaM_col), 3),
                                  ")"),
                           paste0(round(median(sigmaF_col), 3),
                                  " (",
                                  round(min(sigmaF_col), 3),
                                  "-",
                                  round(max(sigmaF_col), 3),
                                  ")"),
                           paste0(round(median(fit_col$pval), 3),
                                  " (",
                                  round(min(fit_col$pval), 3),
                                  "-",
                                  round(max(fit_col$pval), 3),
                                  ")")
)


fit_birds[6, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Galliformes"]), 3),
                                        " (", 
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Galliformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Galliformes"]), 3),
                                        ")"),
                           paste0(round(median(sigmaM_gal), 3),
                                  " (",
                                  round(min(sigmaM_gal), 3),
                                  "-",
                                  round(max(sigmaM_gal), 3),
                                  ")"),
                           paste0(round(median(sigmaF_gal), 3),
                                  " (",
                                  round(min(sigmaF_gal), 3),
                                  "-",
                                  round(max(sigmaF_gal), 3),
                                  ")"),
                           paste0(round(median(fit_gal$pval), 3),
                                  " (",
                                  round(min(fit_gal$pval), 3),
                                  "-",
                                  round(max(fit_gal$pval), 3),
                                  ")")
)


fit_birds[7, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Passeriformes"]), 3),
                                         " (", 
                                         round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Passeriformes"]), 3),
                                         "-",
                                         round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Passeriformes"]), 3),
                                         ")"),
                           paste0(round(median(sigmaM_pas), 3),
                                  " (",
                                  round(min(sigmaM_pas), 3),
                                  "-",
                                  round(max(sigmaM_pas), 3),
                                  ")"),
                           paste0(round(median(sigmaF_pas), 3),
                                  " (",
                                  round(min(sigmaF_pas), 3),
                                  "-",
                                  round(max(sigmaF_pas), 3),
                                  ")"),
                           paste0(round(median(fit_pas$pval), 3),
                                  " (",
                                  round(min(fit_pas$pval), 3),
                                  "-",
                                  round(max(fit_pas$pval), 3),
                                  ")")
)


fit_birds[8, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Piciformes"]), 3),
                                        " (", 
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Piciformes"]), 3), 
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Piciformes"]), 3), 
                                        ")"),
                           paste0(round(median(sigmaM_pic), 3),
                                  " (",
                                  round(min(sigmaM_pic), 3),
                                  "-",
                                  round(max(sigmaM_pic), 3),
                                  ")"),
                           paste0(round(median(sigmaF_pic), 3),
                                  " (",
                                  round(min(sigmaF_pic), 3),
                                  "-",
                                  round(max(sigmaF_pic), 3),
                                  ")"),
                           paste0(round(median(fit_pic$pval), 3),
                                  " (",
                                  round(min(fit_pic$pval), 3),
                                  "-",
                                  round(max(fit_pic$pval), 3),
                                  ")")
)


fit_birds[9, 3:6] <- c(paste0(round(median(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Psittaciformes"]), 3),
                                        " (",
                                        round(min(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Psittaciformes"]), 3),
                                        "-",
                                        round(max(subdat_sdi$sdi[subdat_sdi$Order ==
                                        "Psittaciformes"]), 3), 
                                        ")"),
                           paste0(round(median(sigmaM_psi), 3),
                                  " (",
                                  round(min(sigmaM_psi), 3),
                                  "-",
                                  round(max(sigmaM_psi), 3),
                                  ")"),
                           paste0(round(median(sigmaF_psi), 3),
                                  " (",
                                  round(min(sigmaF_psi), 3),
                                  "-",
                                  round(max(sigmaF_psi), 3),
                                  ")"),
                           paste0(round(median(fit_psi$pval), 3),
                                  " (",
                                  round(min(fit_psi$pval), 3),
                                  "-",
                                  round(max(fit_psi$pval), 3),
                                  ")")
)

write.csv(fit_birds, "tables/Table1_unformatted.csv")

par(mar = c(2, 4, 4, 2))

hist(fit_acc$pval


## Figure 4
pdf("figures/Figure4.pdf", width = 9, height = 21)

layout(matrix(1:18, ncol = 2, byrow = T))

par(mar = c(4, 4, 4, 2))

cols <- "gray"
cols_al <- rgb(190/255, 190/255, 190/255, 0.5)

male <- mako(20)[9]
female <- rocket(20)[9]

male_al <- lighten(male, amount = 0.5)
female_al <- lighten(female, amount = 0.5)

## Sigmas do fit 1 (sigmas diferentes)
sigmaM_acc <- sigmaM_ans <- sigmaM_apo <- sigmaM_cha <- sigmaM_col <-
	sigmaM_gal <- sigmaM_pas <- sigmaM_pic <- sigmaM_psi <- numeric()
sigmaF_acc <- sigmaF_ans <- sigmaF_apo <- sigmaF_cha <- sigmaF_col <-
	sigmaF_gal <- sigmaF_pas <- sigmaF_pic <- sigmaF_psi <- numeric()
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

par(mar = c(2, 4, 4, 2))

hist(fit_acc$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Accipitriformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_acc, main = "", xlab = "", col = male_al, border = male,
     breaks = 20, xlim = c(0.042, 0.068), ylab = "")
hist(sigmaF_acc, add = T, col = female_al, border = female, breaks = 2)
legend("topright", pch = 15, bty = 'n', col = c(male_al, female_al), 
       legend  = c("Male body size", "Female body size"))

par(mar = c(2, 4, 4, 2))

hist(fit_ans$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Anseriformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_ans, main = "", xlab = "", col = male_al, border = male,
     breaks = 20, xlim = c(0.065, 0.31), ylab = "")
hist(sigmaF_ans, add = T, col = female_al, border = female, breaks = 10)

par(mar = c(2, 4, 4, 2))

hist(fit_apo$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Apodiformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_apo, main = "", xlab = "", col = male_al, border = male,
     breaks = 10, xlim = c(0.010, 0.025), ylab = "")
hist(sigmaF_apo, add = T, col = female_al, border = female, breaks = 20)

par(mar = c(2, 4, 4, 2))

hist(fit_cha$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Charadriiformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_cha, main = "", xlab = "", col = male_al, border = male,
     breaks = 5, xlim = c(0.06, 0.08), ylab = "")
hist(sigmaF_cha, add = T, col = female_al, border = female, breaks = 1)

par(mar = c(2, 4, 4, 2))

hist(fit_col$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Columbiformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_col, main = "", xlab = "", col = male_al, border = male,
     breaks = 20, xlim = c(0.02, 0.06), ylab = "")
hist(sigmaF_col, add = T, col = female_al, border = female, breaks = 10)

par(mar = c(2, 4, 4, 2))

hist(fit_gal$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Galliformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_gal, main = "", xlab = "", col = male_al, border = male,
     breaks = 5, xlim = c(0.02, 0.04), ylab = "")
hist(sigmaF_gal, add = T, col = female_al, border = female, breaks = 5)

par(mar = c(2, 4, 4, 2))

hist(fit_pas$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Passeriformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_pas, main = "", xlab = "", col = male_al, border = male,
     breaks = 20, xlim = c(0.01, 0.08), ylab = "")
hist(sigmaF_pas, add = T, col = female_al, border = female, breaks = 20)

par(mar = c(2, 4, 4, 2))

hist(fit_pic$pval, main = "", xlab = "", col = cols_al, border = cols,
     breaks = 20)
title("Piciformes", adj = 0)

par(mar = c(2, 4, 4, 2))

hist(sigmaM_pic, main = "", xlab = "", col = male_al, border = male,
     breaks = 10, xlim = c(0.02, 0.045), ylab = "")
hist(sigmaF_pic, add = T, col = female_al, border = female, breaks = 10)

par(mar = c(4, 4, 4, 2))

hist(fit_psi$pval, main = "", xlab = "p-value", col = cols_al, border = cols,
     breaks = 20)
title("Psittaciformes", adj = 0)

par(mar = c(4, 4, 4, 2))

hist(sigmaM_psi, main = "", xlab = "Sigma", col = male_al, border = male,
     breaks = 20, xlim = c(0.022, 0.035), ylab = "")
hist(sigmaF_psi, add = T, col = female_al, border = female, breaks = 10)

dev.off()

