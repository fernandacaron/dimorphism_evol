rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(mvMORPH)
library(phytools)
library(geiger)

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


## Aves = 4515 spp
fit_ave <- multiModelComp(phy = tr, data = subdat, ntree = 5)
#save(fit_ave, file = "data/aves/fit_ave.RData")

## Columbiformes = 105 spp
fit_col <- multiModelComp(phy = tr, data = subdat, taxon = "Columbiformes",
                          ntree = 100)
#save(fit_col, file = "data/aves/fit_col.RData")

## Psittaciformes = 131 spp
fit_psi <- multiModelComp(phy = tr, data = subdat, taxon = "Psittaciformes",
                          ntree = 100)
#save(fit_psi, file = "data/aves/fit_psi.RData")

## Anseriformes = 150 spp
fit_ans <- multiModelComp(phy = tr, data = subdat, taxon = "Anseriformes",
                          ntree = 100)
#save(fit_ans, file = "data/aves/fit_ans.RData")

## Accipitriformes = 169 spp
fit_acc <- multiModelComp(phy = tr, data = subdat, taxon = "Accipitriformes",
                          ntree = 100)
#save(fit_acc, file = "data/aves/fit_acc.RData")

## Galliformes = 179 spp
fit_gal <- multiModelComp(phy = tr, data = subdat, taxon = "Galliformes",
                          ntree = 100)
#save(fit_gal, file = "data/aves/fit_gal.RData")

## Piciformes = 217 spp
fit_pic <- multiModelComp(phy = tr, data = subdat, taxon = "Piciformes",
                          ntree = 100)
#save(fit_pic, file = "data/aves/fit_pic.RData")

## Apodiformes = 247 spp
fit_apo <- multiModelComp(phy = tr, data = subdat, taxon = "Apodiformes",
                          ntree = 100)
#save(fit_apo, file = "data/aves/fit_apo.RData")

## Charadriiformes = 252 spp
fit_cha <- multiModelComp(phy = tr, data = subdat, taxon = "Charadriiformes",
                          ntree = 100)
#save(fit_cha, file = "data/aves/fit_cha.RData")

## Passeriformes = 2290 spp
fit_pas <- multiModelComp(phy = tr, data = subdat, taxon = "Passeriformes",
                          ntree = 100)
#save(fit_pas, file = "data/aves/fit_pas.RData")

#load(file = "data/aves/fit_ave.RData")
#load(file = "data/aves/fit_col.RData")
#load(file = "data/aves/fit_psi.RData")
#load(file = "data/aves/fit_ans.RData")
#load(file = "data/aves/fit_acc.RData")
#load(file = "data/aves/fit_gal.RData")
#load(file = "data/aves/fit_pic.RData")
#load(file = "data/aves/fit_apo.RData")
#load(file = "data/aves/fit_cha.RData")
#load(file = "data/aves/fit_pas.RData")

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
sigmaM_ave <- sigmaM_acc <- sigmaM_ans <- sigmaM_apo <- sigmaM_cha <- 
  sigmaM_col <- sigmaM_gal <- sigmaM_pas <- sigmaM_pic <- sigmaM_psi <- 
  numeric()
sigmaF_ave <- sigmaF_acc <- sigmaF_ans <- sigmaF_apo <- sigmaF_cha <- 
  sigmaF_col <- sigmaF_gal <- sigmaF_pas <- sigmaF_pic <- sigmaF_psi <- 
  numeric()
for (i in 1:length(fit_acc)) {
  sigmaM_ave[i] <- fit_ave$fit[[i]][1, 3]
  sigmaM_acc[i] <- fit_acc$fit[[i]][1, 3]
	sigmaM_ans[i] <- fit_ans$fit[[i]][1, 3]
	sigmaM_apo[i] <- fit_apo$fit[[i]][1, 3]
	sigmaM_cha[i] <- fit_cha$fit[[i]][1, 3]
	sigmaM_col[i] <- fit_col$fit[[i]][1, 3]
	sigmaM_gal[i] <- fit_gal$fit[[i]][1, 3]
	sigmaM_pas[i] <- fit_pas$fit[[i]][1, 3]
	sigmaM_pic[i] <- fit_pic$fit[[i]][1, 3]
	sigmaM_psi[i] <- fit_psi$fit[[i]][1, 3]

	sigmaF_ave[i] <- fit_ave$fit[[i]][1, 4]
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

fit_birds <- data.frame(nspp = c(4515, 169, 150, 247, 252, 105, 179, 2290, 217, 131),
                        nsim = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
                        SDI = rep(NA, 10),
                        sigmaM = rep(NA, 10), 
                        sigmaF = rep(NA, 10), 
                        pval = rep(NA, 10),
                        row.names = c("Aves", "Accipitriformes", 
                                      "Anseriformes", "Apodiformes", 
                                      "Charadriiformes", "Columbiformes", 
                                      "Galliformes", "Passeriformes", 
                                      "Piciformes", "Psittaciformes"))

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

fit_birds[1, 3:6] <- c(summ(subdat_sdi$sdi, TRUE),
                       summ(sigmaM_ave, FALSE),
                       summ(sigmaF_ave, FALSE),
                       summ(fit_ave$pval, FALSE)
)

fit_birds[2, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order ==
                                             "Accipitriformes"], 
                            TRUE),
                       summ(sigmaM_acc, FALSE),
                       summ(sigmaF_acc, FALSE),
                       summ(fit_acc$pval, FALSE)
)

fit_birds[3, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Anseriformes"], 
                            TRUE),
                       summ(sigmaM_ans, FALSE),
                       summ(sigmaF_ans, FALSE),
                       summ(fit_ans$pval, FALSE)
)

fit_birds[4, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Apodiformes"], 
                            TRUE),
                       summ(sigmaM_apo, FALSE),
                       summ(sigmaF_apo, FALSE),
                       summ(fit_apo$pval, FALSE)
)

fit_birds[5, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == 
                                             "Charadriiformes"], 
                            TRUE),
                       summ(sigmaM_cha, FALSE),
                       summ(sigmaF_cha, FALSE),
                       summ(fit_cha$pval, FALSE)
)

fit_birds[6, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Columbiformes"], 
                            TRUE),
                       summ(sigmaM_col, FALSE),
                       summ(sigmaF_col, FALSE),
                       summ(fit_col$pval, FALSE)
)

fit_birds[7, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Galliformes"], 
                            TRUE),
                       summ(sigmaM_gal, FALSE),
                       summ(sigmaF_gal, FALSE),
                       summ(fit_gal$pval, FALSE)
)

fit_birds[8, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Passeriformes"], 
                            TRUE),
                       summ(sigmaM_pas, FALSE),
                       summ(sigmaF_pas, FALSE),
                       summ(fit_pas$pval, FALSE)
)

fit_birds[9, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == "Piciformes"], 
                            TRUE),
                       summ(sigmaM_pic, FALSE),
                       summ(sigmaF_pic, FALSE),
                       summ(fit_pic$pval, FALSE)
)

fit_birds[10, 3:6] <- c(summ(subdat_sdi$sdi[subdat_sdi$Order == 
                                              "Psittaciformes"], 
                             TRUE),
                        summ(sigmaM_psi, FALSE),
                        summ(sigmaF_psi, FALSE),
                        summ(fit_psi$pval, FALSE)
)

write.csv(fit_birds, "tables/Table1_unformatted.csv")


