rm(list = ls())

library(mvMORPH)
library(phytools)
library(geiger)

tr <- read.tree("tr200.tre")
dat <- read.csv("data_aves.csv", row.names = 1)

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


## Passeriformes = 2290 spp
fit_pas <- multiModelComp(phy = tr, data = subdat, taxon = "Passeriformes",
                          ntree = 100)
save(fit_pas, file = "fit_pas.RData")