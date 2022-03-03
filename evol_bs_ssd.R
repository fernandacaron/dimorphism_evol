rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(phytools)
library(geiger)
library(ratematrix)

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

## Ref Blanckenhorn et al. 2007 (Livro Fairbairn et al. 2007, Chapter 6, p. 65)
## SDI = (female size/male size - 1) when female larger
## SDI = -(male size/female size - 1) when male size. Olhar:
## We use the index of Lovich and Gibbons (1992) to express SSD as (length of 
## larger sex/length of smaller sex)–1, negative by convention when males are 
## the larger sex and positive when females are larger than males.

SDI <- function (male = male, female = female, cutoff = 0.1) {
	if (max(male, female) >= (min(male, female)*cutoff) + min(male, female)){
		if (male <= female) {
				SDI <- (female/male) - 1
			} else {
				if (female < male) {
					SDI <- -((male/female) - 1)
				} 
			}

		return(SDI)

	} else {
		SDI <- 0
	}

	return(SDI)

}

dat_red <- dat[, c(5, 9, 13)]
dat_red <- dat_red[complete.cases(dat_red$Body_mass_g_M_mean) & 
				   complete.cases(dat_red$Body_mass_g_F_mean), ]

SSD <- numeric()
for (i in 1:nrow(dat_red)) {
	SSD[i] <- SDI(male = dat_red$Body_mass_g_M_mean[i],
	              female = dat_red$Body_mass_g_F_mean[i], cutoff = 0)
}
names(SSD) <- dat_red$Scientific_name
SSD <- SSD[complete.cases(SSD)]

body_size <- numeric()
for (i in 1:nrow(dat_red)) {
	body_size[i] <- mean(c(dat_red$Body_mass_g_M_mean[i], 
	                       dat_red$Body_mass_g_F_mean[i]))
}
names(body_size) <- dat_red$Scientific_name 

for (i in 1:1000) tr[[i]] <- treedata(tr[[i]], SSD, warnings = F)$phy

##ratebystate
## Acho que tem que fazer pro grau geral de SSD, grau de male-biased SSD e 
## de female-biased SSD
SSD_g <- abs(SSD)
SSD_g[SSD_g == 0] <- 0.00000000000000001
rate_state_g <- list()
for (i in 1:100) {
	rate_state_g[[i]] <- ratebystate(tr[[i]], log(body_size), log(SSD_g), 
	                                 nsim = 10)
}

SSD_m <- SSD
SSD_m[SSD_m > 0] <- 0
SSD_m <- abs(SSD_m)
rate_state2 <- list()
for (i in 1:1000) {
	rate_state2[[i]] <- ratebystate(tr[[i]], log(body_size), log(SSD_m), 
	                                nsim = 1)
}

SSD_f <- SSD
SSD_f[SSD_f < 0] <- 0
SSD_f <- abs(SSD_f)
rate_state2 <- list()
for (i in 1:1000) {
	rate_state2[[i]] <- ratebystate(tr[[i]], log(body_size), log(SSD_f), 
	                                nsim = 1)
}






#ratematrix
xx <- data.frame(SSD = SSD, Body_size = body_size)
rownames(xx) <- names(SSD)
xx$Body_size <- log(xx$Body_size)
xx$SSD <- abs(xx$SSD)
xx$SSD[xx$SSD == 0] <- 0.00000000000000001
xx$SSD <- log(xx$SSD)

pred <- numeric()
for (i in 1:length(SSD)) {
	pred[i] <- ifelse(SSD[i] > 0, 2, ifelse(SSD[i] < 0, 1, 0))
}
names(pred) <- names(SSD)

phy <- make.simmap(tr[1:100], pred)

handle <- ratematrixMCMC(data = xx, phy = phy, prior = "uniform", 
                         gen = 1000000, dir = "data/aves")
mcmc <- readMCMC(handle, burn = 0.25, thin = 1)
logAnalyzer(handle, burn = 0.25, thin = 1) 

## 0(Ñ) = #002244 (azul); 1(M) = #69BE28 (verde); 2(F) = #A5ACAF (cinza)
plotRatematrix(chain = mcmc) 
plotRootValue(chain = mcmc) 
checkConvergence(anolesPost$chain1, anolesPost$chain2)
testRatematrix(chain = mcmc, par = "correlation")
testRatematrix(chain = mcmc, par = "rates")



