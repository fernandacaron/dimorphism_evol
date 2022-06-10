rm(list = ls())

setwd("Documents/lab/dimorph_evol")

library(phytools)
library(geiger)
library(plotrix)

########## MK + DR ##########

dat <- read.csv("data/aves/BodySizeAves_30may22_edit.csv", row.names = 1)
tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
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

    pruned_tr <- treedata(phy, subsdi, warnings = F)$phy
  	subsdi <- subsdi[names(subsdi) %in% pruned_tr$tip.label]
  	subsdi <- subsdi[match(pruned_tr$tip.label, names(subsdi))]  
	
	fitSYM <- fitMk(pruned_tr, subsdi, model = "SYM", pi = "fitzjohn")
	fitARD <- fitMk(pruned_tr, subsdi, model = "ARD", pi = "fitzjohn")
	aic.res <- AIC(fitSYM, fitARD)

	res <- list()
	res$ratesSYM <- fitSYM$rates
	res$ratesARD <- fitARD$rates
	res$aic <- aic.res$AIC
	names(res$aic) <- c("aicSYM", "aicARD")
	res$aicw <- aic.w(aic.res$AIC)
	names(res$aicw) <- c("weigthSYM", "weightARD")

	if (length(fitSYM$states) == 3) {
		names(res$ratesSYM) <- c("rate_0/-1", "rate_-1/1", "rate_0/1")
		names(res$ratesARD) <- c("rate_0to-1", "rate_1to-1", "rate_-1to0", 
		                         "rate_1to0", "rate_-1to1", "rate_0to1")
	} else {
		names(res$ratesSYM) <-  "rate"
		names(res$ratesARD) <- c("rate_1to-1", "rate_-1to1")
	}

	res$fitSYM <- fitSYM
	res$fitARD <- fitARD

	res
}

# Análises feitas no cluster da UFG

fitMK_col <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Columbiformes")
#save(fitMK_col, file = "data/aves/results/fitMK_col.RData")

fitMK_psi <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Psittaciformes")
#save(fitMK_psi, file = "data/aves/results/fitMK_psi.RData")

fitMK_ans <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Anseriformes")
#save(fitMK_ans, file = "data/aves/results/fitMK_ans.RData")

fitMK_acc <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Accipitriformes")
#save(fitMK_acc, file = "data/aves/results/fitMK_acc.RData")

fitMK_gal <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Galliformes")
#save(fitMK_gal, file = "data/aves/results/fitMK_gal.RData")

fitMK_pic <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Piciformes")
#save(fitMK_pic, file = "data/aves/results/fitMK_pic.RData")

fitMK_apo <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Apodiformes")
#save(fitMK_apo, file = "data/aves/results/fitMK_apo.RData")

fitMK_cha <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Charadriiformes")
#save(fitMK_cha, file = "data/aves/results/fitMK_cha.RData")

fitMK_pas <- lapply(tr, fitMK, sdi = sdi, data = dat, taxon = "Passeriformes")
#save(fitMK_pas, file = "data/aves/results/fitMK_pas.RData")

#load(file = "data/aves/results/fitMK_col.RData")
#load(file = "data/aves/results/fitMK_psi.RData")
#load(file = "data/aves/results/fitMK_ans.RData")
#load(file = "data/aves/results/fitMK_acc.RData")
#load(file = "data/aves/results/fitMK_gal.RData")
#load(file = "data/aves/results/fitMK_pic.RData")
#load(file = "data/aves/results/fitMK_apo.RData")
#load(file = "data/aves/results/fitMK_cha.RData")
#load(file = "data/aves/results/fitMK_pas.RData")

aicw_res <- matrix(nrow = 9, ncol = 2)
colnames(aicw_res) <- c("SYM", "ARD")
rownames(aicw_res) <- c("Accipitriformes", "Anseriformes", "Apodiformes",
                        "Charadriiformes", "Columbiformes", "Galliformes",
                        "Passeriformes", "Piciformes", "Psittaciformes")
stats1 <- function(x) {
	aicw_sym <- aicw_ard <- numeric()
	for (i in 1:100) {
		aicw_sym[i] <- x[[i]]$aicw[[1]]
		aicw_ard[i] <- x[[i]]$aicw[[2]]
	}
	me_sym <- round(mean(aicw_sym), 3)
	me_ard <- round(mean(aicw_ard), 3)
	range_sym <- range(aicw_sym)
	range_ard <- range(aicw_ard)
	return(c(paste0(me_sym, " (", round(range_sym[1], 3), "-", 
	                round(range_sym[2], 3), ")"),
	       paste0(me_ard, " (", round(range_ard[1], 3), "-", 
	              round(range_ard[2], 3), ")")))
}

aicw_res[1, ] <- stats1(fitMK_acc)
aicw_res[2, ] <- stats1(fitMK_ans)
aicw_res[3, ] <- stats1(fitMK_apo)
aicw_res[4, ] <- stats1(fitMK_cha)
aicw_res[5, ] <- stats1(fitMK_col)
aicw_res[6, ] <- stats1(fitMK_gal)
aicw_res[7, ] <- stats1(fitMK_pas)
aicw_res[8, ] <- stats1(fitMK_pic)
aicw_res[9, ] <- stats1(fitMK_psi)

write.csv(aicw_res, "tables/Table1_unformatted.csv")

stats2 <- function(x, taxon) {

	sym <- c("Anseriformes")

	if (taxon %in% sym) {
		stats <- list()
		for (j in 1:length(x[[1]]$ratesSYM)) {	
			stats[[j]] <- numeric()	
			for (i in 1:100) {
				stats[[j]][i] <- x[[i]]$ratesSYM[[j]]
			}
		}
		names(stats) <- c(names(x[[1]]$ratesSYM))

		return(stats)
	} else {
		stats <- list()
		for (j in 1:length(x[[1]]$ratesARD)) {	
			stats[[j]] <- numeric()	
			for (i in 1:100) {
				stats[[j]][i] <- x[[i]]$ratesARD[[j]]
			}
		}
		names(stats) <- c(names(x[[1]]$ratesARD))
	
		return(stats)
	}
}

rates_acc <- stats2(fitMK_acc, "Accipitriformes")
rates_ans <- stats2(fitMK_ans, "Anseriformes")
rates_apo <- stats2(fitMK_apo, "Apodiformes")
rates_cha <- stats2(fitMK_cha, "Charadriiformes")
rates_col <- stats2(fitMK_col, "Columbiformes")
rates_gal <- stats2(fitMK_gal, "Galliformes")
rates_pas <- stats2(fitMK_pas, "Passeriformes")
rates_pic <- stats2(fitMK_pic, "Piciformes")
rates_psi <- stats2(fitMK_psi, "Psittaciformes")

# definir cores
male <- "#9966FF"
monom <- "gray"
female <- "#E69F00"

pdf("figures/Figure4.pdf", height = 9, width = 9)

layout(matrix(1:9, ncol = 3, byrow = T))

par(mar = c(2, 2, 2, 2))
cols2 <- setNames(c(male, female), c("Male-biased SSD", "Female-biased SSD"))

plot_acc <- plot(fitMK_acc[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_acc$x, plot_acc$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_acc$x, plot_acc$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Accipitriformes", adj = 0, line = -1)

plot_ans <- plot(fitMK_ans[[1]]$fitSYM, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_ans$x, plot_ans$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_ans$x, plot_ans$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Anseriformes", adj = 0, line = -1)

plot_apo <- plot(fitMK_apo[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_apo$x, plot_apo$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_apo$x, plot_apo$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Apodiformes", adj = 0, line = -1)

plot_cha <- plot(fitMK_cha[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_cha$x, plot_cha$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_cha$x, plot_cha$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Charadriiformes", adj = 0, line = -1)

plot_col <- plot(fitMK_col[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_col$x, plot_col$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_col$x, plot_col$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Columbiformes", adj = 0, line = -1)

plot_gal <- plot(fitMK_gal[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_gal$x, plot_gal$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_gal$x, plot_gal$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Galliformes", adj = 0, line = -1)

plot_pas <- plot(fitMK_pas[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_pas$x, plot_pas$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_pas$x, plot_pas$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Passeriformes", adj = 0, line = -1)

plot_pic <- plot(fitMK_pic[[1]]$fitARD, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_pic$x, plot_pic$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_pic$x, plot_pic$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Piciformes", adj = 0, line = -1)

plot_psi <- plot(fitMK_psi[[1]]$fitSYM, show.zeros = FALSE, spacer = 0.2, 
                 mar = rep(0.5,4), cex.rates = 1.2)
invisible(mapply(draw.circle, plot_psi$x, plot_psi$y, col = cols2, border = NA,
                 MoreArgs = list(radius = 0.13)))
text(plot_psi$x, plot_psi$y, c("M", "F"), cex = 1.3, col = "white", 
     font = 2)
title("Psittaciformes", adj = 0, line = -1)

dev.off()

##################

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

plot(sdi ~ subset_es[[1]], xlab = "DR statistic", ylab = "SDI", main = "",
     pch = 16, col = cols)
for (i in 2:100) {
	points(sdi ~ subset_es[[i]], pch = 16, col = cols)
}

legend("bottomright", pch = 16, bty = 'n', 
       col = c(male, monom, female), 
       legend  = c("Male-biased SSD", "Monomorphism", "Female-biased SSD"))

dev.off()


