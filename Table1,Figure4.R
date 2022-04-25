rm(list = ls())

setwd("Documents/Lab/dimorph_evol")

library(mvMORPH)
library(phytools)
library(geiger)

tr <- read.nexus("data/aves/aves_Ericson_VertLife_27JUL20.nex")
dat <- read.csv("data/aves/BodySizeAves_18jan22.csv", row.names = 1)

subdat <- dat[complete.cases(dat$Body_mass_g_M_mean) & 
               complete.cases(dat$Body_mass_g_F_mean), 
               c("Body_mass_g_M_mean", "Body_mass_g_F_mean")]
rownames(subdat) <- dat$Scientific_name[complete.cases(dat$Body_mass_g_M_mean) &
										complete.cases(dat$Body_mass_g_F_mean)]

subdat$Body_mass_g_M_mean <- log(subdat$Body_mass_g_M_mean)
subdat$Body_mass_g_F_mean <- log(subdat$Body_mass_g_F_mean)

modelComp<-function(tr, data) {
	xx <- treedata(tr, data, warnings = FALSE)
	tree <- xx$phy
	data <- xx$data
	# BMM multi-rate/multi-selective regimes, BM1 unique rate of evolution/trait
	fit1 <- mvBM(tree, data, model = "BM1", optimization = c("subplex"))
	# Constraint faz que com sigma fique igual para os dois caracteres
	fit2 <- mvBM(tree, data, model = "BM1", param = list(constraint = TRUE),
	             optimization = c("subplex"))
	res <- AIC(fit2) -AIC(fit1)
	res
}

t0<-Sys.time()
res <- numeric()
for (i in 1:10) {
	a <- Sys.time()
	res[i] <- modelComp(tr[[i]], subdat)
	b <- Sys.time()
	print(b-a)
}
t1<-Sys.time()
t1-t0

write.csv(res, "data/aves/rates.csv")

phyloSig<-function(tr, data) {
	xx <- treedata(tr, data, warnings = FALSE)
	tree <- xx$phy
	data <- xx$data
	M <- phylosig(tree, data[,"Body_mass_g_M_mean"], method = "lambda", 
	              test = TRUE)
	F <- phylosig(tree, data[,"Body_mass_g_F_mean"], method = "lambda", 
	              test = TRUE)
	res <- c(M$lambda, M$logL, M$logL0, M$P, F$lambda, F$logL, F$logL0, F$P)
	res
}

t0<-Sys.time()
resPhyl <- matrix(ncol = 8, nrow = 1000)
colnames(resPhyl) <- c("lambdaM", "logLM", "logL0M", "PM", "lambdaF", "logLF", 
                       "logL0F", "PF")
for(i in 1:10) resPhyl[i, ] <- phyloSig(tr[[i]], dat)
t1<-Sys.time()
t1-t0

pdf("figures/Figure4.pdf")

hist(as.data.frame(resPhyl)$lambdaF, breaks = seq(0, 1, 0.05), 
     col = rgb(255, 0, 0, 60, maxColorValue = 255), main = "", xlab = "lambda")
hist(as.data.frame(resPhyl)$lambdaM, breaks = seq(0, 1, 0.05), add = TRUE, 
     col = rgb(0, 0, 255, 60, maxColorValue = 255))

dev.off()

https://www.notion.so/e7299c70ceb8475a8fc9f0db6165e48f?v=588e5fcb320a4b51a41057cf27da0e4c