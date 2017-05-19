setwd("/data/scripts/mining")
library(pls)
library(MASS)

source("load_data.R", chdir = T)
source("model/methods.R", chdir=T)
source("model/reg.pca.f1.R", chdir=T)
source("model/reg.pca.fn.R", chdir=T)


# Test la convergence 
#xbry.fit <- reg.regularization(cookie.app[,-1], cookie.app[, 1], s=0.1221, Niter=1000)

# Test de validation EQM
#set.seed(87)
#cvseg <- cvsegments(nrow(cookie.app), k=4, type="random")
#reg.f1 <- reg.pca.f1.cv(1, cookie.app, cvseg, n.setting = 1000)
#save(reg.f1, file = "xbryf1.RData")         

#plot(reg.f1$s.set ,as.vector(reg.f1$press)/300, ylab="Erreur Moyenne Carree (EMC)",
#    xlab="Parametre de penalisation (s)")

# Test de nombres de composantes
pca <- reg.pca.regularization(cookie.app, nb.comp = 10, n.setting = 1000)

