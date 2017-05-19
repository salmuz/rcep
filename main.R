## Projet stage Master 1
setwd("/Users/salmuz/Documents/git/stage")
## Regression avec Ridge 
library(pls)
library(MASS)
## Regression methods 
source("reg_methods.R", chdir = T)
source("load_data.R", chdir = T)
source("reg_pen_1acp.R", chdir=T)
# Graine du generateur fixee 
set.seed(87)
## Creation des segments aleatoires pour la validation croisée
cvseg <- cvsegments(nrow(cookie.app), k=4, type="random")
n.row <- nrow(cookie.val[,-1])
n.col <- ncol(cookie.val[,-1])

#####################################################################
######## Tests Regression avec Penalisation(Ridge, PCA, LSA) ########
#####################################################################

# Calculer l'erreur de prevision
coeff <- coef(lm.ridge(sucres~.,data = cookie.app,lambda=kappaet))
n.val <- nrow(cookie.val[,-1])
prediction <- as.vector(rep(coeff[1],n.val)) + as.vector(coeff[-1]%*%t(data.matrix(cookie.val[,-1])))
mean((cookie.val[,1]-prediction)^2)

modele.lm=lm(sucres~.,data = cookie.app)
mean((cookie.val[,1]-predict(modele.lm,newdata=cookie.val))^2)

# Plot la convergence de k 
coefflm  <- coef(lm.ridge(sucres~., data=cookie.app, lambda=0))
coeff <- coef(lm.ridge(sucres~.,data = cookie.app,lambda=kappaet))
matplot(t(rbind(coeff, coefflm)), type="l", col=1)

###########################################################################

regression <- reg.regularization(cookie.app[,-1], cookie.app[,1], resp$s_hat)
C <- as.matrix(regression$u) 
response <- reg.cv.smaxacp(1, cvseg, 100, C)
plot(response$set.s, as.vector(response$press), ylab="Erreur Moyenne Carree (EMC)", 
     xlab="Parametre de penalisation (s)")


reg.acp <- reg.acp_reglarization(cookie.app[,-1], cookie.app[,1], response$s_hat, C)
n.val <- nrow(cookie.val[,-1])
Xs <- (cookie.val[,-1] - matrix(as.vector(reg.acp$Xm), nrow(cookie.val), ncol(cookie.val[,-1]), byrow=T)) /
  matrix(as.vector(reg.acp$Xsd), nrow(cookie.val), ncol(cookie.val[,-1]), byrow=T)
Cp <- as.matrix(Xs) %*% as.matrix(reg.acp$u) * reg.acp$delta
prediction <- as.vector(rep(reg.acp$Ym, nrow(cookie.val))) + 
  as.vector(reg.acp$gamma*reg.acp$u %*% t(data.matrix(Xs))) + as.vector(Cp)
mean((cookie.val[,1]-prediction)^2)

save(response, file = "Regf1Varbaise.RData")
