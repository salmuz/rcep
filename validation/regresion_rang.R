################ Charger des données ##############
library(pls)
library(MASS)
source("load_data.R", chdir = T)
source("model/methods.R", chdir=T)
source("model/reg.pca.f1.R", chdir=T)
source("model/reg.pca.fn.R", chdir=T)

source("model/reg.ridge.R", chdir=T)

## Validation of the linear regresion collinear 
## Hypothesis E(e) = 0 et V(e) = sigma^2, e ~ N(0, sigma^2)
## EQM : Error Quadratique Moyenne 
#####################################################################
### Stage de Regression Composantes explicatives par penalisation ###
#####################################################################
set.seed(87)
cvseg <- cvsegments(nrow(cookie.app), k=4, type="random")

fit.eqm <- function(model.coef){
  prediction <- model.coef[1] + model.coef[-1] %*% t(cookie.val[,-1])
  eqm <- mean((cookie.val[, 1] - prediction)^2)
  print(paste("EQM regression:", eqm, sep = " "))
}

rs.fit <- lm(sucres ~. , data = cookie.app)
rs.coef <- coef(rs.fit)
rs.coef[is.na(rs.coef)]  <- 0
fit.eqm(rs.coef)

rid.fit <- lm.ridge(sucres ~. , data = cookie.app, lambda = 0.12)
rid.coef <- coef(rid.fit)
fit.eqm(rid.coef)

xbry.fit <- reg.regularization(cookie.app[,-1], cookie.app[, 1], s=0.2, Niter = 10)
xbry.coef <- reg.coef(xbry.fit)
fit.eqm(xbry.coef)

#####################################################################
rid.cv <- ridge(1, cookie.app , cvseg, 1000)
plot(seq(0,1,length=length(rid.cv$press)), rid.cv$press, ylab="Erreur Moyenne Carree (EMC)", 
     xlab="Parametre de penalisation (s)", type="o")
rid.fit <- lm.ridge(sucres~., data=cookie.app, lambda = 0.1081081)
cv.validation(cookie.val, coef(rid.fit))

reg.f1 <- reg.cv.smax(1, cookie.app, cvseg, 100)
plot(reg.f1$s.set, as.vector(reg.f1$press), ylab="Erreur Moyenne Carree (EMC)", 
     xlab="Parametre de penalisation (s)", type = "o")
reg.fit <- reg.regularization(cookie.app[,-1], cookie.app[,1], 0.002004008)
cv.validation(cookie.val, reg.coef(reg.fit))

#########################################################################
sq <- seq(0, .9, by=.1)
ee <- NULL
for(setting in sq){
  reg.fit <- reg.regularization(cookie.app[,-1], cookie.app[,1], setting)
  eqm <- cv.validation(cookie.val, reg.coef(reg.fit))
  print(eqm)
  ee <- c(ee, eqm)  
}
plot(sq, ee)

cof1 <- NULL
cof1 <- coef(rid.fit)
for(setting in seq(0, .3, by=.1)){
  reg.fit <- reg.regularization(cookie.app[,-1], cookie.app[,1], setting)
  cof1 <- cbind(cof1, reg.coef(reg.fit))
}
cof1 <- as.data.frame(cof1)

plot(cof1[, 1], type='l', lty=2, col="red")
lines(cof1[, 2], lty=3, col="blue")
abline(1,0, lty=2)
lines(cof1[, 3], lty=3, col="black")


# compute erreur quadratique
Xs <- (cookie.val[,-1] - matrix(as.vector(regression$Xm), n.row, n.col, byrow=T)) /
  matrix(as.vector(regression$Xsd), n.row, n.col, byrow=T)
##########################################################################
rg <- reg.pca.cv(cbind(Y,X), vt.comp, cvseg, n.setting=10)

plot(rg$s.set, as.vector(rg$press), ylab="Erreur Moyenne Carree (EMC)", 
     xlab="Parametre de penalisation (s)", type = "o")
for(setting in sq){
  reg.fit <- reg.pca.regression(X, Y, setting, vt.comp)
}

##########################################################################
