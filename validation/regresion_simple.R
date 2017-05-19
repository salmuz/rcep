################ Charger des données ##############
library(pls)
library(MASS)
source("load_data.R", chdir = T)
source("model/methods.R", chdir=T)
source("model/reg.pca.f1.R", chdir=T)
source("model/reg.ridge.R", chdir=T)

data("trees")
## Validation de regression linéaire simple avec les
## deux hypothèse E(e) = 0 et V(e) = sigma^2.
## EQM : Error Quadratique Moyenne 

################ Regression simple ################
rs.fit <- lm(O3 ~ ., data = ozone_long)
rs.coef <- rs.fit$coefficients
rs.eqm <- mean((ozone_long$O3 - rs.coef[1] - 
                  as.matrix(ozone_long[, -1]) %*% as.matrix(rs.coef[-1]))^2)
cat(paste("Simple, EQM:", rs.eqm, seq = " "))

#summary(rs.fit)
#confint(rs.fit,level=0.95) 

################ Regression Ridge ##################
ridge.fit <- function(formule, sample, lambda = 0, idy = 1){
  rid.fit <- lm.ridge(formule , data = sample, lambda = lambda)
  rid.coef <- coef(rid.fit)
  rid.eqm <- mean((sample[, idy] - rid.coef[1] - 
                    as.matrix(sample[, -idy]) %*% as.matrix(rid.coef[-1]))^2)
  print(paste("Ridge, EQM:", rid.eqm, ", lambda:", lambda, seq = " "))
  return(list(fit=rid.fit, eqm=rid.eqm))
}
rfit <- ridge.fit(O3 ~., ozone_long)
coef(rfit$fit)

rfit <- ridge.fit(Volume ~., trees, idy = 3)
coef(rfit$fit)

################ Regression Xbry ##################
xbry.fit <- function(sample, setting = 0, niter = 200, idy = 1){
  X <- as.matrix(sample[, -idy])
  Y <- as.matrix(sample[, idy])
  xb.fit <- reg.regularization(X, Y, s=setting, niter)
  xb.coef <- reg.coef(xb.fit)
  xb.eqm <- mean((sample[, idy] - xb.coef[1] - 
                      as.matrix(sample[, -idy]) %*% as.matrix(xb.coef[-1]))^2)
  print(paste("Xbry, EQM:", xb.eqm, ", setting:", setting, seq = " "))
  return(list(fit=xb.fit, eqm=xb.eqm))
}
xfit <- xbry.fit(ozone_long)
reg.coef(xfit$fit)

xfit <- xbry.fit(trees, setting = 0,  idy = 3)
reg.coef(xfit$fit)

## EQM evolution dans le temps 
pt <- NULL
se <- seq(0, .99, by=.001)
for(st in se){
  xfit <- xbry.fit(ozone_long, setting = st,  idy = 3)
  pt <- c(pt, xfit$eqm)
}
plot(se, pt, type="o")

################ Regression rang(X) > n ###############
ozone_rang <- ozone_long[1:5,]
fit <- xbry.fit(ozone_rang, .1, 2000)
plot(1:length(fit$fit$converg), fit$fit$converg, type = "l")
for(lambda in seq(.1, .9, by=.1)){
  fr <- ridge.fit(O3 ~., ozone_rang, lambda)
  fx <- xbry.fit(ozone_rang, lambda)
  print(paste("EQM Ridge > Xbry", fr$eqm > fx$eqm,  sep = " "))
}


