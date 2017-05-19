# library dependance for computing the first composante
reg.pca.regularization <- function(sample, nb.comp = 1, n.setting = 10, 
                                   cv.seed = 87, cv.k = 4, idx.y = 1, setting.max = 1, 
                                   n.iter = 200 , ...){
  set.seed(cv.seed)
  cv.segments <- cvsegments(nrow(sample), k = cv.k, type="random")
  # matrix de directions de composantes 
  vt.comp <- cum.press <- cum.hat <- cum.coef <- NULL;
  
  # compute the first component
  cv.f1 <- reg.pca.f1.cv(setting.max, sample, cv.segments, n.setting, idx.y)
  cum.press <- cbind(cum.press, cv.f1$press)
  cum.hat <- cbind(cum.hat, cv.f1$s.hat)
  reg.fit <- reg.pca.f1.regression(sample[,-idx.y], sample[, idx.y], cv.f1$s.hat, n.iter)
  vt.comp <- cbind(vt.comp, reg.fit$u)
  cum.coef <- cbind(cum.coef, reg.pca.f1.coef(reg.fit))
  
  # compute the k-1 component
  if(nb.comp > 1)
    for(i in 2:nb.comp){
      cv.pca <- reg.pca.cv(sample, vt.comp, cv.segments, n.setting, setting.max, idx.y)
      cum.press <- cbind(cum.press, cv.pca$press)
      cum.hat <- cbind(cum.hat, cv.pca$s.hat)
      reg.fit <- reg.pca.regression(sample[,-idx.y], sample[, idx.y], cv.pca$s.hat, vt.comp, n.iter)
      vt.comp <- cbind(vt.comp, reg.fit$u)
      cum.coef <- cbind(cum.coef, reg.pca.f1.coef(reg.fit))
    }
  
  return(list(vt.comp = vt.comp, cum.press = cum.press, 
              cum.hat = cum.hat, cum.coef = cum.coef))
}


###############################################################
## cv.segments : des segments aléatoires pour la validation croissé 
## cv.sample : échantillon à utiliser dans la validation croissé 
## idx.y : indice de la variable à expliquer 
###############################################################
reg.pca.cv <- function(cv.sample, vt.comp, cv.segments, n.setting = 10, 
                       setting.max = 1, idx.y = 1, ...){
  set.setting <- seq(0, setting.max, length = n.setting+2)
  set.setting <- set.setting[c(-1,-length(set.setting))]
  press <- rep(0, length(set.setting))
  for(i in 1:length(cv.segments)){
    sample.valid <- cv.sample[unlist(cv.segments[i]), ]
    X <- cv.sample[unlist(cv.segments[-i]), -idx.y]
    Y <- cv.sample[unlist(cv.segments[-i]), idx.y]
    coeff <- NULL
		for(j in 1:length(set.setting)){
			fit <- reg.pca.regression(X, Y, set.setting[j], vt.comp)
			coeff <- rbind(coeff, reg.pca.coef(fit)) 
		}
    # Chaque cellule de la matrice est une prevision + interception
		prediction <- matrix(coeff[,1], nrow(coeff), nrow(sample.valid)) + coeff[,-1] %*% t(sample.valid[,-1])
		press <- press + rowSums((matrix(sample.valid[,1], nrow(coeff), nrow(sample.valid), byrow=T)-prediction)^2)
  }
  s_hat <- set.setting[which.min(press)]
	return(list(s.hat=s_hat, press=press, s.set=set.setting))
}

###############################################################
# Recovery the coefficient of the fitted model
# @param model.fit : model fitting
# @return coeff : the coeff with interception 
###############################################################
reg.pca.coef <- function(model.fit){
  coeff <- model.fit$gamma * model.fit$u
  comps <- model.fit$vt.comp %*% model.fit$delta
  scalecoeff <- t(as.matrix(coeff/model.fit$Xsd))
  scalecomps <- t(as.matrix(comps/model.fit$Xsd))
  intercept <- model.fit$Ym - (scalecoeff + scalecomps) %*% model.fit$Xm
  return(drop(cbind(intercept, scalecoeff + scalecomps)))
}

###############################################################
# @param X : variables explicatives centre et reduites
# @param Y : variables a expliquer centre 
# @param s : paramètre de reglage réel ou un vecteur des réels 
###############################################################
# v.comp est le vector de directions de composantes principales (p,k) {u_1, u_2, ..., u_k}
reg.pca.regression <- function(var.x, var.y, setting = 0.5, vt.comp, Niter = 100, ...){
  scale <- centre_reduit(as.matrix(var.x), var.y)
  X <- scale$X
  Y <- scale$Y
  n <- scale$n 
  p <- scale$p
  vt.comp <- as.matrix(vt.comp)
  
  # initialisation des variables 
  u <- matrix(10^-20, ncol = Niter, nrow = p)
  W <- diag(x = 1/n, n, n)
  N <- t(X) %*% W %*% X
  Tc <- X %*% vt.comp
  tA <- t(vt.comp) %*% N #  Ak'= uk'X'WX 
  iAA <- ginv(tA %*% t(tA)) # A'A
  
  # svd: singular value decomposition 
  Xsvd <- svd(X) #  X <- svd$u %*% diag(svd$d) %*% t(svd$v)  
  d2 <- Xsvd$d ^ 2 
  n.diag <- length(Xsvd$d) 
  Nd2 <- rep(W[1,1], n.diag) * d2 # DW*D , W* = UWU'
  Nvu <- diag(Nd2, n.diag, n.diag) %*% t(Xsvd$v) %*% vt.comp
  du <- diag(Xsvd$d, n.diag, n.diag) %*% t(Xsvd$u)
  duY <- du %*% Y
  
  # one stages 
  one.stage <- function(k = 2){
    Z <- cbind(X %*% u[,k-1], Tc)
    Zsvd <- svd(Z)
    theta <- Zsvd$v %*% ((t(Zsvd$u) %*% Y) / Zsvd$d)
    sigma <-  norm(Y - Z %*% theta, "2")^2 / n
    return(list(gamma = theta[1], delta = theta[-1], sigma=sigma))
  }
  
  two.stage <- function(k = 2, gamma, delta, sigma){
    zeta0 <- (setting-1)*(gamma^2)
    zeta1 <- (setting-1)*gamma
    
    Nu <- N %*% u[,k-1]
    uNu <- t(u[,k-1]) %*% Nu
    
    Xu <- X %*% u[, k-1]
    uX <-  t(Xu)
    TDelta <- Tc %*% delta
    
    #Delta <- uX %*% Y - gamma * (norm(Xu, "2")^2) - uX %*% TDelta
    #print(paste("Delta",Delta,sep = ":"))
    #lambda <- drop(((zeta1 * Delta)/(2*sigma)) + setting*uNu)
    lambda <- drop(setting*uNu)
    
    cst <- (-1*zeta1/sigma)
    first <- cst * t(X) %*% (Y - gamma*Xu - TDelta) 
    tau <- iAA %*% tA %*% (first + 2*setting*Nu)
    
    Sigma <- zeta0 * d2 + 2*sigma*(setting*Nd2 - lambda)
    sRight <- zeta1 * (duY + du %*% TDelta) + sigma * Nvu %*% tau
    uk <- sRight / Sigma
    dim(uk) <- c(n.diag, 1)
    u[,k] <- Xsvd$v %*% uk
    u[,k] <- u[,k]/norm(u[,k], "2")
    return(u[,k])
  }
  
  u[,1] <- u[,1]/norm(u[,1], "2")
  params <- one.stage(2)
  u[,2] <- two.stage(2, params$gamma, params$delta, params$sigma)
  i <- 2
  condition <- (t(u[,i-1])%*%u[,i])^2
  converge <- c(NULL, condition)
  #print(paste("Method avec SVD avec setting", setting, sep = ":"));
  while(condition < (1 - 10^-6) && i < Niter){
    i <- i+1
    params <- one.stage(i)
    u[,i] <- two.stage(i, params$gamma, params$delta, params$sigma)
    condition <- (t(u[,i-1])%*%u[,i])^2
    converge <- c(converge, condition)
  }
  params <- one.stage(i+1)
  return(list(u = u[,i], converge = converge, gamma = params$gamma, sigma=params$sigma, 
              delta=params$delta, vt.comp = vt.comp, 
              Xm=scale$Xm, Ym=scale$Ym, Xsd=scale$Xsd));
}

###############################################################
## Test of the reg.regularization function 
## with cookies data and 3 setting values
###############################################################
reg.pca.test <- function(){
  sets <- c(0.1, 0.5, 0.8)
  reg.fit <- reg.pca.f1.regression(cookie.app[,-1], cookie.app[,1], 0.002004008)
  for(j in 1:length(sets)){
    resp <- reg.pca.regression(cookie.app[,-1], cookie.app[,1], sets[j], reg.fit$u)
    cat(paste("resp[", j, "] : ( gamma=", resp$gamma,", sigma=", resp$sigma, ")", sep = ""), "\n")
  }
}