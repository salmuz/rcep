###############################################################
# @param var.exp : variables explicatives 
# @param var.reg : variables a expliquer
# @param s : paramètre de reglage
###############################################################
reg.regularization <- function(var.exp, var.reg, s = 0.5, Niter = 200, ...){
  scale <- centre_reduit(as.matrix(var.exp), var.reg)
  n <- scale$n 
  p <- scale$p
  u <- matrix(10^-9, ncol=Niter, nrow=p)
  W <- diag(x = 1/n, n, n)
  N <- t(scale$X) %*% W %*% scale$X
  tXX <- t(scale$X) %*% scale$X
  tXY <- t(scale$X) %*% scale$Y
  
  # compute parameters gamma and sigmax
  compute_paremeters <- function(k = 1){
    Xu <- scale$X %*% u[,k-1] 
    uXY <- t(Xu) %*% scale$Y # <Xu|Y>
    Xu.norm <- norm(Xu, "2")^2  # <Xu|Xu>
    uNu <- t(u[,k-1]) %*% N %*% u[,k-1]  
    gamma <- as.numeric(uXY / Xu.norm)
    sigma <- (norm(scale$Y - gamma*Xu, "2")^2)/n
    lambda <- as.numeric(s*uNu)
    return(list(lambda=lambda, sigma=sigma, gamma=gamma));
  }
  
  # compute the u vector with parameters
  compute_u_normed <- function(k, lambda, sigma, gamma){
    Sigma <- ((s-1)*gamma^2* tXX + 2*sigma*(s*N - diag(lambda, p, p)))
    # Compute la pseudo inverse de la matrice Sigma
    u[,k] <- gamma*(s-1)*(ginv(Sigma) %*% tXY)
    u[,k] <- u[,k]/norm(u[,k], "2")
    return(u[,k])
  }
  
  u[,1] <- u[,1]/norm(u[,1], "2")
  params <- compute_paremeters(2)
  u[,2] <- compute_u_normed(2, params$lambda, params$sigma, params$gamma)
  i <- 2
  condition <- (t(u[,i-1])%*%u[,i])^2
  converge <- c(NULL, condition)
  print(paste("Method avec ginv pseudo inverse avec setting", s, sep = ":"));
  while( condition < (1 - 10^-6)  && i < Niter){
    i <- i+1
    params <- compute_paremeters(i)
    u[,i] <- compute_u_normed(i, params$lambda, params$sigma, params$gamma)
    condition <- (t(u[,i-1])%*%u[,i])^2
    converge <- c(converge, condition)
  }
  params <- compute_paremeters(i+1)
  return(list(u = u[,i], converg = converge,
              gamma = params$gamma, sigma=params$sigma, 
              Xm=scale$Xm, Ym=scale$Ym, Xsd=scale$Xsd));
}

###############################################################
# Recovery the coefficient of the fitted model
# @param model.fit : model fitting
# @return coeff : the coeff with interception 
###############################################################
reg.coef <- function(model.fit){
  coeff <- model.fit$gamma * model.fit$u
  scalecoeff <- t(as.matrix(coeff/model.fit$Xsd))
  intercept <- model.fit$Ym - scalecoeff %*% model.fit$Xm
  return(drop(cbind(intercept, scalecoeff)))
}

###############################################################
# Compute the best setting value for the fitted model
# @param s.max : setting parametre
# @param cvseg : crossvalidation segment
# @param nbe : number of settings 
###############################################################
reg.cv.smax <- function(s.max=1, sample, cvseg, nbe=10, idy = 1, ...){
  set_s <- seq(0, s.max, length=nbe)
  set_s <- set_s[c(-1, -length(set_s))] # delete s=(0,1)
  press <- rep(0, length(set_s))
  for(i in 1:length(cvseg)){
    valid <- sample[unlist(cvseg[i]),]
    coeff <- NULL
    for(j in 1:length(set_s)){
      fit <- reg.regularization(sample[unlist(cvseg[-i]),-idy], 
                                sample[unlist(cvseg[-i]),idy], set_s[j])
      coeff <- rbind(coeff, reg.coef(fit)) # 0 .... 1 
    }
    # Chaque cellule de la matrice est une prevision + interception
    prediction <- matrix(coeff[,1], nrow(coeff), nrow(valid)) + coeff[,-1] %*% t(valid[,-1])
    press <- press + rowSums((matrix(valid[,1], nrow(coeff), nrow(valid), byrow=T)-prediction)^2)
  }
  s_hat <- set_s[which.min(press)]
  return(list(s_hat=s_hat, press=press, s.set=set_s))
}	

###############################################################
## Test of the reg.regularization function 
## with cookies data and 3 setting values
###############################################################
test <- function(){
  sets <- c(0.1, 0.5, 0.8)
  for(j in 1:length(sets)){
    resp <- reg.regularization(cookie.app[,-1], cookie.app[,1], sets[j])
    cat(paste("resp[", j, "] : ( gamma=", resp$gamma,", sigma=", resp$sigma, ")", sep = ""), "\n")
  }
}

