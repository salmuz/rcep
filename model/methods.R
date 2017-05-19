
########################################
# Centre et reduite les donnees 
# @param X : variables explicatives 
# @param Y : variables a expliquer
########################################
centre_reduit <- function(X, Y){
  p <- ncol(X) # explicative
  n <- nrow(X) # n sample
  Xm <- colMeans(X)
  Ym <- mean(Y)
  X <- X - rep(Xm, rep(n, p)) # centree
  Xsd <- drop(rep(1/n, n) %*% X^2) ^ 0.5
  X <- X/rep(Xsd, rep(n, p)) # reduite
  Y <- Y - Ym
  return(list(X=X, Y=as.matrix(Y), Xmean=Xm, Ymean=Ym, Xsd=Xsd, p=p, n=n));
}

########################################
# Cross Validation with a validated sample 
# @param X : variables explicatives 
# @param Y : variables a expliquer
########################################
cv.validation <- function(validation, coeff){
  interception <- rep(coeff[1], nrow(validation)) 
  prediction <- interception + drop(coeff[-1]%*%t(data.matrix(validation[,-1])))
  mean((validation[,1]-prediction)^2)
}