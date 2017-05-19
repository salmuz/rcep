## Regression Ridge for data cookies ###
# @param kmax : max setting parameter
# @param sample : data cookies to analyze
# @param cvseg : crossvalidation segments 
# @param nbe : number of lambda
ridge <- function (kmax=1, sample, cvseg, nbe=1000){
  press <- rep(0, nbe)
  for(i in 1:length(cvseg)){
    valid <- sample[unlist(cvseg[i]),]
    modele <- lm.ridge(sucres~., data=sample[unlist(cvseg[-i]),], 
                       lambda=seq(0, kmax, length=nbe))
    coeff <- coef(modele)
    prediction <- matrix(coeff[,1], nrow(coeff), nrow(valid)) + 
      coeff[,-1]%*%t(data.matrix(valid[,-1]))
    press <- press + 
      rowSums((matrix(valid[,1], nrow(coeff), nrow(valid), byrow=T)-prediction)^2)
  }
  kchoice <- seq(0, kmax, length=nbe)[which.min(press)]
  return(list(k=kchoice, press=press))
}

