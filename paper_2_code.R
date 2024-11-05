#### Names of main functions to be called ######################
## calculate_univariate_density_MAR: returns overall density
## calculate_univariate_hazard_MAR: returns overall hazard rate
## estimateCoefVar_conditional_density_MAR: for conditional density
## estimateCoefVar_conditional_hazard_MAR: returns conditional hazard rate
############## CODE FOR DENSITY AND HAZARD ############
### helper function- find cosine basis on [0,1]
tensor_j <- function(y,j){ # on [0,1]
  return(1*(j==0)+sqrt(2)*cos(pi*j*y)*1*(j>0))
}

### estimate univariate density
calculate_univariate_density_MAR <- function(V,delta,A, n,a,cTH=4, knots=100, 
                                             cJ0=4, cJ1 =0.5, cJM = 6,cB = 2){
  ### step 1: find Fourier coef
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  V <- V/a # scale on [0,1]
  est <- estimateCoefVar_univariate_MAR_var(V,delta,A, n,JFourier,cTH) 
  coef <- est$coef
  var <- est$var
  vjn <- est$vjn
  
  ### step 2: find optimal empirical cutoffs J
  d <- coef[1]
  Msummation <- rep(0,Jn)
  Msummation <- cumsum(2*var/n - (coef^2))[0:(Jn+1)] #cutoffs up to Jn
  J <- which(Msummation == min(Msummation), arr.ind = TRUE) 
  thetaCutoff <- coef[0:J] # theta to be used
  thetaCutoff[thetaCutoff^2<=cTH*vjn[0:J]] <-0 # perform thresholding 
  
  # find density over [0,1]
  z1 <- seq(from = 0, to = 1, len = knots)
  if(J==1){
    Basis <- matrix(1,ncol=1,nrow=length(z1))
  } else{ 
    Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                  (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))} 
  func <- Basis%*%thetaCutoff
  
  ### step 3: scale and remove negative
  func[func < 0] <-0 
  denProj<- func/a
  return(denProj)
}

### helper function- calculate Fourier coefficients for density
# pass scaled V 
estimateCoefVar_univariate_MAR_var <-  function(V,delta, A, n,
                                                J_cutoff,cTH, cc=1){
  J <- J_cutoff-1
  CTH <- cTH
  
  # calculate omega for each delta
  omega_1 <- sum(A*delta)/sum(delta)  
  omega_0 <- sum(A*(1-delta))/sum(1-delta) 
  omega_1 <- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  
  # calculate survival of V
  survival_V <- function(V,delta,y,n){
    return(sum(A*V*delta >= y)/(n*omega_1) 
           + sum(A*V*(1-delta) >= y)/(n*omega_0)) 
  }
  surv_V <- sapply(V,survival_V,V=V,delta=delta,n=n)
  
  # calculate survival of C 
  survival_C <- function(V,delta,y,n){ 
    ind1 <-0     
    sum1 <- 0
    for (l in 1:n){
      ind1 <- 0 
      if (A[l]*V[l]<=y) {ind1 <-1 }
      sum1 <- sum1+ (A[l]*(1-delta[l])*ind1)/(omega_0*surv_V[l]) # omega[l]
    }
    return(exp(-sum1/n))
  }
  
  # calculate coefficients
  Mcoef <- rep(0,(J+1))
  Mterm <- matrix(0,ncol=(J+1), nrow=n)
  surv_C <- sapply(V,survival_C,V=V,delta=delta,n=n)
  
  for (j in 0:J){ 
    for (l in 1:n){
      ind <- 0
      if (V[l]>=0. && V[l]<=1 ){ ind <-1 }
      tensorj <- tensor_j(V[l],j)
      Mterm[l,(j+1)] <- (A[l]*delta[l]*ind*tensorj/(omega_1*surv_C[l])) 
    }
  }
  Mcoef <- apply(Mterm, 2, mean)
  Mvar <- apply(Mterm, 2, var)
  vjn <- Mvar/n  # variance estimator
  retlist <- list("coef"=Mcoef,"vjn"=vjn, "var"=Mvar)
  return(retlist)
}

### estimate univariate hazard
calculate_univariate_hazard_MAR <- function(V,delta,A, n,a=0,b=0.75,cTH=4,  
                                            knots=100, cJ0=4, cJ1 =0.5, cJM = 6,cB = 2){
  ### step 1: find Fourier coef
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  
  est <- estimateCoef_hazard_MAR(V,delta,A, n,a,b,JFourier) # missing of V
  coef <- est$coef
  var <- est$var
  vjn <- est$vjn
  
  ### step 2: find optimal empirical cutoffs J
  Msummation <- rep(0,Jn)
  Msummation <- cumsum(2*var/n - (coef^2))[0:(Jn+1)] #cutoffs up to Jn
  J <- which(Msummation == min(Msummation), arr.ind = TRUE) 
  thetaCutoff <- coef[0:J] # theta to be used
  thetaCutoff[thetaCutoff^2<=cTH*vjn[0:J]] <-0 # perform thresholding 
  
  # over [0,1]
  z1 <- seq(from = 0, to = 1, len = knots)
  if(J==1){
    Basis <- matrix(1,ncol=1,nrow=length(z1))
  } else{ 
    Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                  (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))} 
  func <- Basis%*%thetaCutoff
  
  ### step 3: scale and remove negative
  finalHazard <- func/b
  finalHazard[finalHazard < 0] <-0 
  return(finalHazard)
}

### helper function- calculate Fourier coefficients for hazard rate
# pass scaled V, estimate on [a,a+b]
estimateCoef_hazard_MAR <-  function(V,delta, A, n,a=0,b=0.75,J_cutoff){ 
  J <- J_cutoff-1
  
  # calculate omega 
  omega_1 <- sum(A*delta)/sum(delta)   
  omega_0 <- sum(A*(1-delta))/sum(1-delta) 
  # to avoid divison by zero
  omega_1 <- ifelse(omega_1>1/log(n+3), omega_1, 1/log(n+3))
  omega_0 <- ifelse(omega_0>1/log(n+3), omega_0, 1/log(n+3)) 
  
  YSC <- (V-a)/b # scale
  # calculate survival of V
  survival_V <- function(V,delta,y,n){
    return(sum(A*V*delta >= y)/(n*omega_1) 
           + sum(A*V*(1-delta) >= y)/(n*omega_0)) 
  }
  surv_V <- sapply(YSC,survival_V,V=YSC,delta=delta,n=n) 
  
  # calculate coefficients
  Mcoef <- rep(0,(J+1))
  Mterm <- matrix(0,ncol=(J+1), nrow=n)
  indVector <- numeric(n)
  indVector[(YSC>=0 & YSC<=1)]<-1
  for (j in 0:J){ 
    for (l in 1:n){
      tensorj <- tensor_j(YSC[l],j)
      Mterm[l,(j+1)]<-(A[l]*delta[l]*indVector[l]*tensorj/(omega_1*surv_V[l]))
    }
  }
  Mcoef <- apply(Mterm, 2, mean)
  Mvar <- apply(Mterm, 2, var)
  vjn <- Mvar/n  #  variance estimator
  retlist <- list("coef"=Mcoef,"vjn"=vjn, "var"=Mvar)
  return(retlist)
}


############## CODE FOR CONDITIONAL DENSITY AND HAZARD RATE ############
### helper function- find density of X
### pass scaled on [0,1] and returns scaled on [0,1]
density_X <- function(X, n, cTH=5, cJ0=4, cJ1 =0.5, cJM = 6,
                      cB = 2, knots=100){ 
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  J <- JFourier-1
  # calculate coefficients
  qn<-ceiling(log(n+3))
  Mterm <- matrix(0,ncol=(J+1), nrow=n)
  ind <- ifelse( (X>=0 & X<=1), 1, 0)
  tensor_matrix <-outer(X, c(0:J), tensor_j)
  for (j in 0:J){ 
    for (l in 1:n){
      Mterm[l,(j+1)] <- ind[l]*tensor_matrix[l,(j+1)]
    }
  }
  coef <- apply(Mterm, 2, mean)
  var <- apply(Mterm, 2, var) # variance of each Fourier
  vjn <- var/n  # sample variance estimator
  
  # now find optimal cutoff J
  Msummation <- cumsum(2*vjn - (coef^2))[0:(Jn+1)] #cutoffs up to Jn
  J <- which(Msummation == min(Msummation), arr.ind = TRUE) 
  thetaCutoff <- coef[0:J] # theta to be used
  thetaCutoff[thetaCutoff^2<=cTH*vjn[0:J]] <-0 # perform thresholding 
  
  # find density over [0,1]
  z1<-X #[X>=0 & X<=1]
  if(J==1){
    Basis <- matrix(1,ncol=1,nrow=length(z1))
  } else{ 
    Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                  (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))} 
  func <- Basis%*%thetaCutoff
  func[z1<0 | z1>1] <-0  # outside of [0,1], density is 0 for all X's
  denProj <- negden(func, delta = 0.01, FLAGBUMP = 1, cB = cB)  
  denProj <- pmax(1/qn,denProj) # final density
  funcfinal <- cbind(z1,denProj) 
  return(funcfinal) 
}

### helper function- pass scaled on [0,1] and returns scaled on [0,1]
p_nuissance_t <- function(V,delta, A, X, n,cTH=4, cc=1, 
                          cJ0=4, cJ1 =0.5, cJM = 6,cB = 2, knots=100){
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  Jcutoff <- JFourier-1
  # calculate coefficients
  Mterm <- vjn <- matrix(0,ncol=(Jcutoff+1), nrow=n)
  omega_1 <- sum(A*delta)/sum(delta) 
  omega_0 <- sum(A*(1-delta))/sum(1-delta)  
  omega_hat_1<- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_hat_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  if(is.na(omega_hat_0)){omega_hat_0<-cc/log(n)}
  if(is.na(omega_hat_1)){omega_hat_1<-cc/log(n)}
  theta<-numeric(n)
  qn<-ceiling(log(n+3))
  
  tensor_matrix <-  outer(X, c(0:Jcutoff), tensor_j)  
  ind1 <- ifelse( (A*V>=0 & A*V<=1), 1, 0) 
  ind2 <- ifelse( (X>=0 & X<=1), 1, 0) 
  omegahat <- ifelse( (delta==0), omega_hat_0, omega_hat_1)
  ind_fun <- Vectorize(function(x, y){
    ifelse( (x>=y), 1, 0)
  }, vectorize.args=c("x", "y"))
  ind <- outer(X=A*V, Y=A*V, ind_fun) 
  for (k in 1:n){
    for (j in 0:Jcutoff){ 
      value <- numeric(n)
      for (l in 1:n){
        value[l] <- (A[l]*ind[l,k]*ind1[l]*ind2[l]*
                       tensor_matrix[l,(j+1)])/omegahat[l]
      }
      Mterm[k,(j+1)] <- mean(value) 
      vjn[k,(j+1)] <- var(value)/n 
    }
    # for each row: find optimal cutoff J
    Msummation <- cumsum(2*vjn[k,] - (Mterm[k,])^2)[0:(Jn+1)] 
    J <- which(Msummation == min(Msummation), arr.ind = TRUE)[1]
    thetaCutoff <- Mterm[k,0:J] # theta to be used
    thetaCutoff[thetaCutoff^2<=cTH*vjn[k,0:J]] <-0 # perform thresholding 
    
    # find density at each x=X[k]
    z1 <- X[k]
    if(J==1){
      Basis <- matrix(1,ncol=1,nrow=length(z1))
    } else{ 
      Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                    (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))}
    func <- Basis%*%thetaCutoff
    func[func < 0] <-0 # remove negatives
    func <- pmax(1/(2*qn),func) # bound from below
    theta[k] <- func
  }
  func <- cbind(A*V, X, theta) # values for survival function at each point
  return(func)
}

### helper function - estimate S^{V|X}
### pass scaled on [0,1] and return scaled on [0,1]
survival_V_conditional_on_X <- function(V,delta, A, X, n,cTH=10, cTHd=5,  
                                        cc=1,cJ0=4, cJ1 =0.5, cJM = 6,cB = 2, knots=100){ 
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  Jcutoff <- JFourier-1
  Mterm <- matrix(0,ncol=(Jcutoff+1), nrow=n)
  omega_1 <- sum(A*delta)/sum(delta)  
  omega_0 <- sum(A*(1-delta))/sum(1-delta)  
  omega_hat_1<- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_hat_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  if(is.na(omega_hat_0)){omega_hat_0<-cc/log(n)}
  if(is.na(omega_hat_1)){omega_hat_1<-cc/log(n)}
  densityX <- density_X(X, n, cTH=cTHd)
  qn<-ceiling(log(n+3))
  
  theta<-numeric(n)
  # S^{V|X} is evaluated on all V
  ind_fun <- Vectorize(function(x, y){
    ifelse( (x>=y), 1, 0)
  }, vectorize.args=c("x", "y"))
  ind <- outer(X=A*V, Y=A*V, ind_fun) 
  omega_hat <- ifelse( (delta==0), omega_hat_0, omega_hat_1)
  tensor_matrix <-  outer(X, c(0:Jcutoff), tensor_j)
  # calculate coefficients
  for (k in 1:n){
    for (j in 0:Jcutoff){ 
      for (l in 1:n){
        Mterm[l,(j+1)] <- (A[l]*ind[l,k]*
                             tensor_matrix[l,(j+1)])/(omega_hat[l]*densityX[l,2])
      }
    }
    coef <- apply(Mterm, 2, mean)
    var <- apply(Mterm, 2, var) 
    vjn <- var/n  # sample variance estimator
    
    # now find optimal cutoff J
    Msummation <- cumsum(2*vjn - (coef^2))[0:(Jn+1)] #cutoffs up to Jn
    J <- which(Msummation == min(Msummation), arr.ind = TRUE) 
    thetaCutoff <- coef[0:J] # theta to be used
    thetaCutoff[thetaCutoff^2<=cTH*vjn[0:J]] <-0 # perform thresholding 
    
    # multiply by tensor at x=X[k]
    z1 <- X[k]
    if(J==1){
      Basis <- matrix(1,ncol=1,nrow=length(z1))
    } else{ 
      Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                    (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))}
    func <- Basis%*%thetaCutoff
    func <- pmax(0.02,func)
    theta[k] <- func
  }
  func <- cbind(A*V, X, theta) # values for survival function at each point
  return(func)
}

### helper function- estimate S^{C|X}
### pass scaled on [0,1] and return scaled on [0,1]
survival_C_conditional_on_X <- function(V,delta, A, X, n, cTHd=5, cTHv=9, 
                                        cTH=15, cc=1, cJ0=4, cJ1 =0.5, cJM = 6,cB = 2, knots=100){
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  Jcutoff <- JFourier-1
  omega_1 <- sum(A*delta)/sum(delta) 
  omega_0 <- sum(A*(1-delta))/sum(1-delta) 
  omega_hat_1<- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_hat_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  if (is.na(omega_hat_0)){omega_hat_0 <- cc/log(n)}
  densityX <- density_X(X, n, cTH=cTHd)
  survV_X<-survival_V_conditional_on_X(V,delta, A, X, n, cTH=cTHv, cTHd=cTHd) 
  theta<-numeric(n)
  Mterm <- vjn <- matrix(0,ncol=(Jcutoff+1), nrow=n)
  ind_fun <- Vectorize(function(x, y){
    ifelse( (x<=y), 1, 0)
  }, vectorize.args=c("x", "y"))
  ind <- outer(X=A*V, Y=A*V, ind_fun)
  tensor_matrix <-  outer(X, c(0:Jcutoff), tensor_j)
  # calculate coefficients
  for (k in 1:n){
    for (j in 0:Jcutoff){ 
      value<- numeric(n)
      for (l in 1:n){
        denX <- as.numeric(densityX[l,2])
        survVX <- as.numeric(survV_X[l,3])
        value[l]<-(A[l]*(1-delta[l])*ind[l,k]*
                     tensor_matrix[l,(j+1)])/(omega_hat_0*denX*survVX)
      }
      Mterm[k,(j+1)] <- mean(value) # coef for each k and j
      vjn[k,(j+1)] <- var(value)/n 
    }
    
    # now find optimal cutoff J
    Msummation <- cumsum(2*vjn[k,] - (Mterm[k,]^2))[0:(Jn+1)] 
    J <- which(Msummation == min(Msummation), arr.ind = TRUE)[1]
    thetaCutoff <- Mterm[k,0:J] # theta to be used
    thetaCutoff[thetaCutoff^2<=cTH*vjn[k,0:J]] <-0 # perform thresholding 
    
    # multiply by tensor at x=X[k]
    z1 <- X[k]
    if(J==1){
      Basis <- matrix(1,ncol=1,nrow=length(z1))
    } else{ 
      Basis <-cbind(matrix(1,ncol=1,nrow=length(z1)),
                    (2^(1/2)) * cos(outer(z1, pi * (1:(J-1)))))}
    func <- Basis%*%thetaCutoff
    theta[k] <- as.numeric(exp(-func))
  }
  func <- cbind(A*V, X, theta) # values for survival function at each point
  return(func)
}

### helper function- evaluate function using Fourier coefficients
calculateValue <- function(x,y, Mcoef, J1, J2){ # over [0,1]
  value<-0
  for (j in 1:J2){
    for (i in 1:J1){
      tensor_ji <- tensor_j(x,(j-1))*tensor_j(y,(i-1))
      value <- value+ (Mcoef[j,i]*tensor_ji)
    }
  }
  return(value)
}

### conditional hazard function
estimateCoefVar_conditional_hazard_MAR <-  function(V, delta, A, X, n, a, b, 
                                                    cTH=10, cTHp=25, cc=1, cJ0=4, cJ1 =0.5, cJM = 6, cB = 2, knots=100){ 
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  J <- JFourier-1
  
  # calculate omega for each delta
  omega_1 <- sum(A*delta)/sum(delta) 
  omega_0 <- sum(A*(1-delta))/sum(1-delta) 
  omega_1 <- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  
  if (is.na(omega_0)){omega_0 <- cc/log(n)}
  if (is.na(omega_1)){omega_1 <- cc/log(n)}
  
  VSC <- V/b # scale V
  XX <- X/a # scale X
  p_nuissance <- p_nuissance_t(VSC,delta,A,XX,n,cTH=cTHp) 
  terms <-  array(dim=c(J+1,J+1,n))
  ind <- ifelse( (A*VSC>=0 & A*VSC<=1), 1, 0) 
  for (j in 0:J){
    for (i in 0:J){
      for (l in 1:n){
        tensor_ji <- tensor_j(A[l]*VSC[l],j)*tensor_j(XX[l],i)
        pp<-p_nuissance[l,3] 
        terms[(j+1),(i+1),l] <- (A[l]*delta[l]*ind[l]*tensor_ji)/(omega_1*pp)
      }
    }
  }
  
  # calculate coefficients
  arr <- array( unlist(terms) , c(J+1,J+1,n) )
  Mcoef <- apply( arr , 1:2 , mean ) 
  Mvar<- apply( arr , 1:2 , var )
  v_ijn <- Mvar/n 
  
  Msummation <- matrix(0,(Jn+1),(Jn+1))
  Mterm <- matrix(0,(Jn+1),(Jn+1))
  
  # find optimal cutoffs
  for (j in 0:Jn){
    for (i in 0:Jn){
      Mterm[(j+1),(i+1)] <- 2*v_ijn[(j+1),(i+1)] - (Mcoef[(j+1),(i+1)])^2
      Msummation[(j+1),(i+1)] <- sum(Mterm[1:(j+1),1:(i+1)])
    }
  }
  minElement <- which(Msummation == min(Msummation), arr.ind = TRUE)
  J2 <- minElement[1] #j row
  J1 <- minElement[2] #i column
  Mcoef[Mcoef^2<=cTH*v_ijn] <- 0 # 15
  coef <- Mcoef[1:J2,1:J1, drop=FALSE]
  
  # find density over [0,1] 
  z_T <- seq(from = 0, to = 1, len = knots) # TT
  z_X <- seq(from = 0, to = 1, len = knots) # X
  den.est <- outer(z_T,z_X,calculateValue,Mcoef=coef, J1=J1, J2=J2)
  finalHazard<-den.est/(a*b) # to scale back to [0,a]*[0,b]
  finalHazard[finalHazard<0]<-0 # remove negatives
  return(finalHazard)
}

### conditional density function
estimateCoefVar_conditional_density_MAR <-  function(V, delta, A, X, n, a, b, 
                                                     cTH=10, cTHd=5, cTHs=5, cTHv=9, cc=1, cJ0=4, cJ1 =0.5, 
                                                     cJM = 6, cB = 2, knots=100){
  Jn <- floor(cJ0 + cJ1*log(n)) 
  JFourier <- cJM*Jn
  J <- JFourier-1
  
  # calculate omega for each delta
  omega_1 <- sum(A*delta)/sum(delta)
  omega_0 <- sum(A*(1-delta))/sum(1-delta) 
  omega_1 <- ifelse(omega_1>cc/log(n), omega_1, cc/log(n))
  omega_0 <- ifelse(omega_0>cc/log(n), omega_0, cc/log(n)) 
  
  if(is.na(omega_0)){omega_0 <- 1/log(n+3)}
  if(is.na(omega_1)){omega_1 <- 1/log(n+3)}
  
  VSC <- V/b # scale V
  XX <- X/a # scale X
  
  densityX <- density_X(XX, n, cTH=cTHd) 
  survC_X <- survival_C_conditional_on_X(VSC, delta, A, XX, n, 
                                         cTHd =cTHd, cTHv =cTHv, cTH=cTHs) 
  # find coefficients
  terms <-  array(dim=c(J+1,J+1,n))
  ind <- ifelse( (A*VSC>=0 & A*VSC<=1), 1, 0) 
  for (j in 0:J){
    for (i in 0:J){
      for (l in 1:n){
        tensor_ji <- tensor_j(A[l]*VSC[l],j)*tensor_j(XX[l],i)
        denX <- as.numeric(densityX[l,2])
        survCX <- as.numeric(survC_X[l,3])
        terms[(j+1),(i+1),l] <- (A[l]*delta[l]*ind[l]*
                                   tensor_ji)/(omega_1*denX*survCX)
      }
    }
  }
  arr <- array( unlist(terms) , c(J+1,J+1,n) )
  Mcoef <- apply( arr , 1:2 , mean )
  Mvar<- apply( arr , 1:2 , var )
  v_ijn <- Mvar/n
  
  Msummation <- matrix(0,(Jn+1),(Jn+1))
  Mterm <- matrix(0,(Jn+1),(Jn+1))
  
  # find optimal cutoffs
  for (j in 0:Jn){
    for (i in 0:Jn){
      Mterm[(j+1),(i+1)] <- 2*v_ijn[(j+1),(i+1)] - (Mcoef[(j+1),(i+1)])^2
      Msummation[(j+1),(i+1)] <- sum(Mterm[1:(j+1),1:(i+1)])
    }
  }
  minElement <- which(Msummation == min(Msummation), arr.ind = TRUE)
  J2 <- minElement[1] #j 
  J1 <- minElement[2] #i
  Mcoef[Mcoef^2<=cTH*v_ijn] <- 0
  coef <- Mcoef[1:J2,1:J1, drop=FALSE] 
  
  # find density over [0,1] 
  z_T <- seq(from = 0, to = 1, len = knots) # TT
  z_X <- seq(from = 0, to = 1, len = knots) # X
  den.est <- outer(z_T,z_X,calculateValue,Mcoef=coef, J1=J1, J2=J2)
  finalDensity<-den.est/(a*b) # to scale back to [0,a]*[0,b]
  finalDensity[finalDensity<0]<-0 # remove negatives
  return(finalDensity)
}

# helper functions from Efromovich(2018) book
negden<-function(f = NA, delta = 0.01, FLAGBUMP = 1, cB = 1) {
  #this is negden in book1fig that finds nonnegative projection
  #FLAGBUM =1 then the program removes bumps whose int f^2 dx
  # less than cof*\int (f - f.neg)^2 dx
  flag <- 0
  f1 <- f
  k <- length(f)
  AREA <- (k/(k - 1)) * mean(f) - (f[1] + f[k])/(2 * (k - 1))
  if(all(f >= 0)) {
    flag <- 1
  }
  if(all(f <= 2 * delta) | (AREA <= 2 * delta)) {
    flag <- 2
  }
  while(flag == 0) {
    f <- f - delta
    f[f < 0] <- 0
    int <- (k/(k - 1)) * mean(f) - (f[1] + f[k])/(2 * (k - 1))
    if(int <= AREA) {
      if(int > (10 * delta)) {
        f <- f * (AREA/int)
      }
      flag <- 1
    }
  }
  if(FLAGBUMP == 1) {
    AREASQ <- mean((f - f1)^2) # area removed so far
    f <- rem.bump1(f = f, AREASQ = AREASQ, coef = cB)
  }
  if(flag == 1) {
    if(mean(f) > (10 * delta)) {
      f <- f * (AREA/mean(f))
    }
  }
  f[f < 0] <- 0
  f
}

rem.bump1<-function(f = NA, AREASQ = NA, coef = 1) {
  ##this is rem.bump1 in book1fig that removes any bump with L_2 area
  ## larger than coef*AREASQ, 
  n <- length(f)
  vec <- abvec(f) # same for f and den.est
  if(length(vec) > 2) {
    vec <- vec[ - c(1, 2)]
    k <- length(vec)/2
    for(s in 1:k) {
      if(sum((f[vec[2 * s - 1]:vec[2 * s]])^2)/n <= coef * AREASQ) {
        f[vec[2 * s - 1]:vec[2 * s]] <- 0
      }
    }
  }
  f
}

abvec<-function(f = NA) {
  ##ab.vec (abvec)
  ##in book1fPC that returns intervals where the function f is positive
  ###!!!!!!! the first two elements of calculated vec are 00. This
  ####!!!!!!!  indicates that all entries of f are nonpositive.
  n <- length(f) + 1
  f <- c(f, 0)
  vec <- c(0, 0)
  if(all(f > 0)) {
    vec <- c(vec, 1, n)
  }
  else {
    seq.pos <- (1:n)[f > 0]
    seq.neg <- (1:n)[f <= 0]
    seq.neg <- seq.neg[seq.neg > 1]
    a <- 1
    while(length(seq.pos) * length(seq.neg) > 0) {
      if(f[a + 1] > 0) {
        b <- min(seq.neg) - 1
      }
      else {
        a <- min(seq.pos)
        seq.neg <- seq.neg[seq.neg > a]
        b <- min(seq.neg) - 1
      }
      vec <- c(vec, a, b)
      seq.pos <- seq.pos[seq.pos > b]
      seq.neg <- seq.neg[seq.neg > b]
      a <- b
    }
  }
  if(vec[length(vec)] == n) {
    vec[length(vec)] <- n - 1
  }
  vec
}