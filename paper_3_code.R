### applyMethod_bivariate: call for bivariate density estimation under MCAR
applyMethod_bivariate <- function(V,W,delta,gamma,A,n,a,b, cJ0=4, cJ1=0.5, cJM=2, cTH=4, cB=2, knots=100){
  ### step 1: find Fourier coef
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  VV<-V/a # rescale
  WW<-W/b # rescale
  func <- estimateCoefVar_MCAR_V(VV,WW,A,delta, gamma, n,JFourier)
  mcoef <- func$coef
  mvar <- func$var
  v_ijn <- func$v_ijn
  ### step 2: find optimal empirical cutoffs J1 and J2
  d <- mcoef[1,1]
  Msummation <- matrix(0,(Jn+1),(Jn+1))
  Mterm <- matrix(0,(Jn+1),(Jn+1))
  # find J1 and J2 
  for (i in 0:Jn){
    for (j in 0:Jn){
      Mterm[(i+1),(j+1)] <- 2*v_ijn[(i+1),(j+1)] - (mcoef[(i+1),(j+1)])^2
      Msummation[(i+1),(j+1)] <- sum(Mterm[1:(i+1),1:(j+1)])
    }
  }
  minElement <- which(Msummation == min(Msummation), arr.ind = TRUE)
  J1 <- minElement[1] 
  J2 <- minElement[2] 
  # perform thresholding 
  mcoef[mcoef^2<=cTH*v_ijn] <- 0
  z1 <- seq(from = 0, to = 1, len = knots)
  z2 <- seq(from = 0, to = 1, len = knots)
  ### step 3: smoothing weights
  Msmoothing <- matrix(0,JFourier,JFourier)
  Msmoothing <- 1 - d/(n*(mcoef)^2)
  Msmoothing[Msmoothing <0] <- 0
  Msmoothing[1,1] <- 1
  den.est <- outer(z1,z2,calculateCutoffSmoothing,Mcoef=mcoef, Msmoothing=Msmoothing, J1=J1, J2=J2)
  ### step 4: bona fide projection
  df <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  df$z <- as.vector(den.est)
  V <- getVolume(df) 
  # iteratively find delta
  dt <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  delta <- 0.05
  flag<-0
  f<-den.est
  while (flag==0){
    f <- f - delta
    f[f<0]<-0
    dt$z <- as.vector(f)
    Volume <- getVolume(dt)
    if (Volume<=V){
      f <- f * (V / Volume)
      flag<-1}
  }
  dt <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  dt$z <- as.vector((den.est-f)^2)
  Vremoved <- getVolume(dt)
  den1.est<-f
  # find bumps
  # make a matrix with 0 and 1 where density is above 0
  Mbinary <- matrix(0,length(z1),length(z2))
  Mbinary[den1.est > 0.] <- 1
  # label each component of Mbinary
  Mlabel <- cc.label(Mbinary, connect=4)
  img <- Mlabel$image
  cluster_names <- Mlabel$summary$label # find each unique clusters
  # find corresponding i and j for each cluster and area of each cluster
  row_markers <- c()
  col_markers <- c()
  z_axis <- c()
  volume <- numeric(length(cluster_names)) 
  for (k in cluster_names){
    for (i in 1:nrow(img)){
      for (j in 1:ncol(img)){
        if (img[i,j]==k){
          row_markers<- append(row_markers,z1[i])
          col_markers<-append(col_markers,z2[j])
          z_axis<-append(z_axis,den1.est[i,j])
        }
      }
    } 
    # hold i and j indices and corresponding density measurement
    df <- data.frame(row_markers,col_markers,(z_axis^2))
    names(df)[1] <- "x"
    names(df)[2] <- "y"
    names(df)[3] <- "z"
    volume[k] <- getVolume(df)
    f (volume[k] < cB*Vremoved){
      for (i in 1:nrow(img)){
        for (j in 1:ncol(img)){
          if (img[i,j]==k){
            den1.est[i,j] <-0
          }
        }
      } 
    }
    row_markers <- c()
    col_markers <- c()
    z_axis <- c()
  }
  finaldensity <- den1.est/(a*b) # rescale back; don't ensure integration to 1 because may not be the full support
  return(finaldensity)
}

library(geometry)

getVolume <- function(df) {
  #find triangular tesselation of (x,y) grid
  if (nrow(as.matrix(df))>2 & length(unique(df$x))>1) {
    res=delaunayn(as.matrix(df[,-3]),full=TRUE,options="Qz")
    # calulates sum of truncated prism volumes
    if (length(res$areas)!=0){
      sum(mapply(function(triPoints,A) A/3*sum(df[triPoints,"z"]),
                 split.data.frame(res$tri,seq_along(res$areas)),
                 res$areas))
    }
    else{
      return(0)
    }
  }
  else{
    return(0)
  }
}

readw <-function(fun=NA,z1=NA,z2=NA,dL=0,dU=1,NN=10000,FLAGD=0){
  if(FLAGD ==0){
    if(dL <0){stop("dL is negative")}
    if(dU > 1){stop("dU is larger 1")}}
  if((min(z1) < 0) |(max(z1) > 1)){
    stop("The support is beyond [0,1]")}
  if((min(z2) < 0) |(max(z2) > 1)){
    stop("The support is beyond [0,1]")}
  x <- seq(0,1,len=NN)
  y <-x; v <- x; u <- x
  eval(parse(text=paste("f <- ", fun)))
  if(is.na(f[2])){stop("the function must depend on x or y or v")}
  if(FLAGD ==0){
    f[f<dL]<- dL
    f[f>dU] <-dU 
  }
  f
}


# from wvtool library
cc.label <- function(x, connect=8, inv=FALSE, img.show=FALSE,text.size=0.3){
  if (length(x[which(x==0 | x==1)])!=length(x)){
    print("x should be a binary image.")
  } else {
    if(inv){
      x <- abs(x-1)
    }
    img <- img.lb <- matrix(0,nrow(x)+2,ncol(x)+2)
    img[2:(nrow(img)-1),2:(ncol(img)-1)] <- x
    if(connect==8){
      conn <- rbind(c(-1,-1),c(-1,0),c(-1,1),c(0,-1)) #8-connected area
    } else if(connect==4){
      conn <- rbind(c(-1,0),c(0,-1)) #4-connected area
    } else {
      print("connect should be 4 or 8")
    }	
    
    k <- 0
    for(i in 2:(nrow(img.lb)-1)){
      for(j in 2:(ncol(img.lb)-1)){
        z <- t(c(i,j)+t(conn))
        if(img[i,j]==0){
        }else if(img[i,j]==1 && sum(img[z]==0)==nrow(conn)){
          k <- k+1
          img.lb[i,j] <- k
        }else {
          zl <- z[which (img.lb[z] != 0),,drop=F]
          z.lb <- img.lb[zl]
          z.lb <- sort(z.lb[!duplicated(z.lb)])
          img.lb[i,j] <- z.lb[1]
          if(length(z.lb)>1){
            for(l in 2:length(z.lb)){
              img.lb[which(img.lb==z.lb[l])] <- z.lb[1]
            }
          }
        }
      }
    }
    lb <- sort(img.lb[!duplicated(array(img.lb))])
    lb <- lb[-1]
    for(i in lb){
      img.lb[which(img.lb==i)] <- which(lb==i)
    }
    img.lb <- img.lb[2:(nrow(img.lb)-1),2:(ncol(img.lb)-1)]
    n.lb <- length(lb)
    sum.lab <- data.frame(matrix(0,n.lb,7))
    colnames(sum.lab) <- c("label","area","aveX","aveY","dX","dY","edge")
    sum.lab[,1] <- 1:n.lb
    for (i in 1:n.lb){
      cdn <- which(img.lb==i,arr.ind=T)
      if (length(cdn[which(cdn==1 | cdn[,1]== nrow(img.lb) | cdn[,2]==ncol(img.lb))]) ==0){edg <- 0
      }else {edg <- 1}
      ranX <- range(cdn[,1])
      ranY <- range(cdn[,2])
      sum.lab[i,2:7] <- c(nrow(cdn),mean(cdn[,1]),mean(cdn[,2]),ranX[2]-ranX[1]+1,ranY[2]-ranY[1]+1,edg)
    }
    if(img.show){
      image(rot90c(img),col=(255:0)/255,ann=F,axes=F,useRaster=TRUE)
      par(new=T)
      plot(sum.lab$aveY,sum.lab$aveX,type="n",ylim=c(nrow(img),1),xlim=c(1,ncol(img)),xaxs="i",yaxs="i",ann=F,axes=F)
      text(sum.lab$aveY,sum.lab$aveX,sum.lab$label,cex=text.size,ylim=c(nrow(img),1),xlim=c(1,ncol(img)),col="red",xaxs="i",yaxs="i",ann=F)
    }
    return(list(image=img.lb,summary=sum.lab))
  }
}

# matrix rotation 90 degree clockwise
rot90c <- function(x){
  img.rot <- t(apply(x,2,rev))
  attr(img.rot, "bits.per.sample") <- attr(x,"bits.per.sample")
  attr(img.rot, "samples.per.pixel") <- attr(x, "samples.per.pixel")
  return(img.rot)
}

tensor_j <- function(y,j){ # on [0,1]
  return(1*(j==0)+sqrt(2)*cospi(j*y)*1*(j>0))
}

# calculates coefficients for MCAR
estimateCoefVar_MCAR_V <- function(V,W,A,delta, gamma, n,J_cutoff){ # pass the real number of J
  J <- J_cutoff-1
  khat <- mean(A)
  survival <- function(x,y,V,W, n,A) {
    ind <-0
    sum1 <- 0
    for (k in 1:n){ # first summation
      ind <-0
      if (A[k]*V[k] <= x){ ind <- 1 }
      count1 <- sum(A*V >= A[k]*V[k])
      sum1 <- sum1 + A[k]*(1-delta[k])*ind/count1
    }
    sum2 <-0
    for (k in 1:n){ # second summation
      ind1 <- 0
      ind2 <- 0
      if (A[k]*V[k] > x){ ind1 <- 1 }
      if (A[k]*W[k] <= y) { ind2 <- 1 }
      ind3 <- 0
      ind4 <- 0
      sum3 <- 0
      for (r in 1:n){
        ind3 <- 0
        ind4 <- 0
        if (A[r]*V[r]>=x) { ind3 <- 1 }
        if (A[r]*W[r]>=A[k]*W[k]) { ind4 <- 1 }
        sum3 <- sum3 + (ind3*ind4)
      }
      sum2 <- sum2 + ( A[k]*(1-gamma[k])*ind1*ind2 ) /(1+sum3)
    }
    return( (1/n) + exp(-sum1-sum2))
  }
  terms <-  array(dim=c(J+1,J+1,n)) 
  surv <- mapply(survival, V, W, MoreArgs = list(V=V,W=W,n=n, A=A))
  # estimate each Fourier coef 
  for (i in 0:J){
    for (j in 0:J){
      val <-0
      for (l in 1:n){
        ind <- 0
        if (V[l]<=1 & W[l]<=1){
          ind <-1
        }
        tensor_ij <- tensor_j(V[l],i)*tensor_j(W[l],j)
        terms[(i+1),(j+1),l] <-A[l]*delta[l]*gamma[l]*ind*tensor_ij/(surv[l]*khat)
      }
    }
  }
  # now we need to calculate the variance and cutoffs 
  arr <- array( unlist(terms) , c(J+1,J+1,n) )
  Mcoef <- apply( arr , 1:2 , mean )
  Mvar<- apply( arr , 1:2 , var )
  v_ijn <- Mvar/n
  retlist <- list("coef"=Mcoef,"v_ijn"=v_ijn, "var"=Mvar)
  return(retlist)
}

calculateCutoffSmoothing <- function(x,y, Mcoef, Msmoothing, J1, J2){
  value<-0
  for (i in 1:J1){
    for (j in 1:J2){
      tensor_ij <- tensor_j(x,(i-1))*tensor_j(y,(j-1))
      value <- value+ (Mcoef[i,j]*tensor_ij*Msmoothing[i,j])
    }
  }
  return(value)
} 

## applyMethod_bivariateMAR: call for bivariate density estimation under MAR missing
applyMethod_bivariateMAR<- function(V,W,delta,gamma,A,n,a,b, cJ0=4, cJ1=0.5, cJM=2, cTH=4, cB=2, knots=100){ # cTH=1
  ### step 1: find Fourier coef
  Jn <- floor(cJ0 + cJ1*log(n))
  JFourier <- cJM*Jn
  VV<-V/a # rescale
  WW<-W/b # rescale
  func <- estimateCoefVar_MAR_V(VV,WW,A,delta, gamma, n,JFourier)
  mcoef <- func$coef
  mvar <- func$var
  v_ijn <- func$v_ijn
  ### step 2: find optimal empirical cutoffs J1 and J2
  d <- mcoef[1,1]
  Msummation <- matrix(0,(Jn+1),(Jn+1))
  Mterm <- matrix(0,(Jn+1),(Jn+1))
  # find J1 and J2 
  for (i in 0:Jn){
    for (j in 0:Jn){
      Mterm[(i+1),(j+1)] <- 2*v_ijn[(i+1),(j+1)] - (mcoef[(i+1),(j+1)])^2
      Msummation[(i+1),(j+1)] <- sum(Mterm[1:(i+1),1:(j+1)])
    }
  }
  minElement <- which(Msummation == min(Msummation), arr.ind = TRUE)
  J1 <- minElement[1] 
  J2 <- minElement[2]
  # perform thresholding 
  mcoef[mcoef^2<=cTH*v_ijn] <- 0 
  z1 <- seq(from = 0, to = 1, len = knots)
  z2 <- seq(from = 0, to = 1, len = knots)
  ### step 3: smoothing weights
  Msmoothing <- matrix(0,JFourier,JFourier)
  Msmoothing <- 1 - d/(n*(mcoef)^2)
  Msmoothing[Msmoothing <0] <- 0
  Msmoothing[1,1] <- 1
  den.est <- outer(z1,z2,calculateCutoffSmoothing,Mcoef=mcoef, Msmoothing=Msmoothing, J1=J1, J2=J2)
  ### step 4: bona fide projection
  df <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  df$z <- as.vector(den.est)
  V <- getVolume(df)
  # iteratively find delta
  dt <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  delta <- 0.05
  flag<-0
  f<-den.est
  while (flag==0){
    f <- f - delta
    f[f<0]<-0
    dt$z <- as.vector(f)
    Volume <- getVolume(dt)
    if (Volume<=V){
      f <- f * (V / Volume)
      flag<-1}
  }
  dt <- expand.grid(x = as.numeric(z1), y = as.numeric(z2))
  dt$z <- as.vector((den.est-f)^2)
  Vremoved <- getVolume(dt)
  den1.est<-f
  # find bumps
  # make a matrix with 0 and 1 where density is above 0
  Mbinary <- matrix(0,length(z1),length(z2))
  Mbinary[den1.est > 0.] <- 1
  # label each component of Mbinary
  Mlabel <- cc.label(Mbinary, connect=4)
  img <- Mlabel$image
  cluster_names <- Mlabel$summary$label # find each unique clusters
  # find corresponding i and j for each cluster and area of each cluster
  row_markers <- c()
  col_markers <- c()
  z_axis <- c()
  volume <- numeric(length(cluster_names)) #7
  for (k in cluster_names){
    for (i in 1:nrow(img)){
      for (j in 1:ncol(img)){
        if (img[i,j]==k){
          row_markers<- append(row_markers,z1[i])
          col_markers<-append(col_markers,z2[j])
          z_axis<-append(z_axis,den1.est[i,j])
        }
      }
    } 
    # hold i and j indices and corresponding density measurement
    df <- data.frame(row_markers,col_markers,(z_axis^2))
    names(df)[1] <- "x"
    names(df)[2] <- "y"
    names(df)[3] <- "z"
    volume[k] <- getVolume(df)
    # if volume is smaller than cB*volume removed, then remove bump
    if (volume[k] < cB*Vremoved){
      for (i in 1:nrow(img)){
        for (j in 1:ncol(img)){
          if (img[i,j]==k){
            den1.est[i,j] <-0
          }
        }
      } 
    }
    row_markers <- c()
    col_markers <- c()
    z_axis <- c()
  }
  finaldensity <- den1.est/(a*b) # rescale back; don't ensure integration to 1 because may not be the full support
  return(finaldensity)
}

estimateCoefVar_MAR_V <- function(V,W,A,delta, gamma, n,J_cutoff){ # pass the real number of J
  J <- J_cutoff-1
  qn<-ceiling(log(n+3))
  omega_11 <- sum(A*delta*gamma)/(1+sum(delta*gamma))
  omega_11 <- pmax(1/qn,omega_11)
  omega_00 <- sum(A*(1-delta)*(1-gamma))/(1+sum((1-delta)*(1-gamma)))
  omega_00 <- pmax(1/qn,omega_00)
  P_gamma_0 <- sum(1-gamma)/n
  P_gamma_0 <- pmax(1/qn,P_gamma_0)
  P_gamma_1 <- sum(gamma)/n
  P_gamma_1 <- pmax(1/qn,P_gamma_1)
  omega_01 <- sum(A*(1-delta)*gamma)/(1+sum((1-delta)*gamma))
  omega_01 <- pmax(1/qn,omega_01)
  omega_10 <- sum(A*delta*(1-gamma))/(1+sum(delta*(1-gamma)))
  omega_10 <- pmax(1/qn,omega_10)
  
  survival_V <- function(x){ 
    num1 <- sum(A*delta*gamma*V >= x)/n
    num2 <- sum(A*(1-delta)*(1-gamma)*V >= x)/n
    return(1/n+(num1/(omega_11*P_gamma_1))+(num2/(omega_00*P_gamma_0))) # bound by 1/n from below
  }
  survival_VW <- function(x,y){ # joint survival of V and W
    term_11 <- sum((A*delta*gamma*V>=x)*(A*delta*gamma*W>=y))/(n*omega_11)
    term_01 <- sum((A*(1-delta)*gamma*V>=x)*(A*(1-delta)*gamma*W>=y))/(n*omega_01)
    term_00 <- sum((A*(1-delta)*(1-gamma)*V>=x)*(A*(1-delta)*(1-gamma)*W>=y))/(n*omega_00)
    term_10 <- sum((A*delta*(1-gamma)*V>=x)*(A*delta*(1-gamma)*W>=y))/(n*omega_10)
    return( (1/n)+term_11+term_01+term_00+term_10) # bound from below by 1
  }
  survival_CD <- function(x,y,V,W, n,A) {
    # H1(x)
    ind <-0
    sum1 <- 0
    for (k in 1:n){ # first summation
      ind <-0
      if (A[k]*V[k] <= x){ ind <- 1 }
      surv_V <- survival_V(V[k])
      sum1 <- sum1 + (A[k]*(1-delta[k])*(1-gamma[k])*ind/(omega_00*P_gamma_0*surv_V))
    } 
    sum1 <- sum1/n
    # H2(x,y)
    sum2 <-0
    ind1 <- 0
    ind2 <- 0
    val<-0
    for (k in 1:n){ # second summation
      ind1 <- 0
      ind2 <- 0
      if (A[k]*V[k] > x){ ind1 <- 1 }
      if (A[k]*W[k] <= y) { ind2 <- 1 }
      surv_VW <- survival_VW(x,W[k])
      val[k] <- (A[k]*(1-delta[k])*ind1*ind2)/(omega_00*surv_VW)
      sum2 <- sum2 + (A[k]*(1-delta[k])*(1-gamma[k])*ind1*ind2)/(omega_00*surv_VW) # issue is in (1-gamma)
    }
    sum2<- sum2/n
    return( (1/n) + exp(-sum1-sum2))
  }
  terms <-  array(dim=c(J+1,J+1,n)) 
  surv <- mapply(survival_CD, V, W, MoreArgs = list(V=V,W=W,n=n, A=A))
  # estimate each Fourier coef
  for (i in 0:J){
    for (j in 0:J){
      for (l in 1:n){
        ind <- 0
        if (V[l]<=1 & W[l]<=1){ # rescaled to [0,1]
          ind <-1
        }
        tensor_ij <- tensor_j(V[l],i)*tensor_j(W[l],j)
        terms[(i+1),(j+1),l] <-A[l]*delta[l]*gamma[l]*ind*tensor_ij/(surv[l]*omega_11)
      }
    }
  }
  # now we need to calculate the variance and cutoffs 
  arr <- array( unlist(terms) , c(J+1,J+1,n) )
  Mcoef <- apply( arr , 1:2 , mean )
  Mvar<- apply( arr , 1:2 , var )
  v_ijn <- Mvar/n
  retlist <- list("coef"=Mcoef,"v_ijn"=v_ijn, "var"=Mvar)
  return(retlist)
}
