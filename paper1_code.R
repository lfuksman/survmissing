E_Estimate <- function(V, delta, A, n, a = 1, cTH = 4, knots = 100,
                       cJ0 = 4, cJ1 = 0.5, cJM = 6, cB = 2) {
  ### step 1: find Fourier coefficients
  Jn <- floor(cJ0 + cJ1 * log(n))
  JFourier <- cJM * Jn
  est <- estimateCoefVar(V, delta, A, n, a, JFourier) # missing of V
  coef <- est$coef
  vjn <- est$vjn
  
  ### step 2: find optimal empirical cutoffs J
  Msummation <- cumsum(2 * vjn - (coef ^ 2))[0:(Jn + 1)] #cutoffs up to Jn
  J <- which(Msummation == min(Msummation), arr.ind = TRUE) # optimal empirical cutoff
  thetaCutoff <- coef[0:J] # theta to be used
  thetaCutoff[thetaCutoff ^ 2 <= cTH * vjn[0:J]] <- 0 # perform thresholding
  
  # calculate density over [0,1]
  z1 <- seq(from = 0, to = 1, len = knots)
  if (J == 1) {
    Basis <- matrix(1, ncol = 1, nrow = length(z1))
  } else{
    Basis <- cbind(matrix(1, ncol = 1, nrow = length(z1)),
                   (2 ^ (1 / 2)) * cos(outer(z1, pi * (1:(J - 1)))))
  }
  func <- Basis %*% thetaCutoff
  
  ### step 3: bona fide projection - calling functions from Efromovich (1999)
  denProj <- negden(func, delta = 0.01, FLAGBUMP = 1, cB = cB)
  return(denProj) # final density to return
}

tensor_j <- function(y,j){ # on [0,1]
  return(1*(j==0)+sqrt(2)*cos(pi*j*y)*1*(j>0))
}

estimateCoefVar <- function(V, delta, A, n, a = 1, J_cutoff) {
  # returns Fourier coefficients, its variance and sample-mean estimator of variance
  J <- J_cutoff - 1
  khat <- mean(A)
  
  survival <- function(V, delta, y, n) {
    sum1 <- 0
    for (l in 1:n) {
      ind1 <- ifelse((A[l] * V[l] <= y), 1, 0)
      sum2 <- sum(A * V >= A[l] * V[l])
      sum1 <- sum1 + A[l] * (1 - delta[l]) * ind1 / sum2
    }
    return(exp(-sum1))
  }
  
  Mcoef <- rep(0, (J + 1))
  Mterm <- matrix(0, ncol = (J + 1), nrow = n)
  surv <- sapply(V, survival, V = V, delta = delta, n = n)
  
  for (j in 0:J) {
    for (l in 1:n) {
      ind <- ifelse(V[l] >= 0. && V[l] <= a, 1, 0)
      tensorj <- tensor_j(V[l], j)
      Mterm[l, (j + 1)] <- (A[l] * delta[l] * ind * tensorj / (khat * surv[l]))
    }
  }
  Mcoef <- apply(Mterm, 2, mean) # Fourier coefficient
  Mvar <- apply(Mterm, 2, var) # variance of each Fourier
  vjn <- Mvar / n  # variance of each Fourier coefficient
  retlist <- list("coef" = Mcoef, "vjn" = vjn,"var" = Mvar)
  return(retlist)
}

negden <- function(f = NA, delta = 0.01, FLAGBUMP = 1, cB = 1) {
  # finds nonnegative projection
  # FLAGBUM =1 then the program removes bumps whose int f^2 dx
  # less than cof*\int (f - f.neg)^2 dx
  flag <- 0
  f1 <- f
  k <- length(f)
  AREA <- (k / (k - 1)) * mean(f) - (f[1] + f[k]) / (2 * (k - 1))
  if (all(f >= 0)) {
    flag <- 1
  }
  if (all(f <= 2 * delta) | (AREA <= 2 * delta)) {
    flag <- 2
  }
  while (flag == 0) {
    f <- f - delta
    f[f < 0] <- 0
    int <- (k / (k - 1)) * mean(f) - (f[1] + f[k]) / (2 * (k - 1))
    if (int <= AREA) {
      if (int > (10 * delta)) {
        f <- f * (AREA / int)
      }
      flag <- 1
    }
  }
  if (FLAGBUMP == 1) {
    AREASQ <- mean((f - f1) ^ 2) # area removed so far
    f <- rem.bump1(f = f, AREASQ = AREASQ, coef = cB)
  }
  if (flag == 1) {
    if (mean(f) > (10 * delta)) {
      f <- f * (AREA / mean(f))
    }
  }
  f[f < 0] <- 0
  f
}

rem.bump1 <- function(f = NA, AREASQ = NA, coef = 1) {
  ## removes any bump with L_2 area larger than coef*AREASQ
  n <- length(f)
  vec <- abvec(f) 
  if (length(vec) > 2) {
    vec <- vec[-c(1, 2)]
    k <- length(vec) / 2
    for (s in 1:k) {
      if (sum((f[vec[2 * s - 1]:vec[2 * s]]) ^ 2) / n <= coef * AREASQ) {
        f[vec[2 * s - 1]:vec[2 * s]] <- 0
      }
    }
  }
  f
}

abvec <- function(f = NA) {
  ## returns intervals where the function f is positive
  ### The first two elements of calculated vec are 00. 
  #### This indicates that all entries of f are nonpositive.
  n <- length(f) + 1
  f <- c(f, 0)
  vec <- c(0, 0)
  if (all(f > 0)) {
    vec <- c(vec, 1, n)
  }
  else {
    seq.pos <- (1:n)[f > 0]
    seq.neg <- (1:n)[f <= 0]
    seq.neg <- seq.neg[seq.neg > 1]
    a <- 1
    while (length(seq.pos) * length(seq.neg) > 0) {
      if (f[a + 1] > 0) {
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
  if (vec[length(vec)] == n) {
    vec[length(vec)] <- n - 1
  }
  vec
}