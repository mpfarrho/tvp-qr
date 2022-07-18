ffbs_R <- function(y, Z, Ht, Qtt, m, p, t, B0, V0){
  require(MASS)
  bp <- B0 #the prediction at time t=0 is the initial state
  Vp <- V0 #Same for the variance
  bt <- matrix(0,t,m) #Create vector that stores bt conditional on information up to time t
  Vt <- matrix(0,m^2,t) #Same for variances
  
  for (i in 1:t){
    R <- Ht[i,] 
    if(p==1) Qt <- Qtt[i,]*diag(p) else Qt <- diag(Qtt[i,])
    H <- Z[i,,drop = F]
    
    cfe <- y[i] - H%*%bp   # conditional forecast error
    f <- H%*%Vp%*%t(H) + R    # variance of the conditional forecast error
    inv_f <- try(t(H)%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- t(H)%*%ginv(f)
    btt <- bp + Vp%*%inv_f%*%cfe  #updated mean estimate for btt Vp * inv_F is the Kalman gain
    Vtt <- Vp - Vp%*%inv_f%*%H%*%Vp #updated variance estimate for btt
    if (i < t){
      bp <- btt
      Vp <- Vtt + Qt
    }
    bt[i,] <- t(btt)
    Vt[,i] <- matrix(Vtt,m^2,1)
  }
  
  # draw the final value of the latent states using the moments obtained from the KF filters' terminal state
  bdraw <- matrix(0,t,m)
  
  bdraw.temp <- try(btt+t(chol(Vtt))%*%rnorm(nrow(Vtt)), silent=T)
  if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, btt, Vtt+diag(1e-6,m))
  bdraw[t,] <- bdraw.temp
  
  #Now do the backward recurssions
  for (i in 1:(t-1)){
    if(p==1) Qt <- Qtt[t-1,]*diag(p) else Qt <- diag(Qtt[t-1,])
    bf <- t(bdraw[t-i+1,])
    btt <- t(bt[t-i,])
    Vtt <- matrix(Vt[,t-i,drop=FALSE],m,m)
    f <- Vtt + Qt
    
    inv_f <- try(Vtt%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- Vtt%*%ginv(f)
    
    cfe <- bf - btt
    bmean <- t(btt) + inv_f%*%t(cfe)
    bvar <- Vtt - inv_f%*%Vtt
    
    bdraw.temp <- try(bmean+t(chol(bvar))%*%rnorm(nrow(bvar)), silent=T)
    if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, bmean, bvar+diag(1e-6,m))
    
    bdraw[t-i,] <- bdraw.temp
  }
  
  return(bdraw)
}

jpr_R <- function(eps,htt,h0m,h0V,h_mu,h_rho,h_sig,theta,tau2,v,T){
  for(tt in 0:T){
    if(tt==0){
      ht_sig <- 1/(h0m/h0V + h_rho^2/h_sig)
      ht_mu <- ht_sig*(h0m/h0V + h_rho*(htt[1]-h_mu)/h_sig)
      h0d <- ht_mu + sqrt(ht_sig) * rnorm(1)
    }else{
      if(tt==1){
        ht_mu <- ((1-h_rho)*h_mu + h_rho*(h0d + htt[2]))/(1+h_rho^2)
        ht_sig <- h_sig/(1+h_rho^2)
      }else if(tt==T){
        ht_mu <- h_mu + h_rho*htt[T-1]
        ht_sig <- h_sig
      }else{
        ht_mu <- ((1-h_rho)*h_mu + h_rho*(htt[tt-1] + htt[tt+1]))/(1+h_rho^2)
        ht_sig <- h_sig/(1+h_rho^2)
      }
      ht_prop <- ht_mu + sqrt(ht_sig)*rnorm(1)
      
      liki_old <- dnorm(eps[tt],theta*z[tt]*exp(htt[tt]),exp(htt[tt])*sqrt(tau2*z[tt]),log=TRUE)
      liki_prop <- dnorm(eps[tt],theta*z[tt]*exp(ht_prop),exp(ht_prop)*sqrt(tau2*z[tt]),log=TRUE)
      
      alpha <- liki_prop-liki_old
      if(is.nan(alpha)) alpha <- -Inf
      if(log(runif(1))<alpha) htt[tt] <- ht_prop
    }
  }
  return(htt)
}

# -----------------------------------------------------------------------------------------------
# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 2 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# -----------------------------------------------------------------------------------------------
# horseshoe prior -- shrinkage
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# -----------------------------------------------------------------------------------------------
# function to lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

# -----------------------------------------------------------------------------------------------
# linear regression
get.reg <- function(y,x,sig,V0){
  k <- ncol(x)
  Vpo <- solve(crossprod(x)/sig + diag(k)*(1/V0))
  bpo <- Vpo %*% (crossprod(x,y)/sig)
  
  bdraw <- bpo + t(chol(Vpo)) %*% rnorm(k)
  return(bdraw=bdraw)
}

