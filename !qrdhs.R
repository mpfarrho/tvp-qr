# other required packages
require(Rcpp)
require(Matrix)
require(MASS)
require(spam)
require(pgdraw)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # mac

# C++ functions
sourceCpp("ffbs.cpp")   # C++ implementation of FFBS
sourceCpp("jpr_qr.cpp") # C++ implementation of the JPR algorithm for QRs
source("aux.R")         # some auxiliary functions
source("dsp_aux.R")     # functions for dynamic horseshoe: https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12325

# -----------------------------------------------------------------------------------------
# quantile specific regression model
tvpqr <- function(Y,X,Xout=NULL,bt_init=NULL,mode="qr",p=0.5,tvp="dhs",sv=FALSE,cons.cons=FALSE,fhorz=0,
                  nburn=1000,nsave=1000,thinfac=1,quiet=FALSE){
  message("\nStarting with quantile p=",p,".")
  
  # sampler specification
  ntot <- nburn + nsave*thinfac
  thin.set <- floor(seq(nburn+1,ntot,length.out=nsave))
  in.thin <- 0
  pred.dsp <- FALSE
  pred.sv <- FALSE
  
  # some checks
  if(!is.null(Xout)){
    if(nrow(Xout)!=fhorz) stop("Check size of 'Xout' and desired forecast horizon 'fhorz'.")
  }
  if(!(tvp %in% c("dhs","shs","iG"))) stop("Setting 'tvp=",tvp,"' not available. Choose from 'iG', 'shs' or 'dhs'.")
  if(mode=="reg" & p!=0.5 & !sv) stop("Gaussian error term cannot be used to estimate quantile specific models, and SV is required.")
  
  # data and dimensions
  K <- NCOL(X)
  M <- NCOL(Y)
  T <- NROW(Y)
  
  if(cons.cons){
    cons <- TRUE
    alpha <- 0
    iota <- matrix(1,T,1)
  }else{
    cons <- FALSE
    alpha <- 0
    iota <- matrix(1,T,1)
  }

  # horseshoe prior setup
  D <- 1
  XtX <- build_XtX(X)
  sigma_e <- sd(Y, na.rm=TRUE); sigma_et = rep(sigma_e,T)
  chol0 <- initCholReg.spam(obs_sigma_t2 = abs(rnorm(T)),
                            evol_sigma_t2 = matrix(abs(rnorm(T*K)), nrow = T),
                            XtX = XtX, D = D)
  draw <- sampleBTF_reg(Y, X, obs_sigma_t2 = sigma_et^2, 
                        evol_sigma_t2 = matrix(0.01*sigma_et^2, nr = T, nc = K), 
                        XtX = XtX, D = D, chol0 = chol0)
  omega <- diff(draw, differences = D)
  beta0 <- matrix(draw[1:D,], nr = D)
  evolParams <- initEvolParams(omega, evol_error = "DHS")
  evolParams0 <- initEvol0(beta0, commonSD = FALSE)
  tau0 <- rep(1,K)
  dhs_gl <- sqrt(T*K)
  steq_var <- matrix(1,nrow=T,ncol=K)
  
  # in case of static SHS
  lambda.shs <- nu.shs <- matrix(1,T,K)
  tau.shs <- zeta.shs <- rep(1,K)
  
  # prior for scale parameter
  mm_n0 <- 1 # mean for the iG prior
  vv_s0 <- 1 # variance for the iG prior
  n0 <- 2*(2+((mm_n0^2)/vv_s0))
  s0 <- 2*(mm_n0 + ((mm_n0^3)/vv_s0))
  
  if(sv){
    h0m <- 0 # prior mean of initial state
    h0V <- 1 # prior variance of initial state
    
    # state equation parameters
    h_mu <- 0
    h_rho <- 1
    h_sig <- 0.01
  }
  
  # scale parameters for AL_p
  theta <- (1-2*p)/(p*(1-p))
  tau2 <- 2/(p*(1-p))
  
  if(mode=="reg"){
    theta <- 0
    tau2 <- 1
    
    sv_priors <- specify_priors(
      mu = sv_normal(mean = 0, sd = 100),
      phi = sv_constant(1-1e-12),
      sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
      nu = sv_infinity(),
      rho = sv_constant(0)
    )
    
    sv_draw <- list(mu = 0, phi = 0.99999, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    htt <- sv_latent <- matrix(0,T,1)
  }
  
  # MCMC objects
  b <- bv <- b_pr <- matrix(0,K)
  
  if(is.null(bt_init)){
    bt <- array(0,dim=c(T,K))
  }else{
    bt <- bt_init
  }
  
  omega <- omegat <- array(10,dim=c(K))
  dhs_se_para <- matrix(0,nrow=K,ncol=2)
  sig2 <- matrix(1,T)
  v <- matrix(1,T)
  tau2sig2v <- tau2*sig2*v
  
  # storage
  bt_store <- array(NA,dim=c(nsave,T,K))
  shrink_store <- array(NA,dim=c(nsave,T,K))
  sig2_store <- array(NA,dim=c(nsave,T))
  fcst_store <- fcstsig2_store <- array(NA,dim=c(nsave,fhorz))
  
  # main sampling loop
  if(!quiet) pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
  for(irep in 1:ntot){
    # Step 1: Sample quantile-specific coefficients
    if(cons){
      norm <- as.numeric(1/sqrt(tau2sig2v))
      X_new <- iota*norm
      y_new <- (Y-(theta*v)-apply(Xt*bt,1,sum))*norm
      
      alph_V <- solve(crossprod(X_new) + 1/10)
      alph_a <- alph_V %*% crossprod(X_new,y_new)
      alpha <- as.numeric(alph_a + t(chol(alph_V))%*%rnorm(1))
    
      X_new <- X*norm
      y_new <- (Y-theta*v-iota*alpha)*norm
    }else{
      norm <- as.numeric(1/sqrt(tau2sig2v))
      X_new <- X*norm
      y_new <- (Y-theta*v)*norm
    }
    
    if(tvp=="dhs"){
      steq_var <- rbind(matrix(apply(evolParams$sigma_wt^2,2,median),ncol=K),
                        evolParams$sigma_wt^2)
      # steq_var[steq_var>1] <- 1
      # steq_var[steq_var<1e-12] <- 1e-12
      tau0 <- as.numeric(evolParams0$sigma_w0^2)
      tau0[tau0<1] <- 1
      
      draw <- try(t(ffbs(t(as.matrix(y_new)),X_new,matrix(1,T,1),steq_var,K,1,T,matrix(0,K,1),diag(K)*tau0)),silent=TRUE)
      if(is(draw,"try-error")){
        draw <- ffbs_R(as.matrix(y_new),X_new,matrix(1,T,1),steq_var,K,1,tt,matrix(0,K,1),diag(K)*tau0)
      }
      omega <- diff(draw, differences = D)
      beta0 <- matrix(draw[1:D,], nr = D)
      evolParams0 <- sampleEvol0(beta0, evolParams0, A = 1, commonSD = FALSE)
      evolParams <- sampleEvolParams(omega, evolParams, 1/dhs_gl, "DHS")
      
      bt <- draw
      bv <- (evolParams$sigma_wt[T-1,])
      dhs_se_para <- cbind(evolParams$dhs_mean,evolParams$dhs_phi)
    }else if(tvp=="shs"){
      draw <- try(t(ffbs(t(as.matrix(y_new)),X_new,matrix(1,T,1),steq_var,K,1,T,matrix(0,K,1),10*diag(K))),silent=TRUE)
      if(is(draw,"try-error")){
        draw <- ffbs_R(as.matrix(y_new),X_new,matrix(1,T,1),steq_var,K,1,tt,matrix(0,K,1),10*diag(K))
      }
      omega <- rbind(apply(draw,2,function(x){mean(diff(x))}),diff(draw, differences = D))
      for(k in 1:K){
        shs.draw <- get.hs(omega[,k],lambda.hs=lambda.shs[,k],nu.hs=nu.shs[,k],tau.hs=tau.shs[k],zeta.hs=zeta.shs[k])
        lambda.shs[,k] <- shs.draw$lambda
        nu.shs[,k] <- shs.draw$nu
        tau.shs[k] <- shs.draw$tau
        zeta.shs[k] <- shs.draw$zeta
        steq_var[,k] <- shs.draw$psi
        # steq_var[steq_var>1] <- 1
        # steq_var[steq_var<1e-12] <- 1e-12
      }
      bt <- draw
    }else if(tvp=="iG"){
      draw <- try(t(ffbs(t(as.matrix(y_new)),X_new,matrix(1,T,1),steq_var,K,1,T,matrix(0,K,1),10*diag(K))),silent=TRUE)
      if(is(draw,"try-error")){
        draw <- ffbs_R(as.matrix(y_new),X_new,matrix(1,T,1),steq_var,K,1,tt,matrix(0,K,1),10*diag(K))
      }
      omega <- diff(draw, differences = D)
      for(k in 1:K){
        steq_var[,k] <- 1/rgamma(1,3 + (T-1)/2,0.3 + sum(omega[,k]^2)/2)
      }
      # steq_var[steq_var>1] <- 1
      # steq_var[steq_var<1e-8] <- 1e-8
      bt <- draw
    }
    
    # Step 2a: Sample auxiliary variable v
    if(mode=="qr"){
      for(tt in 1:T){
        delta2 <- as.numeric((Y[tt,] - alpha - X[tt,] %*% as.numeric(bt[tt,])))^2 / (tau2 * sig2[tt,])
        gamma2 <- (2/sig2[tt,]) + ((theta^2) / (tau2 * sig2[tt,]))
        v[tt,] <- rgig(n=1,1/2,delta2,gamma2)
      }
    }
    
    # Step 2b: Sample scaling parameter sigma^2
    if(mode=="qr"){
      if(sv){
        htt <- log(sig2)
        z <- v/sig2
        eps <- Y*NA
        for(tt in 1:T){
          eps[tt,] <- Y[tt,] - X[tt,] %*% as.numeric(bt[tt,]) - alpha
        }
        htt <- jpr(eps=eps,htt=htt,h0m=h0m,h0V=h0V,h_mu=h_mu,h_rho=h_rho,h_sig=h_sig,theta=theta,tau2=tau2,z=z,T=T,c=0.1)
        h_sig <- 1/rgamma(1,5 + (T-1)/2, 0.05 + sum((htt[2:T]-htt[1:(T-1)])^2)/2)
        htt[htt<1e-4] <- 1e-4
        
        sig2 <- exp(htt)
        v <- z*sig2
        tau2sig2v <- tau2 * sig2 * v
      }else{
        eps <- Y*NA
        for(tt in 1:T){
          eps[tt,] <- Y[tt,] - X[tt,] %*% as.numeric(bt[tt,]) - theta*v[tt,] - alpha
        }
        
        ntilde <- n0 + 3*T
        stilde <- s0 + 2*sum(v) + sum((((eps)^2)/(tau2 * v)))
        sig2[,] <- 1/rgamma(1,ntilde/2,stilde/2)
        tau2sig2v <- tau2 * sig2 * v
      }
    }else if(mode=="reg"){
      eps <- Y*NA
      for(tt in 1:T){
        eps[tt,] <- Y[tt,] - X[tt,] %*% as.numeric(bt[tt,]) - alpha
      }
      
      sv_draw <- svsample_fast_cpp(eps, startpara = sv_draw, startlatent = sv_latent, priorspec = sv_priors)
      sv_draw[c("mu", "phi", "sigma", "nu", "rho")] <- as.list(sv_draw$para[, c("mu", "phi", "sigma", "nu", "rho")])
      sv_latent <- sv_draw$latent
      htt <- as.numeric(sv_latent)
      sig2 <- exp(htt)
      tau2sig2v <- tau2 * sig2 * v
    }
    
    # storage
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      bt_store[in.thin,,] <- bt
      shrink_store[in.thin,,] <- steq_var
      sig2_store[in.thin,] <- sig2
      if(fhorz>0){
        bfc <- matrix(NA,fhorz,K)
        yfc <- rep(NA,fhorz)
        
        if(tvp=="dhs"){
          bvt <- log(bv^2)
        }else if(tvp=="shs"){
          bvt <- (1/rgamma(K,0.5,1/(1/rgamma(K,0.5,1)))) * tau.shs # sample from the prior
          bvt[bvt>1] <- 1
          bvt <- log(bvt)
        }else if(tvp=="iG"){
          bvt <- log(steq_var[T,])
        }
        
        sig2fc <- sig2[T,]
        if(sv){
          if(pred.sv){
            fcstsig2_store[in.thin,1] <- sig2fc
            for(hh in 2:fhorz){
              sig2fc <- exp(h_mu + h_rho*log(sig2fc) + sqrt(h_sig)*rnorm(1))
              fcstsig2_store[in.thin,hh] <- sig2fc
            }
          }else{
            fcstsig2_store[in.thin,] <- sig2fc
          }
        }else{
          sig2fc <- sig2[T,]
          fcstsig2_store[in.thin,] <- sig2fc
        }
      
        bfc[1,] <- bt[T,] + exp(bvt/2)*rnorm(K,0,1)
        yfc[1] <- Xout[1,]%*%bfc[1,]
        for(hh in 2:fhorz){
          if(pred.dsp){
            if(tvp=="dhs"){
              bvt <- dhs_se_para[,1] + dhs_se_para[,2]*(bvt-dhs_se_para[,1]) + log(LaplacesDemon::rinvbeta(K,1/2,1/2))
            }
            if(tvp=="shs"){
              bvt <- (1/rgamma(K,0.5,1/(1/rgamma(K,0.5,1)))) * tau.shs # sample from the prior
              bvt[bvt>1] <- 1
              bvt <- log(bvt)
            } 
          } 
          bfc[hh,] <- bfc[hh-1,] + exp(bvt/2)*rnorm(K,0,1)
          yfc[hh] <- Xout[hh,] %*% bfc[hh,] + alpha
        }
        fcst_store[in.thin,] <- yfc
      }
    }
    if(!quiet) setTxtProgressBar(pb, irep)
  }
  
  return(list("y"=Y,"x"=X,
              "bt"=bt_store,"sig2"=sig2_store,"shrink"=shrink_store,
              "fcst"=fcst_store,"fcstsig2"=fcstsig2_store,
              "p"=p))
}

tvpqr.grid <- function(Y,X,Xout=NULL,p=seq(0.05,0.95,by=0.05),cpu=1,tvp="dhs",sv=TRUE,cons.cons=FALSE,fhorz=0,
                      nburn=1000,nsave=1000,thinfac=1,out="mcmc"){
  K <- NCOL(X)
  T <- NROW(Y)
  
  p.grid <- p
  P <- length(p.grid)
  if(cpu>1){
    require(parallel)
    require(doParallel)
    require(foreach)
    if(detectCores()<=cpu){
      cpu <- detectCores()-1
      message("CPU argument equal/exceeds available cores. Using ",cpu," cores.")
    }
    registerDoParallel(cores=cpu)
  }
  
  # parallelization of estimation
  message("Sampling independent TVP-QR models.")
  if(cpu==1){
    p.list <- list()
    for(pp in 1:P){
      p.list[[pp]] <- tvpqr(Y=Y,X=X,Xout=Xout,mode="qr",p=p.grid[pp],tvp=tvp,sv=sv,cons.cons=cons.cons,
                           nburn=nburn,nsave=nsave,thinfac=thinfac,quiet=FALSE)
    }
  }else{
    p.list <- foreach(pp = 1:P) %dopar% {
      tvpqr(Y=Y,X=X,Xout=Xout,mode="qr",p=p.grid[pp],tvp=tvp,sv=sv,
           nburn=nburn,nsave=nsave,thinfac=thinfac,quiet=FALSE)
    }
  }
  message("Finished estimation. Starting post-processing.")
  
  # post processing
  bt_store <- array(NA,dim=c(nsave,P,T,K))
  sig2_store <- array(NA,dim=c(nsave,P,T))
  dimnames(bt_store) <- list(paste0("mcmc",1:nsave),paste0("p",p*100),NULL,NULL)
  dimnames(sig2_store) <- list(paste0("mcmc",1:nsave),paste0("p",p*100),NULL)
  
  for(pp in 1:P){
    bt_store[,paste0("p",p.grid[pp]*100),,] <- p.list[[pp]]$bt
    sig2_store[,paste0("p",p.grid[pp]*100),] <- p.list[[pp]]$sig2
  }
  return(list("bt"=bt_store,"sig2"=sig2_store))
}
