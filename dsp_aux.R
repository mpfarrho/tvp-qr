# Functions from replication code for Dynamic horseshoe
# https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12325

#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
sampleBTF = function(y, obs_sigma_t2, evol_sigma_t2, D = 1, chol0 = NULL){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  T = length(y)
  
  # Linear term:
  linht = y/obs_sigma_t2
  
  # Quadratic terms and solutions are computed differently, depending on D:
  
  if(D == 0){
    # Special case: no differencing
    
    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2)
    postMean = (linht)*postSD^2
    
    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)
    
  } else {
    # All other cases: positive integer differencing (D = 1 or D = 2)
    
    # Quadratic term (D = 1 or D = 2)
    QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2, evol_sigma_t2 = evol_sigma_t2, D = D)
    
    if(!is.null(chol0)){
      # New sampler, based on spam package:
      
      # Sample the states:
      mu = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
                                    Rstruct = chol0))
    } else {
      # Original sampler, based on Matrix package:
      
      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)
      
      # Sample the states:
      mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
      
    }
  }
  
  # And return the states:
  mu
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
sampleBTF_sparse = function(y,
                            obs_sigma_t2,
                            evol_sigma_t2,
                            zero_sigma_t2,
                            D = 1, chol0 = NULL){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  T = length(y)
  
  # Linear term:
  linht = y/obs_sigma_t2
  
  # Quadratic terms and solutions are computed differently, depending on D:
  
  if(D == 0){
    # Special case: no differencing
    
    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2 + 1/zero_sigma_t2)
    postMean = (linht)*postSD^2
    
    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)
    
  } else {
    # All other cases: positive integer differencing (D = 1 or D = 2)
    
    # Quadratic term (D = 1 or D = 2)
    #QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2,
    QHt_Matrix = build_Q(obs_sigma_t2 = 1/(1/obs_sigma_t2 + 1/zero_sigma_t2),
                         evol_sigma_t2 = evol_sigma_t2,
                         D = D)
    
    if(!is.null(chol0)){
      # New sampler, based on spam package:
      
      # Sample the states:
      mu = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
                                    Rstruct = chol0))
    } else {
      # Original sampler, based on Matrix package:
      
      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)
      
      # Sample the states:
      mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
      
    }
  }
  
  # And return the states:
  mu
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
sampleBTF_reg = function(y, X, obs_sigma_t2, evol_sigma_t2, XtX, D = 1, chol0 = NULL){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  # Dimensions of X:
  T = nrow(X); p = ncol(X)
  
  if(D == 1){
    # Lagged version of transposed precision matrix, with zeros as appropriate (needed below)
    t_evol_prec_lag_mat = matrix(0, nr = p, nc = T);
    t_evol_prec_lag_mat[,1:(T-1)] = t(1/evol_sigma_t2[-1,])
    
    # Diagonal of quadratic term:
    Q_diag = matrix(t(1/evol_sigma_t2) + t_evol_prec_lag_mat)
    
    # Off-diagonal of quadratic term:
    Q_off = matrix(-t_evol_prec_lag_mat)[-(T*p)]
    
    # Quadratic term:
    Qevol = bandSparse(T*p, k = c(0,p), diag = list(Q_diag, Q_off), symm = TRUE)
    
    # For checking via direct computation:
    # H1 = bandSparse(T, k = c(0,-1), diag = list(rep(1, T), rep(-1, T)), symm = FALSE)
    # IH = kronecker(as.matrix(H1), diag(p));
    # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
    # print(sum((Qevol - Q0)^2))
    
  } else {
    if(D == 2){
      
      # Lagged x2 version of transposed precision matrix (recurring term)
      t_evol_prec_lag2 = t(1/evol_sigma_t2[-(1:2),])
      
      # Diagonal of quadratic term:
      Q_diag = t(1/evol_sigma_t2)
      Q_diag[,2:(T-1)] = Q_diag[,2:(T-1)] + 4*t_evol_prec_lag2
      Q_diag[,1:(T-2)] = Q_diag[,1:(T-2)] + t_evol_prec_lag2
      Q_diag = matrix(Q_diag)
      
      # Off-diagonal (1) of quadratic term:
      Q_off_1 = matrix(0, nr = p, nc = T);
      Q_off_1[,1] = -2/evol_sigma_t2[3,]
      Q_off_1[,2:(T-1)] = Q_off_1[,2:(T-1)] + -2*t_evol_prec_lag2
      Q_off_1[,2:(T-2)] = Q_off_1[,2:(T-2)] + -2*t_evol_prec_lag2[,-1]
      Q_off_1 = matrix(Q_off_1)
      
      # Off-diagonal (2) of quadratic term:
      Q_off_2 =  matrix(0, nr = p, nc = T); Q_off_2[,1:(T-2)] = t_evol_prec_lag2
      Q_off_2 = matrix(Q_off_2)
      
      # Quadratic term:
      Qevol = bandSparse(T*p, k = c(0, p, 2*p), diag = list(Q_diag, Q_off_1, Q_off_2), symm = TRUE)
      
      # For checking via direct computation:
      # H2 = bandSparse(T, k = c(0,-1, -2), diag = list(rep(1, T), c(0, rep(-2, T-1)), rep(1, T)), symm = FALSE)
      # IH = kronecker(as.matrix(H2), diag(p));
      # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
      # print(sum((Qevol - Q0)^2))
      
    } else stop('sampleBTF_reg() requires D=1 or D=2')
  }
  
  # Quadratic term:
  Qobs = 1/rep(obs_sigma_t2, each = p)*XtX
  Qpost = Qobs + Qevol
  
  # Linear term:
  linht = matrix(t(X*as.numeric(y/obs_sigma_t2))) #matrix(t(X*tcrossprod(y/obs_sigma_t2, rep(1,p))))
  
  if(!is.null(chol0)){
    # Use spam sampler (new version)
    
    # Convert to spam object:
    QHt_Matrix = as.spam.dgCMatrix(as(Qpost, "dgCMatrix"))
    
    # NOTE: reorder (opposite of log-vol!)
    beta = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = QHt_Matrix,
                                    Rstruct = chol0),
                  nrow = T, byrow = TRUE)
  } else {
    # Use original sampler:
    
    # Cholesky:
    chQht_Matrix = Matrix::chol(Qpost)
    
    # NOTE: reorder (opposite of log-vol!)
    beta = matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T*p)), nr = T, byrow = TRUE)
  }
  beta
}
#----------------------------------------------------------------------------
#' (Backfitting) Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
sampleBTF_reg_backfit = function(y, X, beta, obs_sigma_t2, evol_sigma_t2, D = 1){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  # Dimensions of X:
  T = nrow(X); p = ncol(X)
  
  # Sample each predictor curve via backfitting:
  for(j in sample(1:p,p)){
    #for(j in 1:p){
    # Subtract off non-j terms:
    y_nj = y - rowSums(X[,-j]*beta[,-j])
    
    # Linear term:
    linht = y_nj*X[,j]/obs_sigma_t2
    
    # Quadratic terms and solutions are computed differently, depending on D:
    
    if(D == 0){
      # Special case: no differencing
      
      # Posterior SDs and posterior means:
      postSD = 1/sqrt((X[,j]^2)/obs_sigma_t2 + 1/evol_sigma_t2)
      postMean = (linht)*postSD^2
      
      # Sample the states:
      beta[,j] = rnorm(n = T, mean = postMean, sd = postSD)
      
    } else {
      # Quadratic term (D = 1 or D = 2)
      # The likelihood precision term is simply X[,j]^2/obs_sigma_t2, so invert here:
      QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2/X[,j]^2,
                           evol_sigma_t2 = evol_sigma_t2[,j],
                           D = D)
      
      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)
      
      # Sample the states:
      beta[,j] = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
    }
  }
  beta
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
sampleBTF_bspline = function(y, X, obs_sigma2, evol_sigma_t2, XtX_bands, Xty = NULL, D = 1){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  # Dimensions of X:
  T = nrow(X); p = ncol(X)
  
  # Linear term:
  if(is.null(Xty)) Xty = crossprod(X, y)
  linht = 1/obs_sigma2*Xty
  
  # Quadratic terms and solutions are computed differently, depending on D:
  if(D == 0){
    # Special case: no differencing
    QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                            diag = list(XtX_bands$XtX_0/obs_sigma2,
                                        XtX_bands$XtX_1/obs_sigma2,
                                        XtX_bands$XtX_2/obs_sigma2,
                                        XtX_bands$XtX_3/obs_sigma2),
                            symm = TRUE)
  } else {
    # Prior/evoluation quadratic term: can construct directly for D = 1 or D = 2
    if(D == 1){
      QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                              diag = list(XtX_bands$XtX_0/obs_sigma2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0),
                                          XtX_bands$XtX_1/obs_sigma2 + -1/evol_sigma_t2[-1],
                                          XtX_bands$XtX_2/obs_sigma2,
                                          XtX_bands$XtX_3/obs_sigma2),
                              symm = TRUE)
    } else {
      if(D == 2){
        QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                                diag = list(XtX_bands$XtX_0/obs_sigma2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0),
                                            XtX_bands$XtX_1/obs_sigma2 + c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0))),
                                            XtX_bands$XtX_2/obs_sigma2 + 1/evol_sigma_t2[-(1:2)],
                                            XtX_bands$XtX_3/obs_sigma2),
                                symm = TRUE)
      } else stop('sampleBTF_bspline() requires D=0, D=1, or D=2')
    }
  }
  
  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)
  
  # And sample the basis coefficients:
  beta = Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(p))
  
  # Return the sampled basis coefficients:
  beta
}
#----------------------------------------------------------------------------
#' Sample the latent log-volatilities
sampleLogVols = function(h_y, h_prev, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0){
  
  # Compute dimensions:
  h_prev = as.matrix(h_prev) # Just to be sure (T x p)
  n = nrow(h_prev); p = ncol(h_prev)
  
  # Mixture params: mean, variance, and weights
  # Kim, Shephard, Chib (1998) 7-component mixture:
  #m_st  = c(-11.40039, -5.24321, -9.83726, 1.50746,  -0.65098, 0.52478,  -2.35859)
  #v_st2 = c(5.795960,  2.613690, 5.179500, 0.167350, 0.640090, 0.340230, 1.262610)
  #q     = c(0.007300,  0.105560, 0.000020, 0.043950, 0.340010, 0.245660, 0.257500)
  
  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  
  # Add an offset: common for all times, but distict for each j=1,...,p
  yoffset = tcrossprod(rep(1,n),
                       apply(as.matrix(h_y), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  
  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)
  
  # Sample the mixture components
  #z = draw.indicators(res = ystar-h_prev, nmix = list(m = m_st, v = v_st2, p = q))
  z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  
  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = matrix(m_st[z], nr=n); v_st2_all = matrix(v_st2[z], nr=n)
  
  # Joint AWOL sampler for j=1,...,p:
  
  # Constant (but j-specific) mean
  h_mu_all = tcrossprod(rep(1,n), h_mu)
  
  # Constant (but j-specific) AR(1) coef
  h_phi_all = tcrossprod(rep(1,n), h_phi)
  
  # Linear term:
  linht = matrix((ystar - m_st_all - h_mu_all)/v_st2_all)
  
  # Evolution precision matrix (n x p)
  evol_prec_mat = matrix(0, nr = n, nc = p);
  evol_prec_mat[1,] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1,] = 1/h_sigma_eta_t^2;
  
  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = matrix(0, nr = n, nc = p);
  evol_prec_lag_mat[1:(n-1),] = evol_prec_mat[-1,]
  
  # Diagonal of quadratic term:
  Q_diag = matrix(1/v_st2_all +  evol_prec_mat + h_phi_all^2*evol_prec_lag_mat)
  
  # Off-diagonal of quadratic term:
  Q_off = matrix(-h_phi_all*evol_prec_lag_mat)[-(n*p)]
  
  # Quadratic term:
  QHt_Matrix = bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE)
  #QHt_Matrix = as.spam.dgCMatrix(as(bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE),"dgCMatrix"))
  
  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)
  
  # Sample the log-vols:
  hsamp = h_mu_all + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(length(linht))), nr = n)
  #hsamp = h_mu_all +matrix(rmvnorm.canonical(n = 1, b = linht, Q = QHt_Matrix, Rstruct = cholDSP0))
  
  # Return the (uncentered) log-vols
  hsamp
}
#----------------------------------------------------------------------------
#' Sampler evolution error variance parameters
sampleEvolParams = function(omega, evolParams,  sigma_e = 1, evol_error = "DHS"){
  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')
  
  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)
  
  if(evol_error == "DHS") return(sampleDSP(omega, evolParams, sigma_e))
  
  if(evol_error == "HS"){
    
    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset
    
    # Local scale params:
    evolParams$tauLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$xiLambdaj + hsInput2/2), nr = n)
    evolParams$xiLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$tauLambdaj + tcrossprod(rep(1,n), evolParams$tauLambda)), nr = n)
    
    # Global scale params:
    evolParams$tauLambda = rgamma(n = p, shape = 0.5 + n/2, colSums(evolParams$xiLambdaj) + evolParams$xiLambda)
    #evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1/sigma_e^2)
    evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1)
    
    evolParams$sigma_wt = 1/sqrt(evolParams$tauLambdaj)
    
    return(evolParams)
  }
  if(evol_error == "BL"){
    
    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset
    
    # 1/tau_j^2 is inverse-gaussian (NOTE: this is very slow!)
    evolParams$tau_j = matrix(sapply(matrix(hsInput2), function(x){1/sqrt(rig(n = 1,
                                                                              mean = sqrt(evolParams$lambda2*sigma_e^2/x), # already square the input
                                                                              scale = 1/evolParams$lambda2))}), nr = n)
    # Note: should be better priors for lambda2
    evolParams$lambda2 = rgamma(n = 1,
                                shape = 1 + n*p,
                                rate = 2 + sum(evolParams$tau_j^2)/2)
    
    # For Bayesian lasso, scale by sigma_e:
    evolParams$sigma_wt = sigma_e*evolParams$tau_j
    
    return(evolParams)
  }
  if(evol_error == "SV") return(sampleSVparams(omega = omega, svParams = evolParams))
  #if(evol_error == "SV") return(sampleSVparams0(omega = omega, svParams = evolParams))
  if(evol_error == "NIG") {
    evolParams = list(sigma_wt = tcrossprod(rep(1,n),
                                            apply(omega, 2,
                                                  function(x) 1/sqrt(rgamma(n = 1, shape = n/2 + 0.01, rate = sum(x^2)/2 + 0.01)))))
    return(evolParams)
  }
}

#----------------------------------------------------------------------------
#' Sample the dynamic shrinkage process parameters
sampleDSP = function(omega, evolParams, sigma_e = 1, prior_dhs_phi = c(10,2), alphaPlusBeta = 1){
  
  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; dhs_phi = evolParams$dhs_phi; sigma_eta_t = evolParams$sigma_eta_t; sigma_eta_0 = evolParams$sigma_eta_0; dhs_mean0 = evolParams$dhs_mean0
  
  # "Local" number of time points
  ht = as.matrix(ht)
  n = nrow(ht); p = ncol(ht)
  
  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0)
  
  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), dhs_mean)
  
  # Sample AR(1) parameters
  # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  if(!all(dhs_phi == 0) && !all(dhs_phi == 1)) dhs_phi = sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, prior_dhs_phi = prior_dhs_phi)
  
  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), dhs_phi)*ht_tilde[-n, ]       # Residuals
  # sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1)*p, h = alphaPlusBeta, z = eta_t)), nc = p) # Sample
  sigma_eta_t <- matrix(1/sqrt(pgdraw::pgdraw(b=alphaPlusBeta,c=eta_t)), nc=p)
  # sigma_eta_0 = 1/sqrt(rpg(num = p, h = 1, z = ht_tilde[1,]))                # Sample the inital
  sigma_eta_0 <- 1/sqrt(pgdraw::pgdraw(b=1,c=ht_tilde[1,]))
  
  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  if(!all(dhs_phi == 1)){
    if(p > 1){
      # Assume a hierarchy of the global shrinkage params across j=1,...,p
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = dhs_mean0);
      dhs_mean = muSample$dhs_mean
      dhs_mean0 = sampleLogVolMu0(h_mu = dhs_mean, h_mu0 = dhs_mean0, dhs_mean_prec_j = muSample$dhs_mean_prec_j, h_log_scale = log(sigma_e^2))
    } else {
      # p = 1
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = log(sigma_e^2));
      dhs_mean = dhs_mean0 = muSample$dhs_mean # save dhs_mean0 = dhs_mean for coding convenience later
    }
  } else {dhs_mean = rep(0, p); dhs_mean0 = 0} # When RW for log-vols, fix unconditional mean for identifiability
  
  # Evolution error SD:
  sigma_wt = exp(ht/2)
  
  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}

#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters
sampleSVparams = function(omega, svParams){
  
  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)
  
  for(j in 1:p){
    # First, check for numerical issues:
    svInput = omega[,j]; #if(all(svInput==0)) {svInput = 10^-8} else svInput = svInput + sd(svInput)/10^8
    #hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    
    # Sample the SV parameters:
    svsamp = stochvol::svsample2(svInput,
                                 startpara = list(
                                   mu = svParams$svParams[1,j],
                                   phi = svParams$svParams[2,j],
                                   sigma = svParams$svParams[3,j]),
                                 startlatent = svParams$ht[,j])# ,priorphi = c(10^4, 10^4));
    # Update the parameters:
    svParams$svParams[,j] = svsamp$para;
    svParams$ht[,j] = svsamp$latent
  }
  # Finally, up the evolution error SD:
  svParams$sigma_wt = exp(svParams$ht/2)
  
  # Check for numerically large values:
  svParams$sigma_wt[which(svParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3
  
  return(svParams)
}

#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters using same functions as DHS prior
sampleSVparams0 = function(omega, svParams){
  
  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)
  
  # Store the parameters locally:
  ht = as.matrix(svParams$ht); sv_mean = svParams$svParams[1,]; sv_phi = svParams$svParams[2,]; sv_sigma = svParams$svParams[3,]
  
  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = sv_mean, h_phi = sv_phi,
                     h_sigma_eta_t = matrix(rep(sv_sigma, each = n-1), nrow = n-1), h_sigma_eta_0 = sv_sigma) # New part
  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), sv_mean)
  
  # Sample the AR(1) parameters:
  sv_phi = sampleAR1(h_yc = ht_tilde, h_phi = sv_phi,
                     h_sigma_eta_t = matrix(rep(sv_sigma, each = n-1), nrow = n-1),
                     prior_dhs_phi = c(10, 2))
  # Sample the evolution error SD of log-vol
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), sv_phi)*ht_tilde[-n, ]       # Residuals
  sv_sigma = apply(eta_t, 2, function(x)
    1/sqrt(rgamma(n = 1, shape = length(x)/2 + 0.01, rate = sum(x^2)/2 + 0.01)))
  
  # Sample the mean parameters:
  y_mu = (ht[-1,] - tcrossprod(rep(1,n-1), sv_phi)*ht[-n,])/matrix(rep(sv_sigma, each = n-1), nrow = n-1);
  x_mu = tcrossprod(rep(1,n-1), 1 - sv_phi)/matrix(rep(sv_sigma, each = n-1), nrow = n-1)
  postSD = 1/sqrt(colSums(x_mu^2) + 1/10^2)
  postMean = (colSums(x_mu*y_mu))*postSD^2
  sv_mean = rnorm(n = p, mean = postMean, sd = postSD)
  
  # Evolution error SD:
  sigma_wt = exp(ht/2)
  
  # Update:
  svParams$sigma_wt = sigma_wt; svParams$ht = ht;
  svParams$svParams[1,] = sv_mean;
  svParams$svParams[2,] = sv_phi;
  svParams$svParams[3,] = sv_sigma
  
  # And return:
  return(svParams)
}

#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.
  
  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }
  
  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1
  
  # Find the log density at the initial point, if not already known.
  
  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }
  
  # Determine the slice level, in log terms.
  
  logy <- gx0 - rexp(1)
  
  # Find the initial interval to sample from.
  
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  {
    x1 <- runif(1,L,R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)
    
    if (gx1>=logy) break
    
    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.
  
  attr(x1,"log.density") <- gx1
  return (x1)
  
}

#----------------------------------------------------------------------------
#' Sample the AR(1) coefficient(s)
sampleAR1 = function(h_yc, h_phi, h_sigma_eta_t, prior_dhs_phi = NULL){
  
  # Compute dimensions:
  n = nrow(h_yc); p = ncol(h_yc)
  
  # Loop over the j=1:p
  for(j in 1:p){
    
    # Compute "regression" terms for dhs_phi_j:
    y_ar = h_yc[-1,j]/h_sigma_eta_t[,j] # Standardized "response"
    x_ar = h_yc[-n,j]/h_sigma_eta_t[,j] # Standardized "predictor"
    
    # Using Beta distribution:
    if(!is.null(prior_dhs_phi)){
      
      # Check to make sure the prior params make sense
      if(length(prior_dhs_phi) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')
      
      dhs_phi01 = (h_phi[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])
      
      # Slice sampler when using Beta prior:
      dhs_phi01 = uni.slice(dhs_phi01, g = function(x){
        -0.5*sum((y_ar - (2*x - 1)*x_ar)^2) +
          dbeta(x, shape1 = prior_dhs_phi[1], shape2 = prior_dhs_phi[2], log = TRUE)
      }, lower = 0, upper = 1)[1]#}, lower = 0.005, upper = 0.995)[1] #
      
      h_phi[j] = 2*dhs_phi01 - 1
      
    } else {
      # For h_phi ~ Unif(-1, 1), the posterior is truncated normal
      h_phi[j] = rtrunc(n = 1, spec = 'norm',
                        a = -1, b = 1,
                        mean = sum(y_ar*x_ar)/sum(x_ar^2),
                        sd = 1/sqrt(sum(x_ar^2)))
    }
  }
  h_phi
}

#----------------------------------------------------------------------------
#' Sample the AR(1) unconditional means
sampleLogVolMu = function(h, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0, h_log_scale = 0){
  
  # Compute "local" dimensions:
  n = nrow(h); p = ncol(h)
  
  # Sample the precision term(s)
  # dhs_mean_prec_j = rpg(num = p, h = 1, z = h_mu - h_log_scale)
  dhs_mean_prec_j <- pgdraw::pgdraw(b=1,c=h_mu-h_log_scale)
  
  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (h[-1,] - tcrossprod(rep(1,n-1), h_phi)*h[-n,])/h_sigma_eta_t;
  x_mu = tcrossprod(rep(1,n-1), 1 - h_phi)/h_sigma_eta_t
  
  # Include the initial sd?
  #if(!is.null(h_sigma_eta_0)){y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu); x_mu = rbind(1/h_sigma_eta_0, x_mu)}
  y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu);
  x_mu = rbind(1/h_sigma_eta_0, x_mu)
  
  # Posterior SD and mean:
  postSD = 1/sqrt(colSums(x_mu^2) + dhs_mean_prec_j)
  postMean = (colSums(x_mu*y_mu) + h_log_scale*dhs_mean_prec_j)*postSD^2
  dhs_mean = rnorm(n = p, mean = postMean, sd = postSD)
  
  list(dhs_mean = dhs_mean, dhs_mean_prec_j = dhs_mean_prec_j)
}

#----------------------------------------------------------------------------
#' Sample the mean of AR(1) unconditional means
sampleLogVolMu0 = function(h_mu, h_mu0, dhs_mean_prec_j, h_log_scale = 0){
  
  # dhs_mean_prec_0 = rpg(num = 1, h = 1, z = h_mu0 - h_log_scale)
  dhs_mean_prec_0 <- pgdraw::pgdraw(b=1,c=h_mu0-h_log_scale)
  
  # Sample the common mean parameter:
  postSD = 1/sqrt(sum(dhs_mean_prec_j) + dhs_mean_prec_j)
  postMean = (sum(dhs_mean_prec_j*h_mu) + dhs_mean_prec_j*h_log_scale)*postSD^2
  rnorm(n = 1, mean = postMean, sd = postSD)
}

#----------------------------------------------------------------------------
#' Sample the parameters for the initial state variance
sampleEvol0 = function(mu0, evolParams0, commonSD = FALSE, A = 1){
  
  # Store length locally:
  p = length(mu0)
  
  # For numerical stability:
  mu02offset = any(mu0^2 < 10^-16)*max(10^-8, mad(mu0)/10^6)
  mu02 = mu0^2 + mu02offset
  
  if(commonSD){
    # (Common) standard deviations:
    evolParams0$sigma_w0 = rep(1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(mu02)/2 + evolParams0$px_sigma_w0[1])), p)
    
    # (Common) paramater expansion:
    evolParams0$px_sigma_w0 = rep(rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0[1]^2 + 1/A^2), p)
    
  } else {
    # (Distinct) standard deviations:
    evolParams0$sigma_w0 = 1/sqrt(rgamma(n = p, shape = 1/2 + 1/2, rate = mu02/2 + evolParams0$px_sigma_w0))
    
    # (Distict) paramater expansion:
    #evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/A^2)
    evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/evolParams0$sigma_00^2)
    
    # Global standard deviations:
    evolParams0$sigma_00 = 1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(evolParams0$px_sigma_w0) + evolParams0$px_sigma_00))
    
    # (Global) parameter expansion:
    evolParams0$px_sigma_00 = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_00^2 + 1/A^2)
  }
  
  # And return the list:
  evolParams0
}

#----------------------------------------------------------------------------
#' Sample a Gaussian vector using the fast sampler of BHATTACHARYA et al.

sampleFastGaussian = function(Phi, Ddiag, alpha){
  
  # Dimensions:
  Phi = as.matrix(Phi); n = nrow(Phi); p = ncol(Phi)
  
  # Step 1:
  u = rnorm(n = p, mean = 0, sd = sqrt(Ddiag))
  delta = rnorm(n = n, mean = 0, sd = 1)
  
  # Step 2:
  v = Phi%*%u + delta
  
  # Step 3:
  w = solve(crossprod(sqrt(Ddiag)*t(Phi)) + diag(n), #Phi%*%diag(Ddiag)%*%t(Phi) + diag(n)
            alpha - v)
  
  # Step 4:
  theta =  u + Ddiag*crossprod(Phi, w)
  
  # Return theta:
  theta
}

#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
initEvolParams = function(omega, evol_error = "DHS"){
  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") ||(evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')
  
  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)
  
  if(evol_error == "DHS") return(initDHS(omega))
  
  if(evol_error == "HS"){
    tauLambdaj = 1/omega^2;
    xiLambdaj = 1/(2*tauLambdaj); tauLambda = 1/(2*colMeans(xiLambdaj)); xiLambda = 1/(tauLambda + 1)
    
    # Parameters to store/return:
    return(list(sigma_wt = 1/sqrt(tauLambdaj), tauLambdaj = tauLambdaj, xiLambdaj = xiLambdaj, tauLambda = tauLambda, xiLambda = xiLambda))
  }
  if(evol_error == "BL"){
    tau_j = abs(omega); lambda2 = mean(tau_j)
    return(list(sigma_wt = tau_j, tau_j = tau_j, lambda2 = lambda2))
  }
  if(evol_error == "SV") return(initSV(omega))
  if(evol_error == "NIG") return(list(sigma_wt = tcrossprod(rep(1,n), apply(omega, 2, function(x) sd(x, na.rm=TRUE)))))
}

#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
initDHS = function(omega){
  
  # "Local" number of time points
  omega = as.matrix(omega)
  n = nrow(omega); p = ncol(omega)
  
  # Initialize the log-volatilities:
  ht = log(omega^2 + 0.0001)
  
  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = apply(ht, 2, function(x){
    params = try(arima(x, c(1,0,0))$coef, silent = TRUE); if(class(params) == "try-error") params = c(0.8, mean(x)/(1 - 0.8))
    params
  })
  dhs_mean = arCoefs[2,]; dhs_phi = arCoefs[1,]; dhs_mean0 = mean(dhs_mean)
  
  # Initialize the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nr = n-1, nc = p)
  sigma_eta_0 = rep(pi, p) # Initial value
  
  # Evolution error SD:
  sigma_wt = exp(ht/2)
  
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}

#----------------------------------------------------------------------------
#' Initialize the stochastic volatility parameters
initSV = function(omega){
  
  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)
  
  # log-volatility:
  ht = log(omega^2 + 0.0001)
  
  # AR(1) pararmeters: check for error in initialization too
  svParams = apply(ht, 2, function(x){
    ar_fit = try(arima(x, c(1,0,0)), silent = TRUE)
    if(class(ar_fit) != "try-error") {
      params = c(ar_fit$coef[2], ar_fit$coef[1], sqrt(ar_fit$sigma2))
    } else params = c(mean(x)/(1 - 0.8),0.8, 1)
    params
  }); rownames(svParams) = c("intercept", "ar1", "sig")
  
  # SDs, log-vols, and other parameters:
  return(list(sigma_wt = exp(ht/2), ht = ht, svParams = svParams))
}

#----------------------------------------------------------------------------
#' Initialize the parameters for the initial state variance
initEvol0 = function(mu0, commonSD = TRUE){
  p = length(mu0)
  
  # Common or distict:
  if(commonSD) {
    sigma_w0 = rep(mean(abs(mu0)), p)
  } else  sigma_w0 = abs(mu0)
  
  # Initialize at 1 for simplicity:
  px_sigma_w0 = rep(1, p)
  
  sigma_00 = px_sigma_00 = 1
  
  list(sigma_w0 = sigma_w0, px_sigma_w0 = px_sigma_w0, sigma_00 = sigma_00, px_sigma_00 = px_sigma_00)
}

#----------------------------------------------------------------------------
#' Compute X'X
build_XtX = function(X){
  # Store the dimensions:
  T = nrow(X); p = ncol(X)
  
  # Store the matrix
  XtX = bandSparse(T*p, k = 0, diag = list(rep(1,T*p)), symm = TRUE)
  
  t.seq.p = seq(1, T*(p+1), by = p)
  
  for(t in 1:T){
    t.ind = t.seq.p[t]:(t.seq.p[t+1]-1)
    XtX[t.ind, t.ind] = tcrossprod(matrix(X[t,]))
  }
  XtX
}

#----------------------------------------------------------------------------
#' Compute initial Cholesky decomposition for Bayesian Trend Filtering
initChol.spam = function(T, D = 1){
  
  # Random initialization
  QHt_Matrix = build_Q(obs_sigma_t2 = abs(rnorm(T)),
                       evol_sigma_t2 = abs(rnorm(T)),
                       D = D)
  
  # And return the Cholesky piece:
  chQht_Matrix0 = chol.spam(as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")))
  
  chQht_Matrix0
}

#----------------------------------------------------------------------------
#' Compute initial Cholesky decomposition for TVP Regression
initCholReg.spam = function(obs_sigma_t2, evol_sigma_t2, XtX, D = 1){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  # Dimensions of X:
  T = nrow(evol_sigma_t2); p = ncol(evol_sigma_t2)
  
  if(D == 1){
    # Lagged version of transposed precision matrix, with zeros as appropriate (needed below)
    t_evol_prec_lag_mat = matrix(0, nr = p, nc = T);
    t_evol_prec_lag_mat[,1:(T-1)] = t(1/evol_sigma_t2[-1,])
    
    # Diagonal of quadratic term:
    Q_diag = matrix(t(1/evol_sigma_t2) + t_evol_prec_lag_mat)
    
    # Off-diagonal of quadratic term:
    Q_off = matrix(-t_evol_prec_lag_mat)[-(T*p)]
    
    # Quadratic term:
    Qevol = bandSparse(T*p, k = c(0,p), diag = list(Q_diag, Q_off), symm = TRUE)
    
    # For checking via direct computation:
    # H1 = bandSparse(T, k = c(0,-1), diag = list(rep(1, T), rep(-1, T)), symm = FALSE)
    # IH = kronecker(as.matrix(H1), diag(p));
    # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
    # print(sum((Qevol - Q0)^2))
    
  } else {
    if(D == 2){
      
      # Lagged x2 version of transposed precision matrix (recurring term)
      t_evol_prec_lag2 = t(1/evol_sigma_t2[-(1:2),])
      
      # Diagonal of quadratic term:
      Q_diag = t(1/evol_sigma_t2)
      Q_diag[,2:(T-1)] = Q_diag[,2:(T-1)] + 4*t_evol_prec_lag2
      Q_diag[,1:(T-2)] = Q_diag[,1:(T-2)] + t_evol_prec_lag2
      Q_diag = matrix(Q_diag)
      
      # Off-diagonal (1) of quadratic term:
      Q_off_1 = matrix(0, nr = p, nc = T);
      Q_off_1[,1] = -2/evol_sigma_t2[3,]
      Q_off_1[,2:(T-1)] = Q_off_1[,2:(T-1)] + -2*t_evol_prec_lag2
      Q_off_1[,2:(T-2)] = Q_off_1[,2:(T-2)] + -2*t_evol_prec_lag2[,-1]
      Q_off_1 = matrix(Q_off_1)
      
      # Off-diagonal (2) of quadratic term:
      Q_off_2 =  matrix(0, nr = p, nc = T); Q_off_2[,1:(T-2)] = t_evol_prec_lag2
      Q_off_2 = matrix(Q_off_2)
      
      # Quadratic term:
      Qevol = bandSparse(T*p, k = c(0, p, 2*p), diag = list(Q_diag, Q_off_1, Q_off_2), symm = TRUE)
      
      # For checking via direct computation:
      # H2 = bandSparse(T, k = c(0,-1, -2), diag = list(rep(1, T), c(0, rep(-2, T-1)), rep(1, T)), symm = FALSE)
      # IH = kronecker(as.matrix(H2), diag(p));
      # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
      # print(sum((Qevol - Q0)^2))
      
    } else stop('Requires D=1 or D=2')
  }
  
  Qobs = 1/rep(obs_sigma_t2, each = p)*XtX
  Qpost = Qobs + Qevol
  
  # New version (NOTE: reorder; opposite of log-vol!)
  QHt_Matrix = as.spam.dgCMatrix(as(Qpost, "dgCMatrix"))
  
  # And return the Cholesky piece:
  chQht_Matrix0 = chol.spam(QHt_Matrix)
  
  chQht_Matrix0
}

#----------------------------------------------------------------------------
#' Compute the quadratic term in Bayesian trend filtering
build_Q = function(obs_sigma_t2, evol_sigma_t2, D = 1){
  
  if(!(D == 1 || D == 2)) stop('build_Q requires D = 1 or D = 2')
  
  T = length(evol_sigma_t2)
  
  # For reference: first and second order difference matrices (not needed below)
  #H1 = bandSparse(T, k = c(0,-1), diag = list(rep(1, T), rep(-1, T)), symm = FALSE)
  #H2 = bandSparse(T, k = c(0,-1, -2), diag = list(rep(1, T), c(0, rep(-2, T-1)), rep(1, T)), symm = FALSE)
  
  # Quadratic term: can construct directly for D = 1 or D = 2 using [diag(1/obs_sigma_t2, T) + (t(HD)%*%diag(1/evol_sigma_t2, T))%*%HD]
  if(D == 1){
    # D = 1 case:
    Q = bandSparse(T, k = c(0,1),
                   diag = list(1/obs_sigma_t2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0),
                               -1/evol_sigma_t2[-1]),
                   symm = TRUE)
  } else {
    # D = 2 case:
    Q = bandSparse(T, k = c(0,1,2),
                   diag = list(1/obs_sigma_t2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0),
                               c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0))),
                               1/evol_sigma_t2[-(1:2)]),
                   symm = TRUE)
  }
  Q
}

#----------------------------------------------------------------------------
#' Compute Non-Zeros (Signals)
getNonZeros = function(post_evol_sigma_t2, post_obs_sigma_t2 = NULL){
  
  # Posterior distribution of shrinkage parameters in (0,1)
  if(is.null(post_obs_sigma_t2)){
    post_kappa = 1/(1 + post_evol_sigma_t2)
  } else {
    
    # Check: if p > 1, then adjust the dimension of post_obs_sigma_t2
    if(length(dim(post_evol_sigma_t2)) > 2) post_obs_sigma_t2 = array(rep(post_obs_sigma_t2, times = dim(post_evol_sigma_t2)[3]), dim(post_evol_sigma_t2))
    
    post_kappa = 1/(1 + post_evol_sigma_t2/post_obs_sigma_t2)
  }
  
  # Indices of non-zeros:
  non_zero = which(colMeans(post_kappa) < 1/2, arr.ind = TRUE)
  
  # Return:
  non_zero
}

#----------------------------------------------------------------------------
#' Sample components from a discrete mixture of normals
ncind = function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}


