### This function is needed if we want to have an unknown rho/phi

# expcov <- nimbleFunction(     
#   run = function(dists = double(2), rho = double(0), sigma = double(0)) {
#     returnType(double(2))
#     n <- dim(dists)[1]
#     result <- matrix(nrow = n, ncol = n, init = FALSE)
#     sigma2 <- sigma*sigma
#     
#     for(i in 1:(n-1)){
#       for(j in (i+1):n){
#         temp <- sigma2*exp(-dists[i,j]/rho)
#         result[i, j] <- temp
#         result[j, i] <- temp
#       }
#     }
#     for(i in 1:(n)){
#       result[i, i] <- sigma2
#     }
#     return(result)
#   })
# 
# cExpcov <- compileNimble(expcov)

### Here, beta represents beta* discussed in the supplement, the product of alpha_k and \beta_{k,j}

nimble_code1 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)

  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  # for(i in 1:p_sigma){
  #   beta_sigma[i] ~ dnorm(0, var = 100)
  # }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  sigma2 ~ dinvgamma(1, 1)
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] 
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = sigma2)
  }
  
})

nimble_code2 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  var_out[1:n] <- exp((X_sigma[1:n,1:p_sigma] %*% beta_sigma[1:p_sigma])[1:n,1])
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] 
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})

nimble_code3 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] 
    var_out[i] <- exp(beta_sigma[1] + mu[i] * beta_sigma[2] + mu[i]^2 * beta_sigma[3] + mu[i]^3 * beta_sigma[4])
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})


nimble_code4 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  # for(i in 1:p_sigma){
  #   beta_sigma[i] ~ dnorm(0, var = 100)
  # }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  sigma2 ~ dinvgamma(1, 1)
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + abs(psi[row_ind[i]] - psi[col_ind[i]])
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = sigma2)
  }
  
})

nimble_code5 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  var_out[1:n] <- exp((X_sigma[1:n,1:p_sigma] %*% beta_sigma[1:p_sigma])[1:n,1])
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + abs(psi[row_ind[i]] - psi[col_ind[i]])
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})


nimble_code6 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + abs(psi[row_ind[i]] - psi[col_ind[i]])
    var_out[i] <- exp(beta_sigma[1] + mu[i] * beta_sigma[2] + mu[i]^2 * beta_sigma[3] + mu[i]^3 * beta_sigma[4])
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})


nimble_code7 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  # for(i in 1:p_sigma){
  #   beta_sigma[i] ~ dnorm(0, var = 100)
  # }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  sigma2 ~ dinvgamma(1, 1)
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + (psi[row_ind[i]] - psi[col_ind[i]])^2
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = sigma2)
  }
  
})

nimble_code8 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  var_out[1:n] <- exp((X_sigma[1:n,1:p_sigma] %*% beta_sigma[1:p_sigma])[1:n,1])
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + (psi[row_ind[i]] - psi[col_ind[i]])^2
    
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})




nimble_code9 <- nimbleCode({
  
  beta_0 ~ dnorm(0, var = 100)
  sig2_psi ~ dinvgamma(1, 1)
  prec_use[1:n_loc, 1:n_loc] <- R_inv[1:n_loc, 1:n_loc] / sig2_psi
  psi[1:n_loc] ~ dmnorm(zeros[1:n_loc], prec = prec_use[1:n_loc, 1:n_loc])
  
  for(i in 1:p){
    log(beta[i]) ~ dnorm(0, var = 100)
  }
  
  for(i in 1:p_sigma){
    beta_sigma[i] ~ dnorm(0, var = 100)
  }
  
  linpred[1:n] <- (x[1:n, 1:p] %*% beta[1:p])[1:n,1]
  
  for(i in 1:n){
    mu[i] <- beta_0 + linpred[i] + (psi[row_ind[i]] - psi[col_ind[i]])^2
    var_out[i] <- exp(beta_sigma[1] + mu[i] * beta_sigma[2] + mu[i]^2 * beta_sigma[3] + mu[i]^3 * beta_sigma[4])
    censored[i] ~ dinterval(log_V[i], c[i])
    log_V[i] ~ dnorm(mu[i], var = var_out[i])
  }
  
})