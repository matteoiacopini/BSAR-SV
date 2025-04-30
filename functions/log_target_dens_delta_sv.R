log_target_dens_delta_sv <- function(delta, beta, h_t,
                                           rho, X_t, nu, y,
                                           min_network, plus_network){
  T = nrow(y)
  num_countries = ncol(y)
  A = list()
  
  summation_term = 0
  for(t in 1:T){
    A[[t]] = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network[[t]] 
                                                   + delta[2]*plus_network[[t]]))
    e_temp = A[[t]] %*% unlist(y[t,])- X_t[,,t] %*% beta
    #e_temp = A[[t]] %*% y[t,] - alpha - x[[t]] %*% beta_tilde
    Sigma_inv = diag(exp(-h_t[t,]))
    summation_term = (summation_term + log(det(A[[t]])) 
                      - (t(e_temp) %*% Sigma_inv %*% e_temp)/2)
  }
  log_llik = summation_term
  log_prior = ddirichlet_R(x = delta, alpha = nu, log = TRUE)
  #print(log_llik)
  #print(log_prior)
  log_posterior = log_prior + log_llik
  #print(log_posterior)
  return(log_posterior)
}


