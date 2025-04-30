sample_rho_sv <- function(rho, beta, h_t, delta, y, X_t,
                           min_network, plus_network,
                           rho_params, prior_rho="beta_-11"){
   # Compute dimensions:
   n = nrow(y)
   p = ncol(y)
   
   # Loop over the j=1:p
   for(j in 1:p){
      
      # Using Beta distribution:
      if(!is.null(rho_params)){
         # Check to make sure the prior params make sense
         if(length(rho_params) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')
         
         if (prior_rho == "beta_-11"){
            u = (rho[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])
         }else if(prior_rho == "beta_01"){
            u = rho[j]         # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])
         }
         
         ##################  SLICE SAMPLER WHEN USING BETA PRIOR ##################
         u = uni.slice(x0 = u, g = function(x){
            
            T = nrow(y); num_countries = ncol(y)
            rho_star = rho
            if(prior_rho == "beta_-11"){
               rho_star[j] = 2*x - 1
            }else if(prior_rho == "beta_01"){
               rho_star[j] = x
            }
            
            A = list()
            summation_term = 0
            for(t in 1:T){
               A[[t]] = diag(num_countries) - (diag(rho_star) %*% (delta[1]*min_network[[t]] 
                                                                   + delta[2]*plus_network[[t]]))
               e_temp = A[[t]] %*% unlist(y[t,]) - X_t[,,t] %*% beta
               Sigma_inv = diag(exp(-h_t[t,]))
               summation_term = (summation_term + log(det(A[[t]])) - (t(e_temp) %*% Sigma_inv %*% e_temp)/2)
            }
            log_llik = summation_term
            log_prior = dbeta(x, shape1 = rho_params[1], shape2 = rho_params[2], log = TRUE)
            log_dens = log_llik + log_prior
            #cat("log llik: \n", log_llik, "\n")
            #cat("log prior: \n", log_prior, "\n")
            #print(log_dens)
            return(log_dens)
            #cat("log density: \n", log_dens, "\n")
         }, w=1, m=Inf, lower = 0, upper = 1, gx0 = NULL)[1]
         
         if (prior_rho == "beta_-11"){
            rho[j] = 2*u - 1
         }else if (prior_rho == "beta_01"){
            rho[j] = u
         }
      }
   }
   return(rho)
}
