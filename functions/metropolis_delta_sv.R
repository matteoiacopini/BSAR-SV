metropolis_delta_sv <- function(start, beta, rho, h_t, nu, X_t, y, min_network, plus_network, print_out=FALSE){
  
  theta_ini = matrix(data = start, nrow = 1, ncol = 2)
  repeat {
    theta_star <- rdirichlet(1, alpha = nu)
    if(max(theta_star) < 0.9) break
  }
  #print(theta_star)
  
  ############ TESTS FOR THETA_STAR ##############
  # print(theta_ini)
  # print(nrow(theta_ini))
  # print(ncol(theta_ini))
  # cat("typeof theta_ini:", typeof(theta_ini), "\n")
  # 
  # print(theta_star)
  # print(nrow(theta_star))
  # print(ncol(theta_star))
  # cat("typeof theta_star:", typeof(theta_star), "\n")
  #theta_star = c(0.4,0.3,0.2,0.1)
  # cat("Theta_ini: \n", theta_ini, "\n")
  # cat("Theta_star: \n", theta_star, "\n")
  
  # repeat {
  #   theta_star <- rdirichlet(1, alpha = nu)
  #   #theta_star <- rdirichlet(1, alpha = theta_ini*10)
  #   #print(theta_star)
  #   if(min(theta_star > 0)) break
  # }
  #print(theta_star)
  
  ############# COMPUTE A AND B ################
  
  b = (log_target_dens_delta_sv(delta = theta_star, beta, h_t,
                                      rho, X_t, nu, y,
                                      min_network, plus_network) 
       + log_proposal_delta(theta_ini, nu))
  
  a = (log_target_dens_delta_sv(delta = theta_ini, beta, h_t,
                                      rho, X_t, nu, y,
                                      min_network, plus_network)
       + log_proposal_delta(theta_star, nu))
  
  H = min(1,exp(b-a)) # Take ratio like this since using logs
  #cat("H:", round(H), "\n")
  runif = runif(1)
  if(runif < H){ 
    theta <- theta_star
    if (print_out){
       cat("Accepted! \n old theta: \n", theta_ini, "\n", "new theta: \n", theta_star, "\n")
       cat("The old posterior density was:", a, "\n")
       cat("And the new posterior density is:", b, "\n")
       cat("H = ", H, "\n")
       cat("The uniform draw was: \n", runif, "\n")
    }
  } else {
    theta <- theta_ini
  }
  return(theta)
}