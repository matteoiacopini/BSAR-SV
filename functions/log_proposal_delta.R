log_proposal_delta <- function(delta, nu){ #not smapling, evaluating
  log_dens = ddirichlet_R(x = matrix(data = delta, ncol = length(delta)), alpha = nu, log=TRUE)
  return(log_dens)
}