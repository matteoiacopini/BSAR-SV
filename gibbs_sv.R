# MCMC for Bayesian SAR-SV model with:
#  1) Time-varying (observed) network
#  2) Multi-layer networks
#  3) Multiple rho
#  4) SV for structural variance

gibbs_sv = function(nsave=1000, nburn=100, nshow=100,
                     ft, nu, rho_params, y,
                     min_network, plus_network, prior_rho,
                     fixed.effects = FALSE,
                     forec.horiz = 0,
                     min_network.pred = NA, plus_network.pred = NA, ft.pred = NA){
   
   niter = nsave + nburn
   
   num_countries = ncol(y)
   ft = as.matrix(ft)
   ft.pred = as.matrix(ft.pred)
   T   = nrow(ft)
   num_factors = ncol(ft)
   
   
   sv_priors <- list()
   for(j in seq_len(num_countries)){
      sv_priors[[j]] <- specify_priors(
         phi = sv_beta(shape1 = 50, shape2 = 1.5),   # mean = 0.94,  sd = 0.05
         mu  = sv_normal(mean = 0.0, sd = 1.0),
         sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
         nu  = sv_infinity(),
         rho = sv_constant(0)
      ) 
   }
   
   #%%%%%%%%%%%% Allocate variables %%%%%%%%%%%%
   if (fixed.effects){
      num_beta = (num_countries * num_factors  +  num_countries)
   }else{
      num_beta = (num_countries * num_factors  +  1)
   }
   mu_beta    = rep(0,num_beta)
   Sigma_beta = diag(num_beta)*2
   A0_beta = solve(Sigma_beta)
   b0_beta = solve(Sigma_beta) %*% mu_beta
   muh  = matrix(NA, num_countries,1)
   phih = matrix(NA, num_countries,1)
   sigh = matrix(NA, num_countries,1)
   
   y_star = matrix(nrow = T, ncol = num_countries)
   X_t = array(0,c(num_countries,num_beta,T))
   h_t = array(0,c(T,num_countries))
   for(t in 1:T){
      if (fixed.effects){
         X_t[,,t] = cbind(diag(num_countries), kronecker(t(as.matrix(ft[t,])), diag(num_countries)))
      }else{
         X_t[,,t] = cbind(matrix(1,num_countries,1), kronecker(t(as.matrix(ft[t,])), diag(num_countries)))
      }
   }
   
   
   #%%%%%%%%%%%% Set up storage matrices %%%%%%%%%%%%
   beta.store  = array(NA, c(nsave, num_beta))
   ht.store    = array(NA, c(nsave, T,num_countries))
   delta.store = array(NA, c(nsave, 2))
   rho.store   = array(NA, c(nsave, num_countries))
   muh.store   = array(NA, c(nsave, num_countries))
   phih.store  = array(NA, c(nsave, num_countries))
   sigh.store  = array(NA, c(nsave, num_countries))
   accepts.store = array(NA, c(nsave, 1))
   y.forec.store = NA
   if (forec.horiz > 0){
      y.forec.store = array(NA, c(nsave, num_countries, forec.horiz))
   }
   
   
   #%%%%%%%%%%%% Set INITIAL VALUES %%%%%%%%%%%%
   sv.idio.draw <- list()
   for (j in seq_len(num_countries)){
      h_t[,j] = log(var(y[,j]))
      sv.idio.draw[[j]] <- list(mu = -10, phi = 0.9, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
   }
   
   beta  = matrix( mvrnorm(n = 1, mu = rep(0,num_beta), Sigma = diag(num_beta)),  nrow = num_beta, ncol = 1 )
   delta = rep(0.5, 2)
   rho   = rep(0.5, num_countries)
   
   A = list()
   for(t in 1:T){
      A[[t]] = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network[[t]] + delta[2]*plus_network[[t]]))
      y_star[t,] = A[[t]] %*% unlist(y[t,])
   }
   accepts = 0
   
   #%%%%%%%%%%%% MAIN loop %%%%%%%%%%%%
   start.time = Sys.time()
   for(irep in 1:niter){
      
      ##### BETA #####
      A1_beta = 0; b1_beta = 0
      for (t in 1:T){
         H_t_inv = diag(exp(-h_t[t,]))
         XSi = crossprod(X_t[,,t],H_t_inv)
         A1_beta = A1_beta + XSi %*% X_t[,,t] 
         b1_beta = b1_beta + XSi %*% y_star[t,]
      }
      An_beta = A0_beta + A1_beta    # posterior precision matrix
      bn_beta = b0_beta + b1_beta
      # beta = matrix( mvrnorm(n=1, mu = solve(An_beta) %*% bn_beta, Sigma = solve(An_beta)) )
      betahat = solve(An_beta,bn_beta)
      beta = betahat + solve(chol(An_beta), matrix(rnorm(num_beta)))
      
      # G_Sigma = G' * kron_I_Sigma_eps_inv;
      # Kh = H'*(S\H) + G_Sigma*G;                   % posterior precision matrix
      # xhat = Kh \ (H'*(S\a0) + G_Sigma*stackY);    % posterior mean
      # thetat_all = xhat  +  chol(Kh,'lower')'\randn(ntheta*T,1);
      
      
      ##### h_t #####
      e_temp = matrix(0, T, num_countries)
      for(t in 1:T){
         e_temp[t,] = A[[t]] %*% unlist(y[t,]) - X_t[,,t] %*% beta
      }
      # sample the idiosyncratic volatility (H)
      shocks = e_temp
      for (j in seq_len(num_countries)){
         sv.idio.draw.i <- svsample_fast_cpp(shocks[ ,j], startpara= sv.idio.draw[[j]], startlatent= h_t[,j], priorspec= sv_priors[[j]])
         sv.idio.draw[[j]][c("mu", "phi", "sigma")] <- as.list(sv.idio.draw.i$para[, c("mu", "phi", "sigma")])
         h_t[,j] <- sv.idio.draw.i$latent   # sigma2.idio.draw
         # Threshold
         h_t[exp(h_t[,j]) < 1e-6,j] <- log(1e-6)
         h_t[exp(h_t[,j]) > 10,j]   <- log(10)
         # hyperparameters
         muh[j]  = as.numeric(sv.idio.draw[[j]][c("mu")])
         phih[j] = as.numeric(sv.idio.draw[[j]][c("phi")])
         sigh[j] = as.numeric(sv.idio.draw[[j]][c("sigma")])
      }
      
      ##### DELTA #####
      delta.old = delta
      delta = metropolis_delta_sv(delta.old, beta, rho, h_t,
                                  nu, X_t, y,
                                  min_network, plus_network, print_out=FALSE)
      if(!setequal(delta,delta.old)){
         accepts = 1
      }
      
      ##### RHO ##### 
      # Network A will be constructed within 'sample_rho()', so no  need to update A already here
      rho.old = rho
      rho = sample_rho_sv(rho.old, beta, h_t, delta, y, X_t,
                          min_network, plus_network,
                          rho_params, prior_rho = prior_rho)
      #paste0("rho = ",rho)
      
      ##### UPDATE MATRIX A #####
      # Since, delta and rho are connected to A, now update A with newest values
      for(t in 1:T){
         A[[t]] = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network[[t]]
                                                        + delta[2]*plus_network[[t]]))
         y_star[t,] = A[[t]] %*% unlist(y[t,])
      }
      
      ## Store posterior quantities of interest 
      if (irep > nburn){
         beta.store[irep-nburn, ]  = beta
         ht.store[irep-nburn, , ]  = h_t
         delta.store[irep-nburn, ] = delta
         rho.store[irep-nburn, ]   = rho
         muh.store[irep-nburn, ]   = muh
         phih.store[irep-nburn, ]  = phih
         sigh.store[irep-nburn, ]  = sigh
         accepts.store[irep-nburn] = accepts
         if (forec.horiz > 0){
            h_t.pred = matrix(0, num_countries,1)
            for (j in 1:num_countries){
               muh  = as.numeric(sv.idio.draw[[j]]["mu"])
               phih = as.numeric(sv.idio.draw[[j]]["phi"])
               sigh = as.numeric(sv.idio.draw[[j]]["sigma"])
               h_t.pred[j] = muh + phih*(h_t[j]-muh) + sigh*rnorm(1)
            }
            Anew = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network.pred[[1]] + delta[2]*plus_network.pred[[1]]))
            Xtnew = cbind(matrix(1,num_countries,1), kronecker(t(ft.pred[1,]), diag(num_countries)))
            y.forec.store[irep-nburn, ,1] = solve(Anew,  Xtnew%*%beta + 
                                                         diag(as.vector(exp(h_t.pred/2)))%*%matrix(rnorm(num_countries),num_countries,1))
         }
         if (forec.horiz > 1){
            for (h in 2:forec.horiz){
               for (j in 1:num_countries){
                  muh  = as.numeric(sv.idio.draw[[j]]["mu"])
                  phih = as.numeric(sv.idio.draw[[j]]["phi"])
                  sigh = as.numeric(sv.idio.draw[[j]]["sigma"])
                  h_t.pred[j] = muh + phih*(h_t.pred[j]-muh) + sigh*rnorm(1)
               }
               Anew = diag(num_countries) - (diag(rho) %*% (delta[1]*min_network.pred[[h]] + delta[2]*plus_network.pred[[h]]))
               Xtnew = cbind(matrix(1,num_countries,1), kronecker(t(ft.pred[h,]), diag(num_countries)))
               y.forec.store[irep-nburn, ,h] = solve(Anew,  Xtnew%*%beta +
                                                            diag(as.vector(exp(h_t.pred/2)))%*%matrix(rnorm(num_countries),num_countries,1))
            }
         }
      }
      
      #%%% Display iteration and time %%%%
      if (irep %% nshow == 0){
         time.taken = round(as.numeric(Sys.time() - start.time, units="secs"),2)
         print( paste0("Iter ",irep," of ",niter,". Computing time = ",time.taken," secs") )
      }
   }
   acceptance_rate = mean(accepts.store)
   
   
   results = list(beta.store    = beta.store,
                  delta.store   = delta.store,
                  rho.store     = rho.store,
                  ht.store      = ht.store,
                  accepts.store = accepts.store,
                  y.forec.store = y.forec.store,
                  muh.store     = muh.store,
                  phih.store    = phih.store,
                  sigh.store    = sigh.store,
                  A = A,
                  X_t = X_t)
   
   return(results)
}
