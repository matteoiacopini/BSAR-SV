sample_eta <- function(rho, zeta, alpha, arho, lambda, upper.bound){
   
   T = dim(rho)[2]
   # transform:   rho in (-1,1) --> (rho+1)/2 in (0,1)
   # as zeta governs the time-variation of the latter
   rho.star = (rho+1)/2
   rho.star[rho.star==0] = 0 + 10^(-10)
   rho.star[rho.star==1] = 1 - 10^(-10)
   
   zetastar   = max(zeta)
   candidates = zetastar:upper.bound
   unnormpdf  = matrix(NA, 1, length(candidates))
   
   jj = 1
   for (mm in zetastar:upper.bound){
      term1 = sum( lgamma(mm+1) + lgamma(1+alpha+mm) - lgamma(arho+mm-zeta) - lgamma(alpha+mm-zeta) )
      term2 = mm * sum( log(1-rho.star[,1:(T-1)]) + log(1-rho.star[,2:T]) )
      logpdf.prior = dpoitrunc(mm, lambda, upper.bound)
      
      unnormpdf[jj] = term1 + term2 + logpdf.prior
      jj = jj+1
   }
   # compute "normalised log-pdf"
   logpdf = unnormpdf - logSumExp(unnormpdf)
   u = runif(1)
   pos = which.max(u <= cumsum(exp(logpdf)))
   nnn = candidates[pos]
   # browser(expr = {pos==1})
   
   return(nnn)
}
