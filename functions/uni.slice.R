# Slice Sampler

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

uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL){
  # Here, I added that all necessary variables are also parameterized
  # Check the validity of the arguments.
  if (!is.numeric(x0) || length(x0)!=1 #<-2nd can be removed when working with multiple rho's 
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