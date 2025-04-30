# BSAR-SV

This README provides instructions on how to run the Bayesian SAR model with stochastic volatility and multiple time-varying weights using an example dataset.


Reference:
-

Costola, M., Iacopini, M., & Wichers, C. (2025). Bayesian SAR model with stochastic volatility and multiple time-varying weights. Journal of Financial Econometrics (23):3, nbae035 (https://doi.org/10.1093/jjfinec/nbae035)


-------------------------------------------------
File Descriptions:
-------------------------------------------------

1. ExampleData.RData
   - A synthetic dataset 

2. main_example.R
   - Main script file to be run in order to estimate the model 

-------------------------------------------------
List of auxiliary files 
-------------------------------------------------

1. gibbs_sv.R
   - Implements the Gibbs sampling algorithm for stochastic volatility models.

2. log_target_dens_delta_sv.R
   - Defines the target density function for delta in the stochastic volatility model.

3. max_row_norm.R
   - Provides a function to calculate the maximum row norm of a matrix.

4. metropolis_delta_sv.R
   - Implements the Metropolis-Hastings algorithm for delta in the stochastic volatility model.

5. norm_rows.R
   - Includes functions to normalize rows of a matrix.

6. sample_eta.R
   - Contains the function to sample the parameter eta.

7. sample_rho_sv.R
   - Contains the function to sample the parameter rho in the stochastic volatility model.

8. spginv.R
   - Provides a function for computing the generalized inverse of a matrix.

9. uni.slice.R
   - Implements a univariate slice sampling algorithm.

-----------------------------------------------------
Reference:
-----------------------------------------------------

Costola, M., Iacopini, M., & Wichers, C. (2025). Bayesian SAR model with stochastic volatility and multiple time-varying weights. Journal of Financial Econometrics (23):3, nbae035 (https://doi.org/10.1093/jjfinec/nbae035)
