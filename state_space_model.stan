//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
data {
  int<lower=1> T;                       // number of time steps
  int<lower=1> K;                       // number of observed variables per time point
  matrix[T, K] y_obs;                // observed values
}

parameters {
  matrix[K, T] X;                // latent states
  matrix[K, K] theta;            // system matrix
  matrix<lower=0>[K] A;            // intercept term since model is not centered.  Should be positive since all observations are positive.
  vector<lower=0>[K] sigma_proc;               // standard deviations
  cholesky_factor_corr[K] L_corr;              // Cholesky factor of correlation matrix
  real<lower=0> sigma_obs;            // observation noise std dev
}

transformed parameters {
  matrix[K, K] L_proc = diag_pre_multiply(sigma_proc, L_corr); // Full Cholesky of covariance
}

model {
  // Priors for noise parameters
  sigma_proc ~ normal(0, .25);         // Half-normal(0,sd) has mean sd * sqrt(2/pi).  with sd = .25 the mean is 0.199.  Note this is for the SD of the STATE process
  L_corr     ~ lkj_corr_cholesky(2);   // Mildly informative prior that shrinks towards the identity matrix
  sigma_obs ~ normal(0, .008);         // Half-normal(0,sd) has mean sd * sqrt(2/pi).  with sd = .008 the mean is 0.0064.  Note this is for the SD of the OBS process
  A ~ normal(0,.25)                    // originally N(0,1), but the observed values are all near zero and all positive, so now it's half normal prior.
  // prior for the model parameters
  to_vector(theta) ~ normal(0, .5);     // Weakly informative prior on each entry, originally N(0,1)
  
  // initial state
  X[,1] ~ multi_normal(A, diag_matrix(rep_vector(1, K))); 
  
  // State evolution
  for (t in 2:T)
    X[, t] ~ multi_normal_cholesky(A+theta * X[, t-1], L_proc);
  
  // Observation model.  All have shared noise process.
  for (t in 1:T) {
    int obs_idx[K];
    int count = 0;
    for (k in 1:K) {
      if (!is_nan(Y[t][k])) {
        count += 1;
        obs_idx[count] = k;
      }
    }
    if (count > 0) {
      vector[count] Y_obs;
      vector[count] X_obs;
      for (i in 1:count) {
        Y_obs[i] = Y[t][obs_idx[i]];
        X_obs[i] = X[t][obs_idx[i]];
      }
      Y_obs ~ multi_normal(X_obs, diag_matrix(rep_vector(sigma_obs[1,1], count)));
    }
  }
}
