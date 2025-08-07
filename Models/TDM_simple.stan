functions {
  
  vector GP_periodic(array[] real X, real rho, real alpha, vector eta, real period) {
    return(cholesky_decompose(add_diag(gp_periodic_cov(X, alpha, rho, period), 1e-9)) * eta);
  }
  
}

data {
  
  int<lower=0> N_samples;
  int<lower=0> N_covariates;
  int<lower=0> N_days;
  int<lower=0> N_bf;
  
  array[N_samples] int<lower=0> Y;
  
  matrix[N_samples,N_covariates] X;
  
  matrix[N_days,N_bf] bf_d0;
  matrix[N_days,N_bf] bf_d1;
  matrix[N_days,N_bf] bf_d2;
  vector[N_samples] c;
  int m;
  real log_m;
  array[m,N_samples] int sampling_intervals;
  array[N_samples] int midpoint;
  array[N_bf] real bf_range;
  
  int<lower=1,upper=4> model_id;
  int<lower=1,upper=2> model_likelihood;
  
}

transformed data {
  
  int<lower=0,upper=1> negbin = model_likelihood-1;
  
}

parameters {
  
  real beta_0;
  vector[N_covariates] Beta;
  
  vector[N_bf] phenology_weights_std;
  real<lower=0> phenology_weights_sd;
  real<lower=0> phenology_weights_ls;
  
  array[negbin] real<lower=0> phi;
  
}

transformed parameters {
  
  vector[N_bf] phenology_weights = GP_periodic(bf_range, phenology_weights_ls, phenology_weights_sd, phenology_weights_std, 2);
  vector[N_days] phenology_d0 = bf_d0 * phenology_weights;
  vector[N_days] phenology_d1 = bf_d1 * phenology_weights;
  vector[N_days] phenology_d2 = bf_d2 * phenology_weights;
  
}

model {
  
  // PRIORS
  target += normal_lpdf(beta_0 | 0, 3);
  target += normal_lpdf(Beta | 0, 3);
  
  target += std_normal_lpdf(phenology_weights_std);
  target += normal_lpdf(phenology_weights_sd | 0, 3);
  target += inv_gamma_lpdf(phenology_weights_ls | 5, 5);
  
  if (model_likelihood == 2) {target += normal_lpdf(phi[1] | 0, 3);}
  
  // LIKELIHOOD
  if (model_likelihood == 1) {
    if (model_id == 1) {
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + log_m);
    }
    if (model_id == 2) {
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + log_m + phenology_d0[midpoint]);
    }
    if (model_id == 3) {
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + log_m + phenology_d0[midpoint] + log1p(fmax(-0.99, c .* (phenology_d1[midpoint]^2 + phenology_d2[midpoint]))));
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = log_sum_exp(phenology_d0[sampling_intervals[,sample]]);}
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + phenology);
    }
  } else if (model_likelihood == 2) {
    if (model_id == 1) {
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + log_m, phi[1]);
    }
    if (model_id == 2) {
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + log_m + phenology_d0[midpoint], phi[1]);
    }
    if (model_id == 3) {
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + log_m + phenology_d0[midpoint] + log1p(fmax(-0.99, c .* (phenology_d1[midpoint]^2 + phenology_d2[midpoint]))), phi[1]);
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = log_sum_exp(phenology_d0[sampling_intervals[,sample]]);}
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + phenology, phi[1]);
    }
  }
}
