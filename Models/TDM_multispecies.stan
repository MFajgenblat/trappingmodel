functions{
  
  vector GP_1D(array[] real x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }
  
  vector GP_nD(array[] vector x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }
  
  vector GP_periodic(array[] real x, real rho, real alpha, vector eta, real period) {
    return(cholesky_decompose(add_diag(gp_periodic_cov(x, alpha, rho, period), 1e-9)) * eta);
  }
  
  real partial_sum_negbin_Taylor1(array[,] int Y_sliced,
    int start, int end, int N_species, array[] int year, array[] int site,
    row_vector beta_0, vector log_m, array[] int midpoint, matrix bf_d0, matrix phenology_weights, matrix X, matrix Beta, matrix trend, matrix species_loadings, matrix site_loadings,
    vector phi) {
    real ll = 0;
    for (sample in start:end) {
      int i = sample - start + 1;
      ll += neg_binomial_2_log_lpmf(Y_sliced[i,] | beta_0 + rep_row_vector(log_m[sample], N_species) + bf_d0[midpoint[sample],]*phenology_weights + X[sample,] * Beta + trend[year[sample],] + site_loadings[site[sample],]*species_loadings, phi);
    }
    return(ll);
  }
  
}

data {
  
  int<lower=0> N_samples;
  int<lower=0> N_species;
  int<lower=0> N_covariates;
  int<lower=0> N_sites;
  int<lower=0> N_spatial_bf;
  int<lower=0> N_days;
  int<lower=0> N_years;
  int<lower=0> N_bf;
  int<lower=0> N_dims;
  
  array[N_samples,N_species] int<lower=0> Y;
  
  matrix[N_samples,N_covariates] X;
  
  array[N_samples] int<lower=1,upper=N_sites> site;
  array[N_samples] int<lower=1,upper=N_years> year;
  
  array[N_years] real<lower=-1,upper=1> year_range;
  
  matrix[N_sites,N_spatial_bf] spatial_bf;
  array[N_spatial_bf] vector<lower=-1,upper=1>[2] spatial_bf_range;
  
  matrix[N_days,N_bf] bf_d0;
  vector[N_samples] log_m;
  array[N_samples] int<lower=1,upper=N_days> midpoint;
  array[N_bf] real bf_range;
  
}

transformed data {
  
  int grainsize = 1;
  
}

parameters {
  
  real beta_0_0;
  row_vector[N_species] beta_0_std;
  real<lower=0> beta_0_sd;
  
  vector[N_covariates] Beta_0;
  matrix[N_covariates,N_species] Beta_std;
  vector<lower=0>[N_covariates] Beta_sd;
  
  array[N_species] real<lower=0> trend_species_ls;
  array[N_species] real<lower=0> trend_species_sd;
  array[N_species] vector[N_years] trend_species_std;
  
  real<lower=0> phenology_weights_overall_sd;
  real<lower=0> phenology_weights_overall_ls;
  vector[N_bf] phenology_weights_overall_std;
  
  array[N_species] real<lower=0> phenology_weights_species_sd;
  array[N_species] real<lower=0> phenology_weights_species_ls;
  array[N_species] vector[N_bf] phenology_weights_species_std;
  
  matrix[N_dims,N_species] species_loadings_std;
  vector<lower=0>[N_dims] species_loadings_delta;
  row_vector<lower=0,upper=1>[N_dims] rho_space;
  matrix[N_sites,N_dims] site_loadings_iid;
  array[N_dims] real<lower=0> site_loadings_spatial_ls;
  matrix[N_spatial_bf,N_dims] site_loadings_spatial_std;
  
  vector<lower=0>[N_species] phi;
  
}

transformed parameters {
  
  row_vector[N_species] beta_0 = beta_0_0 + beta_0_std * beta_0_sd;
  matrix[N_covariates,N_species] Beta = rep_matrix(Beta_0, N_species) + Beta_std .* rep_matrix(Beta_sd, N_species);
  matrix[N_years,N_species] trend;
  vector[N_bf] phenology_weights_overall;
  matrix[N_bf,N_species] phenology_weights;
  matrix[N_dims,N_species] species_loadings;
  matrix[N_sites,N_dims] site_loadings;
  
  for (i in 1:N_species) {trend[,i] = GP_1D(year_range, trend_species_ls[i], trend_species_sd[i], trend_species_std[i]);}
  
  phenology_weights_overall = GP_periodic(bf_range, phenology_weights_overall_ls, phenology_weights_overall_sd, phenology_weights_overall_std, 2);
  for (i in 1:N_species) {phenology_weights[,i] = phenology_weights_overall + GP_periodic(bf_range, phenology_weights_species_ls[i], phenology_weights_species_sd[i], phenology_weights_species_std[i], 2);}

  species_loadings = species_loadings_std ./ rep_matrix(sqrt(exp(cumulative_sum(log(species_loadings_delta)))), N_species);
  
  for (i in 1:N_dims) {site_loadings[,i] = (1 - rho_space[i]) * site_loadings_iid[,i] + rho_space[i] * (spatial_bf * GP_nD(spatial_bf_range, site_loadings_spatial_ls[i], 1.0, site_loadings_spatial_std[,i]));}

}

model {
  
  target += normal_lpdf(beta_0_0 | 0, 3);
  target += std_normal_lpdf(beta_0_std);
  target += normal_lpdf(beta_0_sd | 0, 3);
  
  target += normal_lpdf(Beta_0 | 0, 3);
  target += std_normal_lpdf(to_vector(Beta_std));
  target += normal_lpdf(Beta_sd | 0, 3);
    
  target += normal_lpdf(phenology_weights_overall_sd | 0, 3);
  target += inv_gamma_lpdf(phenology_weights_overall_ls | 5, 5);
  target += std_normal_lpdf(phenology_weights_overall_std);
  
  for (i in 1:N_species) {
    
    target += normal_lpdf(phi[i] | 0, 3);
    
    target += inv_gamma_lpdf(trend_species_ls[i] | 5, 5);
    target += normal_lpdf(trend_species_sd[i] | 0, 3);
    target += std_normal_lpdf(trend_species_std[i]);
  
    target += normal_lpdf(phenology_weights_species_sd[i] | 0, 3);
    target += inv_gamma_lpdf(phenology_weights_species_ls[i] | 5, 5);
    target += std_normal_lpdf(phenology_weights_species_std[i]);
  }
    
  for (j in 1:N_dims) {
    if (j == 1) {target += gamma_lpdf(species_loadings_delta[j] | 2, 1);}
    if (j > 1) {target += gamma_lpdf(species_loadings_delta[j] | 6, 1);}
    target += std_normal_lpdf(species_loadings_std[j,]);
    target += uniform_lpdf(rho_space[j] | 0, 1);
    target += std_normal_lpdf(site_loadings_iid[,j]);
    target += std_normal_lpdf(site_loadings_spatial_std[,j]);
    target += inv_gamma_lpdf(site_loadings_spatial_ls[j] | 5, 5);
  }
  
  target += reduce_sum(partial_sum_negbin_Taylor1,
      Y, grainsize, N_species, year, site,
      beta_0, log_m, midpoint, bf_d0, phenology_weights, X, Beta, trend, species_loadings, site_loadings,
      phi);

}
