functions {
  
  vector GP_1D(array[] real x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }
  
  vector GP_nD(array[] vector x, real rho, real alpha, vector eta) {
    return(cholesky_decompose(add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-9)) * eta);
  }
  
  vector GP_periodic(array[] real X, real rho, real alpha, vector eta, real period) {
    return(cholesky_decompose(add_diag(gp_periodic_cov(X, alpha, rho, period), 1e-9)) * eta);
  }
  
}

data {
  
  int<lower=0> N_samples;
  int<lower=0> N_covariates;
  int<lower=0> N_sites;
  int<lower=0> N_spatial_bf;
  int<lower=0> N_days;
  int<lower=0> N_years;
  int<lower=0> N_bf;
  
  array[N_samples] int<lower=0> Y;
  
  matrix[N_samples,N_covariates] X;
  
  array[N_samples] int<lower=1,upper=N_sites> site;
  array[N_samples] int<lower=1,upper=N_years> year;
  
  array[N_years] real<lower=-1,upper=1> year_range;
  
  matrix[N_sites,N_spatial_bf] spatial_bf;
  array[N_spatial_bf] vector<lower=-1,upper=1>[2] spatial_bf_range;
  
  matrix[N_days,N_bf] bf_d0;
  matrix[N_days,N_bf] bf_d1;
  matrix[N_days,N_bf] bf_d2;
  vector[N_samples] c;
  vector[N_samples] log_m;
  array[N_samples] int m;
  array[N_samples] int<lower=1,upper=N_days> start_day;
  array[N_samples] int<lower=1,upper=N_days> end_day;
  array[N_samples] int<lower=1,upper=N_days> midpoint;
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
  
  real<lower=0> year_effect_ls;
  real<lower=0> year_effect_sd;
  vector[N_years] year_effect_std;
  
  vector[N_bf] phenology_weights_std;
  real<lower=0> phenology_weights_sd;
  real<lower=0> phenology_weights_ls;
  
  real<lower=0> phenology_weights_change_sd;
  vector[N_bf] phenology_weights_change_std;
  
  real<lower=0> phenology_weights_yearly_sd;
  matrix[N_bf,N_years] phenology_weights_yearly_std;
  
  real<lower=0> site_effect_sd;
  real<lower=0,upper=1> site_effect_alpha;
  vector[N_sites] site_effect_iid;
  real<lower=0> site_effect_spatial_ls;
  vector[N_spatial_bf] site_effect_spatial_std;
  
  array[negbin] real<lower=0> phi;
  
}

transformed parameters {
  
  vector[N_years] year_effect = GP_1D(year_range, year_effect_ls, year_effect_sd, year_effect_std);
  vector[N_bf] phenology_weights_overall = GP_periodic(bf_range, phenology_weights_ls, phenology_weights_sd, phenology_weights_std, 2);
  matrix[N_bf,N_years] phenology_weights;
  vector[N_sites] site_effect = site_effect_sd * ((1 - site_effect_alpha) * site_effect_iid + site_effect_alpha * (spatial_bf * GP_nD(spatial_bf_range, site_effect_spatial_ls, 1.0, site_effect_spatial_std)));

  for (i in 1:N_years) {
    phenology_weights[,i] = phenology_weights_overall + ((i-1.0)/(N_years-1.0)) * phenology_weights_change_std * phenology_weights_change_sd + phenology_weights_yearly_std[,i] * phenology_weights_yearly_sd;
  }
  
}

model {
  
  // PRIORS
  
  target += normal_lpdf(beta_0 | 0, 3);
  target += normal_lpdf(Beta | 0, 3);
    
  target += inv_gamma_lpdf(year_effect_ls | 5, 5);
  target += normal_lpdf(year_effect_sd | 0, 3);
  target += std_normal_lpdf(year_effect_std);
  
  target += normal_lpdf(site_effect_sd | 0, 3);
  target += uniform_lpdf(site_effect_alpha | 0, 1);
  target += std_normal_lpdf(site_effect_iid);
  target += inv_gamma_lpdf(site_effect_spatial_ls | 5, 5);
  target += std_normal_lpdf(site_effect_spatial_std);
    
  target += normal_lpdf(phenology_weights_sd | 0, 3);
  target += inv_gamma_lpdf(phenology_weights_ls | 5, 5);
  target += std_normal_lpdf(phenology_weights_std);
  
  target += normal_lpdf(phenology_weights_change_sd | 0, 1);
  target += std_normal_lpdf(phenology_weights_change_std);
  
  target += normal_lpdf(phenology_weights_yearly_sd | 0, 1);
  target += std_normal_lpdf(to_vector(phenology_weights_yearly_std));
  
  if (model_likelihood == 2) {target += normal_lpdf(phi[1] | 0, 3);}
  
  // LIKELIHOOD
  if (model_likelihood == 1) {
    if (model_id == 1) {
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m);
    }
    if (model_id == 2) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]];}
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m + phenology);
    }
    if (model_id == 3) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]] + log1p(fmax(-0.99, c[sample] .* ((bf_d1[midpoint[sample],] * phenology_weights[,year[sample]])^2 + bf_d2[midpoint[sample],] * phenology_weights[,year[sample]])));}
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m + phenology);
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        if (start_day[sample] <= end_day[sample]) {
          phenology[sample] = log_sum_exp(bf_d0[start_day[sample]:end_day[sample],] * phenology_weights[,year[sample]]);
        } else {
          phenology[sample] = log_sum_exp(append_row(bf_d0[start_day[sample]:N_days,] * phenology_weights[,year[sample]], bf_d0[1:end_day[sample],] * phenology_weights[,year[sample]+1]));
        }
      }
      target += poisson_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + phenology);
    }
  } else if (model_likelihood == 2) {
    if (model_id == 1) {
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m, phi[1]);
    }
    if (model_id == 2) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]];}
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m + phenology, phi[1]);
    }
    if (model_id == 3) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]] + log1p(fmax(-0.99, c[sample] .* ((bf_d1[midpoint[sample],] * phenology_weights[,year[sample]])^2 + bf_d2[midpoint[sample],] * phenology_weights[,year[sample]])));}
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + log_m + phenology, phi[1]);
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        if (start_day[sample] <= end_day[sample]) {
          phenology[sample] = log_sum_exp(bf_d0[start_day[sample]:end_day[sample],] * phenology_weights[,year[sample]]);
        } else {
          phenology[sample] = log_sum_exp(append_row(bf_d0[start_day[sample]:N_days,] * phenology_weights[,year[sample]], bf_d0[1:end_day[sample],] * phenology_weights[,year[sample]+1]));
        }
      }
      target += neg_binomial_2_log_lpmf(Y | beta_0 + X*Beta + year_effect[year] + site_effect[site] + phenology, phi[1]);
    }
  }
  
}

generated quantities {
  
  array[N_samples] real ll;
  matrix[N_days,N_years] yearly_phenology = bf_d0 * phenology_weights;
  
  if (model_likelihood == 1) {
    if (model_id == 1) {
      for (sample in 1:N_samples) {
        ll[sample] = poisson_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample]);
      }
    }
    if (model_id == 2) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]];
        ll[sample] = poisson_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample] + phenology[sample]);
      }
    }
    if (model_id == 3) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]] + log1p(fmax(-0.99, c[sample] .* ((bf_d1[midpoint[sample],] * phenology_weights[,year[sample]])^2 + bf_d2[midpoint[sample],] * phenology_weights[,year[sample]])));
        ll[sample] = poisson_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample] + phenology[sample]);
      }
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        if (start_day[sample] <= end_day[sample]) {
          phenology[sample] = log_sum_exp(bf_d0[start_day[sample]:end_day[sample],] * phenology_weights[,year[sample]]);
        } else {
          phenology[sample] = log_sum_exp(append_row(bf_d0[start_day[sample]:N_days,] * phenology_weights[,year[sample]], bf_d0[1:end_day[sample],] * phenology_weights[,year[sample]+1]));
        }
        ll[sample] = poisson_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + phenology[sample]);
      }
    }
  } else if (model_likelihood == 2) {
    if (model_id == 1) {
      for (sample in 1:N_samples) {
        ll[sample] = neg_binomial_2_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample], phi[1]);
      }
    }
    if (model_id == 2) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]];
        ll[sample] = neg_binomial_2_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample] + phenology[sample], phi[1]);
      }
    }
    if (model_id == 3) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        phenology[sample] = bf_d0[midpoint[sample],] * phenology_weights[,year[sample]] + log1p(fmax(-0.99, c[sample] .* ((bf_d1[midpoint[sample],] * phenology_weights[,year[sample]])^2 + bf_d2[midpoint[sample],] * phenology_weights[,year[sample]])));
        ll[sample] = neg_binomial_2_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + log_m[sample] + phenology[sample], phi[1]);
      }
    }
    if (model_id == 4) {
      vector[N_samples] phenology;
      for (sample in 1:N_samples) {
        if (start_day[sample] <= end_day[sample]) {
          phenology[sample] = log_sum_exp(bf_d0[start_day[sample]:end_day[sample],] * phenology_weights[,year[sample]]);
        } else {
          phenology[sample] = log_sum_exp(append_row(bf_d0[start_day[sample]:N_days,] * phenology_weights[,year[sample]], bf_d0[1:end_day[sample],] * phenology_weights[,year[sample]+1]));
        }
        ll[sample] = neg_binomial_2_log_lpmf(Y[sample] | beta_0 + X[sample,]*Beta + year_effect[year[sample]] + site_effect[site[sample]] + phenology[sample], phi[1]);
      }
    }
  }
  
}
