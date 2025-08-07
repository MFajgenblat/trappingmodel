################################################################################
# TEMPORAL DISAGGREGATION MODEL
# Script A - Simulation study
################################################################################

library(mgcv)
library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(ggtext)
library(ggh4x)

#-------------------------------------------------------------------------------
# Set working directory
#-------------------------------------------------------------------------------

# Set working directory to the folder where you cloned/downloaded the repository
# Example (update this to your own path):

# setwd("~/path/to/repository-folder")

#-------------------------------------------------------------------------------
# Specify simulation settings
#-------------------------------------------------------------------------------

# Parameters
N_samples <- 1000
N_covariates <- 1
N_days <- 365
beta_0 <- -4
Beta = c(-1)
gamma_sd <- 5
N_bf <- 30

#-------------------------------------------------------------------------------
# Loop through all simulation runs
#-------------------------------------------------------------------------------

for (run_id in 1:50) {
  
  for (N_bf_datagen in c(8, 16)) {
    
    for (N_duration in c(15, 30)) {
      
      # Basis function design matrix (incl. derivatives)
      bf_d0 <- cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf_datagen + 1), ord = 4, derivs = 0)
      bf_d1 <- cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf_datagen + 1), ord = 4, derivs = 1)
      bf_d2 <- cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf_datagen + 1), ord = 4, derivs = 2)
      
      # Randomly generate zero-centered normally distributed loadings for basis functions
      gamma <- rnorm(N_bf_datagen, mean = 0, sd = gamma_sd)
      gamma <- gamma - mean(gamma)
      
      # Compute the spline over the entire year
      f <- bf_d0 %*% gamma
      plot(f)
      
      # Generate random 30-day sampling_intervals
      start_days <- sample(1:(N_days - N_duration + 1), N_samples, replace = TRUE)
      sampling_intervals <- sapply(start_days, function(s) seq(s, s + N_duration - 1) %% N_days + 1)
      
      # Randomly generate predictor values for the interannual linear trend
      X <- replicate(N_covariates, runif(N_samples, 0, 1))
      
      # Evaluate spline, exponentiate, and sum for each interval
      Y <- sapply(1:N_samples, function(i) {sum(rpois(N_duration, exp(beta_0 + c(X[i,] %*% Beta) + f[sampling_intervals[,i]])))})
      hist(Y)
      
      # Compute the necessary metadata for the Taylor series approximations
      midpoint <- round(apply(sampling_intervals, 2, mean))
      c <- 0.5*sapply(1:N_samples, function(i) sum((sampling_intervals[,i] - midpoint[i])^2)/30)
      
      # Compile data list
      datalist <- list(N_samples = N_samples,
                       N_covariates = N_covariates,
                       N_days = N_days,
                       N_bf = N_bf,
                       X = X,
                       Y = Y,
                       bf_d0 = cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf + 1), ord = 4, derivs = 0),
                       bf_d1 = cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf + 1), ord = 4, derivs = 1),
                       bf_d2 = cSplineDes(1:N_days, seq(1, N_days, length.out = N_bf + 1), ord = 4, derivs = 2),
                       c = c,
                       m = N_duration,
                       log_m = log(N_duration),
                       sampling_intervals = sampling_intervals,
                       midpoint = midpoint,
                       bf_range = seq(-1, 1, length.out = N_bf + 1)[1:N_bf])
      
      model <- cmdstan_model("Models/TDM_singlespecies_simple.stan")
      
      fit1 <- model$sample(data = c(datalist, model_id = 1, model_likelihood = 1), chains = 4, iter_warmup = 100, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit2 <- model$sample(data = c(datalist, model_id = 2, model_likelihood = 1), chains = 4, iter_warmup = 100, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit3 <- model$sample(data = c(datalist, model_id = 3, model_likelihood = 1), chains = 4, iter_warmup = 100, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit4 <- model$sample(data = c(datalist, model_id = 4, model_likelihood = 1), chains = 4, iter_warmup = 100, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit5 <- model$sample(data = c(datalist, model_id = 1, model_likelihood = 2), chains = 4, iter_warmup = 100, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit6 <- model$sample(data = c(datalist, model_id = 2, model_likelihood = 2), chains = 4, iter_warmup = 1, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit7 <- model$sample(data = c(datalist, model_id = 3, model_likelihood = 2), chains = 4, iter_warmup = 1, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      fit8 <- model$sample(data = c(datalist, model_id = 4, model_likelihood = 2), chains = 4, iter_warmup = 1, iter_sampling = 1, refresh = 500, parallel_chains = 4)
      
      stanfit1 <- read_cmdstan_csv(fit1$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit2 <- read_cmdstan_csv(fit2$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit3 <- read_cmdstan_csv(fit3$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit4 <- read_cmdstan_csv(fit4$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit5 <- read_cmdstan_csv(fit5$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit6 <- read_cmdstan_csv(fit6$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit7 <- read_cmdstan_csv(fit7$output_files(), format = "draws_matrix")$post_warmup_draws
      stanfit8 <- read_cmdstan_csv(fit8$output_files(), format = "draws_matrix")$post_warmup_draws
      
      phenology1 <- stanfit2[,substr(colnames(stanfit1), 1, 13) == "phenology_d0["] - rowMeans(stanfit1[,substr(colnames(stanfit1), 1, 13) == "phenology_d0["])
      phenology2 <- stanfit2[,substr(colnames(stanfit2), 1, 13) == "phenology_d0["] - rowMeans(stanfit2[,substr(colnames(stanfit2), 1, 13) == "phenology_d0["])
      phenology3 <- stanfit3[,substr(colnames(stanfit3), 1, 13) == "phenology_d0["] - rowMeans(stanfit3[,substr(colnames(stanfit3), 1, 13) == "phenology_d0["])
      phenology4 <- stanfit4[,substr(colnames(stanfit4), 1, 13) == "phenology_d0["] - rowMeans(stanfit4[,substr(colnames(stanfit4), 1, 13) == "phenology_d0["])
      phenology5 <- stanfit5[,substr(colnames(stanfit5), 1, 13) == "phenology_d0["] - rowMeans(stanfit5[,substr(colnames(stanfit5), 1, 13) == "phenology_d0["])
      phenology6 <- stanfit6[,substr(colnames(stanfit6), 1, 13) == "phenology_d0["] - rowMeans(stanfit6[,substr(colnames(stanfit6), 1, 13) == "phenology_d0["])
      phenology7 <- stanfit7[,substr(colnames(stanfit7), 1, 13) == "phenology_d0["] - rowMeans(stanfit7[,substr(colnames(stanfit7), 1, 13) == "phenology_d0["])
      phenology8 <- stanfit8[,substr(colnames(stanfit8), 1, 13) == "phenology_d0["] - rowMeans(stanfit8[,substr(colnames(stanfit8), 1, 13) == "phenology_d0["])
      
      beta1 <- stanfit1[,colnames(stanfit1) == "Beta[1]"]
      beta2 <- stanfit2[,colnames(stanfit2) == "Beta[1]"]
      beta3 <- stanfit3[,colnames(stanfit3) == "Beta[1]"]
      beta4 <- stanfit4[,colnames(stanfit4) == "Beta[1]"]
      beta5 <- stanfit5[,colnames(stanfit5) == "Beta[1]"]
      beta6 <- stanfit6[,colnames(stanfit6) == "Beta[1]"]
      beta7 <- stanfit7[,colnames(stanfit7) == "Beta[1]"]
      beta8 <- stanfit8[,colnames(stanfit8) == "Beta[1]"]
      
      write.table(data.frame(run_id = run_id,
                             sampling_interval = N_duration,
                             bf = N_bf_datagen,
                             model_id = rep(c("Naive", "Taylor 1", "Taylor 2", "Exact"), 2),
                             model_likelihood = rep(c("Poisson", "Negative binomial"), each = 4),
                             phenology_coverage = c(mean(sapply(1:365, function(i) f[i] > quantile(phenology1[,i], 0.025) & f[i] < quantile(phenology1[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology2[,i], 0.025) & f[i] < quantile(phenology2[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology3[,i], 0.025) & f[i] < quantile(phenology3[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology4[,i], 0.025) & f[i] < quantile(phenology4[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology5[,i], 0.025) & f[i] < quantile(phenology5[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology6[,i], 0.025) & f[i] < quantile(phenology6[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology7[,i], 0.025) & f[i] < quantile(phenology7[,i], 0.975))),
                                                    mean(sapply(1:365, function(i) f[i] > quantile(phenology8[,i], 0.025) & f[i] < quantile(phenology8[,i], 0.975)))),
                             phenology_mae = c(mean(sapply(1:365, function(i) abs(f[i] - mean(phenology1[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology2[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology3[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology4[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology5[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology6[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology7[,i])))),
                                               mean(sapply(1:365, function(i) abs(f[i] - mean(phenology8[,i]))))),
                             trend_mae = c(abs(Beta[1] - mean(beta1)),
                                           abs(Beta[1] - mean(beta2)),
                                           abs(Beta[1] - mean(beta3)),
                                           abs(Beta[1] - mean(beta4)),
                                           abs(Beta[1] - mean(beta5)),
                                           abs(Beta[1] - mean(beta6)),
                                           abs(Beta[1] - mean(beta7)),
                                           abs(Beta[1] - mean(beta8))),
                             trend_width = c(quantile(beta1, 0.975) - quantile(beta1, 0.025),
                                             quantile(beta2, 0.975) - quantile(beta2, 0.025),
                                             quantile(beta3, 0.975) - quantile(beta3, 0.025),
                                             quantile(beta4, 0.975) - quantile(beta4, 0.025),
                                             quantile(beta5, 0.975) - quantile(beta5, 0.025),
                                             quantile(beta6, 0.975) - quantile(beta6, 0.025),
                                             quantile(beta7, 0.975) - quantile(beta7, 0.025),
                                             quantile(beta8, 0.975) - quantile(beta8, 0.025)),
                             runtime = c(fit1$time()$total,
                                         fit2$time()$total,
                                         fit3$time()$total,
                                         fit4$time()$total,
                                         fit5$time()$total,
                                         fit6$time()$total,
                                         fit7$time()$total,
                                         fit8$time()$total)),
                  paste0("Simulation/simulation_output_", N_duration, "_", N_bf_datagen, "_", run_id, ".csv"), row.names = F, col.names = T, sep = ",")
      
    }
    
  }
  
}

#-------------------------------------------------------------------------------
# Visualize simulation results
#-------------------------------------------------------------------------------

simulation_data <- do.call(rbind, lapply(list.files("Simulation", "^simulation_output_.*\\.csv$", full.names = T), read.csv))
simulation_data %>%
  pivot_longer(!c(run_id, sampling_interval, bf, model_id, model_likelihood)) %>%
  mutate(sampling_interval = case_when(sampling_interval == 15 ~ "<b>Short sampling interval (15 days)</b>",
                                       sampling_interval == 30 ~ "<b>Short sampling interval (30 days)</b>"),
         bf = case_when(bf == 8 ~ "<b>Slowly changing<br>phenology (8 bf)</b>",
                        bf == 16 ~ "<b>Rapidly changing<br>phenology (16 bf)</b>"),
         name = factor(case_when(name == "phenology_coverage" ~ "<b>Phenology<br>CrI coverage</b><br>Higher is better",
                                 name == "phenology_mae" ~ "<b>Phenology mean<br>absolute error</b><br>Lower is better",
                                 name == "runtime" ~ "<b>Run time</b><br>Lower is better",
                                 name == "trend_mae" ~ "<b>Trend mean<br>absolute error</b><br><i>log10-transformed</i><br>Lower is better",
                                 name == "trend_width" ~ "<b>Trend<br>CrI width</b><br><i>log10-transformed</i>"),
                       levels = c("<b>Phenology<br>CrI coverage</b><br>Higher is better", "<b>Phenology mean<br>absolute error</b><br>Lower is better", "<b>Trend mean<br>absolute error</b><br><i>log10-transformed</i><br>Lower is better", "<b>Trend<br>CrI width</b><br><i>log10-transformed</i>", "<b>Run time</b><br>Lower is better")),
         model_likelihood = factor(model_likelihood,
                                   levels = c("Poisson", "Negative binomial"))) %>%
  ggplot(aes(x = model_id, y = value, color = model_likelihood, group = run_id), position = position_jitter(width = 0.2, height = 0)) +
  geom_rect(data = data.frame(xmin = c(1.5, 3.5), xmax = c(2.5, 4.5)), aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = +Inf, fill = "black", alpha = 0.05, inherit.aes = F) +
  geom_point(shape = 16, size = 0.05) +
  geom_path(linewidth = 0.1) +
  stat_summary(aes(x = model_id, y = value, color = model_likelihood),
               fun.data = mean_cl_normal, fun.args = list(mult = 1),
               position = position_jitter(width = 0.2, height = 0), inherit.aes = F, shape = 16, size = 0.25, linewidth = 0.5) +
  scale_color_brewer("Model likelihood", palette = "Dark2") +
  scale_x_discrete("Model") +
  facet_nested(name ~ sampling_interval + bf, scales = "free_y", nest_line = element_line(color = "black"), switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_markdown(size = 7),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 6),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(4, "mm"),
        legend.position = "bottom")
ggsave("Plots/Simulation_study.png", width = 16, height = 18, units = "cm", dpi = 600)
