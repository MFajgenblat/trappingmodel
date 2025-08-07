################################################################################
# TEMPORAL DISAGGREGATION MODEL
# Script C - Case study: inferring varying phenology
################################################################################

#-------------------------------------------------------------------------------
# Set working directory
#-------------------------------------------------------------------------------

# Set working directory to the folder where you cloned/downloaded the repository
# Example (update this to your own path):

# setwd("~/path/to/repository-folder")

#-------------------------------------------------------------------------------
# Load packages
#-------------------------------------------------------------------------------

library(cmdstanr)
library(reshape2)
library(tidyverse)
library(patchwork)
library(vroom)
library(tidybayes)
library(ggtext)

#-------------------------------------------------------------------------------
# Compiling the Stan model
#-------------------------------------------------------------------------------

set_cmdstan_path("C:/Users/maxim/Documents/.cmdstan/cmdstan-2.36.0")
model <- cmdstan_model("Models/TDM_yearlyphen.stan", cpp_options = list(stan_threads = TRUE))

for (i in c(1:50)) {
  
  #-------------------------------------------------------------------------------
  # Reading data object
  #-------------------------------------------------------------------------------
  
  species_id <- i
  datalist <- readRDS("Data/datalist_2025_07_17.rds")
  datalist$Y <- datalist$Y[,species_id]
  datalist$model_likelihood = 2
  
  #-------------------------------------------------------------------------------
  # Fit the model
  #-------------------------------------------------------------------------------
  
  datalist$model_id <- 1
  pathfit <- model$pathfinder(data = datalist, refresh = 500, num_threads = 4, num_paths = 4, draws = 4, max_lbfgs_iters = 2000, psis_resample = T)
  fit1 <- model$sample(data = datalist,
                       iter_warmup = 500, iter_sampling = 500, refresh = 50,
                       chains = 4, parallel_chains = 4, threads_per_chain = 1,
                       init = pathfit$draws(format = "list"))
  
  datalist$model_id <- 2
  pathfit <- model$pathfinder(data = datalist, refresh = 500, num_threads = 4, num_paths = 4, draws = 4, max_lbfgs_iters = 2000, psis_resample = T)
  fit2 <- model$sample(data = datalist,
                       iter_warmup = 500, iter_sampling = 500, refresh = 50,
                       chains = 4, parallel_chains = 4, threads_per_chain = 1,
                       init = pathfit$draws(format = "list"))
  
  datalist$model_id <- 3
  pathfit <- model$pathfinder(data = datalist, refresh = 500, num_threads = 4, num_paths = 4, draws = 4, max_lbfgs_iters = 2000, psis_resample = T)
  fit3 <- model$sample(data = datalist,
                       iter_warmup = 500, iter_sampling = 500, refresh = 50,
                       chains = 4, parallel_chains = 4, threads_per_chain = 1,
                       init = pathfit$draws(format = "list"))
  
  datalist$model_id <- 4
  pathfit <- model$pathfinder(data = datalist, refresh = 500, num_threads = 4, num_paths = 4, draws = 4, max_lbfgs_iters = 2000, psis_resample = T)
  fit4 <- model$sample(data = datalist,
                       iter_warmup = 500, iter_sampling = 500, refresh = 50,
                       chains = 4, parallel_chains = 4, threads_per_chain = 1,
                       init = pathfit$draws(format = "list"))
  
  #-------------------------------------------------------------------------------
  # Extract the posterior samples
  #-------------------------------------------------------------------------------
  
  stanfit1 <- read_cmdstan_csv(fit1$output_files(), format = "draws_matrix")$post_warmup_draws
  stanfit2 <- read_cmdstan_csv(fit2$output_files(), format = "draws_matrix")$post_warmup_draws
  stanfit3 <- read_cmdstan_csv(fit3$output_files(), format = "draws_matrix")$post_warmup_draws
  stanfit4 <- read_cmdstan_csv(fit4$output_files(), format = "draws_matrix")$post_warmup_draws
  
  #-------------------------------------------------------------------------------
  # Save the loo-values
  #-------------------------------------------------------------------------------
  
  loo1 <- fit1$loo(variables = "ll")
  loo2 <- fit2$loo(variables = "ll")
  loo3 <- fit3$loo(variables = "ll")
  loo4 <- fit4$loo(variables = "ll")
  
  loo <- rbind(data.frame(species_id = species_id, model_id = 1, metric = c("elpd_loo", "p_loo", "looic"), loo1$estimates),
               data.frame(species_id = species_id, model_id = 2, metric = c("elpd_loo", "p_loo", "looic"), loo2$estimates),
               data.frame(species_id = species_id, model_id = 3, metric = c("elpd_loo", "p_loo", "looic"), loo3$estimates),
               data.frame(species_id = species_id, model_id = 4, metric = c("elpd_loo", "p_loo", "looic"), loo4$estimates))
  
  write.table(loo, paste0("Output/loo_", species_id, ".csv"), row.names = F, col.names = T, sep = ",")
  
  #-------------------------------------------------------------------------------
  # Save the posterior samples on phenology
  #-------------------------------------------------------------------------------
  
  write.table(stanfit2[,grep("^yearly_phenology\\[", colnames(stanfit2), value = TRUE)], paste0("Output/yearly_phenology_2_", species_id, ".csv"), row.names = F, col.names = T, sep = ",")
  write.table(stanfit3[,grep("^yearly_phenology\\[", colnames(stanfit3), value = TRUE)], paste0("Output/yearly_phenology_3_", species_id, ".csv"), row.names = F, col.names = T, sep = ",")
  write.table(stanfit4[,grep("^yearly_phenology\\[", colnames(stanfit4), value = TRUE)], paste0("Output/yearly_phenology_4_", species_id, ".csv"), row.names = F, col.names = T, sep = ",")
  
}

#-------------------------------------------------------------------------------
# Visualize cross-validation and phenological inferential quality results
#-------------------------------------------------------------------------------

# Compute Wasserstein distances
wasserstein_1d <- function(x, y) {
  # Sort the samples
  x_sorted <- sort(x)
  y_sorted <- sort(y)
  
  # Match sample sizes by interpolation if needed
  n <- max(length(x_sorted), length(y_sorted))
  x_ecdf <- quantile(x_sorted, probs = seq(0, 1, length.out = n), type = 1)
  y_ecdf <- quantile(y_sorted, probs = seq(0, 1, length.out = n), type = 1)
  
  # Compute the 1-Wasserstein distance (L1 norm of difference)
  mean(abs(x_ecdf - y_ecdf))
}

WD <- do.call(rbind, lapply(list.files("Outpur", "slopes", full.names = T), function(i) data.frame(species_id = as.numeric(sub(".*_(\\d+)\\..*", "\\1", i)), M2 = wasserstein_1d(subset(read.csv(i), model == 2)$slope, subset(read.csv(i), model == 4)$slope), M3 = wasserstein_1d(subset(read.csv(i), model == 3)$slope, subset(read.csv(i), model == 4)$slope))))

# Visualize
do.call(rbind, lapply(list.files("Output", "loo_", full.names = T), read.csv)) %>%
  filter(metric == "looic") %>%
  ggplot(aes(color = factor(model_id), x = Estimate, y = gsub("\\.", " ", gsub("Sp\\.\\.", "", colnames(datalist$Y)[species_id])))) +
  geom_rect(data = data.frame(ymin = seq(1.5, 49.5, by = 2), ymax = seq(2.5, 50.5, by = 2)), aes(ymin = ymin, ymax = ymax), xmin = -Inf, xmax = +Inf, fill = "black", alpha = 0.1, inherit.aes = F) +
  geom_errorbar(aes(xmin = Estimate - SE, xmax = Estimate + SE, y = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[species_id]))), width = 0, position = position_dodge(width = 0.8), linewidth = 0.2) +
  geom_point(position = position_dodge(width = 0.8), size = 0.25, shape = 16) +
  scale_x_continuous("<b>Leave-one-out cross-validation information criterion</b><br>(LOOIC; lower is better)") +
  scale_y_discrete("Species", limits = c(gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[c(52,50:15,13:1)])))) +
  scale_color_brewer("Model", palette = "Dark2", limits = c("1", "2", "3", "4"), labels = c("Naïve", "Taylor 1", "Taylor 2", "Exact")) +
  ggtitle("(a) Predictive performance") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey93"),
        panel.grid.minor.x = element_line(color = "grey93"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 7),
        axis.title.x = element_markdown(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(face = "italic", size = 6),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        legend.position = "bottom") +
  WD %>%
  mutate(species_id = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[species_id]))) %>%
  ggplot() +
  geom_rect(data = data.frame(ymin = seq(1.5, 49.5, by = 2), ymax = seq(2.5, 50.5, by = 2)), aes(ymin = ymin, ymax = ymax), xmin = -Inf, xmax = +Inf, fill = "black", alpha = 0.1, inherit.aes = F) +
  geom_point(aes(y = species_id, x = M2, color = "2"), position = position_nudge(y = -0.2), size = 1, shape = 16) +
  geom_point(aes(y = species_id, x = M3, color = "3"), position = position_nudge(y = 0.2), size = 1, shape = 16) +
  geom_segment(aes(y = species_id, yend = species_id, x = M2, xend = 0, color = "2"), position = position_nudge(y = -0.2), linewidth = 0.3) +
  geom_segment(aes(y = species_id, yend = species_id, x = M3, xend = 0, color = "3"), position = position_nudge(y = 0.2), linewidth = 0.3) +
  geom_vline(aes(xintercept = mean(M2), color = "2"), linetype = "dashed") +
  geom_vline(aes(xintercept = mean(M3), color = "3"), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(clip = "off") +
  scale_y_discrete("Species", limits = c(gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[c(52,50:15,13:1)]))), guide = "none") +
  scale_x_continuous("<b>1-Wasserstein distance (days)</b><br>(lower is better)", expand = c(0,0), limits = c(0,9)) +
  scale_color_brewer("Model", palette = "Dark2", limits = c("1", "2", "3", "4"), labels = c("Naïve", "Taylor 1", "Taylor 2", "Exact"), guide = "none") +
  ggtitle("(b) Ability to infer phenological change") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey93"),
        panel.grid.minor.x = element_line(color = "grey93"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 7),
        axis.title.x = element_markdown(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(face = "italic", size = 6),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        legend.position = "bottom") + plot_layout(guides = "collect") & theme(legend.position = "bottom",
                                                                              legend.margin = margin(0,0,0,0),
                                                                              legend.box.margin = margin(-7.5,0,0,0),
                                                                              plot.margin = margin(5,5,0,5))
ggsave("Plots/Performance_output.png", width = 16, height = 16, dpi = 600, units = "cm")

# Formal comparison
summary(WD[,c(2,3)])
wilcox.test(WD$M2, WD$M3, paired = TRUE)

#-------------------------------------------------------------------------------
# Visualize yearly phenology and phenological change
#-------------------------------------------------------------------------------

year_range <- seq(-0.5, 0.5, length.out = 49)

# Create empty lists to compile all species-specific estimates
posteriors_list <- list()
peaks_list <- list()
intercepts_list <- list()
slopes_list <- list()
interceptslopes_list <- list()

# Loop through species to compute posterior mean activity functions, peak day intercepts and slopes
for (i in 1:4) {
  
  posterior <- vroom(paste0("Output/yearly_phenology_4_", i, ".csv")) %>%
    melt() %>%
    mutate(day = as.numeric(sub(".*yearly_phenology\\[(\\d+),.*", "\\1", variable)),
           year = as.numeric(sub(".*yearly_phenology\\[\\d+,([0-9]+)\\]", "\\1", variable)))
  
  peaks <- posterior %>%
    group_by(day, year) %>%
    mutate(.draw = row_number()) %>%
    group_by(.draw, year) %>%
    filter(value == max(value)) %>%
    filter(row_number() == 1) %>%
    data.frame(Species = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[i])))
  
  intercepts <- peaks %>% 
    select(day, year, .draw) %>%
    pivot_wider(names_from = year, values_from = day, names_sort = T) %>%
    ungroup() %>%
    select(-.draw) %>%
    apply(., 1, function(x) summary(lm(x ~ year_range))$coef[1]) %>%
    data.frame(model = 4, intercept = .) %>%
    data.frame(Species = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[i])))
  
  slopes <- peaks %>% 
    select(day, year, .draw) %>%
    pivot_wider(names_from = year, values_from = day, names_sort = T) %>%
    ungroup() %>%
    select(-.draw) %>%
    apply(., 1, function(x) summary(lm(x ~ year_range))$coef[2]) %>%
    data.frame(model = 4, slope = .) %>%
    data.frame(Species = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[i])))
  
  posteriors_list <- append(posteriors_list, list(posterior %>%
                                                    group_by(day, year) %>%
                                                    summarize(value = mean(value)) %>%
                                                    ungroup() %>%
                                                    group_by(year) %>%
                                                    mutate(value = (value - min(value))/(max(value) - min(value))) %>%
                                                    data.frame(Species = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[i])))))
  peaks_list <- append(peaks_list, list(peaks))
  intercepts_list <- append(intercepts_list, list(intercepts))
  slopes_list <- append(slopes_list, list(slopes))
  interceptslopes_list <- append(interceptslopes_list, list(left_join(intercepts, slopes)))
  
}

# Visualization
slopes_list %>%
  do.call(rbind, .) %>%
  ggplot() +
  geom_density(aes(x = slope), fill = "black", alpha = 0.2, bw = 4, linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_x_continuous("Phenological change (days)") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ factor(Species, levels = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[1:4]))), nrow = 1) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.minor.x = element_line(color = "grey93", linewidth = 0.15),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 7, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold.italic", hjust = 0.5),
        legend.position = "bottom") +
  ggplot() +
  geom_tile(data = do.call(rbind, posteriors_list), aes(x = day, y = 1975 + year, fill = value, color = value)) +
  geom_polygon(data = do.call(rbind, lapply(interceptslopes_list, function(i)
    data.frame(y = c(1976, 1976, 2024, 2024),
               x = c(quantile(i$intercept, 0.025) -  0.5*quantile(i$slope, 0.975),
                     quantile(i$intercept, 0.975) - 0.5*quantile(i$slope, 0.025),
                     quantile(i$intercept, 0.975) + 0.5*quantile(i$slope, 0.975),
                     quantile(i$intercept, 0.025) + 0.5*quantile(i$slope, 0.025)),
               Species = unique(i$Species)))),
    aes(x = x, y = y), fill = "white", alpha = 0.75) +
  geom_segment(data = do.call(rbind, interceptslopes_list) %>%
                 group_by(Species) %>%
                 summarize(intercept = mean(intercept),
                           slope = mean(slope)),
               aes(x = intercept - 0.5*slope, y = 1976, xend = intercept + 0.5*slope, yend = 2024),
               linewidth = 0.2) +
  geom_segment(data = do.call(rbind, peaks_list) %>%
                 group_by(year, Species) %>%
                 summarize(L = quantile(day, 0.025),
                           U = quantile(day, 0.975),
                           M = mean(day)),
               aes(x = L, xend = U, y = 1975 + year, yend = 1975 + year), linewidth = 0.1) +
  geom_point(data = do.call(rbind, peaks_list) %>%
               group_by(year, Species) %>%
               summarize(L = quantile(day, 0.025),
                         U = quantile(day, 0.975),
                         M = mean(day)),
             aes(x = M, y = 1975 + year), size = 0.2, shape = 16) +
  scale_x_continuous("Julian date", breaks = seq(0, 366, by = 100), expand = c(0,0)) +
  scale_y_continuous("Year", breaks = seq(1975, 2050, by = 5), expand = c(0,0)) +
  scale_color_viridis_c("Posterior mean phenological intensity", 
                        guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = 0)) +
  scale_fill_viridis_c("Posterior mean phenological intensity", 
                       guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = 0)) +
  facet_wrap(~ factor(Species, levels = gsub("\\.", " ", gsub("Sp\\.\\.","", colnames(datalist$Y)[1:4]))), nrow = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-2.5,0,0,0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6),
        legend.position = "bottom") +
  plot_layout(nrow = 2, heights = c(0.2, 1))
ggsave("Plots/Varyingphenology_output.png", width = 16, height = 10, dpi = 600, units = "cm")

# Posterior mean and 95% CrI phenological change for the for most abundantly caught spider species
interceptslopes_list %>%
  do.call(rbind, .) %>%
  group_by(Species) %>%
  summarize(Mean = mean(slope),
            LCrI = quantile(slope, 0.025),
            UCrI = quantile(slope, 0.975))
