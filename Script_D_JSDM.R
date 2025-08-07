################################################################################
# TEMPORAL DISAGGREGATION MODEL
# Script D - Case study: Joint Species Distribution Modelling
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
library(tidybayes)
library(abind)
library(patchwork)
library(sf)

#-------------------------------------------------------------------------------
# Compiling the Stan model
#-------------------------------------------------------------------------------

set_cmdstan_path("C:/Users/maxim/Documents/.cmdstan/cmdstan-2.36.0")
model <- cmdstan_model("Models/TDM_multispecies.stan", cpp_options = list(stan_threads = TRUE))

datalist <- readRDS("Data/datalist_2025_07_17.rds")
datalist$N_dims <- 10

pathfit <- model$pathfinder(data = datalist, refresh = 50, num_threads = 80, num_paths = 8, draws = 1, max_lbfgs_iters = 2000, psis_resample = T)
fit1 <- model$sample(data = datalist,
                     iter_warmup = 500, iter_sampling = 500, refresh = 10,
                     chains = 8, parallel_chains = 8, threads_per_chain = 10,
                     init = pathfit$draws(format = "list"))
stanfit1 <- read_cmdstan_csv(fit1$output_files(), format = "draws_matrix")$post_warmup_draws

saveRDS(stanfit1, "Output/stanfit.rds")

#-------------------------------------------------------------------------------
# Phenology
#-------------------------------------------------------------------------------

phenology_matrix <- abind(lapply(1:nrow(stanfit), function(i) datalist$bf_d0 %*% matrix(stanfit[i, grep("^phenology_weights\\[", colnames(stanfit))], ncol = datalist$N_species, nrow = datalist$N_bf, byrow = F)), along = 3) %>%
  melt(varnames = c("Day", "Species", ".draw")) %>%
  group_by(Species, Day) %>%
  summarize(value = mean(value)) %>%
  group_by(Species) %>%
  mutate(value = (value - min(value))/(max(value) - min(value)))

phenology_order <- phenology_matrix %>%
  group_by(Species) %>%
  filter(value == max(value)) %>%
  slice_min(Day) %>%
  ungroup() %>%
  arrange(Day) %>%
  pull(Species)

A <- phenology_matrix %>%
  mutate(Species = factor(Species, levels = rev(phenology_order))) %>%
  ggplot() +
  geom_tile(aes(x = Day, y = Species, fill = value), color = NA) +
  scale_x_continuous("Julian date", expand = c(0,0), breaks = seq(0, 366, by = 50)) +
  scale_y_discrete("Species", expand = c(0,0)) +
  scale_fill_viridis_c("Posterior mean phenological intensity", 
                       option = "D", limits = c(0,1),
                       guide = guide_colorbar(barwidth = unit(5.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  ggtitle("(a) Phenological activity patterns") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7, vjust = 2),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        aspect.ratio = 1,
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
A

#-------------------------------------------------------------------------------
# Inter-annual trends
#-------------------------------------------------------------------------------

trend_matrix <- abind(lapply(1:nrow(stanfit), function(i) matrix(stanfit[i, grep("^trend\\[", colnames(stanfit))], ncol = datalist$N_species, nrow = datalist$N_years, byrow = F)), along = 3)
year_range <- seq(0, 1, length.out = 49)
trend_slopes <- do.call(rbind, lapply(1:datalist$N_species, function(j) sapply(1:nrow(stanfit), function(i) summary(lm(trend_matrix[,j,i] ~ year_range))$coef[2])))

B <- melt(trend_slopes, varnames = c("Species", ".draw")) %>%
  group_by(Species) %>%
  mutate(Trend_prob = mean(value > 0)*2 - 1) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = phenology_order)) %>%
  ggplot() +
  stat_interval(aes(x = value, y = reorder(Species, value, mean), color = Trend_prob), size = 0.5, .width = 0.99, alpha = 1/8) +
  stat_interval(aes(x = value, y = reorder(Species, value, mean), color = Trend_prob), size = 0.5, .width = 0.95, alpha = 1/4) +
  stat_interval(aes(x = value, y = reorder(Species, value, mean), color = Trend_prob), size = 0.5, .width = 0.8, alpha = 1/2) +
  stat_interval(aes(x = value, y = reorder(Species, value, mean), color = Trend_prob), size = 0.5, .width = 0.5, alpha = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous("Linearised trend slope", expand = c(0,0), breaks = seq(-15, 15, by = 5)) +
  scale_y_discrete("Species", expand = c(0,0)) +
  scale_color_gradient2("Posterior probability of a negative (red) or positive (blue) trend",
                        breaks = seq(-1, 1, by = 0.5), labels = c(1, 0.5, 0, 0.5, 1),
                        low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                        guide = guide_colorbar(barheight = unit(0.2, "cm"), barwidth = unit(5.5, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  coord_cartesian(xlim = c(-16, 16)) +
  ggtitle("(b) Inter-annual linearised trend slopes") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.minor.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7, vjust = 2),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        aspect.ratio = 1,
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
B

# Tallying the number of declining and increasing species
melt(trend_slopes, varnames = c("Species", ".draw")) %>%
  group_by(Species) %>%
  summarize(Trend_class = factor(case_when(mean(value > 0) > 0.95 ~ "Increasing",
                                           mean(value < 0) > 0.95 ~ "Decreasing",
                                           T ~ "Uncertain"),
                                 levels = c("Decreasing", "Uncertain", "Increasing"))) %>%
  ungroup() %>%
  group_by(Trend_class) %>%
  tally()

#-------------------------------------------------------------------------------
# Predictors
#-------------------------------------------------------------------------------

predictor_matrix <- abind(lapply(1:nrow(stanfit), function(i) matrix(stanfit[i, grep("^Beta\\[", colnames(stanfit))], ncol = datalist$N_covariates, nrow = datalist$N_species, byrow = T)), along = 3)

C <- melt(predictor_matrix, varnames = c("Species", "Covariate", ".draw")) %>%
  group_by(Species, Covariate) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  mutate(value = case_when(value > 5 ~ 5,
                           value < -5 ~ -5,
                           T ~ value),
         Covariate = colnames(datalist$X)[Covariate],
         Covariate = case_when(Covariate == "MethodHolteval" ~ "Crevace",
                               Covariate == "MethodOndergrondse bodemval" ~ "Underground",
                               Covariate == "MethodPyramideval" ~ "Pyramidal",
                               Covariate == "MethodRamptrap" ~ "Ramp",
                               Covariate == "MethodWaterbodemval" ~ "Aquatic",
                               Covariate == "Built.up" ~ "Built-up",
                               Covariate == "Decidious_forest" ~ "Deciduous forest",
                               Covariate == "Grassland" ~ "Semi-nat. grassland",
                               Covariate == "Grassland_agri" ~ "Agric. grassland",
                               Covariate == "Pine_forest" ~ "Pine forest",
                               Covariate == "Valley_forest" ~ "Valley forest",
                               T ~ Covariate),
         Group = case_when(Covariate %in% c("Crevace", "Underground", "Pyramidal", "Ramp", "Aquatic") ~ "Pitfall trapping\nvariant",
                           T ~ "Land cover\ntype")) %>%
  ggplot() +
  geom_tile(aes(x = Covariate, y = factor(Species), fill = value)) +
  scale_x_discrete("Covariate", expand = c(0,0), position = "top") +
  scale_y_discrete("Species", expand = c(0,0)) +
  scale_fill_gradient2("Posterior\nmean\neffect\n",
                       breaks = seq(-5, 5, by = 2.5), labels = c("<-5", '-2.5', "0", "2.5", ">5.0"),
                       limits = c(-5, 5),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(5.5, "cm"), title.position = "top", title.hjust = 0, title.vjust = -2)) +
  facet_grid(~ Group, scales = "free_x", space = "free_x", switch = "x") +
  ggtitle("(c) Covariate influences") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        plot.title = element_text(face = "bold", size = 7),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        strip.clip = "off",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        legend.position = "right",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
C

#-------------------------------------------------------------------------------
# Site loadings
#-------------------------------------------------------------------------------

grid <- st_read("Data/limburg_grid.shp")

site_loadings <- abind(lapply(1:nrow(stanfit), function(i) matrix(stanfit[i, grep("^site_loadings\\[", colnames(stanfit))], ncol = 10, nrow = datalist$N_sites, byrow = F)), along = 3)
site_loadings_singledim <- site_loadings[,,1:(nrow(stanfit)/8)]

D <- melt(site_loadings_singledim, varnames = c("site", "dim", ".draw")) %>%
  group_by(site, dim) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  filter(dim <= 2) %>%
  mutate(dim = paste0("LV", dim),
         value = case_when(value > 1 ~ 1,
                           value < -1 ~ -1,
                           T ~ value)) %>%
  left_join(grid) %>%
  st_as_sf(sf_column_name = "geometry") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = NA) +
  scale_fill_gradient2("Posterior mean site loading",
                       limits = c(-1,1), breaks = seq(-1, 1, by = 0.5), labels = c("<-1", "-0.5", "0", "0.5", ">1"),
                       low = "#c23662", mid = "#FEF4DE", high = "#0e9296", midpoint = 0,
                       guide = guide_colorbar(barheight = unit(0.2, "cm"), barwidth = unit(5.5, "cm"), title.position = "top", title.hjust = 0.5, title.vjust = -2)) +
  facet_wrap(~ dim, ncol = 2, nrow = 1, strip.position = "top") +
  ggtitle("(d) Constrained ordination loadings") +
  theme_void() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-7.5,0,0,0),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 6),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 6, vjust = -10, hjust = 0.85),
        strip.clip = "off",
        plot.title = element_text(face = "bold", size = 7),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
D

#-------------------------------------------------------------------------------
# Species loadings
#-------------------------------------------------------------------------------

species_loadings <- abind(lapply(1:nrow(stanfit), function(i) matrix(stanfit[i, grep("^species_loadings\\[", colnames(stanfit))], nrow = 10, ncol = datalist$N_species, byrow = F)), along = 3)
species_loadings_singledim <- species_loadings[,,1:(nrow(stanfit)/8)]

E <- melt(species_loadings_singledim, varnames = c("dim", "species", ".draw")) %>%
  group_by(species, dim) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(names_from = dim, values_from = value) %>%
  ggplot() +
  geom_point(aes(x = `1`, y = `2`), size = 0.2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous("Species loading (LV1)") +
  scale_y_continuous("Species loading (LV2)") +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        plot.title = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 6),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
E

#-------------------------------------------------------------------------------
# Composite plot
#-------------------------------------------------------------------------------

(A | B) / ((C | free((D / free(E)) + plot_layout(heights = c(0.55, 1)))) + plot_layout(widths = c(0.8, 1))) + plot_layout(heights = c(0.8, 1))
ggsave("Plots/JSDM_output.png", width = 16, height = 22, dpi = 600, units = "cm")
