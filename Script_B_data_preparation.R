################################################################################
# TEMPORAL DISAGGREGATION MODEL
# Script B - Data preparation for the case study
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

library(tidyverse)
library(lubridate)
library(mgcv)
library(splines)
library(sf)
library(ggspatial)

#-------------------------------------------------------------------------------
# Reading data
#-------------------------------------------------------------------------------

# Main pitfall trapping data
data <- read.csv("Data/spider_data_limburg.csv", sep = ",") %>%
  mutate(Startdate = as.Date(Startdate),
         Enddate = as.Date(Enddate),
         Midpoint = as.Date(Midpoint)) %>%
  filter(year(Midpoint) >= 1900)

# Shapefile of the province of Limburg
limburg <- read_sf("Data/Refprv.shp") %>%
  filter(NAAM == "Limburg") %>%
  st_transform(crs = 31370)

#-------------------------------------------------------------------------------
# Visualizing sampling locations and sampling intervals
#-------------------------------------------------------------------------------

# Sampling effort across space
ggplot() +
  geom_sf(data = limburg) +
  geom_sf(data = data %>%
            st_as_sf(coords = c("X", "Y"), crs = 31370), size = 0.25, shape = 16, alpha = 1) +
  annotation_scale(location = "br", style = "bar", width_hint = 0.5/2, height = unit(0.6/2, "lines"), tick_height = 0.55/2, text_cex = 1.5/2, line_width = 0.6/2, text_pad = unit(0.5/2, "lines")) + #, pad_y = unit(0.2, "lines"), text_pad = unit(-0.7, "lines")
  annotation_north_arrow(location = "tr", which_north = "true", width = unit(3/2, "lines"), height = unit(3/2, "lines"), style = north_arrow_orienteering(line_width = 0.1/2, text_size = 14/2)) +
  theme_void()
ggsave("Plots/Sampling_effort_space.pdf", width = 16, height = 16, units = "cm", dpi = 600)

# Sampling effort across years
data %>%
  select(Startdate, Enddate) %>%
  arrange(Startdate) %>%
  rowid_to_column() %>%
  ggplot() +
  geom_segment(aes(y = rowid, yend = rowid, x = Startdate, xend = Enddate), color = "black", linewidth = 0.1) +
  scale_x_date("Date", expand = c(0,0)) +
  scale_y_reverse("Sampling event (ordered by date of deployment)", expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey93"),
        panel.grid.minor.x = element_line(color = "grey93"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7))
ggsave("Plots/Sampling_effort_years.pdf", width = 16, height = 16, units = "cm", dpi = 600)

# Sampling effort across days
plotdata <- data %>%
  select(Startdate, Enddate) %>%
  mutate(Startdate = yday(Startdate),
         Enddate = yday(Enddate)) %>%
  arrange(Startdate) %>%
  rowid_to_column()
ggplot() +
  geom_segment(data = plotdata %>%
                 filter(Enddate > Startdate),
               aes(y = rowid, yend = rowid, x = Startdate, xend = Enddate), color = "black", linewidth = 0.1) +
  geom_segment(data = plotdata %>%
                 filter(Enddate < Startdate),
               aes(y = rowid, yend = rowid, x = 1, xend = Enddate), color = "black", linewidth = 0.1) +
  geom_segment(data = plotdata %>%
                 filter(Enddate < Startdate),
               aes(y = rowid, yend = rowid, x = Startdate, xend = 365), color = "black", linewidth = 0.1) +
  scale_x_continuous("Julian date", expand = c(0,0)) +
  scale_y_reverse("Sampling event (ordered by Julian date of deployment)", expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey93"),
        panel.grid.minor.x = element_line(color = "grey93"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7))
ggsave("Plots/Sampling_effort_days.pdf", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Preparing phenological basis functions
#-------------------------------------------------------------------------------

N_bf <- 15
bf_d0 <- cSplineDes(1:366, seq(1, 366, length.out = N_bf + 1), ord = 4, derivs = 0)
bf_d1 <- cSplineDes(1:366, seq(1, 366, length.out = N_bf + 1), ord = 4, derivs = 1)
bf_d2 <- cSplineDes(1:366, seq(1, 366, length.out = N_bf + 1), ord = 4, derivs = 2)

#-------------------------------------------------------------------------------
# Preparing spatial basis functions
#-------------------------------------------------------------------------------

limburg_grid <- limburg %>%
  st_make_grid(1000) %>%
  st_intersection(st_geometry(limburg)) %>%
  data.frame() %>%
  rowid_to_column(var = "site") %>%
  st_as_sf() %>%
  cbind(st_coordinates(st_centroid(.)))
#st_write(limburg_grid, "Data/limburg_grid.shp")

site_ids <- data %>%
  select(X, Y) %>%
  st_as_sf(coords = c("X","Y"), crs = 31370) %>%
  st_join(limburg_grid) %>%
  pull(site)

xy <- limburg_grid %>%
  st_drop_geometry() %>%
  select(site, X, Y) %>%
  mutate(site = factor(site)) %>%
  arrange(site) %>%
  select(X,Y) %>%
  mutate(X = X - min(X),
         Y = Y - min(Y),
         X_rescaled = 2 * (X / max(c(X,Y))),
         Y_rescaled = 2 * (Y / max(c(X,Y))),
         X_rescaled = X_rescaled - max(X_rescaled)/2,
         Y_rescaled = Y_rescaled - max(Y_rescaled)/2) %>%
  dplyr::select(X_rescaled, Y_rescaled) %>%
  as.matrix()

N_spatial_bf_1D <- 15
spatial_bf_X <- bs(xy[,1], Boundary.knots = c(-1.2,1.2), knots = seq(-1.2, 1.2, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)]
spatial_bf_Y <- bs(xy[,2], Boundary.knots = c(-1.2,1.2), knots = seq(-1.2, 1.2, length.out = N_spatial_bf_1D - 2))[,-c(N_spatial_bf_1D + 1)]
spatial_bf <- do.call(rbind, lapply(1:nrow(xy), function(i) kronecker(spatial_bf_X[i,], spatial_bf_Y[i,])))
spatial_bf_range <- expand.grid(x = seq(-1, 1, length.out = N_spatial_bf_1D),
                                y = seq(-1, 1, length.out = N_spatial_bf_1D))[colSums(spatial_bf) > 0,]
spatial_bf <- spatial_bf[,colSums(spatial_bf) > 0]

# Visualize a single random realization to assess the maximum spatial grain that can be modelled
ggplot(data.frame(xy, value = spatial_bf %*% rnorm(ncol(spatial_bf)))) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()
# Visualize standard deviation of a large number of B-spline realizations to ensure no strong artefacts are present
ggplot(data.frame(xy, value = apply(replicate(1000, spatial_bf %*% rnorm(ncol(spatial_bf)), simplify = T), 1, sd))) + geom_point(aes(x = X_rescaled, y = Y_rescaled, color = value)) + scale_color_viridis_c() + coord_equal()

#-------------------------------------------------------------------------------
# Prepare the community matrix Y (sample x species)
#-------------------------------------------------------------------------------

Y <- as.matrix(select(data, starts_with("Sp..")))
# Only species detected in at least 10 samples are retained
Y <- Y[,colSums(apply(apply(Y, 1:2, as.logical), 1:2, as.numeric)) >= 10]
Y <- Y[,order(-colSums(Y))]

#-------------------------------------------------------------------------------
# Prepare the trapping method and habitat predictor matrix
#-------------------------------------------------------------------------------

# Read and preprocess land cover data
Habitat <- read.csv("Data/Limburg_grid.csv") %>%
  right_join(limburg_grid) %>%
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
  select(-site, -X, -Y, -geometry)
Habitat <- Habitat[,order(colnames(Habitat))]

# Combine trapping method data with land cover data 
X <- model.matrix(~ 0 + Method, data = data)[,-1] %>%
  cbind(Habitat[site_ids,-1])

# Prepare matrix for predictive purposes
X_pred <- cbind(matrix(0, nrow = nrow(Habitat), ncol = 5), Habitat[,-1])

#-------------------------------------------------------------------------------
# Prepare the data list for Stan
#-------------------------------------------------------------------------------

datalist <- list(N_samples = nrow(Y),
                 N_species = ncol(Y),
                 N_covariates = ncol(X),
                 N_sites = nrow(limburg_grid),
                 N_spatial_bf = ncol(spatial_bf),
                 N_days = 366,
                 N_years = max(data$Year) - min(data$Year) + 1,
                 N_bf = N_bf,
                 Y = Y,
                 X = X,
                 X_pred = X_pred,
                 site = site_ids,
                 year = data$Year - min(data$Year) + 1,
                 year_range = seq(-1, 1, length.out = max(data$Year) - min(data$Year) + 1),
                 spatial_bf = spatial_bf,
                 spatial_bf_range = spatial_bf_range,
                 bf_d0 = bf_d0,
                 bf_d1 = bf_d1,
                 bf_d2 = bf_d2,
                 c = 0.5*sapply(1:nrow(data), function(i) sum((as.numeric(seq(data$Startdate[i], data$Enddate[i], by = "1 day") - data$Midpoint[i]))^2)/data$Duration[i]),
                 log_m = log(data$Duration),
                 m = data$Duration,
                 start_day = yday(data$Startdate),
                 end_day = yday(data$Enddate),
                 midpoint = yday(data$Midpoint),
                 bf_range = seq(-1, 1, length.out = N_bf + 1)[1:N_bf])
saveRDS(datalist, "Data/datalist_2025_07_17.rds")

