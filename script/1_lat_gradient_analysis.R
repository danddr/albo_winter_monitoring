# CONTINENTAL-SCALE ANALYSIS: AEDES ALBOPICTUS WINTER ACTIVITY
# Author: Daniele Da Re
# Date: April 2026
# Purpose:  European Latitudinal Gradient (2024-2025) GAM analysis of winter ovitrap surveillance data
setwd("/home/dared/GitHub/albo_winter_monitoring/")

# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(lubridate)
library(ISOweek)
library(sf)
library(ggrepel)
library(viridis)
library(geodata)
library(terra)
library(mgcv)
library(metR)
library(patchwork)
library(yardstick)
library(pROC)
library(corrplot)

outdir <- "outputs/latgrad/"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 1. Study locations ------------------------------------------------------
city_coords <- tibble(
  City = c("Palermo (ITA)", "Bari (ITA)", "Fondi (ITA)", "Roma (ITA)", 
           "Siena (ITA)", "Modena (ITA)", "Trento (ITA)", "Vienna (AUT)", 
           "Graz (AUT)", "Allschwil (CHE)", "Strasbourg (FRA)", "Bordeaux (FRA)", 
           "Montpellier (FRA)", "Bayonne (FRA)"),
  Latitude = c(38.1157, 41.1171, 41.3565, 41.9028, 43.3188, 44.6471, 46.0667, 
               48.2082, 47.0707, 47.5596, 48.5734, 44.8378, 43.6119, 43.4934),
  Longitude = c(13.3615, 16.8719, 13.4435, 12.4964, 11.3316, 10.9252, 11.1256, 
                16.3738, 15.4395, 7.5886, 7.7521, -0.5792, 3.8772, -1.4757)
) %>%
  arrange(desc(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)


# Country boundaries
mask <- geodata::gadm(
  country = c("ITA", "BEL", "ESP", "ALB", "NL", "CH", "LUX", "FRA", "SVN", 
              "AUT", "LIE", "DEU", "SRB", "SVK", "ROU", "CZE", "HUN", "TUN", 
              "DZA", "AND", "POL", "HRV", "MNE", "BIH", "MKD", "GRC", "XK"),
  level = 0,
  path = "data/"
)

ER_mask <- geodata::gadm(
  country = c("ITA"),
  level = 1,
  path = "data/"
) 
ER_mask <- subset(ER_mask, ER_mask$NAME_1 =="Emilia-Romagna")

# Color palette
myCols <- viridis(nrow(city_coords), option = "C")

# Map
city_df <- city_coords %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )

p_map <- ggplot() +
  geom_sf(data = st_as_sf(mask), linewidth = 0.2, fill = "grey90", color = "grey60") +
  geom_sf(data = st_as_sf(ER_mask), linewidth = 0.2, fill = "grey50", color =  "grey60") +
  geom_sf(data = city_coords, aes(fill = City), size = 4.5, shape = 21, 
          color = "black", stroke = 0.6) +
  geom_text_repel(
    data = city_df,
    aes(x = lon, y = lat, label = City), color = "black",
    size = 4.5, 
    box.padding = 0.7,
    fontface = "bold", max.overlaps = Inf
  ) +
  labs(x="Longitude", y="Latitude")+
  ggspatial::annotation_scale(width_hint = 0.25)+ # proportion of the plot
  scale_fill_manual(values = myCols, breaks = city_coords$City) +
  coord_sf(xlim = c(-5, 20), ylim = c(36, 51.1), expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

p_map

ggsave(paste0(outdir, "map_locations.png"), p_map, 
       width = 16, height = 10, dpi = 320)

# 2. Load and prepare data ------------------------------------------------

# Entomological observations
df <- read_xlsx("data/Albo_WinterActivity_2024_2025_monitoring.xlsx", sheet = 1) %>% 
  pivot_longer(!c(year, week)) %>% 
  rename(loc = name, pos = value) %>% 
  mutate(
    date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")),
    week = as.numeric(week),
    year = as.factor(year),
    loc = as.factor(loc),
    pos = as.integer(pos)
  ) %>% 
  left_join(
    city_coords %>% 
      bind_cols(st_coordinates(city_coords)) %>% 
      st_drop_geometry() %>% 
      rename(lat = Y) %>% 
      select(loc = City, lat),
    by = "loc"
  )

# summary stats by location
sum.stat <- df %>% 
  group_by(loc, year) %>% 
  reframe(minw=min(week),
          maxw=max(week))

# Climate data
mypath <- "/home/dared/Documents/PoD/marieCurie2022/ZanZemap/ML/Download_Process_ClimData/2025/data/weekly/"
tas <- rast(paste0(mypath,"t2m_min_201001-202512_C.nc"))
tp <- rast(paste0(mypath,"tp_weekly_cumulative_201001-202512.nc"))

# Set dates
myDates <- seq.Date(as.Date("2010-01-07"), length.out = nlyr(tas), by = "week")
myNames <- paste0("w_", year(myDates), "_", week(myDates))
time(tas) <- myDates
time(tp) <- myDates
names(tas) <- myNames
names(tp) <- myNames

# Subset from 2023 onwards
tas <- tas[[679:nlyr(tas)]]
tp <- tp[[679:nlyr(tp)]]

# Gap-fill coastal NAs
w <- matrix(1, 5, 5)
tas <- terra::focal(tas, w, mean, na.rm = TRUE, NAonly = TRUE, pad = TRUE)
tp <- terra::focal(tp, w, mean, na.rm = TRUE, NAonly = TRUE, pad = TRUE)

# Extract climate data
loc <- bind_cols(
  Location = city_coords$City, 
  city_coords %>% st_coordinates()
) %>% 
  rename(Longitude = X, Latitude = Y)

loc.v <- vect(loc, geom = c("Longitude", "Latitude"), crs = "EPSG:4326") 

tas.df <- extract(tas, loc.v)
names(tas.df) <- c("ID", as.character(time(tas)))
tas.df <- tas.df %>%
  as_tibble() %>%
  mutate(Location = loc$Location) %>%
  select(-ID) %>%
  pivot_longer(-Location, names_to = "date", values_to = "tas") %>%
  mutate(
    date = as.Date(date),  # Convert BEFORE any other operations
    year = year(date), 
    week = week(date), 
    tid = paste(Location, year, week, sep = "_")
  )

# Extract precipitation data
tp.df <- extract(tp, loc.v)
names(tp.df) <- c("ID", as.character(time(tp)))
tp.df <- tp.df %>%
  as_tibble() %>%
  mutate(Location = loc$Location) %>%
  select(-ID) %>%
  pivot_longer(-Location, names_to = "date", values_to = "tp") %>%
  mutate(
    date = as.Date(date),
    year = year(date), 
    week = week(date), 
    tid = paste(Location, year, week, sep = "_")
  )

# Combine
tas.df$tp <- tp.df$tp

# Compute photoperiod
ph <- loc %>% 
  select(Location, Latitude) %>%
  rename(lat = Latitude) %>% 
  distinct()

myDates <- data.frame(date = seq.Date(as.Date('2023-01-01'), 
                                      as.Date('2026-01-01'), by = "day"))
myDates$week <- week(myDates$date)
myDates$year <- year(myDates$date)
phOut <- list()
for(i in 1:nrow(ph)){
  tmpDat <- myDates
  phOut[[i]] <- tmpDat %>%
    mutate(
      photop = geosphere::daylength(lat = ph$lat[i], doy = tmpDat$date)
    ) %>%
    group_by(year, week) %>%
    summarise(meanPhotoweek = mean(photop, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Location = ph$Location[i],
      tid = paste(Location, year, week, sep = "_")
    ) %>%
    filter(year >= 2023, !(year == 2025 & week > 25))
}
phOut <- bind_rows(phOut)

tas.df <- tas.df %>% 
  left_join(phOut %>% select(tid, meanPhotoweek), by = "tid")

# Combine observations with climate
df <- df %>% 
  mutate(tid = paste(loc, year, week, sep = "_")) %>%
  select(-date) %>% 
  left_join(tas.df %>% select(tid, tas, tp, meanPhotoweek, date), by = "tid")

# 3. Create lagged predictors ---------------------------------------------
df.m <- df %>%
  group_by(loc) %>%
  arrange(year, week, .by_group = TRUE) %>%
  mutate(
    tas_lag1 = lag(tas, n = 1),
    tas_lag2 = lag(tas, n = 2),
    tas_lag3 = lag(tas, n = 3),
    tp_lag1 = lag(tp, n = 1),
    tp_lag2 = lag(tp, n = 2),
    tp_lag3 = lag(tp, n = 3),
    loc = factor(loc)
  ) %>%
  rename(ph = meanPhotoweek) %>%
  ungroup() %>%
  drop_na() %>%
  filter(!loc %in% c("Siena (ITA)", "Vienna (AUT)")) # because they are always zero

# 4. Model selection ------------------------------------------------------
 

# COLLINEARITY ANALYSIS
clim_cors <- df.m %>%
  select(tas, tas_lag1, tas_lag2, tas_lag3, tp, tp_lag1, tp_lag2, tp_lag3, ph) %>%
  drop_na() %>%
  cor()

corrplot(clim_cors, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45, tl.cex = 1.2,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Climate Variable Correlations",
         mar = c(0, 0, 2, 0))

clim_cors
# I keep only tas_lag3 because is the least correlated to photoperiod (<0.7) 
# but still highly correlated with the other temperatures

# Full model with all predictors
m_full_env <- gam(
  pos ~ 
    s(tas_lag3, k = 10) +
    tp + tp_lag1 + tp_lag2 + tp_lag3 +
    s(ph, k = 10) +
    s(loc, bs = "re"),
  family = binomial,
  method = "REML",
  data = df.m
)
# Check overdispersion
overdisp <- sum(residuals(m_full_env, type = "pearson")^2) / m_full_env$df.residual
overdisp
# mild  undersdispersion

summary(m_full_env)

# Remove precipitation terms sequentially
m_no_precip1 <- update(m_full_env, . ~ . - tp)
AIC(m_full_env, m_no_precip1)
anova(m_full_env, m_no_precip1, test = "Chisq")

m_no_precip2 <- update(m_no_precip1, . ~ . - tp_lag1)
AIC(m_full_env, m_no_precip2)
anova(m_full_env, m_no_precip2, test = "Chisq")

m_no_precip3 <- update(m_no_precip2, . ~ . - tp_lag2)
AIC(m_full_env, m_no_precip3)

m_no_precip4 <- update(m_no_precip3, . ~ . - tp_lag3)
AIC(m_full_env, m_no_precip4)
anova(m_full_env, m_no_precip4, test = "Chisq")
summary(m_no_precip4)
gam.check(m_no_precip4)
concurvity(m_no_precip4)

 # Simplify to linear temperature terms (EDF = 1.0 indicates linear)
m_final_simple <- gam(
  pos ~ 
    tas_lag3+
    s(ph, k = 10) +
    s(loc, bs = "re"),
  family = binomial,
  method = "REML",
  data = df.m
)

summary(m_final_simple)
concurvity(m_final_simple)
gam.check(m_final_simple)
overdisp <- sum(residuals(m_final_simple, type = "pearson")^2) / m_final_simple$df.residual

AIC(m_no_precip4, m_final_simple)
anova(m_no_precip4, m_final_simple, test = "Chisq")

# Final model comparison table
model_comparison <- data.frame(
  Model = c("Full", "No precip (tp, lag1, lag2, lag3)",  "Linear temp (final)"),
  AIC = c(AIC(m_full_env), AIC(m_no_precip4), AIC(m_final_simple)),
  R2 = c(summary(m_full_env)$r.sq, summary(m_no_precip4)$r.sq, 
          summary(m_final_simple)$r.sq),
  Dev_Expl = c(summary(m_full_env)$dev.expl, summary(m_no_precip4)$dev.expl,
                summary(m_final_simple)$dev.expl)
)

print(model_comparison)
write.csv(model_comparison, paste0(outdir, "model_selection.csv"), row.names = FALSE)

# 5. Partial dependence plots ---------------------------------------------
ref_loc <- levels(df.m$loc)[1]

# Function to create partial plots
create_partial_plot <- function(var, var_name, color, data, model) {
  
  pred_data <- data.frame(x = seq(min(data[[var]], na.rm = TRUE), 
                                  max(data[[var]], na.rm = TRUE), 
                                  length.out = 200))
  names(pred_data) <- var
  
  # Hold other predictors at median
  # if(var != "tas") pred_data$tas <- median(data$tas, na.rm = TRUE)
  if(var != "tas_lag3") pred_data$tas_lag3 <- median(data$tas_lag3, na.rm = TRUE)
  # if(var != "tp_lag3") pred_data$tp_lag3 <- median(data$tp_lag3, na.rm = TRUE)
  if(var != "ph") pred_data$ph <- median(data$ph, na.rm = TRUE)
  pred_data$loc <- factor(ref_loc, levels = levels(data$loc))
  
  preds <- predict(model, newdata = pred_data, 
                   type = "link", se.fit = TRUE, exclude = "s(loc)")
  
  pred_data$fit <- plogis(preds$fit)
  pred_data$lwr <- plogis(preds$fit - 1.96 * preds$se.fit)
  pred_data$upr <- plogis(preds$fit + 1.96 * preds$se.fit)
  
  ggplot(pred_data, aes(x = .data[[var]], y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = color, alpha = 0.3) +
    geom_line(color = color, linewidth = 1.2) +
    geom_rug(data = data, aes(x = .data[[var]], y = NULL), 
             sides = "b", alpha = 0.3, length = unit(0.02, "npc")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = var_name, y = "Predicted probability") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p1 <- create_partial_plot("tas_lag3", "Average minimum temperature 3 weeks before (°C)", "#0072B2", df.m, m_final_simple) +
  labs(title = "A.")

p2 <- create_partial_plot("ph", "Photoperiod (hours)", "#E69F00" , df.m, m_final_simple) +
  labs(title = "B.")

combined_plot <- (p1 | p2) 
  # plot_annotation(
  # subtitle = paste0("Model R² = 0.664 | n = ", nrow(df.m), " trap-weeks"))
combined_plot
ggsave(paste0(outdir, "gam_partial_dependence.png"), combined_plot, 
       width = 12, height = 10, dpi = 300)

# 6. Heatmap: Photoperiod × Temperature lag-2 ----------------------------

# Prediction grid (for both marginal and location-specific plots)
ph_seq <- seq(min(df.m$ph, na.rm = TRUE), max(df.m$ph, na.rm = TRUE), length.out = 100)
tas_lag3_seq <- seq(min(df.m$tas_lag3, na.rm = TRUE), max(df.m$tas_lag3, na.rm = TRUE), length.out = 100)

## ---- 6.1 Overall heatmap (marginal effect, location RE excluded) ----
grid_2d <- expand.grid(ph = ph_seq, tas_lag3 = tas_lag3_seq) %>%
  mutate(
    tas = median(df.m$tas, na.rm = TRUE),
    tp_lag3 = median(df.m$tp_lag3, na.rm = TRUE),
    loc = factor(ref_loc, levels = levels(df.m$loc))
  )

pred_2d <- predict(m_final_simple, newdata = grid_2d, 
                   type = "link", se.fit = TRUE, exclude = "s(loc)")
grid_2d$prob <- plogis(pred_2d$fit)

obs_df <- df.m %>%
  filter(!is.na(tas_lag3), !is.na(ph)) %>%
  mutate(state = factor(pos, levels = c(0, 1), labels = c("Absence", "Presence")))

p_heat_overall <- ggplot() +
  geom_tile(data = grid_2d, aes(x = tas_lag3, y = ph, fill = prob)) +
  
  geom_point(data = obs_df, aes(x = tas_lag3, y = ph, shape = state, colour = state),
             size = 3, alpha = 0.7) +
  scale_fill_viridis_c(name = "Predicted\nProbability", limits = c(0, 1), option = "plasma") +
  scale_shape_manual(values = c(Absence = 1, Presence = 16), name = "Observed") +
  scale_colour_manual(values = c(Absence = "white", Presence = "red3"), name = "Observed") +
  geom_contour(data = grid_2d, aes(x = tas_lag3, y = ph, z = prob),
               breaks = c(0.1, 0.2, 0.5, 0.8), color = "white", linewidth = 0.6) +
  geom_text_contour(data = grid_2d, aes(x = tas_lag3, y = ph, z = prob),
                    breaks = c(0.1, 0.2, 0.5, 0.8), color = "white", size = 6, skip = 0) +
  labs(
    x = expression(paste("Average minimum temperature 3 weeks before (°C)")),
    y = "Photoperiod (hours)"
    # title = "A. Marginal Effect (Location Random Effects Excluded)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )
p_heat_overall
ggsave(paste0(outdir, "heatmap_overall.png"), p_heat_overall, 
       width = 12, height = 8, dpi = 300)

## ---- 6.2 Location-specific heatmap (location RE included) ----

# Create prediction grid for ALL locations
grid_loc <- expand.grid(
  ph = ph_seq,
  tas_lag3 = tas_lag3_seq,
  loc = levels(df.m$loc)
) %>% 
  filter(!loc %in% c("Siena (ITA)", "Vienna (AUT)")) # because they are always zero
 
# Predict (INCLUDING location random effect)
pred_loc <- predict(m_final_simple, newdata = grid_loc, type = "response")
grid_loc$prob <- pred_loc

# Prepare observed data
obs_df_facet <- df.m %>%
  filter(!is.na(tas_lag3), !is.na(ph)) %>%
  mutate(
    state = factor(pos, levels = c(0, 1), labels = c("Absence", "Presence")),
    loc = factor(loc, levels = levels(grid_loc$loc))
  )

# Get city coordinates with latitude
city_coords_ordered <- city_coords %>%
  bind_cols(st_coordinates(city_coords)) %>% 
  st_drop_geometry() %>% 
  rename(lon = X, lat = Y) %>% 
  arrange(desc(lat)) %>%  # North to south (change to arrange(lat) for south to north)
  mutate(
    Color = myCols,
    # Create ordered factor for faceting
    City_ordered = factor(City, levels = City)
  )

# Add latitude and ordered factor to prediction grid
grid_loc <- grid_loc %>% 
  left_join(
    city_coords_ordered %>% select(City, lat, City_ordered), 
    by = c("loc" = "City")
  ) %>%
  # Convert loc to ordered factor matching latitude order
  mutate(
    loc = factor(loc, levels = levels(city_coords_ordered$City_ordered))
  )

# Add ordered factor to observations
obs_df_facet <- obs_df_facet %>%
  mutate(
    loc = factor(loc, levels = levels(city_coords_ordered$City_ordered))
  )

# Create faceted heatmap with ordered facets
p_heat_facets <- ggplot() +
  geom_tile(data = grid_loc, aes(x = tas_lag3, y = ph, fill = prob)) +
  geom_contour(data = grid_loc, aes(x = tas_lag3, y = ph, z = prob),
               breaks = c(0.1, 0.2, 0.5, 0.8), color = "white", linewidth = 0.6) +
  geom_text_contour(data = grid_loc, aes(x = tas_lag3, y = ph, z = prob),
                    breaks = c(0.1, 0.2, 0.5, 0.8), color = "white", size = 4, skip = 0) +
  geom_point(data = obs_df_facet, aes(x = tas_lag3, y = ph, shape = state, colour = state),
             size = 2, alpha = 0.8) +
  # Key: facet_wrap with ordered factor
  facet_wrap(~loc, ncol = 3) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Predicted\nProbability", option = "plasma") +
  scale_shape_manual(values = c(Absence = 1, Presence = 16), name = "Observed") +
  scale_colour_manual(values = c(Absence = "grey95", Presence = "red3"), name = "Observed") +
  labs(
    x = expression( "Average minimum temperature 3 weeks before (°C)"),
    y = "Photoperiod (hours)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

p_heat_facets

ggsave(paste0(outdir, "heatmap_by_location.png"), p_heat_facets, 
       width = 14, height = 12, dpi = 300)

# 7. Model evaluation -----------------------------------------------------
pred_vec <- predict(m_final_simple, newdata = df.m, type = "response")

df_eval <- df.m %>%
  mutate(
    pred = pred_vec,
    class_pred = factor(ifelse(pred > 0.5, 1, 0), levels = c(0, 1)),
    truth = factor(pos, levels = c(0, 1))
  )

## ---- 7.1 Overall performance metrics ----
roc_obj <- roc(df_eval$pos, df_eval$pred)

perf_table <- tibble(
  Metric = c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1 Score", "ROC AUC"),
  Value = c(
    accuracy(df_eval, truth, class_pred)$.estimate,
    sensitivity(df_eval, truth, class_pred)$.estimate,
    specificity(df_eval, truth, class_pred)$.estimate,
    precision(df_eval, truth, class_pred)$.estimate,
    f_meas(df_eval, truth, class_pred)$.estimate,
    as.numeric(auc(roc_obj))
  )
)
print(perf_table)

## ---- 7.2 Optimal threshold (Youden's index) ----
coords_optimal <- coords(roc_obj, x = "best", 
                         best.method = "closest.topleft",
                         ret = c("threshold", "sensitivity", "specificity"))

optimal_threshold <- coords_optimal$threshold
round(optimal_threshold,3)

# Recalculate metrics with optimal threshold
df_eval_optimal <- df_eval %>%
  mutate(
    class_pred_optimal = factor(ifelse(pred > optimal_threshold, 1, 0), 
                                levels = c(0, 1))
  )

perf_table_optimal <- tibble(
  Metric = c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1 Score", "ROC AUC", "Threshold"),
  Default_0.5 = c(
    accuracy(df_eval, truth, class_pred)$.estimate,
    sensitivity(df_eval, truth, class_pred)$.estimate,
    specificity(df_eval, truth, class_pred)$.estimate,
    precision(df_eval, truth, class_pred)$.estimate,
    f_meas(df_eval, truth, class_pred)$.estimate,
    as.numeric(auc(roc_obj)),
    0.5
  ),
  Optimal = c(
    accuracy(df_eval_optimal, truth, class_pred_optimal)$.estimate,
    sensitivity(df_eval_optimal, truth, class_pred_optimal)$.estimate,
    specificity(df_eval_optimal, truth, class_pred_optimal)$.estimate,
    precision(df_eval_optimal, truth, class_pred_optimal)$.estimate,
    f_meas(df_eval_optimal, truth, class_pred_optimal)$.estimate,
    as.numeric(auc(roc_obj)),
    optimal_threshold
  )
)

print(perf_table_optimal)
write_csv(perf_table_optimal, paste0(outdir, "model_performance_thresholds.csv"))

## ---- 7.3 Per-location performance ----

# Helper function for safe AUC calculation
safe_auc <- function(truth, prob) {
  truth_num <- as.numeric(as.character(truth))
  
  if (length(unique(truth_num)) < 2) return(NA_real_)
  
  roc_obj <- pROC::roc(truth_num, prob, quiet = TRUE)
  as.numeric(pROC::auc(roc_obj))
}

df_loc_perf <- df_eval %>%
  group_by(loc) %>%
  summarise(
    n_obs = n(),
    n_positive = sum(pos),
    prevalence = mean(pos),
    sensitivity = tryCatch(sensitivity_vec(truth, class_pred), 
                           warning = function(w) NA_real_),
    specificity = tryCatch(specificity_vec(truth, class_pred), 
                           warning = function(w) NA_real_),
    precision = tryCatch(precision_vec(truth, class_pred), 
                         warning = function(w) NA_real_),
    f1 = tryCatch(f_meas_vec(truth, class_pred), 
                  warning = function(w) NA_real_),
    accuracy = tryCatch(accuracy_vec(truth, class_pred), 
                        warning = function(w) NA_real_),
    auc = safe_auc(truth, pred),
    .groups = "drop"
  ) %>%
  arrange(desc(auc))

print(df_loc_perf)
write_csv(df_loc_perf, paste0(outdir, "model_performance_by_location.csv"))

## ---- 7.4 ROC curve plot ----
roc_data <- tibble(
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities
)
p_roc <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#0072B2", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(data = coords_optimal, 
             aes(x = 1 - specificity, y = sensitivity),
             color = "red", size = 4, shape = 18) +
  annotate("text", 
           x = 1 - coords_optimal$specificity + 0.1, 
           y = coords_optimal$sensitivity,
           label = paste0("Optimal threshold = ", round(optimal_threshold, 3)),
           hjust = 0, size = 4) +
  labs(
    title = "ROC Curve",
    subtitle = paste0("AUC = ", round(auc(roc_obj), 3)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_roc
ggsave(paste0(outdir, "roc_curve.png"), p_roc, 
       width = 8, height = 6, dpi = 300)

## ---- 7.5 Performance visualization by location ----
p_loc_perf <- df_loc_perf %>%
  filter(!is.na(auc)) %>%
  pivot_longer(cols = c(sensitivity, specificity, precision, accuracy, auc),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = reorder(loc, -value), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric, scales = "free_x") +
  coord_flip() +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Model Performance by Location",
    x = "Location",
    y = "Metric Value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )
p_loc_perf
# ggsave(paste0(outdir, "performance_by_location.png"), p_loc_perf, 
#        width = 12, height = 10, dpi = 300)

# 8. Seasonal predictions -------------------------------------------------
# Align city colors
city_coords_ordered <- city_coords %>%
  bind_cols(st_coordinates(city_coords)) %>% 
  st_drop_geometry() %>% 
  rename(lon = X, lat = Y) %>% 
  arrange(desc(lat)) %>%
  mutate(Color = myCols)

city_pal <- setNames(city_coords_ordered$Color, city_coords_ordered$City)

# Prepare prediction data
tas_pred <- tas.df %>%
  arrange(Location, date) %>% 
  group_by(Location) %>%
  mutate(
    tas_lag3 = lag(tas, 3),   
    # tp_lag3 = lag(tp, 3)
  ) %>%
  ungroup() %>% 
  filter(!is.na(tas_lag3)) %>%
  rename(loc = Location, ph = meanPhotoweek) %>%
  mutate(loc = factor(loc, levels = levels(df.m$loc))) %>% 
  filter(date >= "2024-08-31", date <= "2025-06-15")  # Limit until July 2025

pred_link <- predict(m_final_simple, tas_pred, type = "link", se.fit = TRUE)

tas_pred <- tas_pred %>%
  mutate(
    prob = plogis(pred_link$fit),
    lower = plogis(pred_link$fit - 1.96 * pred_link$se.fit),
    upper = plogis(pred_link$fit + 1.96 * pred_link$se.fit)
  ) %>%
  left_join(city_coords_ordered %>% select(City, Color), by = c("loc" = "City"))

# Filter observed data to same date range
df <- df %>%
  left_join(city_coords_ordered %>% select(City, Color), by = c("loc" = "City")) %>%
  filter(date >= "2024-08-31", date <= "2025-06-15")

# Ensure factor order
ordered_locs <- city_coords_ordered$City
df$loc <- factor(df$loc, levels = ordered_locs)
tas_pred$loc <- factor(tas_pred$loc, levels = ordered_locs)

p_seasonal <- ggplot() +
  geom_point(data = df %>%   filter(!loc %in% c("Siena (ITA)", "Vienna (AUT)")), 
             aes(x = date, y = pos),
             position = position_jitter(width = 0.3, height = 0.02),
             alpha = 0.4, color = "grey40", size = 3) +
  geom_ribbon(data = tas_pred %>% filter(!loc %in% c("Siena (ITA)", "Vienna (AUT)")), aes(x = date, ymin = lower, ymax = upper, fill = loc),
              alpha = 0.15, color = NA) +
  geom_line(data = tas_pred %>%  filter(!loc %in% c("Siena (ITA)", "Vienna (AUT)")), aes(x = date, y = prob, color = loc), linewidth = 1.1) +
  facet_wrap(~loc, scales = "free_y", ncol = 2) +
  # scale_x_date(
  #   limits = as.Date(c("2024-08-31", "2025-06-15")),
  #   date_breaks = "2 months",
  #   date_labels = "%b %Y"
  # ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = city_pal) +
  scale_fill_manual(values = city_pal) +
  labs(
    x = "Date (Oct 2024 - Jul 2025)",
    y = "Predicted probability of ovitrap positivity",
    # title = "Observed and Predicted Seasonal Ovitrap Positivity by Location"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p_seasonal
ggsave(paste0(outdir, "seasonal_predictions.png"), p_seasonal, 
       width = 16, height = 10, dpi = 320)
