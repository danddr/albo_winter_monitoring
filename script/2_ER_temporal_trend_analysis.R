# EMILIA-ROMAGNA REGIONAL ANALYSIS: Winter Aedes albopictus Activity (2011-2024)
# Author: Daniele Da Re
# Date: April 2026
# Purpose: Multi-year GAM analysis of winter ovitrap surveillance data

# Load packages -----------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(ISOweek)
library(sf)
library(terra)
library(mgcv)
library(ggplot2)
library(gratia)
library(patchwork)
library(viridis)
library(corrplot)
library(readr)
library(ggrepel)
# Setup -------------------------------------------------------------------
setwd("/home/dared/GitHub/albo_winter_monitoring/")
outdir <- "outputs/emilia_romagna/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. LOAD AND PREPARE OVITRAP DATA ----------------------------------------
# Load observations
df_raw <- read_xlsx(
  "data/dataset_ovitrappole_invernaliER_2011_2025.xlsx", 
  sheet = "dati"
) %>% 
  rename(
    istat = `istat comune`, 
    prov = provincia, 
    year = anno, 
    winter = inverno, 
    date = datamonitoraggio,  
    id = `cod esteso ovitrappola`
  ) %>% 
  select(id, istat, prov, year, winter, date, n_uova) %>% 
  mutate(
    week = week(date), 
    month = month(date), 
    winter = ifelse(
      month %in% 10:12, 
      paste0(year, "-", year + 1), 
      paste0(year - 1, "-", year)
    )
  )

# Load municipality codes
mun_codes <- read_xlsx(
  "data/dataset_ovitrappole_invernaliER_2011_2025.xlsx", 
  sheet = "codici_istat"
) %>% 
  select(CodiceIstat, Comune, lon, lat) %>% 
  rename(istat = CodiceIstat, mun = Comune)

# Join and create spatial object
df_obs <- df_raw %>% 
  left_join(mun_codes, by = "istat")

locations <- df_obs %>% 
  select(mun, lon, lat) %>% 
  distinct() %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Country boundaries
mask <- geodata::gadm(
  country = c("ITA"),
  level = 1,
  path = "data/"
)

smr <- geodata::gadm(
  country = c("SMR"),
  level = 0,
  path = "data/"
)
mask <- rbind(mask, smr)

p_map <- ggplot() +
  geom_sf(data = st_as_sf(mask), linewidth = 0.2, fill = "grey90", color = "grey60") +
  geom_sf(data = locations, aes(fill = mun), size = 4.5, shape = 21, 
          color = "black", stroke = 0.6) +
  geom_text_repel(
    data = df_obs %>% 
      select(mun, lon, lat) %>% 
      distinct(),
    aes(x = lon, y = lat, label = mun), color = "black",
    size = 4.5, 
    box.padding = 0.7,
    fontface = "bold", max.overlaps = Inf
  ) +
  labs(x="Longitude", y="Latitude")+
  ggspatial::annotation_scale(width_hint = 0.25)+ # proportion of the plot
  # scale_fill_manual(values = , breaks = city_coords$City) +
  coord_sf(xlim = c(9, 13), ylim = c(43.5, 45.5), expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

p_map

ggsave(paste0(outdir, "er_map_locations.png"), p_map, 
       width = 16, height = 10, dpi = 320)

# df_obs %>% 
#   group_by(mun, month, winter) %>% 
#   tally() %>% 
#   # filter(mun=="Bologna") %>% 
#   ggplot(aes(x=winter, y=n, col=mun))+
#   geom_point()+
#   facet_wrap(~as.factor(month))


# Aggregate by municipality and week
df_aggregated <- df_obs %>% 
  group_by(mun, year, winter, week, month) %>% 
  summarise(
    n_ovi = n(),
    n_pos.ovi = sum(n_uova > 1),
    prop.ovi.posi = n_pos.ovi / n_ovi,
    .groups = "drop"
  ) %>% 
  mutate(
    iso_week_string = paste0(year, "-W", sprintf("%02d", week), "-1"), 
    fake_date = ISOweek2date(iso_week_string)
  ) 

table(df_aggregated$winter, df_aggregated$mun)

# Remove years with less samplings and anomalous observation (Parma week 43, 2022-2023)
df_aggregated <- df_aggregated %>%
  filter(!winter %in% c("2010-2011", "2025-2026")) %>% 
  filter(!(mun == "Parma" & winter == "2022-2023" & week == 43)) %>% 
  filter(!(mun == "Bologna" & winter == "2011-2012")) %>% 
  filter(!(mun == "Bologna" & winter == "2014-2015"))

table(df_aggregated$winter, df_aggregated$mun)

df_aggregated %>% 
  ggplot(aes(x=fake_date, y=prop.ovi.posi, col =mun))+
  geom_line()+
  theme_minimal(base_size = 16) +
  facet_wrap(~winter, scales="free_x")+
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

str(df_aggregated)

df_aggregated %>% 
  ggplot(aes(x=fake_date, y=prop.ovi.posi, col =mun))+
  geom_line()+
  theme_minimal(base_size = 16) +
  facet_wrap(~winter, scales="free_x")+
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Calculate mean and SE by winter season and municipality
season_summary <- df_aggregated %>%
  group_by(winter, mun) %>%
  summarise(
    n_weeks = n(),
    mean_pos = mean(prop.ovi.posi, na.rm = TRUE),
    se_pos = sd(prop.ovi.posi, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  # Extract start year for ordering
  mutate(
    year_start = as.numeric(substr(winter, 1, 4)),
    winter_label = winter
  ) %>%
  arrange(year_start)

# Create color palette for municipalities
n_mun <- n_distinct(season_summary$mun)
mun_colors <- setNames(
  scales::hue_pal()(n_mun),
  sort(unique(season_summary$mun))
)

# Plot with lines and error bars
p_season <- ggplot(season_summary, aes(x = winter_label, y = mean_pos, color = mun, group = mun)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_pos - 1.96 * se_pos, 
        ymax = mean_pos + 1.96 * se_pos),
    width = 0.3,
    alpha = 0.6
  ) +
  scale_color_manual(values = mun_colors, name = "Municipality") +
  labs(
    title = "Winter Ovitrap Positivity by Season and Municipality",
    subtitle = "Mean ± 95% CI across all weeks within each winter season",
    x = "Winter Season",
    y = "Proportion Positive Ovitraps",
    caption = paste0("n = ", nrow(df_aggregated), " trap-weeks across ", 
                     n_distinct(df_aggregated$winter), " winter seasons")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )
p_season

p_season_facet <- ggplot(season_summary, aes(x = winter_label, y = mean_pos, group=mun)) +
  geom_ribbon(
    aes(ymin = mean_pos - 1.96 * se_pos, 
        ymax = mean_pos + 1.96 * se_pos),
    fill = "steelblue",
    alpha = 0.3
  ) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~mun, ncol = 2) +
  labs(
    title = "Winter Ovitrap Positivity by Municipality (2011-2024)",
    x = "Winter Season",
    y = "Proportion Positive (mean ± 95% CI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

p_season_facet

ggsave(
  paste0(outdir, "fig_seasonal_positivity_faceted.png"), 
  p_season_facet, 
  width = 12, 
  height = 10, 
  dpi = 300
)

# ALTERNATIVE: Heatmap style
p_heatmap <- ggplot(season_summary, aes(x = winter_label, y = mun, fill = mean_pos)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = round(mean_pos, 2)), 
    size = 3, 
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_viridis_c(
    name = "Mean\nPositivity",
    option = "plasma"
  ) +
  labs(
    title = "Winter Ovitrap Positivity Heatmap",
    subtitle = "Mean proportion positive by municipality and season",
    x = "Winter Season",
    y = "Municipality"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
p_heatmap

# Calculate overall mean and SE by winter season
season_overall <- df_aggregated %>%
  group_by(winter) %>%
  summarise(
    n_weeks = n(),
    n_municipalities = n_distinct(mun),
    mean_pos = mean(prop.ovi.posi, na.rm = TRUE),
    se_pos = sd(prop.ovi.posi, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    year_start = as.numeric(substr(winter, 1, 4)),
    winter_label = winter
  ) %>%
  arrange(year_start)

# Plot overall trend
p_season_overall <- ggplot(season_overall, aes(x = winter_label, y = mean_pos)) +
  geom_ribbon(
    aes(ymin = mean_pos - 1.96 * se_pos,
        ymax = mean_pos + 1.96 * se_pos),
    fill = "steelblue",
    alpha = 0.5,
    group=1
  ) +
  geom_errorbar(
    aes(ymin = mean_pos - 1.96 * se_pos, 
        ymax = mean_pos + 1.96 * se_pos),
    width = 0.4,
    color = "steelblue",
    linewidth = 1,
    alpha = 0.7
  ) +
  geom_line(color = "steelblue", linewidth = 1.5, group = 1) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(
    aes(label = round(mean_pos, 2)),
    vjust = -1.5,
    size = 3,
    fontface = "bold"
  ) +
  # scale_y_continuous(limits = c(0, max(season_overall$mean_pos) * 1.15)) +
  labs(
    title = "Overall Winter Ovitrap Positivity by Season (2011-2024)",
    subtitle = "Mean across all municipalities and weeks within each winter season",
    x = "Winter Season",
    y = "Proportion Positive Ovitraps (mean ± 95% CI)",
    caption = paste0(
      "Total: ", sum(season_overall$n_weeks), " trap-weeks | ",
      n_distinct(df_aggregated$mun), " municipalities | ",
      nrow(season_overall), " winter seasons"
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
p_season_overall

season_overall
summary(season_overall)
m <- lm(mean_pos~year_start, data=season_overall)
summary(m) # apparently no significant linear trend over the years


# 2. LOAD ENVIRONMENTAL COVARIATES ----------------------------------------

# Define spatial extent
bbox <- st_bbox(locations)
extent_buffer <- ext(c(bbox[1], bbox[3], bbox[2], bbox[4])) + 0.2

# Load climate rasters
climate_path <- "/home/dared/Documents/PoD/marieCurie2022/ZanZemap/ML/Download_Process_ClimData/2025/data/weekly/"
tas <- rast(paste0(climate_path, "t2m_min_201001-202512_C.nc"))
tp <- rast(paste0(climate_path, "tp_weekly_cumulative_201001-202512.nc"))

# Crop to study area
tas <- crop(tas, extent_buffer)
tp <- crop(tp, extent_buffer)

# Set temporal reference
date_seq <- seq.Date(as.Date("2010-01-07"), length.out = nlyr(tas), by = "week")
week_names <- paste0("w_", year(date_seq), "_", week(date_seq))

time(tas) <- date_seq
time(tp) <- date_seq
names(tas) <- week_names
names(tp) <- week_names

# Fill coastal NAs using focal mean
focal_window <- matrix(1, 5, 5)
tas <- focal(tas, focal_window, mean, na.rm = TRUE, NAonly = TRUE, pad = TRUE)
tp <- focal(tp, focal_window, mean, na.rm = TRUE, NAonly = TRUE, pad = TRUE)

# 3. EXTRACT CLIMATE DATA AT LOCATIONS ------------------------------------

loc_coords <- bind_cols(
  Location = locations$mun, 
  st_coordinates(locations)
) %>% 
  rename(Longitude = X, Latitude = Y)

loc_vect <- vect(
  loc_coords, 
  geom = c("Longitude", "Latitude"), 
  crs = "EPSG:4326"
)

# Extract temperature
tas_df <- extract(tas, loc_vect)  
names(tas_df) <- c("ID", as.character(time(tas)))  
tas_df <- tas_df %>% 
  mutate(Location = loc_coords$Location) %>%
  select(-ID) %>%
  pivot_longer(-Location, names_to = "date", values_to = "tas") %>%
  mutate(
    date = as.Date(date),
    year = year(date), 
    week = week(date), 
    tid = paste(Location, year, week, sep = "_")
  )

# Extract precipitation
tp_df <- extract(tp, loc_vect)
names(tp_df) <- c("ID", as.character(time(tp)))  
tp_df <- tp_df %>% 
  mutate(Location = loc_coords$Location) %>%
  select(-ID) %>%
  pivot_longer(-Location, names_to = "date", values_to = "tp") %>%
  mutate(
    date = as.Date(date),
    year = year(date), 
    week = week(date), 
    tid = paste(Location, year, week, sep = "_")
  )

# Combine climate data
climate_df <- tas_df %>%
  left_join(tp_df %>% select(tid, tp), by = "tid")

# 4. CREATE ANALYSIS DATASET ----------------------------------------------

# Join observations with climate
df_combined <- df_aggregated %>%
  rename(loc = mun) %>% 
  mutate(tid = paste(loc, year, week, sep = "_")) %>% 
  left_join(climate_df %>% select(tid, tas, tp, date), by = "tid")

# Create lagged predictors and temporal variables
df_model <- df_combined %>%
  group_by(loc) %>%
  arrange(year, week, .by_group = TRUE) %>%
  mutate(
    tas_lag1 = lag(tas, 1),
    tas_lag2 = lag(tas, 2),
    tas_lag3 = lag(tas, 3),
    tp_lag1 = lag(tp, 1),
    tp_lag2 = lag(tp, 2),
    tp_lag3 = lag(tp, 3)
  ) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(
    winter_week = ifelse(week >= 42, week - 42, week + 11),
    year_fac = factor(year),
    year_numeric = as.numeric(as.character(year)),
    year_centered = year_numeric - mean(year_numeric),
    period = ifelse(month %in% 10:12, "Autumn", "Spring"),
    period = factor(period, levels = c("Autumn", "Spring")),
    loc = factor(loc)
  )
summary(df_model)
write.csv(df_model, "outputs/er_gam_calibration_dataset.csv")
# 5. EXPLORATORY DATA ANALYSIS --------------------------------------------

# Summary statistics
cat("\n=== Dataset Summary ===\n")
cat("Sample size:", nrow(df_model), "trap-weeks\n")
cat("Time span:", paste(range(df_model$year), collapse = "-"), "\n")
cat("Municipalities:", n_distinct(df_model$loc), "\n")
cat("Overall positivity:", 
    round(mean(df_model$n_pos.ovi > 0) * 100, 1), "%\n\n")

# Correlation analysis
climate_cors <- df_model %>%
  select(tas, tas_lag1, tas_lag2, tas_lag3, 
         tp, tp_lag1, tp_lag2, tp_lag3) %>%
  drop_na() %>%
  cor()
climate_cors
# png(paste0(outdir, "climate_correlations.png"), 
#     width = 10, height = 10, units = "in", res = 300)
corrplot(
  climate_cors, 
  method = "color", 
  type = "upper", 
  addCoef.col = "black", 
  number.cex = 0.8,
  tl.col = "black", 
  tl.srt = 45, 
  tl.cex = 1.2,
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
  title = "Climate Variable Correlations",
  mar = c(0, 0, 2, 0)
)
# dev.off()

#I keep only tas_lag3 because is the least correlated to the other temperatures (<0.7) 
# no photoperiod used because the region is the same, we use winter week as a proxy for seasonality
# 6. MODEL SELECTION ------------------------------------------------------

# less correlated temperature and annual smooth to detect non linear trends
m_full <- gam(
  cbind(n_pos.ovi, n_ovi - n_pos.ovi) ~ 
    s(winter_week, bs = "cc", k = 24) +
    s(tas_lag3, k = 10) +
    tp + tp_lag1 + tp_lag2 + tp_lag3 +
    period +
    s(year, k = 10) +
    s(loc, bs = "re"),
  family = binomial,
  method = "REML",
  data = df_model,
  knots = list(winter_week = c(1, 31))
)

# Check overdispersion
overdisp <- sum(residuals(m_full, type = "pearson")^2) / m_full$df.residual
overdisp

# Refit with quasibinomial if overdispersed
if (overdisp > 1.5) {
  m_full <- update(m_full, family = quasibinomial)
}
summary(m_full)
concurvity(m_full)

# 2. TEMPERATURE LAG SELECTION  
no_m_tas_lag3  <- update(m_full, . ~ . -  s(tas_lag3, k = 10))
anova(m_full, no_m_tas_lag3)
summary(no_m_tas_lag3)
 
# 4. FINAL MODEL DIAGNOSTICS  
m_final <- no_m_tas_lag3
summary(m_final)
concurvity(m_final)
gam.check(m_final)

# Temporal autocorrelation
resid_df <- data.frame(
  residuals = residuals(m_final, type = "response"),
  loc = df_model$loc,
  date = df_model$date
) %>%
  arrange(loc, date)

loc_acf <- resid_df %>%
  group_by(loc) %>%
  group_map(~ acf(.x$residuals, plot = FALSE, lag.max = 10))

avg_acf <- Reduce("+", lapply(loc_acf, function(x) x$acf)) / length(loc_acf)

# png(paste0(outdir, "acf_by_location.png"), 
# width = 8, height = 6, units = "in", res = 300)
plot(0:10, avg_acf, type = "h", xlab = "Lag", ylab = "ACF", 
     main = "Average ACF by Location", lwd = 2)
abline(h = 0)
abline(h = c(-1.96, 1.96) / sqrt(nrow(resid_df) / 10), lty = 2, col = "blue")
# dev.off()
 

# 8. VISUALIZATION: PARTIAL EFFECTS --------------------------------------
# Panel A: Seasonal pattern  
# p_seasonal_draw <- draw(m_final, select = "s(winter_week)")
# p_seasonal <- p_seasonal_draw[[1]] +  # Extract first (and only) plot
#   labs(
#     title = "A. Seasonality",
#     x = "Winter week (1 = early Oct, 31 = late May)",
#     y = "Effect (log-odds)"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(face = "bold"))

# Create prediction grid at reference values for other predictors

# Create prediction grid
seasonal_pred <- expand.grid(
  winter_week = 1:31,
  tp = median(df_model$tp, na.rm = TRUE),
  tp_lag1 = median(df_model$tp_lag1, na.rm = TRUE),
  tp_lag2 = median(df_model$tp_lag2, na.rm = TRUE),
  tp_lag3 = median(df_model$tp_lag3, na.rm = TRUE),
  period = "Autumn",
  year = median(df_model$year),
  loc = levels(df_model$loc)[1]
)

# Predict on LINK scale (log-odds) for proper SEs
preds_link <- predict(m_final, newdata = seasonal_pred, 
                      type = "link", se.fit = TRUE,
                      exclude = c("s(year)", "s(loc)"))

# Transform to probability scale with proper CIs
seasonal_pred <- seasonal_pred %>%
  mutate(
    fit_link = preds_link$fit,
    se_link = preds_link$se.fit,
    prob = plogis(fit_link),
    lower = plogis(fit_link - 1.96 * se_link),
    upper = plogis(fit_link + 1.96 * se_link)
  )

# Plot with month labels and confidence bands
p_seasonal_response <- ggplot(seasonal_pred, aes(x = winter_week, y = prob)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = "grey70", alpha = 0.4) +  # Slightly darker for visibility
  geom_line(linewidth = 1.2, color = "black") +
  scale_x_continuous(
    breaks = c(2, 6, 11, 15, 20, 24, 28),
    labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"),
    name = NULL
  ) +
  scale_y_continuous(
    # labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "A. Seasonal Pattern",
    y = "Predicted probability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 
# # Panel B: Temperature  
# temp_seq <- seq(
#   min(df_model$tas_lag3, na.rm = TRUE), 
#   max(df_model$tas_lag3, na.rm = TRUE), 
#   length = 200
# )
# 
# pred_temp <- expand.grid(
#   winter_week = median(df_model$winter_week),
#   tas_lag3 = temp_seq,
#   tp = median(df_model$tp, na.rm = TRUE),
#   tp_lag1 = median(df_model$tp_lag1, na.rm = TRUE),
#   tp_lag2 = median(df_model$tp_lag2, na.rm = TRUE),
#   tp_lag3 = median(df_model$tp_lag3, na.rm = TRUE),
#   period = "Autumn",
#   year = median(df_model$year),
#   loc = levels(df_model$loc)[1]
# )
# 
# pred_temp$prob <- predict(
#   m_final, 
#   newdata = pred_temp, 
#   type = "response",
#   exclude = c("s(loc)", "s(year)")
# )
# 
# p_temperature <- ggplot(pred_temp, aes(x = tas_lag3, y = prob)) +
#   geom_line(linewidth = 1.2, color =  "#0072B2") +
#   geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", alpha = 0.5) +
#   # annotate("text", x = 6, y = max(pred_temp$prob) * 0.95, 
#   #          label = "~5°C", size = 3, color = "gray30") +
#   labs(
#     # title = "(B) Temperature Effect (lag 3 weeks, p=0.10)",
#     x = "Average minimum temperature 3 weeks ago (°C)",
#     y = "Predicted probability"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(plot.title = element_text(face = "bold"))

# Panel C: Precipitation (ALL 4 lags together)
# Create predictions for each lag
precip_seq <- seq(0, max(c(df_model$tp, df_model$tp_lag1, 
                           df_model$tp_lag2, df_model$tp_lag3), 
                         na.rm = TRUE), 
                  length = 200)

# Base data for predictions
base_pred <- data.frame(
  winter_week = median(df_model$winter_week),
  tas_lag3 = median(df_model$tas_lag3, na.rm = TRUE),
  period = "Autumn",
  year = median(df_model$year),
  loc = levels(df_model$loc)[1]
)

# Vary tp (lag 0)
pred_tp <- base_pred %>%
  slice(rep(1, length(precip_seq))) %>%
  mutate(
    tp = precip_seq,
    tp_lag1 = median(df_model$tp_lag1, na.rm = TRUE),
    tp_lag2 = median(df_model$tp_lag2, na.rm = TRUE),
    tp_lag3 = median(df_model$tp_lag3, na.rm = TRUE),
    lag = "Lag 0 (current week)",
    precip = tp
  )

# Vary tp_lag1
pred_tp_lag1 <- base_pred %>%
  slice(rep(1, length(precip_seq))) %>%
  mutate(
    tp = median(df_model$tp, na.rm = TRUE),
    tp_lag1 = precip_seq,
    tp_lag2 = median(df_model$tp_lag2, na.rm = TRUE),
    tp_lag3 = median(df_model$tp_lag3, na.rm = TRUE),
    lag = "Lag 1 (1 week prior)",
    precip = tp_lag1
  )

# Vary tp_lag2
pred_tp_lag2 <- base_pred %>%
  slice(rep(1, length(precip_seq))) %>%
  mutate(
    tp = median(df_model$tp, na.rm = TRUE),
    tp_lag1 = median(df_model$tp_lag1, na.rm = TRUE),
    tp_lag2 = precip_seq,
    tp_lag3 = median(df_model$tp_lag3, na.rm = TRUE),
    lag = "Lag 2 (2 weeks prior)",
    precip = tp_lag2
  )

# Vary tp_lag3
pred_tp_lag3 <- base_pred %>%
  slice(rep(1, length(precip_seq))) %>%
  mutate(
    tp = median(df_model$tp, na.rm = TRUE),
    tp_lag1 = median(df_model$tp_lag1, na.rm = TRUE),
    tp_lag2 = median(df_model$tp_lag2, na.rm = TRUE),
    tp_lag3 = precip_seq,
    lag = "Lag 3 (3 weeks prior)",
    precip = tp_lag3
  )

# Combine all precipitation predictions
pred_precip_all <- bind_rows(pred_tp, pred_tp_lag1, pred_tp_lag2, pred_tp_lag3)

# Get predictions
pred_precip_all$prob <- predict(
  m_final, 
  newdata = pred_precip_all, 
  type = "response",
  exclude = c("s(loc)", "s(year)")
)

# Order factor for plotting
pred_precip_all$lag <- factor(
  pred_precip_all$lag,
  levels = c("Lag 0 (current week)", "Lag 1 (1 week prior)", 
             "Lag 2 (2 weeks prior)", "Lag 3 (3 weeks prior)")
)

p_precipitation <- ggplot(pred_precip_all, aes(x = precip, y = prob, color = lag)) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis_d(
    option = "mako",
    name = "Precipitation lag",
    begin = 0.1,
    end = 0.9,
    direction=-1
  ) +
  labs(
    title = "B. Precipitation Effects",
    x = "Cumulative weekly precipitation (mm)",
    y = "Predicted probability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.width = unit(1.5, "cm")
  )

# Panel D: Period effect  
period_data <- data.frame(
  period = c("Autumn", "Spring"),
  estimate = c(0, coef(m_final)["periodSpring"]),
  se = c(0, summary(m_final)$p.table["periodSpring", "Std. Error"])
) %>%
  mutate(
    logit = coef(m_final)["(Intercept)"] + estimate,
    prob = plogis(logit),
    lower = plogis(logit - 1.96 * se),
    upper = plogis(logit + 1.96 * se)
  )

p_period <- ggplot(period_data, aes(x = period, y = prob, fill = period)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1) +
  geom_text(
    aes(label = paste0(round(prob * 100, 1), "%")), 
    vjust = -0.5, 
    fontface = "bold", 
    size = 5
  ) +
  scale_fill_manual(values = c("Autumn" = "#E69F00", "Spring" = "#56B4E9")) +
  labs(
    title = "C. Autumn vs Spring",
    x = NULL,
    y = "Predicted probability"
  ) +
  ylim(0, max(period_data$upper) * 1.15) +
  theme_minimal(base_size = 12) +
  theme(
    # plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# Combine panels  
fig1 <- (p_seasonal_response |p_precipitation) / ( p_period) + # p_temperature
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    # title = "Figure 1. Partial Effects on Winter Ovitrap Positivity",
    # caption = "Predictions at median covariate values, excluding municipality and year random effects.\nShaded regions/error bars = 95% confidence intervals.",
    theme = theme(
      plot.title = element_text(size = 15),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )
fig1

ggsave(
  paste0(outdir, "erMod_partial_effects_final.png"), 
  fig1, 
  width = 14, 
  height = 10, 
  dpi = 300
)

 # 9. VISUALIZATION: TEMPORAL TREND ----------------------------------------

# Get fitted values for all observations
df_model$fitted <- predict(m_final, type = "response")

# Summarize by year AND period
annual_by_period <- df_model %>%
  group_by(year, period) %>%
  summarise(
    n_weeks = n(),
    observed_mean = mean(prop.ovi.posi),
    observed_se = sd(prop.ovi.posi) / sqrt(n()),
    fitted_mean = mean(fitted),
    fitted_se = sd(fitted) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    year_numeric = as.numeric(as.character(year)),
    residual = observed_mean - fitted_mean,
    pct_error = (observed_mean - fitted_mean) / observed_mean * 100
  )

# Print summary
print(annual_by_period, n = Inf)

# Check overall fit by period
annual_by_period %>%
  group_by(period) %>%
  summarise(
    mean_observed = mean(observed_mean, na.rm = TRUE),
    mean_fitted = mean(fitted_mean, na.rm = TRUE),
    mean_abs_error = mean(abs(residual), na.rm = TRUE),
    mean_pct_error = mean(abs(pct_error), na.rm = TRUE)
  )


# PLOT: Separate panels for Autumn and Spring

# Autumn only
p_autumn <- annual_by_period %>%
  filter(period == "Autumn") %>%
  ggplot(aes(x = year_numeric)) +
  # Fitted
  geom_ribbon(
    aes(ymin = fitted_mean - 1.96*fitted_se,
        ymax = fitted_mean + 1.96*fitted_se),
    fill = "#E69F00",
    alpha = 0.3
  ) +
  geom_line(aes(y = fitted_mean), color = "#E69F00", linewidth = 1.5) +
  geom_point(aes(y = fitted_mean), color = "#E69F00", size = 3) +
  # Observed
  geom_errorbar(
    aes(ymin = observed_mean - 1.96*observed_se,
        ymax = observed_mean + 1.96*observed_se),
    width = 0.3,
    color = "black"
  ) +
  geom_point(aes(y = observed_mean), size = 4, color = "black") +
  scale_x_continuous(breaks = seq(2012, 2024, by = 2)) +
  ylim(0, 1) +
  labs(
    title = "Autumn (Oct-Dec)",
    x = "Year",
    y = "Proportion positive"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#E69F00"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Spring only
p_spring <- annual_by_period %>%
  filter(period == "Spring") %>%
  ggplot(aes(x = year_numeric)) +
  # Fitted
  geom_ribbon(
    aes(ymin = fitted_mean - 1.96*fitted_se,
        ymax = fitted_mean + 1.96*fitted_se),
    fill = "#56B4E9",
    alpha = 0.3
  ) +
  geom_line(aes(y = fitted_mean), color = "#56B4E9", linewidth = 1.5) +
  geom_point(aes(y = fitted_mean), color = "#56B4E9", size = 3) +
  # Observed
  geom_errorbar(
    aes(ymin = observed_mean - 1.96*observed_se,
        ymax = observed_mean + 1.96*observed_se),
    width = 0.3,
    color = "black"
  ) +
  geom_point(aes(y = observed_mean), size = 4, color = "black") +
  # # Highlight 2021
  # geom_point(
  #   data = annual_by_period %>% filter(period == "Spring", year == 2021),
  #   aes(y = observed_mean),
  #   size = 6,
  #   shape = 21,
  #   fill = NA,
  #   color = "red",
  #   stroke = 1.5
  # ) +
  scale_x_continuous(breaks = seq(2012, 2024, by = 2)) +
  ylim(0, 0.4) +
  labs(
    title = "Spring (Jan-May)",
    x = "Year",
    y = "Proportion positive"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#56B4E9"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine
fig_by_period <- p_autumn / p_spring 
  # plot_annotation(
  #   title = "Model Fit by Period: Observed vs Fitted Values",
  #   subtitle = "Colored lines = model fitted means ± 95% CI | Black points = observed means ± 95% CI",
  #   caption = paste0("Model R² = ", round(summary(m_final)$r.sq, 3))
  # )

fig_by_period

ggsave(
  paste0(outdir, "fig_diagnostic_by_period.png"),
  fig_by_period,
  width = 10,
  height = 10,
  dpi = 300
)

# COMBINED PLOT: Both periods with clear separation
#fix error bar for the observations + dashed trend for spring only

p_combined <- ggplot(annual_by_period %>% filter(year >= 2012), 
                     aes(x = year_numeric)) +
  
  # ── Autumn ────────────────────────────────────────────────────
  geom_ribbon(
    data = . %>% filter(period == "Autumn"),
    aes(ymin = fitted_mean - 1.96 * fitted_se,
        ymax = fitted_mean + 1.96 * fitted_se),
    fill = "#E69F00", alpha = 0.2
  ) +
  geom_line(
    data = . %>% filter(period == "Autumn"),
    aes(y = fitted_mean, color = "Autumn"),
    linewidth = 1.5
  ) +
  geom_errorbar(
    data = . %>% filter(period == "Autumn"),
    aes(ymin = observed_mean - 1.96 * observed_se,
        ymax = observed_mean + 1.96 * observed_se),
    width = 0.3, color = "#E69F00", linewidth = 0.6
  ) +
  geom_point(
    data = . %>% filter(period == "Autumn"),
    aes(y = observed_mean, fill = "Autumn"),
    size = 4, shape = 21, color = "black"
  ) +
  
  # ── Spring ────────────────────────────────────────────────────
  geom_ribbon(
    data = . %>% filter(period == "Spring"),
    aes(ymin = fitted_mean - 1.96 * fitted_se,
        ymax = fitted_mean + 1.96 * fitted_se),
    fill = "#56B4E9", alpha = 0.2
  ) +
  geom_line(
    data = . %>% filter(period == "Spring"),
    aes(y = fitted_mean, color = "Spring"),
    linewidth = 1.5
  ) +
  geom_errorbar(
    data = . %>% filter(period == "Spring"),
    aes(ymin = observed_mean - 1.96 * observed_se,
        ymax = observed_mean + 1.96 * observed_se),
    width = 0.3, color = "#56B4E9", linewidth = 0.6
  ) +
  geom_point(
    data = . %>% filter(period == "Spring"),
    aes(y = observed_mean, fill = "Spring"),
    size = 4, shape = 21, color = "black"
  ) +
  
  # ── Significant linear trend for Spring only ──────────────────
  geom_smooth(
    data = . %>% filter(period == "Spring"),
    aes(y = observed_mean),
    method = "lm", formula = y ~ x,
    se = FALSE,
    color = "#1a6fa8", linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  # ── Scales & labels ───────────────────────────────────────────
  scale_color_manual(
    values = c("Autumn" = "#E69F00", "Spring" = "#56B4E9"),
    name = "Period"
  ) +
  scale_fill_manual(
    values = c("Autumn" = "#E69F00", "Spring" = "#56B4E9"),
    name = "Period"
  ) +
  scale_x_continuous(breaks = seq(2012, 2024, by = 2)) +
  labs(
    # title = "Autumn–Spring Asymmetry in Winter Ovitrap Positivity (2012–2024)",
    subtitle = "Spring trend: slope = 0.007/year, p = 0.010 | Autumn trend non-significant (p = 0.098)",
    x = "Year",
    y = "Proportion positive ovitraps",
    caption = paste0(
      "Points = observed annual means ± 95% CI (1.96 × SE) | Lines = model fitted values ± 95% CI\n",
      "Dashed line = significant linear trend (Spring only) | ",
      "Period effect: β = −1.227, p = 1.92 × 10⁻¹⁰ | ",
      "Model R² = ", round(summary(m_final)$r.sq, 3)
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption  = element_text(hjust = 0, size = 9),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

p_combined


ggsave(
  paste0(outdir, "fig6_autumn_spring_combined.png"),
  p_combined,
  width = 12,
  height = 7,
  dpi = 300
)

# NUMERICAL SUMMARY

# Overall fit statistics
fit_summary <- annual_by_period %>%
  group_by(period) %>%
  summarise(
    n_years = n(),
    mean_observed = mean(observed_mean, na.rm = TRUE),
    mean_fitted = mean(fitted_mean, na.rm = TRUE),
    correlation = cor(observed_mean, fitted_mean, use = "complete.obs"),
    mean_abs_error = mean(abs(residual), na.rm = TRUE),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    .groups = "drop"
  )

print("=== FIT SUMMARY BY PERIOD ===")
print(fit_summary)

# Years with worst fit
worst_fit <- annual_by_period %>%
  arrange(desc(abs(residual))) %>%
  select(year, period, observed_mean, fitted_mean, residual, pct_error)  

print("=== WORST FIT YEARS ===")
print(worst_fit)


# 10 FORMAL TREND TEST: Linear model on annual means by period ----

# Test Autumn trend
autumn_data <- annual_by_period %>% 
  filter(period == "Autumn") #, year >= 2012)  # Exclude 2011 (only 16 obs)

m_autumn_trend <- lm(observed_mean ~ year_numeric, data = autumn_data)
summary(m_autumn_trend)

# Test Spring trend  
spring_data <- annual_by_period %>% 
  filter(period == "Spring")

m_spring_trend <- lm(observed_mean ~ year_numeric, data = spring_data)
summary(m_spring_trend)

# Combined test (interaction)
both_periods <- annual_by_period %>% 
  filter(!(period == "Autumn" & year == 2011))

m_interaction <- lm(observed_mean ~ year_numeric * period, data = both_periods)
summary(m_interaction)

# Summary table
trend_summary <- data.frame(
  Period = c("Autumn", "Spring", "Interaction"),
  Slope = c(
    coef(m_autumn_trend)["year_numeric"],
    coef(m_spring_trend)["year_numeric"],
    coef(m_interaction)["year_numeric:periodSpring"]
  ),
  P_value = c(
    summary(m_autumn_trend)$coefficients["year_numeric", "Pr(>|t|)"],
    summary(m_spring_trend)$coefficients["year_numeric", "Pr(>|t|)"],
    summary(m_interaction)$coefficients["year_numeric:periodSpring", "Pr(>|t|)"]
  ),
  R2 = c(
    summary(m_autumn_trend)$r.squared,
    summary(m_spring_trend)$r.squared,
    summary(m_interaction)$r.squared
  )
)

print(trend_summary)

# Save
write_csv(trend_summary, paste0(outdir, "temporal_trend_tests.csv"))

# 10. SAVE SUMMARY STATISTICS ---------------------------------------------

# Annual summary by period
annual_by_period <- df_model %>%
  mutate(fitted = predict(m_final, type = "response")) %>%
  group_by(year, period) %>%
  summarise(
    n_weeks = n(),
    observed_mean = mean(prop.ovi.posi),
    observed_se = sd(prop.ovi.posi) / sqrt(n()),
    predicted_mean = mean(fitted),
    predicted_se = sd(fitted) / sqrt(n()),
    .groups = "drop"
  )

write_csv(annual_by_period, paste0(outdir, "annual_summary_by_period.csv"))

# Overall annual summary
annual_overall <- df_model %>%
  mutate(fitted = predict(m_final, type = "response")) %>%
  group_by(year) %>%
  summarise(
    n_weeks = n(),
    observed_mean = mean(prop.ovi.posi),
    observed_se = sd(prop.ovi.posi) / sqrt(n()),
    predicted_mean = mean(fitted),
    predicted_se = sd(fitted) / sqrt(n()),
    .groups = "drop"
  )

write_csv(annual_overall, paste0(outdir, "annual_summary_overall.csv"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", outdir, "\n\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================