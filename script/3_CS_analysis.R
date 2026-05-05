setwd("GitHub/albo_winter_monitoring/")
library(lubridate)
library(tidyverse)
library(patchwork)
outdir <- "outputs/"

df <- read_csv("data/mosquito_alert_adults_europe_2020_2025.csv") %>% 
  select(creation_date, lon, lat, country, nuts_2, nuts_3, movelab_annotation_euro.class_name) %>% 
  rename(date=creation_date, species= movelab_annotation_euro.class_name) %>%
  filter(species =="Aedes albopictus") %>% 
  mutate(week=week(date),
         year=year(date))
df

str(df)

# check "sampling effort
df %>%
  count(year) %>%
  ggplot(aes(x = year, y = n)) +
  geom_col() +
  labs(title = "Reports per year — check for effort trend")

#exclude 2018 and 2025
# Custom month order: Oct→Dec first, then Jan→Sep
month_order  <- c(10:12, 1:9)
month_labels <- month.abb[month_order]

df %>%
  filter(year > 2019, year < 2025) %>%
  mutate(
    lat_bin = floor(lat),
    month   = month(date)
  ) %>%
  group_by(year, month, lat_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year, lat_bin) %>%
  mutate(n_rel = n / sum(n)) %>%
  group_by(month, lat_bin) %>%
  summarise(
    n_rel  = mean(n_rel),
    n_mean = mean(n),          # average observations per year for this cell
    .groups = "drop"
  ) %>%
  mutate(month_plot = match(month, month_order)) %>%
  ggplot(aes(x = month_plot, y = lat_bin, size = n_rel, alpha = n_mean)) +
  geom_point(colour = "firebrick") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_size_continuous(name = "Relative frequency", range = c(1, 12)) +
  scale_alpha_continuous(name = "Mean observations/year", range = c(0.1, 0.9)) +
  labs(
    x     = "Month",
    y     = "Latitude (1° bin)"
  ) +
  theme_minimal()


# Shared aggregated data
df_agg <- df %>%
  filter(year > 2019, year < 2025) %>%
  mutate(lat_bin = floor(lat),
         month   = month(date)) %>%
  group_by(year, month, lat_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year, lat_bin) %>%
  mutate(n_rel = n / sum(n)) %>%
  group_by(month, lat_bin) %>%
  summarise(
    n_rel  = mean(n_rel),
    n_mean = mean(n),
    .groups = "drop"
  ) %>%
  mutate(month_plot = match(month, month_order))

# Main plot
p1 <- df_agg %>%
  ggplot(aes(x = month_plot, y = lat_bin, size = n_rel, alpha = n_mean)) +
  geom_point(colour = "firebrick") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_size_continuous(name = "Relative frequency", range = c(1, 12)) +
  scale_alpha_continuous(name = "Mean obs/year", range = c(0.1, 0.9)) +
  labs(x = "Month", y = "Latitude (1° bin)"
       ) +
  theme_minimal()
p1

# Winter subplot disaggregated by country
p2 <- df %>%
  filter(year > 2019, year < 2025,
         month(date) %in% c(12, 1, 2, 3)) %>%
  mutate(lat_bin = floor(lat),
         month   = month(date)) %>%
  group_by(year, month, lat_bin, country) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year, lat_bin, country) %>%
  mutate(n_rel = n / sum(n)) %>%
  group_by(month, lat_bin, country) %>%
  summarise(
    n_rel  = mean(n_rel),
    n_mean = mean(n),
    .groups = "drop"
  ) %>%
  mutate(month_plot = match(month, month_order)) %>%
  ggplot(aes(x = month_plot, y = lat_bin, size = n_rel, alpha = n_mean)) +
  geom_point(colour = "darkorange") +
  scale_x_continuous(breaks = match(c(12,1,2,3), month_order),
                     labels = c("Dec","Jan","Feb","Mar")) +
  scale_size_continuous(name = "Relative frequency", range = c(1, 10)) +
  scale_alpha_continuous(name = "Mean obs/year", range = c(0.1, 0.9)) +
  geom_vline(xintercept = match(c(12,1,2,3), month_order) + 0.5,
             linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  facet_wrap(~country) +
  labs(x = "Month", y = "Latitude (1° bin)"
        ) +
  theme_minimal()+
  guides(
    alpha = guide_legend(order = 1),
    size  = guide_legend(order = 2)
  )
p2

# Combine
p<- p1  / p2  +
  plot_annotation(tag_levels = "A") 
p
ggsave(paste0(outdir, "cs_overview.png"), p, 
       width = 16, height = 10, dpi = 320)
