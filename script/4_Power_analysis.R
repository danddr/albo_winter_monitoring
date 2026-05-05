# TRAP EFFORT OPTIMIZATION ANALYSIS
# Emilia-Romagna Winter Surveillance (2011-2024)
# Daniele Da Re
# April 2026
library(tidyverse)
library(scales)

# Setup
outdir <- "outputs/emilia_romagna/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load data
df_model <- read.csv("outputs/emilia_romagna/er_gam_calibration_dataset.csv")

# ============================================================================
# 1. EMPIRICAL DETECTION RATES
# ============================================================================

# Summarize trap effort and detection rates by location-year-period
effort_summary <- df_model %>%
  group_by(period, year, loc) %>%
  summarise(
    n_weeks = n(),
    mean_traps = mean(n_ovi),
    mean_positivity = mean(prop.ovi.posi),
    .groups = "drop"
  )

# Overall summary by period
period_summary <- effort_summary %>%
  group_by(period) %>%
  summarise(
    n_observations = n(),
    mean_traps = mean(mean_traps),
    median_traps = median(mean_traps),
    mean_positivity = mean(mean_positivity),
    sd_positivity = sd(mean_positivity),
    .groups = "drop"
  )

print("=== EMPIRICAL DETECTION RATES BY PERIOD ===")
print(period_summary)

write_csv(period_summary, paste0(outdir, "trap_effort_summary_by_period.csv"))

# ============================================================================
# 2. REQUIRED TRAP NUMBERS FOR DETECTION
# ============================================================================

# Vectorized function: Calculate traps needed for desired detection probability
calculate_required_traps <- function(p_positive, target_detection = 0.8) {
  result <- ifelse(
    p_positive <= 0 | p_positive >= 1,
    NA_real_,
    ceiling(log(1 - target_detection) / log(1 - p_positive))
  )
  return(result)
}

# Calculate detection probability for given trap number
calculate_detection_prob <- function(p_positive, n_traps) {
  1 - (1 - p_positive)^n_traps
}

# Calculate required traps for different detection thresholds
detection_requirements <- period_summary %>%
  mutate(
    n_traps_50pct = calculate_required_traps(mean_positivity, 0.50),
    n_traps_80pct = calculate_required_traps(mean_positivity, 0.80),
    n_traps_95pct = calculate_required_traps(mean_positivity, 0.95),
    n_traps_99pct = calculate_required_traps(mean_positivity, 0.99)
  )

print("=== REQUIRED TRAP NUMBERS ===")
print(detection_requirements)

# ============================================================================
# 3. ADEQUACY ASSESSMENT
# ============================================================================

# Assess each location-year-period combination
adequacy_assessment <- effort_summary %>%
  mutate(
    actual_detection_prob = calculate_detection_prob(mean_positivity, mean_traps),
    adequate_80pct = actual_detection_prob >= 0.80,
    adequate_95pct = actual_detection_prob >= 0.95,
    recommended_80pct = calculate_required_traps(mean_positivity, 0.80),
    recommended_95pct = calculate_required_traps(mean_positivity, 0.95)
  )

# Summary by period
adequacy_summary <- adequacy_assessment %>%
  group_by(period) %>%
  summarise(
    n_obs = n(),
    mean_actual_traps = mean(mean_traps),
    pct_adequate_80 = mean(adequate_80pct, na.rm = TRUE) * 100,
    pct_adequate_95 = mean(adequate_95pct, na.rm = TRUE) * 100,
    mean_recommended_80 = mean(recommended_80pct, na.rm = TRUE),
    mean_recommended_95 = mean(recommended_95pct, na.rm = TRUE),
    .groups = "drop"
  )

print("=== ADEQUACY ASSESSMENT ===")
print(adequacy_summary)

write_csv(adequacy_summary, paste0(outdir, "trap_adequacy_summary.csv"))

# ============================================================================
# 4. MANUSCRIPT TABLE
# ============================================================================

manuscript_table <- data.frame(
  Detection_threshold = c("50%", "80%", "95%", "99%"),
  Autumn_traps = c(
    detection_requirements$n_traps_50pct[detection_requirements$period == "Autumn"],
    detection_requirements$n_traps_80pct[detection_requirements$period == "Autumn"],
    detection_requirements$n_traps_95pct[detection_requirements$period == "Autumn"],
    detection_requirements$n_traps_99pct[detection_requirements$period == "Autumn"]
  ),
  Spring_traps = c(
    detection_requirements$n_traps_50pct[detection_requirements$period == "Spring"],
    detection_requirements$n_traps_80pct[detection_requirements$period == "Spring"],
    detection_requirements$n_traps_95pct[detection_requirements$period == "Spring"],
    detection_requirements$n_traps_99pct[detection_requirements$period == "Spring"]
  )
)

manuscript_table <- manuscript_table %>%
  mutate(Actual_ER = 10)  # Add actual deployment column

print("=== MANUSCRIPT TABLE ===")
print(manuscript_table)

write_csv(manuscript_table, paste0(outdir, "table_trap_requirements.csv"))

# ============================================================================
# 5. RAREFACTION CURVE VISUALIZATION
# ============================================================================

# Create grid for rarefaction curves
trap_numbers <- seq(1, 40, by = 1)

rarefaction_data <- expand.grid(
  period = c("Autumn", "Spring"),
  n_traps = trap_numbers,
  stringsAsFactors = FALSE
) %>%
  left_join(
    period_summary %>% select(period, mean_positivity),
    by = "period"
  ) %>%
  mutate(
    detection_prob = calculate_detection_prob(mean_positivity, n_traps)
  )

# Plot rarefaction curves
p_rarefaction <- ggplot(rarefaction_data, aes(x = n_traps, y = detection_prob, color = period)) +
  geom_line(linewidth = 1.5) +
  geom_point(size=3)+
  geom_hline(yintercept = c(0.50, 0.80, 0.95), 
             linetype = "dashed", color = "gray40", alpha = 0.5) +
  geom_vline(xintercept = 10, linetype = "dotted", color = "black", linewidth = 1) +
  annotate("text", x = 11, y = 0.1, label = "Actual\n(10 traps)", 
           hjust = 0, size = 3.5) +
  annotate("text", x = 1.5, y = 0.52, label = "50%", size = 3, color = "gray40") +
  annotate("text", x = 1.5, y = 0.82, label = "80%", size = 3, color = "gray40") +
  annotate("text", x = 1.5, y = 0.97, label = "95%", size = 3, color = "gray40") +
  scale_color_manual(
    values = c("Autumn" = "#E69F00", "Spring" = "#56B4E9"),
    name = "Period"
  ) +
  scale_x_continuous(breaks = seq(0, 40, by = 5)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    title = "Detection Probability vs Trap Effort",
    subtitle = "Emilia-Romagna winter surveillance (2011-2025)",
    x = "Number of ovitraps per municipality",
    y = "Probability of detecting ≥1 positive trap"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

print(p_rarefaction)

ggsave(paste0(outdir, "fig_rarefaction_trap_effort.png"), p_rarefaction,
       width = 8, height = 6, dpi = 300)

# ============================================================================
# 6. SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("TRAP EFFORT ANALYSIS SUMMARY\n")
cat("========================================\n\n")

cat("OBSERVED POSITIVITY:\n")
cat("Autumn: ", round(period_summary$mean_positivity[1] * 100, 1), "% (SD=", 
    round(period_summary$sd_positivity[1] * 100, 1), "%)\n", sep = "")
cat("Spring: ", round(period_summary$mean_positivity[2] * 100, 1), "% (SD=", 
    round(period_summary$sd_positivity[2] * 100, 1), "%)\n\n", sep = "")

cat("ACTUAL DEPLOYMENT: ", round(period_summary$median_traps[1]), " traps (median)\n\n", sep = "")

cat("REQUIRED FOR 80% DETECTION:\n")
cat("Autumn: ", detection_requirements$n_traps_80pct[1], " traps\n", sep = "")
cat("Spring: ", detection_requirements$n_traps_80pct[2], " traps\n\n", sep = "")

cat("ADEQUACY:\n")
cat("Autumn: ", round(adequacy_summary$pct_adequate_80[1], 1), "% adequate\n", sep = "")
cat("Spring: ", round(adequacy_summary$pct_adequate_80[2], 1), "% adequate\n\n", sep = "")

cat("========================================\n")


