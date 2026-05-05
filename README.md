# albo_winter_monitoring
This repository hosts the data and code for the study of winter egg-laying activity of *Aedes albopictus* (Asian tiger mosquito). The project integrates multi-scale spatial analysis, long-term temporal monitoring in Northern Italy, and Citizen Science observations to understand better the winter egg-laying activity of this vector in Europe.

## Project Overview

As *Aedes albopictus* continues to expand across Europe, understanding its winter activity is crucial for vector control. This project analyses:
1.  **Latitudinal Gradients:** Presence/absence of winter egg-laying across a broad European transect.
2.  **Longitudinal Trends:** A 15-year temporal analysis of activity within the Emilia-Romagna (ER) region, Italy.
3.  **Citizen Science observations:** Spatial and temporal patterns of volunteer-contributed data for monitoring invasive mosquitoes.
4.  **Sampling Optimisation:** Power analysis to determine the effort required for reliable winter surveillance.

## Repository Structure

* **`data/`**: Contains raw and processed datasets (e.g., European coordinates, 15-year ER monitoring logs, and CS observation data).
* **`script/`**: Sequential R scripts for data processing and statistical modelling.
* **`output/`**: Stores generated results, including plots, maps, and statistical model summaries.

## Analysis Pipeline

### 1. European Latitudinal Gradient
**File:** `1_lat_gradient_analysis.R`  
Analyses how the probability of winter egg-laying changes with latitude across Europe.

### 2. Emilia-Romagna Temporal Trends
**File:** `2_ER_temporal_trend_analysis.R`  
A longitudinal study using a 15-year dataset from Northern Italy. 

### 3. Citizen Science (CS)
**File:** `3_CS_analysis.R`  
Exploratory analysis of spatio-temporal patterns of CS observations

### 4. Power Analysis & Sampling Effort
**File:** `4_Power_analysis.R`  
Determines the statistical power of the monitoring design. It calculates the necessary sampling intensity (number of traps/frequency) required to detect low-density winter egg-laying with high confidence.
