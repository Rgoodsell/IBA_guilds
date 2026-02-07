Readme
================
Rob Goodsell
2026-02-07

# Repository Overview

This repository contains data, scripts, and figures for the preliminary
analyses of invertebrate functional trait data derived from the [Insect
Biome Atlas project](https://www.insectbiomeatlas.org/).

------------------------------------------------------------------------

## Data

### `data/rawdata`

Tidy data (data/tidydata) were derived from raw amplicon sequence data
from a large ‘deep-metabarcoding’ study of Swedish invertebrates. Raw
amplicon sequence data needed to replicate these analyses can be found
here: <https://doi.org/10.17044/scilifelab.25480681>.

### `data/tidydata/`

- **all_data.rds**: Combined environmental and guild-level species
  richness data.
- **all_traits.rds**: Tidied trait data for all OTU’s within the IBA
  data - designates each species as one of 6 functional guilds
  (Phytophage, Saprophage, Predator, or parasitoids of each.)
- **env_pred_layer.rds**: Environmental layers for spatial predictions.
- **hab_pred_layer.rds**: Another prediction layer - habitat cover
  variables only.
- **niche_SR_data.rds**: Weekly species richness observations aggregated
  by Guild.
- **OTU_trait_data.rds**: Guild designations for each individual OTU.
- **split_data.rds**: A list of K folded data used to run a quick CV
  exercise to check model performance.
- **subfamily_classifications.rds**: Sub-family classifications required
  to classify parasitoid guilds. Retrived from NCBI database.
- **trap_habitat.rds**: A summary of habitat variables for each trap -
  used to create environmental prediction layers.
- **weekly_prec.rds**: Weekly precipitation data derived from
  <https://doi.org/10.1002/joc.4436>
- **weekly_temp.rds**: Weekly temperature data derived from
  <https://doi.org/10.1002/joc.4436>
- **mod_pk.rds**: Model fitting script will save a model of this name to
  this directory.

------------------------------------------------------------------------

## Scripts

### `R/`

- **1_tidy_trait_data.R**: Cleans and tidies trait data from raw data
  trait data.
- **2_format_otu_data.R**: Cleans and tidies OTU-level observations from
  raw amplicon data found here:
  <https://doi.org/10.17044/scilifelab.25480681>.
- **3_format_env_data.R**: Formats environmental covariate data.
- **4_run_model.R**: Example model fitting. Models are computationally
  intensive and were fit on the
  [dardel](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system/about-the-dardel-system-1.1053338)
  HPC cluster.
- **5_plot_coeffs.R**: Plots coefficients and habitat level-predictions
  from the fitted model
- **6_plot_spatemp.R**: Plots spatial and temporal trends estimated from
  the fitted model.
- **functions.R**: Contains custom functions.

------------------------------------------------------------------------

## Figures

### `figures/`

- **guild_distributions.tiff**: Figure of spatial patterns of guilds.
- **guild_phenology.tiff**: Figure of temporal patterns of guilds.
- **habitat_coeffs.tiff**: Figure of habitat coefficients.
- **smooth_coefficient_plots.tiff**: Figure of smooth coefficients.
