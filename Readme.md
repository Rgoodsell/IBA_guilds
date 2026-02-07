# Repository Overview

## This repository contains data, scripts, and figures for the analyses of invertebrate functional trait data derived from the Insect Biome Atlas project: <https://www.insectbiomeatlas.org/>.

## Data

### `data/rawdata`

Raw amplicon sequence data needed to replicate these analyses can be
found here: <https://doi.org/10.17044/scilifelab.25480681>.

### `data/tidydata/`

-   **all\_data.rds**: \[Brief description, e.g., “Combined dataset for
    all analyses.”\]
-   **all\_traits.csv**: \[Brief description, e.g., “Trait data for all
    species.”\]
-   **all\_traits.rds**: \[Brief description, e.g., “Trait data in RDS
    format.”\]
-   **env\_pred\_layer.rds**: \[Brief description, e.g., “Environmental
    prediction layers.”\]
-   **hab\_pred\_layer.rds**: \[Brief description, e.g., “Habitat
    prediction layers.”\]
-   **niche\_SR\_data.rds**: \[Brief description, e.g., “Species
    richness and niche data.”\]
-   **OTU\_trait\_data.rds**: \[Brief description, e.g., “OTU-level
    trait data.”\]
-   **split\_data.rds**: \[Brief description, e.g., “Training/testing
    split data.”\]
-   **subfamily\_classifications.rds**: \[Brief description, e.g.,
    “Subfamily classification data.”\]
-   **trap\_habitat.rds**: \[Brief description, e.g., “Trap and habitat
    association data.”\]
-   **weekly\_prec.rds**: \[Brief description, e.g., “Weekly
    precipitation data.”\]
-   **weekly\_temp.rds**: \[Brief description, e.g., “Weekly temperature
    data.”\]
-   **mod\_pk.rds**: Model fitting script will save a model of this name
    to this directory.

------------------------------------------------------------------------

## Scripts

### `R/`

-   **1\_tidy\_trait\_data.R**: \[Brief description, e.g., “Cleans and
    tidies trait data.”\]
-   **2\_format\_otu\_data.R**: \[Brief description, e.g., “Formats OTU
    data for analysis.”\]
-   **3\_format\_env\_data.R**: \[Brief description, e.g., “Prepares
    environmental data.”\]
-   **4\_run\_model.R**: \[Brief description, e.g., “Runs statistical
    models.”\]
-   **5\_plot\_coeffs.R**: \[Brief description, e.g., “Plots model
    coefficients.”\]
-   **6\_plot\_coeffs.R**: \[Brief description, e.g., “Alternative
    coefficient plots.”\]
-   **6\_plot\_spatemp.R**: \[Brief description, e.g., “Plots
    spatiotemporal patterns.”\]
-   **functions.R**: \[Brief description, e.g., “Custom functions for
    analysis.”\]

------------------------------------------------------------------------

## Figures

### `figures/`

-   **guild\_distributions.tiff**: \[Brief description, e.g.,
    “Distribution of guilds.”\]
-   **guild\_phenology.tiff**: \[Brief description, e.g., “Phenology of
    guilds.”\]
-   **habitat\_coeffs.tiff**: \[Brief description, e.g., “Habitat
    coefficient plots.”\]
-   **smooth\_coefficient\_plots.tiff**: \[Brief description, e.g.,
    “Smoothed coefficient plots.”\]

------------------------------------------------------------------------

## Usage

1.  **Data Preparation**: Run scripts 1–3 to prepare data.
2.  **Modeling**: Run `4_run_model.R` to fit models.
3.  **Visualization**: Use scripts 5–6 to generate figures.

------------------------------------------------------------------------

## Dependencies

\`\`\`r
