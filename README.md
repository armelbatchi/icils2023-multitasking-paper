# Paper B — Academic-Media Multitasking and Digital Performance

[![DOI](https://zenodo.org/badge/1201820421.svg)](https://doi.org/10.5281/zenodo.19426981)

R code for the ICILS 2023 analysis pipeline and the TIMSS 2023 secondary analyses used in this project.

## Repository structure

- `run_all.R`: runs the full pipeline in sequence
- `R/00_setup.R`: package loading, local directory resolution, constants, and shared helpers
- `R/01_load_icils.R`: loads and prepares ICILS 2023 data
- `R/02_icils_descriptives_and_figures.R`: descriptives, tables, and initial ICILS figures
- `R/03_icils_country_models_and_forests.R`: ICILS country models, meta-analysis, and main figures
- `R/04_icils_secondary_analyses.R`: moderation, multilevel, quantile, and robustness analyses
- `R/05_load_timss.R`: loads and prepares TIMSS 2023 data
- `R/06_timss_main_analyses.R`: TIMSS secondary analyses, forest plots, and robustness checks
- `R/07_cross_study_summary.R`: cross-study comparison tables and pooled summaries
- `R/08_supplementary_materials.R`: supplementary tables, diagnostics, and derived outputs

## Data availability

The data used in this project are not included in this repository. ICILS data are available from the [IEA ICILS Data Repository](https://www.iea.nl/data-tools/repository/icils), and TIMSS data are available from the [IEA TIMSS Data Repository](https://www.iea.nl/data-tools/repository/timss), subject to the IEA’s terms and conditions of access and use.

## Running the code

This repository is code only. Before running the scripts, download the source data from the official IEA websites and store them locally outside GitHub tracking, or point the scripts to your local folders with environment variables.

The scripts look for these local directories by default:

- `ICILS_DIR = .local_data/icils`
- `TIMSS_DIR = .local_data/timss`
- `OUT_DIR = .local_output`

You can override them in your R session, for example:

```r
Sys.setenv(
  ICILS_DIR = "/path/to/local/icils",
  TIMSS_DIR = "/path/to/local/timss",
  OUT_DIR   = "/path/to/local/output"
)
source("run_all.R")
```

