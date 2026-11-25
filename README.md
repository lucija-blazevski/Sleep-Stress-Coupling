# Individual Differences in Sleep-Related Distress Alleviation Predict Change in Anxiety Severity

Analysis code for the Brief Report of the same title. The scripts reproduce the
descriptive tables, the main-text analyses, Figure 2, and the analyses in
Supplement 1.

## Scripts

| Script | Reproduces |
|---|---|
| `00_helpers.R` | Checks and attaches the required packages; shared helper functions |
| `00_data.R` | Data preparation shared by `01` and `02` (sourced automatically) |
| `00_tables.R` | Table 1, Table S1 and Table S2 |
| `01_main_analysis.R` | Sample description, Analysis 1, Analysis 2 and Figure 2 (`Figure 2.png`) |
| `02_supplement_analysis.R` | Tables S3–S10 and the sensitivity analyses in Supplement 1, plus the analyses reported only in the response to reviewers |
| `marginal.stan` | Joint Bayesian model used for Table S10 |

`02_supplement_analysis.R` runs `01_main_analysis.R` first, which runs
`00_data.R`. Run the scripts with this folder as the working directory, for
example:

```r
source("00_tables.R")
source("02_supplement_analysis.R")   # also runs 01_main_analysis.R
```

## Data

The scripts read `data_all_studies_merged_30072026.csv` from a `data/` folder
next to the scripts. To use another location, set the environment variable
`STRESS_SLEEP_DATA` to the folder that contains the file, for example
`Sys.setenv(STRESS_SLEEP_DATA = "path/to/folder")` before sourcing a script.
The data are not included in this repository.

## Requirements

- R 4.4.2
- CRAN packages: dplyr, stringr, tibble, readr, tidyr, ggplot2, patchwork,
  lme4, lmerTest, posterior
- For Table S10 only: [cmdstanr](https://mc-stan.org/cmdstanr/) with CmdStan
  installed. Sampling takes several minutes; set `run_bayesian_model <- FALSE`
  in `02_supplement_analysis.R` to skip it.
