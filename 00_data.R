# =============================================================================
# Data preparation
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Shared by 01_main_analysis.R and 02_supplement_analysis.R. Creates:
#
#   d             all eligible nights (both distress ratings and REM duration)
#   d_analysis1   nights of participants with at least 2 eligible nights
#   subject_df    one row per Analysis 1 participant with subject-level
#                 variables and effect-coded treatment condition
#
# R version   : 4.4.2
# =============================================================================


# Attach packages and helper functions
source("00_helpers.R")


# Folder with the merged input CSV: data/ next to the scripts by default (the
# data are not part of the repository). Set the environment variable
# STRESS_SLEEP_DATA to use a different location.
root <- Sys.getenv(
  "STRESS_SLEEP_DATA",
  "data"
)


# =============================================================================
# 1. NIGHT-LEVEL DATA
# =============================================================================

# Standardize IDs, derive a numeric session index, and create within-subject
# centered night-level predictors. Centering uses all of a participant's
# recorded nights; eligible nights require both distress ratings and REM
# duration.
#
# NREM = N1 + N2 + N3; REM % = REM / total sleep time x 100.

d <- read_csv(
  file.path(root, "data_all_studies_merged_30072026.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    subject = standardize_nsr_id(subject_id),
    ses_num = as.integer(str_remove(ses, "^N")),
    nrem_min = n1_of_sleep_period_min + n2_of_sleep_period_min +
      n3_of_sleep_period_min,
    rem_perc = 100 * r_of_sleep_period_min / tst_corrected_min
  ) %>%
  group_by(subject) %>%
  mutate(
    E_c = within_center(stress_rating_e),
    R_c = within_center(r_of_sleep_period_min),
    N1_c = within_center(n1_of_sleep_period_min),
    N2_c = within_center(n2_of_sleep_period_min),
    N3_c = within_center(n3_of_sleep_period_min),
    NREM_c = within_center(nrem_min),
    TST_c = within_center(tst_corrected_min),
    RPERC_c = within_center(rem_perc),
    SE_c = within_center(sleep_eff_corrected_perc),
    WASO_c = within_center(
      wake_after_sleep_onset_of_sleep_period_min
    ),
    SOL_c = within_center(sleep_onset_delay_min)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(stress_rating_m),
    !is.na(stress_rating_e),
    !is.na(r_of_sleep_period_min)
  )


# Analysis 1 requires at least 2 eligible nights per participant.
d_analysis1 <- d %>%
  group_by(subject) %>%
  filter(n() >= 2) %>%
  ungroup()


# =============================================================================
# 2. SUBJECT-LEVEL DATA
# =============================================================================

# Subject-level variables carried into the outcome analyses.
subject_level_vars <- c(
  "condition",
  "age", "female", "ISI_T0", "any_clinical_diagnosis",
  "RMT20_GAD_T0", "RMT20_GAD_T1",
  "RMT20_SAD_T0", "RMT20_SAD_T1",
  "RMT20_panic_disorder_T0", "RMT20_panic_disorder_T1",
  "RMT20_PTSD_T0", "RMT20_PTSD_T1",
  "BPD_C_T0", "BPD_C_T1")


# One row per Analysis 1 participant. Treatment is effect-coded:
#   Control = -0.5
#   CBT-I   = +0.5
# With this coding, a predictor's main effect represents its average
# association across the two treatment conditions.

subject_df <- d_analysis1 %>%
  mutate(
    subject = as.character(subject),
    condition = as.character(condition)
  ) %>%
  group_by(subject) %>%
  summarise(
    across(all_of(subject_level_vars), first_nonmissing),
    n_nights = n(),
    mean_rem_min = mean(r_of_sleep_period_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    condition_ec = case_when(
      condition == "1" ~ -0.5,
      condition == "2" ~  0.5,
      TRUE ~ NA_real_
    )
  )
