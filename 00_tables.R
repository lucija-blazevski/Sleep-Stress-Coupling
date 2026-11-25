# =============================================================================
# Descriptive tables
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Builds the descriptive tables reported in the manuscript:
#
#   Table 1     Participant characteristics, overall Analysis 1 sample
#   Table S1    Baseline characteristics by treatment condition
#   Table S2    Baseline sleep characteristics, overall and by treatment arm
#
# R version   : 4.4.2
# =============================================================================


# Attach packages and helper functions
source("00_helpers.R")


# Folder with the merged input CSV (same as 00_data.R): data/ next to the
# scripts by default. Set the environment variable STRESS_SLEEP_DATA to use a
# different location.
root <- Sys.getenv(
  "STRESS_SLEEP_DATA",
  "data"
)


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Standardize IDs, derive a numeric session variable, and add within-subject
# centered night-level predictors. At this stage, retain nights with complete
# evening distress, morning distress, and REM duration.

d <- read_csv(
  file.path(root, "data_all_studies_merged_30072026.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    subject = standardize_nsr_id(subject_id),
    ses_num = as.integer(str_remove(ses, "^N"))
  ) %>%
  group_by(subject) %>%
  mutate(
    E_c = within_center(stress_rating_e),
    R_c = within_center(r_of_sleep_period_min),
    N1_c = within_center(n1_of_sleep_period_min),
    N2_c = within_center(n2_of_sleep_period_min),
    N3_c = within_center(n3_of_sleep_period_min),
    TST_c = within_center(tst_corrected_min),
    SE_c = within_center(sleep_eff_corrected_perc),
    WASO_c = within_center(
      wake_after_sleep_onset_of_sleep_period_min),
    SOL_c = within_center(sleep_onset_delay_min)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(stress_rating_m),
    !is.na(stress_rating_e),
    !is.na(r_of_sleep_period_min))


# -----------------------------------------------------------------------------
# Define Analysis 1 sample
# -----------------------------------------------------------------------------

# Analysis 1 requires at least 2 eligible evening-night-morning observations.
# This same participant sample is used for all baseline descriptive tables.

d_analysis1 <- d %>%
  group_by(subject) %>%
  filter(n() >= 2) %>%
  ungroup()


analysis1_sample_check <- d_analysis1 %>%
  summarise(
    participants = n_distinct(subject),
    nights = n()
  )

print(analysis1_sample_check)


# Expected current analytic sample
stopifnot(
  n_distinct(d_analysis1$subject) == 223,
  nrow(d_analysis1) == 662
)


# -----------------------------------------------------------------------------
# Create participant-level baseline dataset
# -----------------------------------------------------------------------------

diagnosis_vars <- c(
  "GAD_diagnosis",
  "SAD_diagnosis",
  "PD_diagnosis",
  "PTSD_diagnosis"
)


# One row per participant, exclusively from the Analysis 1 sample.
# The earliest eligible session is retained for participant-level variables.
# Diagnosis variables are normalized to lower-case yes/no strings.

bl_all <- d_analysis1 %>%
  arrange(subject, ses_num) %>%
  distinct(subject, .keep_all = TRUE) %>%
  mutate(
    treatment = factor(
      condition,
      levels = c(2, 1),
      labels = c("CBT-I", "Control")
    ),
    across(
      all_of(diagnosis_vars),
      ~ str_to_lower(str_trim(as.character(.x)))
    ),
    n_diagnoses_available = rowSums(
      across(
        all_of(diagnosis_vars),
        ~ .x %in% c("yes", "no")
      )
    ),
    any_clinical_diagnosis = case_when(
      if_any(
        all_of(diagnosis_vars),
        ~ replace_na(.x == "yes", FALSE)
      ) ~ "yes",
      n_diagnoses_available > 0 ~ "no",
      TRUE ~ NA
    )
  )


stopifnot(nrow(bl_all) == 223)

print(table(bl_all$treatment, useNA = "ifany"))


# Treatment-group sample sizes used in supplement tables
n_cbti <- sum(bl_all$treatment == "CBT-I", na.rm = TRUE)
n_control <- sum(bl_all$treatment == "Control", na.rm = TRUE)


# =============================================================================
# 2. TABLE 1 — PARTICIPANT CHARACTERISTICS, OVERALL SAMPLE ONLY
# =============================================================================

# The main-manuscript table contains the overall Analysis 1 sample only.
# Treatment-specific descriptives are reported in Table S1.

table1 <- bind_rows(
  
  section_total_row("Demographics"),
  total_summary_row(
    "Age, mean (SD)",
    "age",
    format_mean_sd
  ),
  total_summary_row(
    "Female, No. (%)",
    "female",
    format_n_percent
  ),
  
  section_total_row("Insomnia severity"),
  total_summary_row(
    "ISI total, mean (SD)",
    "ISI_T0",
    format_mean_sd
  ),
  
  section_total_row("RMT20 subscale, mean (SD)"),
  total_summary_row(
    "Generalized anxiety disorder",
    "RMT20_GAD_T0",
    format_mean_sd
  ),
  total_summary_row(
    "Social anxiety disorder",
    "RMT20_SAD_T0",
    format_mean_sd
  ),
  total_summary_row(
    "Panic disorder",
    "RMT20_panic_disorder_T0",
    format_mean_sd
  ),
  total_summary_row(
    "Posttraumatic stress disorder",
    "RMT20_PTSD_T0",
    format_mean_sd
  ),
  
  section_total_row("BPD-C severity, mean (SD)"),
  total_summary_row(
    "Borderline personality disorder",
    "BPD_C_T0",
    format_mean_sd
  ),
  
  section_total_row("MINI-confirmed clinical diagnoses, No. (%)"),
  total_summary_row(
    "Generalized anxiety disorder",
    "GAD_diagnosis",
    format_diagnosis
  ),
  total_summary_row(
    "Social anxiety disorder",
    "SAD_diagnosis",
    format_diagnosis
  ),
  total_summary_row(
    "Panic disorder",
    "PD_diagnosis",
    format_diagnosis
  ),
  total_summary_row(
    "Posttraumatic stress disorder",
    "PTSD_diagnosis",
    format_diagnosis
  ),
  total_summary_row(
    "Any of the above diagnoses",
    "any_clinical_diagnosis",
    format_diagnosis
  )
)


names(table1) <- c(
  "Characteristic",
  sprintf("Total (N = %d)", nrow(bl_all))
)


print(table1, n = Inf, width = Inf)


# =============================================================================
# 3. TABLE S1 — BASELINE CHARACTERISTICS BY TREATMENT CONDITION
# =============================================================================

# Continuous variables are compared descriptively using Welch two-sample
# t-tests. Binary and diagnosis variables are compared using Pearson
# chi-square tests. P values are unadjusted and descriptive.

etable1 <- bind_rows(
  
  section_group_row("Demographics"),
  continuous_group_row(
    "Age, mean (SD)",
    "age"
  ),
  binary_group_row(
    "Female, No. (%)",
    "female"
  ),
  
  section_group_row("Insomnia severity"),
  continuous_group_row(
    "ISI total, mean (SD)",
    "ISI_T0"
  ),
  
  section_group_row("RMT20 subscale, mean (SD)"),
  continuous_group_row(
    "Generalized anxiety disorder",
    "RMT20_GAD_T0"
  ),
  continuous_group_row(
    "Social anxiety disorder",
    "RMT20_SAD_T0"
  ),
  continuous_group_row(
    "Panic disorder",
    "RMT20_panic_disorder_T0"
  ),
  continuous_group_row(
    "Posttraumatic stress disorder",
    "RMT20_PTSD_T0"
  ),
  
  section_group_row("BPD-C severity, mean (SD)"),
  continuous_group_row(
    "Borderline personality disorder",
    "BPD_C_T0"
  ),
  
  section_group_row("MINI-confirmed clinical diagnoses, No. (%)"),
  diagnosis_group_row(
    "Generalized anxiety disorder",
    "GAD_diagnosis"
  ),
  diagnosis_group_row(
    "Social anxiety disorder",
    "SAD_diagnosis"
  ),
  diagnosis_group_row(
    "Panic disorder",
    "PD_diagnosis"
  ),
  diagnosis_group_row(
    "Posttraumatic stress disorder",
    "PTSD_diagnosis"
  ),
  diagnosis_group_row(
    "Any of the above diagnoses",
    "any_clinical_diagnosis"
  )
)


names(etable1) <- c(
  "Characteristic",
  sprintf("CBT-I (n = %d)", n_cbti),
  sprintf("Control (n = %d)", n_control),
  "P Value"
)


print(etable1, n = Inf, width = Inf)


# =============================================================================
# 4. TABLE S2 — BASELINE SLEEP CHARACTERISTICS
# =============================================================================

# For each participant, average each sleep variable across their eligible
# Analysis 1 nights. Participant-level means are then summarized as mean (SD)
# overall and by treatment condition.

sleep_variables <- c(
  "N1 duration, min" =
    "n1_of_sleep_period_min",
  "N2 duration, min" =
    "n2_of_sleep_period_min",
  "N3 duration, min" =
    "n3_of_sleep_period_min",
  "REM duration, min" =
    "r_of_sleep_period_min",
  "Total sleep time, min" =
    "tst_corrected_min",
  "Sleep efficiency, %" =
    "sleep_eff_corrected_perc",
  "Wake after sleep onset, min" =
    "wake_after_sleep_onset_of_sleep_period_min",
  "Sleep onset latency, min" =
    "sleep_onset_delay_min"
)


sleep_long <- d_analysis1 %>%
  mutate(
    treatment = factor(
      condition,
      levels = c(2, 1),
      labels = c("CBT-I", "Control")
    )
  ) %>%
  group_by(subject, treatment) %>%
  summarise(
    across(
      all_of(unname(sleep_variables)),
      mean_or_na
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(unname(sleep_variables)),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Characteristic = factor(
      variable,
      levels = unname(sleep_variables),
      labels = names(sleep_variables)
    )
  )


stopifnot(
  n_distinct(sleep_long$subject) == 223
)


sleep_overall <- sleep_long %>%
  group_by(Characteristic) %>%
  summarise(
    value = format_mean_sd(value),
    .groups = "drop"
  ) %>%
  mutate(group = "Total")


sleep_by_group <- sleep_long %>%
  filter(!is.na(treatment)) %>%
  group_by(
    Characteristic,
    treatment
  ) %>%
  summarise(
    value = format_mean_sd(value),
    .groups = "drop"
  ) %>%
  rename(group = treatment)


etable2 <- bind_rows(
  sleep_overall,
  sleep_by_group
) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Total",
        "CBT-I",
        "Control"
      )
    )
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = value
  ) %>%
  arrange(Characteristic)


names(etable2) <- c(
  "Sleep variable",
  sprintf("Total (N = %d)", nrow(bl_all)),
  sprintf("CBT-I (n = %d)", n_cbti),
  sprintf("Control (n = %d)", n_control)
)


print(etable2, n = Inf, width = Inf)
