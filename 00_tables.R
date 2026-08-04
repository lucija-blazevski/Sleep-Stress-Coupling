# =============================================================================
# Descriptive tables
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Builds the descriptive tables reported in the manuscript:
#
#   Table 1     Baseline demographic and clinical characteristics, overall and
#               by treatment arm.
#   eTable 1    Baseline sleep-architecture characteristics, overall and by arm.
#
# R version   : 4.4.2
# =============================================================================

# Attach helpers
source("00_helpers.R")

# Folder with the merged input CSV.
root <- "path/to/the/folder/with/csv"


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================
# Standardize IDs, derive a numeric session, and add within-subject
# night-level predictors: mean-centered (`_c`) versions used in the primary
# models and, later, within-subject z-scores for the standardized table.

d2 <- read_csv(
  file.path(root, "name_of_the_df.csv"),
  show_col_types = FALSE) %>%
  mutate(
    subject = standardize_nsr_id(subject_id),
    ses_num = as.integer(str_remove(ses, "^N"))
  ) %>%
  filter(!is.na(ses_num)) %>%
  group_by(subject) %>%
  mutate(
    E_c = within_center(stress_rating_e),
    R_c = within_center(r_of_sleep_period_min),
    N1_c = within_center(n1_of_sleep_period_min),
    N2_c = within_center(n2_of_sleep_period_min),
    N3_c = within_center(n3_of_sleep_period_min),
    TST_c = within_center(tst_corrected_min),
    SE_c = within_center(sleep_eff_corrected_perc),
    WASO_c = within_center(wake_after_sleep_onset_of_sleep_period_min),
    SOL_c = within_center(sleep_onset_delay_min),
    ISI_change = ISI_T1 - ISI_T0
  ) %>%
  ungroup()


# =============================================================================
# 2. TABLE 1 — BASELINE CHARACTERISTICS
# =============================================================================
# One row per subject (first night with a valid REM value). Diagnoses are
# normalized to yes/no; `any_clinical_diagnosis` is "yes" if any dimension is
# "yes", "no" only when at least one dimension is observed, else NA.

diagnosis_vars <- c(
  "GAD_diagnosis",
  "SAD_diagnosis",
  "PD_diagnosis",
  "PTSD_diagnosis"
)

bl_all <- d2 %>%
  filter(!is.na(r_of_sleep_period_min)) %>%
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
      across(all_of(diagnosis_vars), ~ .x %in% c("yes", "no"))
    ),
    any_clinical_diagnosis = case_when(
      if_any(
        all_of(diagnosis_vars),
        ~ replace_na(.x == "yes", FALSE)
      ) ~ "yes",
      n_diagnoses_available > 0 ~ "no",
      TRUE ~ NA_character_
    )
  )

# Randomized subset (those with a treatment assignment) for the by-arm columns.
bl_randomized <- bl_all %>%
  filter(!is.na(treatment))

# `section_row()` inserts a header; `summary_row()` applies a formatter to the
# whole sample and each arm (both read bl_all / bl_randomized; see helpers.R).
table1 <- bind_rows(
  section_row("Demographics"),
  summary_row("Age, mean (SD)", "age", format_mean_sd),
  summary_row("Female, No. (%)", "female", format_n_percent),
  section_row("Insomnia severity"),
  summary_row("ISI total, mean (SD)", "ISI_T0", format_mean_sd),
  section_row("RMT20 anxiety dimensions, mean (SD)"),
  summary_row("Generalized anxiety disorder", "RMT20_GAD_T0", format_mean_sd),
  summary_row("Social anxiety disorder", "RMT20_SAD_T0", format_mean_sd),
  summary_row("Panic disorder", "RMT20_panic_disorder_T0", format_mean_sd),
  summary_row("Posttraumatic stress disorder", "RMT20_PTSD_T0", format_mean_sd),
  section_row("BPD-C severity, mean (SD)"),
  summary_row("Borderline personality disorder", "BPD_C_T0", format_mean_sd),
  section_row("MINI-confirmed clinical diagnoses, No. (%)"),
  summary_row("Generalized anxiety disorder", "GAD_diagnosis", format_diagnosis),
  summary_row("Social anxiety disorder", "SAD_diagnosis", format_diagnosis),
  summary_row("Panic disorder", "PD_diagnosis", format_diagnosis),
  summary_row("Posttraumatic stress disorder", "PTSD_diagnosis", format_diagnosis),
  summary_row("Any of the above diagnoses", "any_clinical_diagnosis", format_diagnosis)
)

# Column headers with sample sizes.
names(table1) <- c(
  "Characteristic",
  sprintf("Total (N = %d)", nrow(bl_all)),
  sprintf("CBT-I (n = %d)", sum(bl_randomized$treatment == "CBT-I")),
  sprintf("Control (n = %d)", sum(bl_randomized$treatment == "Control"))
)

print(table1, n = Inf, width = Inf)


# =============================================================================
# 3. eTABLE 1 — BASELINE SLEEP CHARACTERISTICS
# =============================================================================
# Subject-level mean of each sleep variable, summarized as mean (SD) overall
# and by arm, then reshaped to one row per characteristic.

sleep_variables <- c(
  "N1 duration, min" = "n1_of_sleep_period_min",
  "N2 duration, min" = "n2_of_sleep_period_min",
  "N3 duration, min" = "n3_of_sleep_period_min",
  "REM duration, min" = "r_of_sleep_period_min",
  "Total sleep time, min" = "tst_corrected_min",
  "Sleep efficiency, %" = "sleep_eff_corrected_perc",
  "Wake after sleep onset, min" = "wake_after_sleep_onset_of_sleep_period_min",
  "Sleep onset latency, min" = "sleep_onset_delay_min"
)

sleep_long <- d2 %>%
  mutate(
    treatment = factor(
      condition,
      levels = c(2, 1),
      labels = c("CBT-I", "Control")
    )
  ) %>%
  group_by(subject, treatment) %>%
  summarise(
    across(all_of(unname(sleep_variables)), mean_or_na),
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

sleep_overall <- sleep_long %>%
  group_by(Characteristic) %>%
  summarise(value = format_mean_sd(value), .groups = "drop") %>%
  mutate(group = "Total")

sleep_by_group <- sleep_long %>%
  filter(!is.na(treatment)) %>%
  group_by(Characteristic, treatment) %>%
  summarise(value = format_mean_sd(value), .groups = "drop") %>%
  rename(group = treatment)

etable1 <- bind_rows(sleep_overall, sleep_by_group) %>%
  mutate(group = factor(group, levels = c("Total", "CBT-I", "Control"))) %>%
  pivot_wider(names_from = group, values_from = value) %>%
  arrange(Characteristic)

print(etable1, n = Inf, width = Inf)
