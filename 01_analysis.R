# =============================================================================
# Analysis code
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Reproduces the analyses reported in the manuscript:
#
#   Analysis 1  Within-person REM duration predicts next-morning distress.
#               The same mixed model provides participant-specific REM slopes.
#   eTable 3    Other sleep variables as predictors of overnight distress
#               alleviation, with Holm-Bonferroni correction across 7 secondary measures.
#   Analysis 2  Participant-specific REM–distress slopes predict GAD severity
#               at 2-month follow-up, adjusting for baseline GAD severity.
#   eTable 4    Estimates from the primary GAD outcome model.
#   eTable 5    Secondary symptom dimensions, Holm-corrected across 4 outcomes.
#   Figure 2    Panels A–C summarizing Analyses 1 and 2.
#
# R version   : 4.4.2
# =============================================================================


# Attach packages and helper functions
source("00_helpers.R")


# Folder with the merged input CSV
root <- "path/to/the/folder/with/csv"


# Figure palette
BLUE <- "#312783"
GREEN <- "#0F8A2E"


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Standardize IDs, derive a numeric session index, and create within-subject
# centered night-level predictors. Eligible nights require both distress ratings
# and REM duration.

d <- read_csv(
  file.path(root, "name_of_the_df.csv"),
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


analysis1_sample_check <- d_analysis1 %>%
  summarise(
    participants = n_distinct(subject),
    nights = n()
  )

analysis1_sample_check


# =============================================================================
# 2. ANALYSIS 1 — REM SLEEP AND OVERNIGHT DISTRESS ALLEVIATION
# =============================================================================

# The primary model serves two purposes:
#   1. its fixed REM effect tests the average within-person association;
#   2. its participant-specific REM slopes provide the individualized
#      REM–distress coupling measure used in Analysis 2.

m_primary <- lmer(
  stress_rating_m ~
    R_c + E_c +
    (1 + R_c || subject),
  data = d_analysis1,
  REML = TRUE
)

summary(m_primary)
confint(m_primary)


# =============================================================================
# 3. eTABLE 3 — OTHER SLEEP VARIABLES
# =============================================================================

# Exploratory comparison of REM with other sleep architecture and continuity
# measures. Each sleep predictor is standardized within participant inside
# `fit_sleep_predictor()`.
#
# For each measure, the random-intercept + random-slope specification is
# retained only when it converges and is nonsingular otherwise the function
# falls back to a random-intercept-only model.
#
# REM is the prespecified primary sleep variable and is therefore excluded from
# the Bonferroni-Holm correction across the 7 secondary sleep measures.

sleep_predictors <- c(
  "REM duration" = "r_of_sleep_period_min",
  "N1 duration" = "n1_of_sleep_period_min",
  "N2 duration" = "n2_of_sleep_period_min",
  "N3 duration" = "n3_of_sleep_period_min",
  "Total sleep time" = "tst_corrected_min",
  "Sleep efficiency" = "sleep_eff_corrected_perc",
  "Wake after sleep onset" = "wake_after_sleep_onset_of_sleep_period_min",
  "Sleep onset latency" = "sleep_onset_delay_min")


etable3 <- lapply(
  names(sleep_predictors),
  function(label) {
    fit_sleep_predictor(
      label = label,
      predictor = sleep_predictors[[label]]
    )
  }) %>%
  bind_rows()

etable3$p_holm <- NA_real_

secondary_sleep <- etable3$`Sleep variable` != "REM duration"

etable3$p_holm[secondary_sleep] <- p.adjust(
  etable3$p[secondary_sleep],
  method = "holm")


print(etable3, n = Inf, width = Inf)


# =============================================================================
# 4. ANALYSIS 2 — REM-RELATED DISTRESS ALLEVIATION AND FOLLOW-UP GAD
# =============================================================================


coupling_info <- coef(m_primary)$subject %>%
  as.data.frame() %>%
  rownames_to_column("subject") %>%
  transmute(
    subject = as.character(subject),
    rem_slope = R_c,
    rem_slope_sd = as.numeric(-scale(rem_slope)))
# Higher rem_slope_sd values indicate stronger REM-related overnight distress
# alleviation because the original slopes are reverse-coded after
# standardization.


# Subject-level variables carried into the outcome analyses.
subject_level_vars <- c(
  "condition",
  "RMT20_GAD_T0", "RMT20_GAD_T1",
  "RMT20_SAD_T0", "RMT20_SAD_T1",
  "RMT20_panic_disorder_T0", "RMT20_panic_disorder_T1",
  "RMT20_PTSD_T0", "RMT20_PTSD_T1",
  "BPD_C_T0", "BPD_C_T1")


# One row per Analysis 1 participant, joined to the individualized REM slope.
# Treatment is effect-coded:
#   Control = -0.5
#   CBT-I   = +0.5
#
# With this coding, the rem_slope_sd coefficient represents the average
# association across the two treatment conditions.

prediction_df <- d_analysis1 %>%
  mutate(
    subject = as.character(subject),
    condition = as.character(condition)
  ) %>%
  group_by(subject) %>%
  summarise(
    across(
      all_of(subject_level_vars),
      first_nonmissing
    ),
    .groups = "drop"
  ) %>%
  inner_join(
    coupling_info,
    by = "subject"
  ) %>%
  mutate(
    condition_ec = case_when(
      condition == "1" ~ -0.5,
      condition == "2" ~  0.5,
      TRUE ~ NA
    )
  )


# Primary outcome model: follow-up GAD severity, adjusted for baseline GAD,
# including treatment condition and the coupling x treatment interaction.

m_gad_outcome <- lm(
  RMT20_GAD_T1 ~
    rem_slope_sd * condition_ec +
    RMT20_GAD_T0,
  data = prediction_df
)


summary(m_gad_outcome)
confint(m_gad_outcome)

dat_gad <- model.frame(m_gad_outcome) %>%
  as_tibble()


# =============================================================================
# 5. eTABLE 4 — COMPLETE GAD OUTCOME MODEL
# =============================================================================


gad_coupling <- model_term(m_gad_outcome,"rem_slope_sd")
gad_condition <- model_term(m_gad_outcome,"condition_ec")
gad_baseline <- model_term(m_gad_outcome,"RMT20_GAD_T0")
gad_interaction <- model_term(m_gad_outcome, "rem_slope_sd:condition_ec")

etable4 <- bind_rows(
  tibble(
    Variable = "REM-related overnight distress alleviation, per 1 SD",
    b = gad_coupling$estimate,
    lower_95 = gad_coupling$lower_95,
    upper_95 = gad_coupling$upper_95,
    p = gad_coupling$p
  ),
  tibble(
    Variable = "Treatment condition (CBT-I vs control)",
    b = gad_condition$estimate,
    lower_95 = gad_condition$lower_95,
    upper_95 = gad_condition$upper_95,
    p = gad_condition$p
  ),
  tibble(
    Variable = "Baseline GAD severity",
    b = gad_baseline$estimate,
    lower_95 = gad_baseline$lower_95,
    upper_95 = gad_baseline$upper_95,
    p = gad_baseline$p
  ),
  tibble(
    Variable = "REM-related distress alleviation x treatment condition",
    b = gad_interaction$estimate,
    lower_95 = gad_interaction$lower_95,
    upper_95 = gad_interaction$upper_95,
    p = gad_interaction$p
  )
)


print(etable4, n = Inf, width = Inf)


# =============================================================================
# 6. eTABLE 5 — SECONDARY SYMPTOM DIMENSIONS
# =============================================================================

# Repeat the Analysis 2 model for the 4 prespecified secondary outcomes.
# The main coupling effects are Holm-corrected across these 4 outcomes.
# Interaction tests are corrected as a separate family for the supplementary
# narrative; they are not columns in the manuscript's eTable 5.

secondary_outcomes <- list(
  SAD = c(
    "RMT20_SAD_T1",
    "RMT20_SAD_T0"),
  PD = c(
    "RMT20_panic_disorder_T1",
    "RMT20_panic_disorder_T0"),
  PTSD = c(
    "RMT20_PTSD_T1",
    "RMT20_PTSD_T0"),
  BPD = c(
    "BPD_C_T1",
    "BPD_C_T0"))

etable5_full <- lapply(
  names(secondary_outcomes),
  function(outcome) {
    fit_outcome_model(
      outcome = outcome,
      variables = secondary_outcomes[[outcome]]
    )
  }) %>%
  bind_rows()

# Holm correction across the 4 secondary main effects.
etable5_full$p_holm <- p.adjust(
  etable5_full$p,
  method = "holm")

# Holm correction across the 4 secondary interaction tests.
etable5_full$interaction_p_holm <- p.adjust(
  etable5_full$interaction_p,
  method = "holm")

# Manuscript eTable 5 contains the main coupling association only, the rest is in text.
etable5 <- etable5_full %>%
  select(
    Outcome,
    n,
    b,
    SE,
    lower_95,
    upper_95,
    p,
    p_holm)

print(etable5, n = Inf, width = Inf)


# Keep interaction results available for the supplementary narrative.
secondary_interactions <- etable5_full %>%
  select(
    Outcome,
    n,
    interaction_b,
    interaction_SE,
    interaction_lower_95,
    interaction_upper_95,
    interaction_p,
    interaction_p_holm)

print(secondary_interactions, n = Inf, width = Inf)


# =============================================================================
# 7. FIGURE 2
# =============================================================================

# Panel A: fitted average within-person REM association from Analysis 1.
# Panel B: unadjusted participant-level REM slopes, shown for illustration.
# Panel C: adjusted association between individualized REM-related distress
#          alleviation and GAD severity at follow-up.

primary_rem_label <- format_effect_label(m_primary, "R_c")
analysis2_label <- format_effect_label( m_gad_outcome, "rem_slope_sd")

pred_primary <- ggpredict(m_primary, terms = "R_c [all]")
pred_gad <- ggpredict(m_gad_outcome, terms = "rem_slope_sd [all]")


# -----------------------------------------------------------------------------
# Panel A — REM duration and next-morning distress
# -----------------------------------------------------------------------------

p_A <- ggplot(
  pred_primary,
  aes(
    x = x,
    y = predicted
  )
) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = GREEN,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    colour = GREEN,
    linewidth = 1
  ) +
  richtext_annotation(
    x = 0,
    y = 5.5,
    label = primary_rem_label,
    colour = GREEN
  ) +
  scale_x_continuous(
    breaks = seq(-100, 100, 50)
  ) +
  scale_y_continuous(
    breaks = seq(0, 6, 2)
  ) +
  coord_cartesian(
    xlim = c(-100, 100),
    ylim = c(0, 6)
  ) +
  labs(
    title = "A",
    x = "Within-person REM duration (min)",
    y = "Morning distress"
  ) +
  theme_paper()


# -----------------------------------------------------------------------------
# Panel B — participant-level variability in REM–distress slopes
# -----------------------------------------------------------------------------

d_primary_model <- model.frame(m_primary) %>%
  as_tibble()

analysis1_slopes <- d_primary_model %>%
  group_by(subject) %>%
  summarise(
    slope = coef(
      lm(
        stress_rating_m ~ R_c
      )
    )[["R_c"]],
    .groups = "drop"
  ) %>%
  mutate(
    slope_direction = factor(
      if_else(
        slope < 0,
        "Negative slope",
        "Zero or positive slope"
      ),
      levels = c(
        "Negative slope",
        "Zero or positive slope"
      )
    )
  )


spaghetti_data <- d_primary_model %>%
  inner_join(
    analysis1_slopes,
    by = "subject"
  )


p_B <- ggplot(
  spaghetti_data,
  aes(
    x = R_c,
    y = stress_rating_m,
    group = subject,
    colour = slope_direction
  )
) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.5,
    alpha = 0.55
  ) +
  scale_colour_manual(
    values = c(
      "Negative slope" = GREEN,
      "Zero or positive slope" = "grey"
    ),
    name = NULL
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        linewidth = 1.1
      )
    )
  ) +
  theme_paper() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "vertical",
    legend.text = element_text(size = 17),
    legend.spacing.y = grid::unit(0.1, "cm"),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq(-100, 100, 50)
  ) +
  scale_y_continuous(
    breaks = seq(0, 10, 2)
  ) +
  coord_cartesian(
    xlim = c(-100, 100),
    ylim = c(0, 10)
  ) +
  labs(
    title = "B",
    x = "Within-person REM duration (min)",
    y = "Morning distress"
  )


# -----------------------------------------------------------------------------
# Panel C — REM-related distress alleviation and follow-up GAD
# -----------------------------------------------------------------------------

p_C <- ggplot() +
  geom_point(
    data = dat_gad,
    aes(
      x = rem_slope_sd,
      y = RMT20_GAD_T1
    ),
    alpha = 0.25,
    size = 1.8,
    colour = "darkgrey"
  ) +
  geom_ribbon(
    data = pred_gad,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = BLUE,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    data = pred_gad,
    aes(
      x = x,
      y = predicted
    ),
    colour = BLUE,
    linewidth = 1
  ) +
  richtext_annotation(
    x = -2.5,
    y = 18.5,
    label = analysis2_label,
    colour = BLUE,
    hjust = 0
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 2)
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, 4)
  ) +
  coord_cartesian(
    ylim = c(0, 20)
  ) +
  labs(
    title = "C",
    x = "REM-related distress alleviation (SD)",
    y = "GAD severity at T1"
  ) +
  theme_paper()


# -----------------------------------------------------------------------------
# Assemble and export
# -----------------------------------------------------------------------------

figure2 <- p_A | p_B | p_C
figure2


ggsave(
  filename = "figure2_combined.png",
  plot = figure2,
  width = 14.5,
  height = 6,
  dpi = 300
)
