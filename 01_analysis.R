# =============================================================================
# Analysis code
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Reproduces the coupling analyses and Figure 3 from the manuscript:
#
#   Analysis 1  Within-person REM duration predicts next-morning distress
#               (linear mixed model, random subject intercept).
#   Table       Same specification for other sleep variables, within-subject
#               standardized, with Holm correction.
#   Analysis 2  Subject-level REM–distress coupling predicts the anxiety (GAD)
#               outcome after treatment.
#   Table       Secondary anxiety dimensions, coupling main effect and
#               coupling x condition interaction, Holm-corrected within family.
#   Figure 3    Panels A–C summarizing the above.
#
#
# Data availability: the input CSV contains individual-level clinical data and
# is not distributed with the code. See the manuscript's data-availability
# statement.
#
# R version   : 4.4.2
# =============================================================================

# Attach helpers
source("00_helpers.R")

# Folder with the merged input CSV.
root <- "path/to/the/folder/with/csv"

# Figure palette.
BLUE <- "#312783"
GREEN <- "#0F8A2E"


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================
# Standardize IDs, derive a numeric session index, and add within-subject
# night-level predictors: mean-centered versions used in the primary
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
# 2. ANALYSIS 1 — REM DURATION PREDICTS NEXT-MORNING DISTRESS
# =============================================================================
# Within-person REM duration (R_c) predicting morning distress, adjusting for
# evening distress and baseline covariates, with a random subject intercept.

m_primary <- lmer(
  stress_rating_m ~
    R_c + E_c + age + RMT20_GAD_T0 + ISI_T0 +
    (1 | subject),
  data = d2,
  REML = TRUE)

summary(m_primary)
confint(m_primary)


# =============================================================================
# 3. OTHER SLEEP VARIABLES
# =============================================================================
# Re-run the Analysis 1 specification for each sleep variable, within-subject
# z-standardized so the coefficients are comparable, then Holm-correct across
# the family of sleep predictors.

d_sleep_standardized <- d2 %>%
  group_by(subject) %>%
  mutate(
    REM = within_z(r_of_sleep_period_min),
    N1 = within_z(n1_of_sleep_period_min),
    N2 = within_z(n2_of_sleep_period_min),
    N3 = within_z(n3_of_sleep_period_min),
    TST = within_z(tst_corrected_min),
    SE = within_z(sleep_eff_corrected_perc),
    WASO = within_z(wake_after_sleep_onset_of_sleep_period_min),
    SOL = within_z(sleep_onset_delay_min)) %>%
  ungroup()

sleep_predictors <- c(
  "REM duration" = "REM",
  "N1 duration" = "N1",
  "N2 duration" = "N2",
  "N3 duration" = "N3",
  "Total sleep time" = "TST",
  "Sleep efficiency" = "SE",
  "Wake after sleep onset" = "WASO",
  "Sleep onset latency" = "SOL"
)

etable2 <- lapply(
  names(sleep_predictors),
  function(label) fit_sleep_predictor(label, sleep_predictors[[label]])) %>%
  bind_rows() %>%
  mutate(
    p_holm = p.adjust(p, method = "holm"))

print(etable2, n = Inf, width = Inf)


# =============================================================================
# 4. ANALYSIS 2 — COUPLING PREDICTS ANXIETY OUTCOME
# =============================================================================
# Estimate each subject's REM–distress coupling as a random slope, then use it
# to predict the post-treatment GAD outcome.

# Complete cases for the coupling model, restricted to subjects with >= 2
# nights so a within-subject slope is identifiable.
coupling_model_vars <- c(
  "stress_rating_m",
  "R_c",
  "E_c",
  "age",
  "RMT20_GAD_T0",
  "ISI_T0"
)

d_coupling <- d2 %>%
  filter(if_all(
    all_of(coupling_model_vars),
    ~ !is.na(.x)
  )) %>%
  group_by(subject) %>%
  filter(
    n() >= 2
  ) %>%
  ungroup()

# Random intercept and uncorrelated random slope of R_c per subject.
m_rem_coupling <- lmer(
  stress_rating_m ~
    R_c + E_c + age + RMT20_GAD_T0 + ISI_T0 +
    (1 + R_c || subject),
  data = d_coupling,
  REML = TRUE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)))

# Subject-level coupling slopes. 'rem_slope_sd' is the standardized slope with
# its sign flipped, so that higher values denote stronger (more negative)
# REM–distress coupling — the direction used as the predictor in Analysis 2.
coupling_info <- coef(m_rem_coupling)$subject %>%
  as.data.frame() %>%
  rownames_to_column("subject") %>%
  transmute(
    subject = as.character(subject),
    rem_slope = R_c,
    rem_slope_sd = as.numeric(-scale(rem_slope)))

# Subject-level variables carried into the outcome analysis.
subject_level_vars <- c(
  "condition",
  "RMT20_GAD_T0", "RMT20_GAD_T1",
  "RMT20_SAD_T0", "RMT20_SAD_T1",
  "RMT20_panic_disorder_T0", "RMT20_panic_disorder_T1",
  "RMT20_PTSD_T0", "RMT20_PTSD_T1",
  "RMT20_MDD_T0", "RMT20_MDD_T1",
  "BPD_C_T0", "BPD_C_T1",
  "ISI_T0", "ISI_T1", "ISI_change")

# One row per subject: baseline/outcome values joined to the coupling slope,
# with treatment condition effect-coded (-0.5 / +0.5) so that the
# main coupling-slope coefficient represents its average association
# with follow-up GAD across treatment conditions.
prediction_df <- d2 %>%
  mutate(subject = as.character(subject)) %>%
  group_by(subject) %>%
  mutate (condition = as.character(condition)) %>%
  summarise(
    across(all_of(subject_level_vars), first_nonmissing),
    .groups = "drop"
  ) %>%
  inner_join(coupling_info, by = "subject") %>%
  mutate(
    condition_ec = case_when(
      condition == "1" ~ -0.5,
      condition == "2"   ~  0.5,
      TRUE ~ NA_real_))

# Primary outcome: GAD at T1, adjusting for GAD at T0.
m_gad_outcome <- lm(
  RMT20_GAD_T1 ~
    rem_slope_sd * condition_ec +
    RMT20_GAD_T0,
  data = prediction_df)

summary(m_gad_outcome)
confint(m_gad_outcome)

dat_gad <- model.frame(m_gad_outcome) %>% as_tibble()


# =============================================================================
# 5. SECONDARY ANXIETY DIMENSIONS
# =============================================================================
# Repeat the Analysis 2 model for the secondary outcomes and Holm-correct the
# secondary main effects and interactions as two separate families. GAD is the
# primary outcome and is therefore excluded from both corrections.

outcomes <- list(
  GAD = c("RMT20_GAD_T1", "RMT20_GAD_T0"),
  SAD = c("RMT20_SAD_T1", "RMT20_SAD_T0"),
  PD= c("RMT20_panic_disorder_T1", "RMT20_panic_disorder_T0"),
  PTSD = c("RMT20_PTSD_T1", "RMT20_PTSD_T0"),
  BPD = c("BPD_C_T1", "BPD_C_T0"))


etable3 <- lapply(
  names(outcomes),
  function(outcome) {
    fit_outcome_model(
      outcome = outcome,
      variables = outcomes[[outcome]]
    )}) %>%
  bind_rows()


# GAD is primary, so it is not included in the Holm correction
etable3$p_holm <- NA_real_

secondary_rows <- etable3$Outcome != "GAD"

etable3$p_holm[secondary_rows] <- p.adjust(
  etable3$p[secondary_rows],
  method = "holm")

# Correct the secondary interaction tests as a separate family
etable3$interaction_p_holm <- NA_real_

etable3$interaction_p_holm[secondary_rows] <- p.adjust(
  etable3$interaction_p[secondary_rows],
  method = "holm")


print(etable3, n = Inf, width = Inf)


# =============================================================================
# 6. FIGURE 3
# =============================================================================
# Panel A: Analysis 1 fitted line.
# Panel B: per-subject slopes (spaghetti), coloured by slope direction. 
# Panel C: Analysis 2 fitted line over the subject-level points.

primary_rem_label <- format_effect_label(m_primary, "R_c")
analysis2_label <- format_effect_label(m_gad_outcome, "rem_slope_sd")

pred_primary <- ggpredict(m_primary, terms = "R_c [all]")
pred_gad <- ggpredict(m_gad_outcome, terms = "rem_slope_sd [all]")

# --- Panel A: REM duration -> morning distress --------------------------------
p_A <- ggplot(pred_primary, aes(x = x, y = predicted)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = GREEN,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(colour = GREEN, linewidth = 1) +
  richtext_annotation(
    x = 0,
    y = 5.5,
    label = primary_rem_label,
    colour = GREEN
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 50)) +
  scale_y_continuous(breaks = seq(0, 6, 2)) +
  coord_cartesian(xlim = c(-100, 100), ylim = c(0, 6)) +
  labs(
    title = "A",
    x = "REM duration (min)",
    y = "Morning distress"
  ) +
  theme_paper()

# --- Panel B: per-subject slopes, coloured by direction ----------------------
d_analysis1 <- model.frame(m_primary) %>%
  as_tibble()

analysis1_slopes <- d_analysis1 %>%
  group_by(subject) %>%
  filter(
    n() >= 2
  ) %>%
  summarise(
    slope = coef(lm(stress_rating_m ~ R_c))[["R_c"]],
    .groups = "drop"
  ) %>%
  mutate(
    slope_direction = factor(
      if_else(slope < 0, "Negative slope", "Zero or positive slope"),
      levels = c("Negative slope", "Zero or positive slope")
    )
  )

spaghetti_data <- d_analysis1 %>%
  inner_join(analysis1_slopes, by = "subject")

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
    legend.spacing.y = unit(0.1, "cm"),
    panel.grid.major.y = element_blank()
  )+
  scale_x_continuous(breaks = seq(-100, 100, 50)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(-100, 100), ylim = c(0, 10)) +
  labs(
    title = "B",
    x = "REM duration (min)",
    y = "Morning distress"
  )

# --- Panel C: coupling -> GAD severity at T1 ---------------------------------
p_C <- ggplot() +
  geom_point(
    data = dat_gad,
    aes(x = rem_slope_sd, y = RMT20_GAD_T1),
    alpha = 0.25,
    size = 1.8,
    colour = "darkgrey"
  ) +
  geom_ribbon(
    data = pred_gad,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = BLUE,
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    data = pred_gad,
    aes(x = x, y = predicted),
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
  scale_x_continuous(breaks = seq(-4, 4,2)) +
  scale_y_continuous(breaks = seq(0, 20, 4)) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(
    title = "C",
    x = "REM–distress coupling",
    y = "GAD severity at T1"
  ) +
  theme_paper()

# --- Assemble and export -----------------------------------------------------
figure3 <- p_A | p_B | p_C
figure3

ggsave(
  filename = "figure3_combined.png",
  plot = figure3,
  width = 14.5,
  height = 6,
  dpi = 300
)
