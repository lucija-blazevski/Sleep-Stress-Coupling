# =============================================================================
# Main-text analyses
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Reproduces the analyses reported in the main text of the manuscript:
#
#   Sample      Participants, nights, and number of nights per participant
#               (Methods).
#   Analysis 1  Within-person REM duration predicts next-morning distress.
#               The same mixed model provides participant-specific REM slopes.
#   Analysis 2  Participant-specific REM–distress slopes predict GAD severity
#               at 2-month follow-up, adjusting for baseline GAD severity.
#   Figure 2    Panels A–C summarizing Analyses 1 and 2.
#
# R version   : 4.4.2
# =============================================================================


# Data, packages and helper functions
source("00_data.R")


# Figure palette
BLUE <- "#312783"
GREEN <- "#0F8A2E"


# =============================================================================
# 1. SAMPLE
# =============================================================================

analysis1_sample_check <- d_analysis1 %>%
  summarise(
    participants = n_distinct(subject),
    nights = n()
  )

analysis1_sample_check


# Number of eligible nights per participant (Methods, Procedure)
subject_df %>%
  summarise(
    mean = mean(n_nights),
    sd = sd(n_nights),
    median = median(n_nights),
    min = min(n_nights),
    max = max(n_nights)
  )

subject_df %>%
  count(n_nights) %>%
  mutate(percent = 100 * n / sum(n))


# Participants with at least one MINI-confirmed diagnosis (Methods)
subject_df %>%
  summarise(
    n = sum(any_clinical_diagnosis == "yes", na.rm = TRUE),
    percent = 100 * n / n()
  )


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
# 3. ANALYSIS 2 — REM-RELATED DISTRESS ALLEVIATION AND FOLLOW-UP GAD
# =============================================================================

# Higher rem_slope_sd values indicate stronger REM-related overnight distress
# alleviation because the original slopes are reverse-coded after
# standardization.
coupling_info <- person_slopes(m_primary, "R_c")


# One row per Analysis 1 participant, joined to the individualized REM slope.
prediction_df <- subject_df %>%
  inner_join(coupling_info, by = "subject")


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
# 4. FIGURE 2
# =============================================================================

# Panel A: overnight change in distress (morning distress − evening distress)
#          against the night's deviation from the participant's average REM
#          duration (hours). Points are nights; the line is the fixed REM
#          effect of the primary model (Analysis 1), shown per hour of REM,
#          with its 95% CI.
# Panel B: participant-specific (unadjusted) linear fits of the overnight
#          change in distress on REM deviation, split by slope direction.
# Panel C: change in GAD severity (GAD at T1 − GAD at T0) against
#          REM-related distress alleviation, i.e. the participant-specific
#          REM slope of the primary model reversed in sign and expressed per
#          hour (decrease in next-morning distress per additional hour of REM
#          sleep). Points are participants; the line is the Analysis 2 model
#          prediction at average baseline GAD severity, averaged over
#          treatment conditions, with its 95% CI.

fig_data <- d_analysis1 %>%
  mutate(
    subject = as.character(subject),
    rem_dev_h = R_c / 60,
    distress_change = stress_rating_m - stress_rating_e
  )

zero_line <- geom_hline(
  yintercept = 0,
  linetype = "dashed",
  colour = "grey60",
  linewidth = 0.4
)


# -----------------------------------------------------------------------------
# Panel A — REM duration and overnight change in distress
# -----------------------------------------------------------------------------

b_rem <- fixef(m_primary)[["R_c"]] * 60                 # per hour of REM
se_rem <- sqrt(vcov(m_primary)["R_c", "R_c"]) * 60

line_A <- tibble(
  rem_dev_h = seq(
    min(fig_data$rem_dev_h),
    max(fig_data$rem_dev_h),
    length.out = 200
  )
) %>%
  mutate(
    fit = b_rem * rem_dev_h,
    lower = pmin((b_rem - 1.96 * se_rem) * rem_dev_h,
                 (b_rem + 1.96 * se_rem) * rem_dev_h),
    upper = pmax((b_rem - 1.96 * se_rem) * rem_dev_h,
                 (b_rem + 1.96 * se_rem) * rem_dev_h)
  )

p_A <- ggplot() +
  zero_line +
  geom_point(
    data = fig_data,
    aes(x = rem_dev_h, y = distress_change),
    colour = "grey60",
    alpha = 0.22,
    size = 1.3
  ) +
  geom_ribbon(
    data = line_A,
    aes(x = rem_dev_h, ymin = lower, ymax = upper),
    fill = GREEN,
    alpha = 0.15
  ) +
  geom_line(
    data = line_A,
    aes(x = rem_dev_h, y = fit),
    colour = GREEN,
    linewidth = 1
  ) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  scale_y_continuous(breaks = seq(-8, 8, 4)) +
  expand_limits(y = 8) +
  labs(
    title = "A",
    x = "Deviation from individual average\nREM duration (h)",
    y = "Overnight change in distress\n(morning − evening)"
  ) +
  theme_paper()


# -----------------------------------------------------------------------------
# Panel B — participant-level variability in REM–distress slopes
# -----------------------------------------------------------------------------

participant_slopes <- fig_data %>%
  group_by(subject) %>%
  summarise(
    slope = coef(lm(distress_change ~ rem_dev_h))[["rem_dev_h"]],
    .groups = "drop"
  ) %>%
  mutate(
    slope_direction = factor(
      if_else(slope < 0, "Negative slope", "Zero or positive slope"),
      levels = c("Negative slope", "Zero or positive slope")
    )
  )

spaghetti_data <- fig_data %>%
  inner_join(participant_slopes, by = "subject")

p_B <- ggplot(
  spaghetti_data,
  aes(
    x = rem_dev_h,
    y = distress_change,
    group = subject,
    colour = slope_direction
  )
) +
  zero_line +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.35,
    alpha = 0.6
  ) +
  facet_wrap(~ slope_direction, ncol = 1) +
  scale_colour_manual(
    values = c(
      "Negative slope" = GREEN,
      "Zero or positive slope" = "grey60"
    ),
    guide = "none"
  ) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  scale_y_continuous(breaks = seq(-8, 8, 4)) +
  expand_limits(y = c(-5, 8)) +
  labs(
    title = "B",
    x = "Deviation from individual average\nREM duration (h)",
    y = "Overnight change in distress\n(morning − evening)"
  ) +
  theme_paper() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )


# -----------------------------------------------------------------------------
# Panel C — REM-related distress alleviation and change in GAD severity
# -----------------------------------------------------------------------------

# REM-related distress alleviation per hour of REM sleep: the participant's
# REM slope from the primary model (change in morning distress per minute),
# multiplied by −60. Its standardized version is `rem_slope_sd`, the
# predictor in Analysis 2, so the Analysis 2 prediction can be mapped onto
# this scale.
alleviation_h <- coef(m_primary)$subject %>%
  rownames_to_column("subject") %>%
  transmute(subject = as.character(subject), alleviation_h = -60 * R_c)

alleviation_mean <- mean(alleviation_h$alleviation_h)
alleviation_sd <- sd(alleviation_h$alleviation_h)

dat_C <- dat_gad %>%
  mutate(
    alleviation_h = alleviation_mean + rem_slope_sd * alleviation_sd,
    gad_change = RMT20_GAD_T1 - RMT20_GAD_T0
  )

gad_t0_mean <- mean(dat_C$RMT20_GAD_T0)

line_C <- tibble(
  alleviation_h = seq(
    min(dat_C$alleviation_h),
    max(dat_C$alleviation_h),
    length.out = 200
  )
) %>%
  mutate(
    rem_slope_sd = (alleviation_h - alleviation_mean) / alleviation_sd,
    condition_ec = 0,
    RMT20_GAD_T0 = gad_t0_mean
  )

pred_C <- predict(m_gad_outcome, newdata = line_C, interval = "confidence")

line_C <- line_C %>%
  mutate(
    fit = pred_C[, "fit"] - gad_t0_mean,
    lower = pred_C[, "lwr"] - gad_t0_mean,
    upper = pred_C[, "upr"] - gad_t0_mean
  )

p_C <- ggplot() +
  zero_line +
  geom_point(
    data = dat_C,
    aes(x = alleviation_h, y = gad_change),
    colour = "grey60",
    alpha = 0.22,
    size = 1.5
  ) +
  geom_ribbon(
    data = line_C,
    aes(x = alleviation_h, ymin = lower, ymax = upper),
    fill = BLUE,
    alpha = 0.12
  ) +
  geom_line(
    data = line_C,
    aes(x = alleviation_h, y = fit),
    colour = BLUE,
    linewidth = 1
  ) +
  scale_y_continuous(breaks = seq(-8, 8, 2)) +
  labs(
    title = "C",
    x = "REM-related distress alleviation",
    y = "Change in GAD severity\n(T1 − T0)"
  ) +
  theme_paper() +
  theme(plot.margin = margin(5.5, 25, 5.5, 5.5))


# -----------------------------------------------------------------------------
# Assemble and export
# -----------------------------------------------------------------------------

figure2 <- p_A | p_B | p_C
figure2


ggsave(
  filename = "Figure 2.png",
  plot = figure2,
  width = 16,
  height = 6,
  dpi = 600
)
