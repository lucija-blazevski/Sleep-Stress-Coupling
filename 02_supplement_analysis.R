# =============================================================================
# Supplementary analyses
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Reproduces the analyses reported in Supplement 1, in the order of the
# supplement:
#
#   Table S3    Sleep stages entered simultaneously (standardized), and total
#               sleep time in a separate model.
#   Table S4    Complete Analysis 2 model (GAD severity at follow-up).
#   Table S5    Average REM duration as predictor of GAD severity.
#   Table S6    Correlates of REM-related distress alleviation.
#   Table S7    Secondary symptom dimensions (Holm-corrected).
#   Table S8    Sleep-stage–related distress alleviation as predictor of GAD
#               severity (N1, N2, N3, TST; REM, N1, N2 and N3 jointly).
#   Text        Sensitivity analysis: REM percentage instead of REM duration.
#   Text        Sensitivity analysis: participants with a MINI-confirmed
#               diagnosis.
#   Table S9    Sensitivity to the number of nights per participant.
#   Table S10   Joint Bayesian model (marginal.stan; requires cmdstanr).
#
# Runs 01_main_analysis.R first (data, primary models, Figure 2).
#
# R version   : 4.4.2
# =============================================================================


source("01_main_analysis.R")


# =============================================================================
# TABLE S3 — SLEEP STAGES ENTERED SIMULTANEOUSLY
# =============================================================================

# N1, N2, N3 and REM duration in one model of next-morning distress, adjusted
# for evening distress, with a random intercept per participant. Each
# within-person centered stage is divided by its SD across all nights and
# morning distress is z-scored, so coefficients are standardized and
# comparable across stages.

d_std <- d_analysis1 %>%
  mutate(
    across(
      c(R_c, N1_c, N2_c, N3_c, TST_c, E_c),
      ~ .x / sd(.x),
      .names = "{.col}_z"
    ),
    distress_m_z = as.numeric(scale(stress_rating_m))
  )

m_all_stages_std <- lmer(
  distress_m_z ~ R_c_z + N1_c_z + N2_c_z + N3_c_z + E_c_z + (1 | subject),
  data = d_std,
  REML = TRUE
)

summary(m_all_stages_std)


# Total sleep time is the sum of the four stages, so it cannot be entered in
# the same model; it is tested on its own with the same specification.
stopifnot(isTRUE(all.equal(
  d_analysis1$tst_corrected_min,
  with(
    d_analysis1,
    n1_of_sleep_period_min + n2_of_sleep_period_min +
      n3_of_sleep_period_min + r_of_sleep_period_min
  )
)))

m_tst_std <- lmer(
  distress_m_z ~ TST_c_z + E_c_z + (1 | subject),
  data = d_std,
  REML = TRUE
)

tableS3 <- bind_rows(
  lapply(
    c(
      "REM duration" = "R_c_z",
      "N1 duration" = "N1_c_z",
      "N2 duration" = "N2_c_z",
      "N3 duration" = "N3_c_z"
    ),
    function(term) model_term(m_all_stages_std, term)
  ),
  .id = "Sleep stage"
) %>%
  bind_rows(
    model_term(m_tst_std, "TST_c_z") %>%
      mutate(`Sleep stage` = "Total sleep time (separate model)")
  ) %>%
  select(`Sleep stage`, beta = estimate, lower_95, upper_95, p)

print(tableS3, n = Inf, width = Inf)


# Within-person correlation between REM duration and total sleep time
cor(d_analysis1$R_c, d_analysis1$TST_c)


# =============================================================================
# TABLE S4 — COMPLETE ANALYSIS 2 MODEL
# =============================================================================

tableS4 <- bind_rows(
  lapply(
    c(
      "REM-related overnight distress alleviation, per 1 SD" = "rem_slope_sd",
      "Treatment condition (CBT-I vs control)" = "condition_ec",
      "Baseline GAD severity" = "RMT20_GAD_T0",
      "REM-related distress alleviation x treatment condition" =
        "rem_slope_sd:condition_ec"
    ),
    function(term) model_term(m_gad_outcome, term)
  ),
  .id = "Variable"
) %>%
  select(Variable, b = estimate, lower_95, upper_95, p)

print(tableS4, n = Inf, width = Inf)


# =============================================================================
# TABLE S5 — AVERAGE REM DURATION AS PREDICTOR OF GAD SEVERITY
# =============================================================================

# Same model as Analysis 2, with the participant's average REM duration
# (minutes, across eligible nights) instead of REM-related distress
# alleviation.

m_sens_gad <- lm(
  RMT20_GAD_T1 ~
    mean_rem_min * condition_ec +
    RMT20_GAD_T0,
  data = prediction_df
)

tableS5 <- bind_rows(
  lapply(
    c(
      "Average REM duration at baseline" = "mean_rem_min",
      "Treatment condition (CBT-I vs control)" = "condition_ec",
      "Baseline GAD severity" = "RMT20_GAD_T0",
      "REM duration x treatment condition" = "mean_rem_min:condition_ec"
    ),
    function(term) model_term(m_sens_gad, term)
  ),
  .id = "Variable"
) %>%
  select(Variable, b = estimate, lower_95, upper_95, p)

print(tableS5, n = Inf, width = Inf)


# =============================================================================
# TABLE S6 — CORRELATES OF REM-RELATED DISTRESS ALLEVIATION
# =============================================================================

# Bivariate associations of the participant-specific REM slopes (Analysis 2
# predictor) with age, baseline insomnia severity (ISI) and baseline GAD
# severity, and the difference between female and male participants.

tableS6_correlates <- bind_rows(
  lapply(
    c(
      "Age" = "age",
      "Baseline insomnia severity (ISI)" = "ISI_T0",
      "Baseline GAD severity (RMT20)" = "RMT20_GAD_T0"
    ),
    function(variable) {
      test <- cor.test(prediction_df$rem_slope_sd, prediction_df[[variable]])
      tibble(
        n = sum(complete.cases(prediction_df[, c("rem_slope_sd", variable)])),
        r = unname(test$estimate),
        lower_95 = test$conf.int[1],
        upper_95 = test$conf.int[2],
        p = test$p.value
      )
    }
  ),
  .id = "Characteristic"
)

print(tableS6_correlates, n = Inf, width = Inf)


# Sex: Welch t test (female minus male)
prediction_df %>%
  group_by(female) %>%
  summarise(n = n(), mean_slope_sd = mean(rem_slope_sd))

t.test(rem_slope_sd ~ female, data = prediction_df)


# =============================================================================
# TABLE S7 — SECONDARY SYMPTOM DIMENSIONS
# =============================================================================

# Analysis 2 model repeated for the 4 secondary outcomes. Main effects and
# interactions are Holm-corrected as separate families of 4 tests.

secondary_outcomes <- list(
  SAD = c("RMT20_SAD_T1", "RMT20_SAD_T0"),
  PD = c("RMT20_panic_disorder_T1", "RMT20_panic_disorder_T0"),
  PTSD = c("RMT20_PTSD_T1", "RMT20_PTSD_T0"),
  BPD = c("BPD_C_T1", "BPD_C_T0")
)

tableS7_full <- lapply(
  names(secondary_outcomes),
  function(outcome) {
    fit_outcome_model(
      outcome = outcome,
      variables = secondary_outcomes[[outcome]]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    p_holm = p.adjust(p, method = "holm"),
    interaction_p_holm = p.adjust(interaction_p, method = "holm")
  )

tableS7 <- tableS7_full %>%
  select(Outcome, n, b, lower_95, upper_95, p, p_holm)

print(tableS7, n = Inf, width = Inf)


# Interactions with treatment condition (reported in the supplement text)
tableS7_full %>%
  select(
    Outcome, interaction_b, interaction_lower_95, interaction_upper_95,
    interaction_p, interaction_p_holm
  ) %>%
  print(n = Inf, width = Inf)


# =============================================================================
# TABLE S8 — SLEEP-STAGE–RELATED DISTRESS ALLEVIATION AND GAD SEVERITY
# =============================================================================

# Participant-specific slopes for N1, N2, N3 and total sleep time are derived
# with the same model as the REM slopes, standardized and reverse-coded
# (higher = stronger distress alleviation).

stage_terms <- c(N1 = "N1_c", N2 = "N2_c", N3 = "N3_c", TST = "TST_c")

stage_slopes <- lapply(names(stage_terms), function(label) {
  person_slopes(
    fit_primary_type(stage_terms[[label]]),
    stage_terms[[label]],
    name = paste0(label, "_slope_sd")
  )
}) %>%
  Reduce(function(a, b) full_join(a, b, by = "subject"), .)

stage_df <- prediction_df %>%
  inner_join(stage_slopes, by = "subject")


# Separate models: same specification as Analysis 2
slope_labels <- c(
  "REM-related alleviation" = "rem_slope_sd",
  "N1-related alleviation" = "N1_slope_sd",
  "N2-related alleviation" = "N2_slope_sd",
  "N3-related alleviation" = "N3_slope_sd",
  "TST-related alleviation" = "TST_slope_sd"
)

tableS8_separate <- bind_rows(
  lapply(slope_labels, function(slope) {
    fit <- fit_gad_model(slope, stage_df)
    model_term(fit, slope) %>% mutate(n = nobs(fit))
  }),
  .id = "Predictor"
) %>%
  mutate(Model = "Separate models")


# REM, N1, N2 and N3 entered simultaneously 
m_rem_stages <- lm(
  RMT20_GAD_T1 ~ rem_slope_sd + N1_slope_sd + N2_slope_sd + N3_slope_sd +
    condition_ec + RMT20_GAD_T0,
  data = stage_df
)

# Do N1, N2 and N3 jointly improve the model beyond REM?
anova(
  update(m_rem_stages, . ~ . - N1_slope_sd - N2_slope_sd - N3_slope_sd),
  m_rem_stages
)


joint_terms <- function(model, terms, label) {
  bind_rows(
    lapply(terms, function(term) model_term(model, term)),
    .id = "Predictor"
  ) %>%
    mutate(n = nobs(model), Model = label)
}

tableS8 <- bind_rows(
  tableS8_separate,
  joint_terms(
    m_rem_stages,
    slope_labels[c(
      "REM-related alleviation", "N1-related alleviation",
      "N2-related alleviation", "N3-related alleviation"
    )],
    "REM, N1, N2, and N3 entered simultaneously"
  )
) %>%
  select(Model, Predictor, n, b = estimate, lower_95, upper_95, p)

print(tableS8, n = Inf, width = Inf)


# Correlations between the REM, N1, N2, N3 and TST slopes (supplement text and
# response to reviewers)
round(cor(stage_df[, c("rem_slope_sd", "N1_slope_sd", "N2_slope_sd", "N3_slope_sd", "TST_slope_sd")]), 2)


# =============================================================================
# SENSITIVITY ANALYSIS — REM PERCENTAGE (TEXT)
# =============================================================================

# Both primary analyses with REM sleep as a percentage of total sleep time
# (within-person centered). Participant-specific slopes are re-derived from
# this model.

m_primary_rperc <- fit_primary_type("RPERC_c")

model_term(m_primary_rperc, "RPERC_c") %>%
  mutate(nights = nobs(m_primary_rperc),
         participants = ngrps(m_primary_rperc)[["subject"]])

rperc_df <- subject_df %>%
  inner_join(person_slopes(m_primary_rperc, "RPERC_c"), by = "subject")

m_gad_rperc <- fit_gad_model("rem_slope_sd", rperc_df)

model_term(m_gad_rperc, "rem_slope_sd") %>%
  mutate(n = nobs(m_gad_rperc))


# =============================================================================
# SENSITIVITY ANALYSIS — MINI-CONFIRMED DIAGNOSIS (TEXT)
# =============================================================================

# Eligibility was based on questionnaire scores, not the MINI. Both primary
# analyses are repeated in participants with at least one MINI-confirmed
# diagnosis (GAD, SAD, PD or PTSD). Participant-specific slopes are re-derived
# within this subsample.

mini_subjects <- subject_df %>%
  filter(any_clinical_diagnosis == "yes") %>%
  pull(subject)

d_mini <- d_analysis1 %>%
  filter(as.character(subject) %in% mini_subjects)

m_primary_mini <- fit_primary_type("R_c", data = d_mini)

model_term(m_primary_mini, "R_c") %>%
  mutate(nights = nobs(m_primary_mini),
         participants = ngrps(m_primary_mini)[["subject"]])

mini_df <- subject_df %>%
  inner_join(person_slopes(m_primary_mini, "R_c"), by = "subject")

m_gad_mini <- fit_gad_model("rem_slope_sd", mini_df)

model_term(m_gad_mini, "rem_slope_sd") %>%
  mutate(n = nobs(m_gad_mini))


# =============================================================================
# TABLE S9 — SENSITIVITY TO THE NUMBER OF NIGHTS
# =============================================================================

# Both primary analyses restricted to participants with at least 3 and with
# 4 (the maximum) eligible nights. The slope model is refitted within each
# subset; the procedure is otherwise identical to the main analysis.

nights_sensitivity <- function(min_nights) {
  keep <- subject_df %>% filter(n_nights >= min_nights) %>% pull(subject)
  d_sub <- d_analysis1 %>% filter(as.character(subject) %in% keep)

  m_night <- fit_primary_type("R_c", data = d_sub)
  sub_df <- subject_df %>%
    inner_join(person_slopes(m_night, "R_c"), by = "subject")
  m_gad <- fit_gad_model("rem_slope_sd", sub_df)

  night <- model_term(m_night, "R_c")
  gad <- model_term(m_gad, "rem_slope_sd")

  tibble(
    Sample = paste0(">= ", min_nights, " nights"),
    participants = ngrps(m_night)[["subject"]],
    nights = nobs(m_night),
    night_b = night$estimate,
    night_lower_95 = night$lower_95,
    night_upper_95 = night$upper_95,
    night_p = night$p,
    outcome_n = nobs(m_gad),
    b = gad$estimate,
    lower_95 = gad$lower_95,
    upper_95 = gad$upper_95,
    p = gad$p
  )
}

full_sample <- tibble(
  Sample = "All (2-4 nights)",
  participants = ngrps(m_primary)[["subject"]],
  nights = nobs(m_primary),
  outcome_n = nobs(m_gad_outcome),
  b = model_term(m_gad_outcome, "rem_slope_sd")$estimate,
  lower_95 = model_term(m_gad_outcome, "rem_slope_sd")$lower_95,
  upper_95 = model_term(m_gad_outcome, "rem_slope_sd")$upper_95,
  p = model_term(m_gad_outcome, "rem_slope_sd")$p
)

tableS9 <- bind_rows(
  full_sample,
  nights_sensitivity(3),
  nights_sensitivity(4)
)

print(tableS9, n = Inf, width = Inf)


# =============================================================================
# TABLE S10 — JOINT BAYESIAN MODEL
# =============================================================================

# Robustness check for Analysis 2. Nightly morning distress and GAD severity at
# follow-up are modelled jointly, so each participant's REM–distress slope is a
# latent variable instead of an estimate that is treated as known. The model
# is defined in marginal.stan. Participants without follow-up GAD contribute 
# only their nightly data.
# A negative beta_rem means that stronger REM-related distress alleviation is
# associated with lower follow-up GAD severity.
#
# Requires cmdstanr with CmdStan installed (https://mc-stan.org/cmdstanr/).
# Sampling takes several minutes; set run_bayesian_model to FALSE to skip it.

run_bayesian_model <- TRUE

if (run_bayesian_model) {

  for (pkg in c("cmdstanr", "posterior")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Install the ", pkg, " package to run the joint Bayesian model.")
    }
  }

  # Nightly observations, sorted by participant
  d_bayes <- d_analysis1 %>%
    mutate(subject = as.character(subject)) %>%
    arrange(subject)

  # One row per participant, including those without GAD follow-up
  bayes_df <- subject_df %>%
    arrange(subject) %>%
    mutate(
      has_gad = as.integer(
        complete.cases(RMT20_GAD_T0, RMT20_GAD_T1, condition_ec)
      ),
      # row where each participant's nightly observations begin
      start = cumsum(lag(n_nights, default = 0)) + 1,
      gad0_c = RMT20_GAD_T0 - mean(RMT20_GAD_T0[has_gad == 1]),
      gad0_c = if_else(has_gad == 1, gad0_c, 0),
      gad = if_else(has_gad == 1, RMT20_GAD_T1, 0),
      condition_bayes = if_else(has_gad == 1, condition_ec, 0)
    )

  stan_data <- list(
    N = nrow(d_bayes),
    J = nrow(bayes_df),
    n_nights = bayes_df$n_nights,
    start = bayes_df$start,
    has_gad = bayes_df$has_gad,

    R = d_bayes$R_c,
    E = d_bayes$E_c,
    stress = d_bayes$stress_rating_m,

    gad = bayes_df$gad,
    gad0 = bayes_df$gad0_c,
    condition = bayes_df$condition_bayes,

    mean_stress = mean(d_bayes$stress_rating_m),
    mean_gad = mean(bayes_df$gad[bayes_df$has_gad == 1])
  )

  # Check sample and ordering
  stopifnot(
    sum(stan_data$n_nights) == stan_data$N,
    identical(unique(d_bayes$subject), bayes_df$subject),
    !anyNA(unlist(stan_data))
  )

  c(
    participants = stan_data$J,
    nights = stan_data$N,
    participants_with_GAD = sum(stan_data$has_gad)
  )

  joint_model <- cmdstanr::cmdstan_model("marginal.stan")

  fit_joint <- joint_model$sample(
    data = stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.99,
    max_treedepth = 13
  )

  print(fit_joint$summary(
    variables = c(
      "beta_rem",
      "beta_R",
      "tau_R",
      "sigma_gad",
      "beta_baseline",
      "beta_condition"
    )
  ), n = Inf, width = Inf)

  print(fit_joint$diagnostic_summary())

  # Posterior of the REM association
  beta_draws <- posterior::as_draws_df(fit_joint$draws("beta_rem"))$beta_rem

  tableS10 <- tibble(
    Parameter = "REM-related distress alleviation (beta_rem)",
    median = median(beta_draws),
    lower_95 = unname(quantile(beta_draws, 0.025)),
    upper_95 = unname(quantile(beta_draws, 0.975)),
    prob_negative = mean(beta_draws < 0)
  )

  print(tableS10, width = Inf)
}


# =============================================================================
# EXTRAS — DO THE DISTRESS RATINGS REFLECT INSOMNIA SEVERITY?
# =============================================================================

d_isi <- d_analysis1 %>%
  mutate(subject = as.character(subject)) %>%
  select(-any_of("ISI_T0")) %>%
  inner_join(subject_df %>% select(subject, ISI_T0), by = "subject") %>%
  mutate(ISI_T0_z = as.numeric(scale(ISI_T0)))

m_isi_moderation <- lmer(
  stress_rating_m ~ R_c * ISI_T0_z + E_c + (1 + R_c || subject),
  data = d_isi,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

model_term(m_isi_moderation, "R_c:ISI_T0_z")

d_distress_long <- d_analysis1 %>%
  mutate(subject = as.character(subject)) %>%
  select(subject, ses_num, stress_rating_e, stress_rating_m) %>%
  pivot_longer(
    c(stress_rating_e, stress_rating_m),
    names_to = "time",
    values_to = "distress"
  ) %>%
  mutate(
    time = factor(
      if_else(time == "stress_rating_m", "morning", "evening"),
      levels = c("evening", "morning")
    )
  ) %>%
  left_join(subject_df %>% select(subject, ISI_T0), by = "subject")

contrasts(d_distress_long$time) <- c(-0.5, 0.5)

m_distress_isi <- lmer(
  distress ~ time * ISI_T0 + (1 | subject),
  data = d_distress_long,
  REML = TRUE
)

anova(m_distress_isi)
summary(m_distress_isi)
confint(m_distress_isi, parm = c("ISI_T0", "time1:ISI_T0"), method = "Wald")
