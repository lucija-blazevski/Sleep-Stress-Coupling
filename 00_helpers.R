# =============================================================================
# Helper functions
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Functions sourced by the analysis script. Sourcing this file checks
# that the required packages are installed, attaches them, and then defines
# helpers:
#
#   1. Subject IDs and missing-data reduction
#   2. Within-subject centering and standardization
#   3. Descriptive-table formatting
#   4. Model-term extraction and effect labels
#   5. Descriptive-table builders
#   6. Model-fitting wrappers for the results tables
#   7. Shared ggplot2 theme and annotation
#
# R version   : R version 4.4.2
# =============================================================================

required_packages <- c(
  "dplyr",
  "stringr",
  "tibble",
  "readr",
  "tidyr",
  "ggplot2",
  "ggtext",
  "ggeffects",
  "patchwork",
  "lme4",
  "lmerTest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before using this helper file: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(ggeffects)
  library(patchwork)
  library(lme4)
  library(lmerTest)
})


# =============================================================================
# 1. SUBJECT IDS AND MISSING-DATA REDUCTION
# =============================================================================

#' Standardize NSR subject identifiers.
#'
#' Extracts the numeric part of an identifier and left-pads it to five digits,
#' so that IDs from different source formats line up for merging.
#'
#' @param x Vector coercible to character containing subject IDs.
#' @return Character vector of 5-digit, zero-padded IDs.
standardize_nsr_id <- function(x) {
  x %>%
    as.character() %>%
    str_extract("\\d+") %>%
    str_pad(width = 5, side = "left", pad = "0")
}


#' First non-missing value of a vector.
#'
#' Used to collapse a subject-level constant (e.g. condition, baseline score)
#' that is repeated across a subject's nightly rows.
#'
#' @param x A vector.
#' @return The first non-missing element, or an NA if all are missing.
first_nonmissing <- function(x) {
  observed <- which(!is.na(x))
  
  if (length(observed) == 0L) {
    return(x[NA][1])
  }
  
  x[observed[1]]
}


#' Mean ignoring NAs, returning NA when every value is missing.
#'
#' @param x Numeric vector.
#' @return The mean of the non-missing values, or NA_real_ if none remain.
mean_or_na <- function(x) {
  if (all(is.na(x))) {
    NA
  } else {
    mean(x, na.rm = TRUE)
  }
}


# =============================================================================
# 2. WITHIN-SUBJECT CENTERING AND STANDARDIZATION
# =============================================================================

#' Within-subject mean-centering.
#'
#' Subtracts a subject's own mean, isolating night-to-night deviations from
#' their typical value. Apply inside `group_by(subject)`.
#'
#' @param x Numeric vector for one subject.
#' @return Mean-centered vector, or all-NA if the subject has no valid values.
within_center <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  
  x - mean(x, na.rm = TRUE)
}


#' Within-subject z-standardization.
#'
#' Centers and scales to a subject's own mean and SD. Returns all-NA when the
#' subject has fewer than two valid observations or no within-subject variance,
#' used only for the standardized, supplement models.
#'
#' @param x Numeric vector for one subject.
#' @return Standardized vector, or all-NA if standardization is undefined.
within_z <- function(x) {
  valid <- !is.na(x)
  
  if (sum(valid) < 2 || isTRUE(sd(x[valid]) == 0)) {
    return(rep(NA, length(x)))
  }
  
  as.numeric(scale(x))
}


# =============================================================================
# 3. DESCRIPTIVE-TABLE FORMATTING
# =============================================================================

#' Format a numeric vector as "mean (SD)" with one decimal.
#'
#' @param x Numeric vector.
#' @return A single character string, or NA_character_ if all values missing.
format_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}


#' Format a binary vector as "n (percent)".
#'
#' @param x Vector of category codes.
#' @param positive_value Value counted as positive (default 1).
#' @return A single character string, or NA_character_ if all values missing.
format_n_percent <- function(x, positive_value = 1) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  n_positive <- sum(x == positive_value)
  
  sprintf(
    "%d (%.1f)",
    n_positive,
    100 * n_positive / length(x)
  )
}


#' Format a P-value for reporting.
#'
#' Uses "<.001" below that threshold, three decimals below .01, otherwise two,
#' with the leading zero stripped per style conventions. Returns a blank string
#' when the P-value is unavailable, which is convenient for manuscript tables.
#'
#' @param p Numeric P-value.
#' @return A formatted character string, or "" if `p` is NA.
format_p <- function(p) {
  if (is.na(p)) {
    return("")
  }
  
  if (p < .001) {
    return("<.001")
  }
  
  if (p < .01) {
    return(sub("^0", "", sprintf("%.3f", p)))
  }
  
  sub("^0", "", sprintf("%.2f", p))
}


#' Summarize a diagnosis indicator as "n/N (percent)".
#'
#' Accepts yes and no for present vs. absent, respectively, and ignores
#' anything it cannot classify (which is excluded from the denominator).
#'
#' @param x Vector of diagnosis indicators.
#' @return A "n/N (pct%)" string, or "NA" if nothing is classifiable.
format_diagnosis <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  
  diagnosis_present <- x %in% c("yes")
  
  diagnosis_absent <- x %in% c("no")
  
  available <- diagnosis_present | diagnosis_absent
  denominator <- sum(available)
  numerator <- sum(diagnosis_present)
  
  if (denominator == 0) {
    return("NA")
  }
  
  sprintf(
    "%d/%d (%.1f%%)",
    numerator,
    denominator,
    100 * numerator / denominator
  )
}


# =============================================================================
# 4. MODEL-TERM EXTRACTION AND EFFECT LABELS
# =============================================================================

#' Extract a tidy summary of a single model term.
#'
#' Returns the estimate, SE, df, test statistic, P-value and 95% CI for one
#' coefficient, choosing an appropriate CI method per model class (Wald for
#' mixed models, normal approximation for GLMs, and `confint()` for `lm`).
#'
#' @param model A fitted model (`lm`, `glm`, or `merMod`).
#' @param term Name of the coefficient to extract.
#' @return A one-row tibble: term, estimate, SE, df, statistic, p, lower_95,
#'   upper_95.
model_term <- function(model, term) {
  coefficients <- coef(summary(model))
  
  if (!term %in% rownames(coefficients)) {
    stop(
      "Term '", term, "' was not found in the model. Available terms are: ",
      paste(rownames(coefficients), collapse = ", ")
    )
  }
  
  estimate <- unname(coefficients[term, "Estimate"])
  se <- unname(coefficients[term, "Std. Error"])
  
  ci <- if (inherits(model, "merMod")) {
    confint(model, parm = term, method = "Wald")
  } else if (inherits(model, "glm")) {
    c(
      estimate - 1.96 * se,
      estimate + 1.96 * se
    )
  } else {
    confint(model, parm = term, level = 0.95)
  }
  
  p_column <- grep(
    "^Pr",
    colnames(coefficients),
    value = TRUE
  )[1]
  
  if (is.na(p_column)) {
    stop("No P-value column was found in the model summary.")
  }
  
  tibble(
    term = term,
    estimate = estimate,
    SE = se,
    df = if ("df" %in% colnames(coefficients)) {
      unname(coefficients[term, "df"])
    } else {
      NA_real_
    },
    statistic = if ("t value" %in% colnames(coefficients)) {
      unname(coefficients[term, "t value"])
    } else if ("z value" %in% colnames(coefficients)) {
      unname(coefficients[term, "z value"])
    } else {
      NA_real_
    },
    p = unname(coefficients[term, p_column]),
    lower_95 = unname(ci[1]),
    upper_95 = unname(ci[2])
  )
}


#' Build a nice effect label for a model term.
#'
#' Produces an italicised "b = ...; P ..." string with HTML markup for use in
#' `richtext_annotation()` for plotting.
#'
#' @param model A fitted model.
#' @param term Name of the coefficient to label.
#' @param digits Decimal places for the estimate (default 3).
#' @return A character string containing HTML/markdown.
format_effect_label <- function(model, term, digits = 3) {
  result <- model_term(model, term)
  
  p_text <- if (result$p < .001) {
    "&lt; .001"
  } else {
    paste0("= ", format_p(result$p))
  }
  
  paste0(
    "<i>b</i> = ",
    sprintf(paste0("%.", digits, "f"), result$estimate),
    "; <i>P</i> ",
    p_text
  )
}


# =============================================================================
# 5. DESCRIPTIVE-TABLE BUILDERS
# =============================================================================

# -----------------------------------------------------------------------------
# Overall-sample rows
# -----------------------------------------------------------------------------

#' Blank section-header row for an overall-sample table.
#'
#' @param label Row label to display in the Characteristic column.
#' @return A one-row tibble with an empty Total column.
section_total_row <- function(label) {
  tibble(
    Characteristic = label,
    Total = ""
  )
}


#' Summary row for an overall-sample table.
#'
#' Applies a formatting function to one variable in the participant-level
#' baseline dataset. Expects `bl_all` in the calling environment.
#'
#' @param label Row label.
#' @param variable Name of the variable to summarize.
#' @param formatter Formatting function such as `format_mean_sd`,
#'   `format_n_percent`, or `format_diagnosis`.
#' @return A one-row tibble with Characteristic and Total columns.
total_summary_row <- function(label, variable, formatter) {
  tibble(
    Characteristic = label,
    Total = formatter(bl_all[[variable]])
  )
}


# -----------------------------------------------------------------------------
# Between-treatment tests
# -----------------------------------------------------------------------------

#' Welch two-sample t-test between treatment groups.
#'
#' Uses non-missing observations from participants with a treatment assignment.
#' Expects `bl_all` with a `treatment` column in the calling environment.
#'
#' @param variable Name of a continuous variable.
#' @return Numeric P-value, or NA_real_ if the comparison cannot be performed.
welch_p <- function(variable) {
  
  dat <- bl_all %>%
    select(treatment, all_of(variable)) %>%
    filter(
      !is.na(treatment),
      !is.na(.data[[variable]])
    )
  
  if (
    nrow(dat) == 0L ||
    n_distinct(dat$treatment) < 2L
  ) {
    return(NA_real_)
  }
  
  t.test(
    dat[[variable]] ~ dat$treatment,
    var.equal = FALSE
  )$p.value
}


#' Pearson chi-square test between treatment groups.
#'
#' Uses non-missing observations from participants with a treatment assignment.
#' Optionally restricts the analysis to a specified set of valid category
#' values. Expects `bl_all` with a `treatment` column in the calling environment.
#'
#' @param variable Name of a categorical variable.
#' @param valid_values Optional vector of category values to retain.
#' @return Numeric P-value, or NA_real_ if the comparison cannot be performed.
chisq_p <- function(variable, valid_values = NULL) {
  
  dat <- bl_all %>%
    select(treatment, all_of(variable)) %>%
    filter(
      !is.na(treatment),
      !is.na(.data[[variable]])
    )
  
  if (!is.null(valid_values)) {
    dat <- dat %>%
      filter(.data[[variable]] %in% valid_values)
  }
  
  tab <- table(
    dat[[variable]],
    dat$treatment
  )
  
  if (
    nrow(tab) < 2 ||
    ncol(tab) < 2
  ) {
    return(NA_real_)
  }
  
  suppressWarnings(
    chisq.test(
      tab,
      correct = FALSE
    )$p.value
  )
}


# -----------------------------------------------------------------------------
# Treatment-group rows
# -----------------------------------------------------------------------------

#' Blank section-header row for a treatment-group table.
#'
#' @param label Row label to display in the Characteristic column.
#' @return A one-row tibble with empty treatment and P-value columns.
section_group_row <- function(label) {
  tibble(
    Characteristic = label,
    `CBT-I` = "",
    Control = "",
    `P Value` = ""
  )
}


#' Continuous-variable row by treatment condition.
#'
#' Reports mean (SD) within each treatment group and the P-value from a
#' Welch two-sample t-test. Expects `bl_all` with a `treatment` column in the
#' calling environment.
#'
#' @param label Row label.
#' @param variable Name of the continuous variable.
#' @return A one-row tibble with treatment summaries and P-value.
continuous_group_row <- function(label, variable) {
  
  tibble(
    Characteristic = label,
    `CBT-I` = format_mean_sd(
      bl_all[[variable]][bl_all$treatment == "CBT-I"]
    ),
    Control = format_mean_sd(
      bl_all[[variable]][bl_all$treatment == "Control"]
    ),
    `P Value` = format_p(
      welch_p(variable)
    )
  )
}


#' Binary-variable row by treatment condition.
#'
#' Reports number (percentage) within each treatment group and the P-value
#' from a Pearson chi-square test. Expects `bl_all` with a `treatment` column
#' in the calling environment.
#'
#' @param label Row label.
#' @param variable Name of the binary variable.
#' @return A one-row tibble with treatment summaries and P-value.
binary_group_row <- function(label, variable) {
  
  tibble(
    Characteristic = label,
    `CBT-I` = format_n_percent(
      bl_all[[variable]][bl_all$treatment == "CBT-I"]
    ),
    Control = format_n_percent(
      bl_all[[variable]][bl_all$treatment == "Control"]
    ),
    `P Value` = format_p(
      chisq_p(variable)
    )
  )
}


#' Diagnosis row by treatment condition.
#'
#' Reports diagnosis frequency within each treatment group using
#' `format_diagnosis()` and the P-value from a Pearson chi-square test.
#' Expects `bl_all` with a `treatment` column in the calling environment.
#'
#' @param label Row label.
#' @param variable Name of the diagnosis variable.
#' @return A one-row tibble with treatment summaries and P-value.
diagnosis_group_row <- function(label, variable) {
  
  tibble(
    Characteristic = label,
    `CBT-I` = format_diagnosis(
      bl_all[[variable]][bl_all$treatment == "CBT-I"]
    ),
    Control = format_diagnosis(
      bl_all[[variable]][bl_all$treatment == "Control"]
    ),
    `P Value` = format_p(
      chisq_p(
        variable,
        valid_values = c("yes", "no")
      )
    )
  )
}


# =============================================================================
# 6. MODEL-FITTING WRAPPERS FOR THE RESULTS TABLES
# =============================================================================

#' Check whether a mixed-effects model converged successfully.
#'
#' @param model A fitted `merMod` object or an error returned by `tryCatch()`.
#' @return Logical scalar indicating whether the model converged.
model_converged <- function(model) {
  
  if (!inherits(model, "merMod")) {
    return(FALSE)
  }
  
  opt_code <- model@optinfo$conv$opt
  conv_messages <- model@optinfo$conv$lme4$messages
  
  opt_ok <- is.null(opt_code) ||
    length(opt_code) == 0L ||
    all(opt_code == 0)
  
  messages_ok <- is.null(conv_messages)
  
  opt_ok && messages_ok
}


#' Fit one within-subject sleep-predictor model (Analysis 1 supplement).
#'
#' Fits morning distress as a function of a within-subject standardized sleep
#' predictor and within-subject centered evening distress. The analytic sample
#' is defined separately for each sleep predictor and requires at least two
#' usable observations per participant.
#'
#' A model with an uncorrelated random intercept and random slope for the sleep
#' predictor is attempted first. This model is retained only if it converges
#' and is nonsingular. Otherwise, a random-intercept-only model is fitted.
#'
#' Expects `d_analysis1` in the calling environment and uses `within_z()` and
#' `model_term()` defined above.
#'
#' @param label Human-readable label for the sleep variable.
#' @param predictor Column name of the sleep predictor.
#' @return A one-row tibble containing the sleep variable, participant and
#'   night counts, random-effects structure, random-slope fit status,
#'   coefficient, SE, 95% CI, and P-value.
fit_sleep_predictor <- function(label, predictor) {
  
  # ---------------------------------------------------------------------------
  # Define model-specific analytic sample
  # ---------------------------------------------------------------------------
  
  dat <- d_analysis1 %>%
    filter(
      !is.na(.data[[predictor]]),
      !is.na(stress_rating_m),
      !is.na(E_c)
    ) %>%
    group_by(subject) %>%
    filter(n() >= 2) %>%
    mutate(
      sleep_predictor_z = within_z(.data[[predictor]])
    ) %>%
    ungroup() %>%
    filter(is.finite(sleep_predictor_z))
  
  
  if (nrow(dat) == 0) {
    stop(
      "No usable observations remain for sleep predictor '",
      predictor,
      "'."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Attempt random intercept + uncorrelated random slope
  # ---------------------------------------------------------------------------
  
  random_slope_model <- tryCatch(
    lmer(
      stress_rating_m ~
        sleep_predictor_z + E_c +
        (1 + sleep_predictor_z || subject),
      data = dat,
      REML = TRUE
    ),
    error = function(e) e
  )
  
  
  random_slope_converged <- model_converged(random_slope_model)
  
  random_slope_singular <- if (random_slope_converged) {
    isSingular(random_slope_model, tol = 1e-4)
  } else {
    NA
  }
  
  
  # ---------------------------------------------------------------------------
  # Retain random slope only if converged AND nonsingular
  # ---------------------------------------------------------------------------
  
  if (
    random_slope_converged &&
    !random_slope_singular
  ) {
    
    model <- random_slope_model
    random_effects <- "Random intercept + slope"
    random_slope_status <- "Converged, nonsingular"
    
  } else {
    
    model <- lmer(
      stress_rating_m ~
        sleep_predictor_z + E_c +
        (1 | subject),
      data = dat,
      REML = TRUE
    )
    
    random_effects <- "Random intercept only"
    
    random_slope_status <- case_when(
      !random_slope_converged ~
        "Random-slope model did not converge",
      isTRUE(random_slope_singular) ~
        "Random-slope model singular",
      TRUE ~
        "Random-slope model not retained"
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Extract sleep-predictor effect
  # ---------------------------------------------------------------------------
  
  result <- model_term(
    model,
    "sleep_predictor_z"
  )
  
  
  tibble(
    `Sleep variable` = label,
    participants = n_distinct(model.frame(model)$subject),
    nights = nobs(model),
    `Random-effects structure` = random_effects,
    `Random-slope status` = random_slope_status,
    b = result$estimate,
    SE = result$SE,
    lower_95 = result$lower_95,
    upper_95 = result$upper_95,
    p = result$p
  )
}

#' Fit one outcome model with the coupling x condition interaction.
#'
#' Regresses a T1 outcome on standardized coupling, effect-coded condition,
#' their interaction, and the matching T0 baseline. Returns both the overall
#' coupling association and the between-condition difference (interaction).
#' Expects `prediction_df` in the calling environment.
#'
#' @param outcome Label for the outcome.
#' @param variables Length-2 vector: c(T1 variable, T0 variable).
#' @return A one-row tibble with main-effect and interaction estimates, SEs,
#'   CIs, and P-values.
fit_outcome_model <- function(outcome, variables) {
  
  outcome_t1 <- variables[1]
  outcome_t0 <- variables[2]
  
  model_formula <- as.formula(
    paste0(
      outcome_t1,
      " ~ rem_slope_sd * condition_ec + ",
      outcome_t0
    )
  )
  
  model <- lm(
    model_formula,
    data = prediction_df
  )
  
  model_summary <- summary(model)$coefficients
  model_ci <- confint(model)
  
  tibble(
    Outcome = outcome,
    n = nobs(model),
    
    # Overall coupling association across conditions
    b = model_summary["rem_slope_sd", "Estimate"],
    SE = model_summary["rem_slope_sd", "Std. Error"],
    lower_95 = model_ci["rem_slope_sd", 1],
    upper_95 = model_ci["rem_slope_sd", 2],
    p = model_summary["rem_slope_sd", "Pr(>|t|)"],
    
    # Difference in coupling association between conditions
    interaction_b =
      model_summary["rem_slope_sd:condition_ec", "Estimate"],
    interaction_SE =
      model_summary["rem_slope_sd:condition_ec", "Std. Error"],
    interaction_lower_95 =
      model_ci["rem_slope_sd:condition_ec", 1],
    interaction_upper_95 =
      model_ci["rem_slope_sd:condition_ec", 2],
    interaction_p =
      model_summary["rem_slope_sd:condition_ec", "Pr(>|t|)"]
  )
}


# =============================================================================
# 7. SHARED GGPLOT2 THEME AND ANNOTATION
# =============================================================================

#' Boxed pretty annotation for effect labels on a plot.
#'
#' @param x,y Data coordinates for the label.
#' @param label Rich-text (HTML/markdown) string, e.g. from
#'   `format_effect_label()`.
#' @param colour Text and border colour.
#' @param hjust Horizontal justification (default 0.5).
#' @return A ggplot2 annotation layer.
richtext_annotation <- function(
    x,
    y,
    label,
    colour,
    hjust = 0.5
) {
  annotate(
    "richtext",
    x = x,
    y = y,
    label = label,
    fill = "white",
    colour = colour,
    label.size = 0.5,
    label.r = grid::unit(0.15, "lines"),
    label.padding = grid::unit(0.25, "lines"),
    size = 5.5,
    fontface = "bold",
    hjust = hjust
  )
}


#' Shared minimal theme for the manuscript figures.
#'
#' @param base_size Base font size (default 18).
#' @return A ggplot2 theme object.
theme_paper <- function(base_size = 18) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        colour = "grey88",
        linewidth = 0.4
      ),
      axis.line = element_line(
        colour = "grey70",
        linewidth = 0.4
      ),
      axis.ticks = element_line(
        colour = "grey70",
        linewidth = 0.4
      ),
      axis.text = element_text(
        colour = "grey25",
        size = base_size
      ),
      axis.title = element_text(
        colour = "black",
        size = base_size + 2
      ),
      plot.title = element_text(
        hjust = -0.1,
        size = base_size + 2,
        face = "bold",
        margin = margin(b = 20)
      ),
      legend.text = element_text(size = base_size),
      legend.key.width = grid::unit(1.1, "cm"),
      plot.margin = margin(10, 18, 10, 18)
    )
}