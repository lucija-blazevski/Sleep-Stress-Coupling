# =============================================================================
# Helper functions
# -----------------------------------------------------------------------------
# Project     : Within-person REM sleep–distress coupling and anxiety outcomes
# Authors     : Lucija Blazevski, Lara Rosler
# Contact     : l.blazevski@nin.knaw.nl
# -----------------------------------------------------------------------------
# Functions sourced by the analysis script. Sourcing this file checks
# that the required packages are installed, attaches them, and then defines
# helpers grouped as:
#
#   1. Subject IDs and missing-data reduction
#   2. Within-subject centering and standardization
#   3. Descriptive-table formatting
#   4. Model-term extraction and effect labels
#   5. Table 1 (baseline characteristics) builders
#   6. Model-fitting wrappers for the results tables
#   7. Shared ggplot2 theme and annotation
#
# R version   : R version 4.4.2
# =============================================================================

required_packages <- c(
  "dplyr",
  "stringr",
  "tibble",
  "ggplot2",
  "ggtext",
  "lme4",      
  "lmerTest",
  "tidyr",
  "patchwork",
  "tidyr",
  "ggeffects",
  "readr"
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
#' that is repeated across a subject's nightly rows but may contain NAs.
#'
#' @param x A vector.
#' @return The first non-missing element, or a typed NA if all are missing.
first_nonmissing <- function(x) {
  observed <- which(!is.na(x))
  
  if (length(observed) == 0L) {
    return(x[NA_integer_][1])
  }
  
  x[observed[1]]
}


#' Mean ignoring NAs, returning NA when every value is missing.
#'
#' @param x Numeric vector.
#' @return The mean of the non-missing values, or NA_real_ if none remain.
mean_or_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
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
    return(rep(NA_real_, length(x)))
  }
  
  x - mean(x, na.rm = TRUE)
}


#' Within-subject z-standardization.
#'
#' Centers and scales to a subject's own mean and SD. Returns all-NA when the
#' subject has fewer than two valid observations or no within-subject variance,
#' which keeps degenerate subjects out of the standardized models.
#'
#' @param x Numeric vector for one subject.
#' @return Standardized vector, or all-NA if standardization is undefined.
within_z <- function(x) {
  valid <- !is.na(x)
  
  if (sum(valid) < 2L || isTRUE(sd(x[valid]) == 0)) {
    return(rep(NA_real_, length(x)))
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
#' with the leading zero stripped per style conventions.
#'
#' @param p Numeric P-value.
#' @return A character string, or NA_character_ if `p` is NA.
format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < .001 ~ "<.001",
    p < .01  ~ sub("^0", "", sprintf("%.3f", p)),
    TRUE     ~ sub("^0", "", sprintf("%.2f", p))
  )
}


#' Summarize a diagnosis indicator as "n/N (percent)".
#'
#' Accepts Yes and No for present vs. absent respectively, ignore
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
#' mixed models, normal approximation for GLMs, profile for `lm`).
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
#' `richtext_annotation()`.
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
# 5. TABLE 1 (BASELINE CHARACTERISTICS) BUILDERS
# =============================================================================

#' Blank section header row for the baseline-characteristics table.
#'
#' @param label Row label to display in the Characteristic column.
#' @return A one-row tibble with empty value columns.
section_row <- function(label) {
  tibble(
    Characteristic = label,
    Total = "",
    `CBT-I` = "",
    Control = ""
  )
}

#' Summary row for the baseline-characteristics table.
#'
#' Applies `formatter` to the whole sample and to each treatment arm.
#' Expects `bl_all` and `bl_randomized` (with a `treatment` column) in the
#' calling environment.
#'
#' @param label Row label.
#' @param variable Name of the variable to summarize.
#' @param formatter A formatting function such as `format_mean_sd`.
#' @return A one-row tibble with Total, CBT-I, and Control columns.
summary_row <- function(label, variable, formatter) {
  tibble(
    Characteristic = label,
    Total = formatter(bl_all[[variable]]),
    `CBT-I` = formatter(
      bl_randomized[[variable]][bl_randomized$treatment == "CBT-I"]
    ),
    Control = formatter(
      bl_randomized[[variable]][bl_randomized$treatment == "Control"]
    )
  )
}


# =============================================================================
# 6. MODEL-FITTING WRAPPERS FOR THE RESULTS TABLES
# =============================================================================

#' Fit one within-subject sleep-predictor model (Analysis 1, supplement table).
#'
#' Regresses morning distress on a within-subject standardized sleep variable
#' plus the fixed covariates, with a random subject intercept. Fitted with
#' REML = FALSE so that models differing only in the sleep predictor are
#' comparable on the same (ML) footing.
#' Expects `d_sleep_standardized` in the calling environment.
#'
#' @param label Human-readable label for the sleep variable.
#' @param predictor Column name of the (standardized) predictor.
#' @return A one-row tibble: Sleep variable, n, b, SE, lower_95, upper_95, p.
fit_sleep_predictor <- function(label, predictor) {
  fit <- lmer(
    as.formula(
      sprintf(
        "stress_rating_m ~ %s + E_c + age + RMT20_GAD_T0 + ISI_T0 + (1 | subject)",
        predictor
      )
    ),
    data = d_sleep_standardized,
    REML = FALSE
  )
  
  result <- model_term(fit, predictor)
  
  tibble(
    `Sleep variable` = label,
    n = nrow(model.frame(fit)),
    b = result$estimate,
    SE = result$SE,
    lower_95 = result$lower_95,
    upper_95 = result$upper_95,
    p = result$p
  )
}

#' Fit one secondary-outcome model without the condition interaction.
#'
#' Regresses a T1 outcome on the standardized coupling slope, adjusting for the
#' matching T0 baseline, using only complete cases for the terms involved.
#' Expects `prediction_df` in the calling environment.
#'
#' @param outcome Label for the outcome.
#' @param variables Length-2 vector: c(T1 variable, T0 variable).
#' @return A one-row tibble: Outcome, n, b, SE, lower_95, upper_95, p.
fit_secondary_outcome <- function(outcome, variables) {
  t1_var <- variables[[1]]
  t0_var <- variables[[2]]
  
  d_outcome <- prediction_df %>%
    filter(
      if_all(
        all_of(c(t1_var, t0_var, "rem_slope_sd")),
        ~ !is.na(.x)
      )
    )
  
  fit <- lm(
    as.formula(sprintf("%s ~ rem_slope_sd + %s", t1_var, t0_var)),
    data = d_outcome
  )
  
  result <- model_term(fit, "rem_slope_sd")
  
  tibble(
    Outcome = outcome,
    n = nrow(d_outcome),
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