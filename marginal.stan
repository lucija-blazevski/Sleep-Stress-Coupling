data {
  // Total number of nightly observations
  int<lower=1> N;
  // Number of participants
  int<lower=1> J;

  // Number of usable nights for each participant
  array[J] int<lower=1> n_nights;
  
  // Row where each participant's nightly observations begin
  // in the dataset, which is sorted by subject
  array[J] int<lower=1> start;
  
  // Whether the participant has complete anxiety outcome data
  // 1 = yes; 0 = no
  array[J] int<lower=0, upper=1> has_gad;

 // Nightly variables
  // R: within-person-centered REM duration
  // E: within-person-centered evening distress
  // stress: morning distress
  vector[N] R;
  vector[N] E;
  vector[N] stress;

// Participant-level outcome variables
  // gad: follow-up GAD severity
  // gad0: baseline GAD severity, centered
  // condition: treatment condition, coded -0.5 / +0.5
  vector[J] gad;
  vector[J] gad0;
  vector[J] condition;

  // Observed means used to center the intercept priors
  real mean_stress;
  real mean_gad;
}

parameters {
// ============================================================
  // 2. POPULATION-LEVEL PARAMETERS
  // ============================================================

  // ---------- Nightly distress model ----------

  // Population intercept for morning distress
  real alpha_stress;
  // Average within-person REM-morning distress association
  real beta_R;
  // Average within-person evening-morning distress association
  real beta_E;
  
  // Unexplained nightly variation in morning distress
  real<lower=0> sigma_stress;
  // Between-participant variability in distress intercepts
  real<lower=0> sigma_intercept;
   // Between-participant variability in REM slopes
  // This parameter is particularly important because it
  // represents heterogeneity in individual REM coupling.
  real<lower=0> tau_R;


  // ---------- Follow-up anxiety model ----------

  // Population intercept for follow-up GAD
  real alpha_gad;
  // Association of baseline GAD with follow-up GAD
  real beta_baseline;
  // Treatment-condition effect on follow-up GAD
  real beta_condition;
   // Association between latent REM slope deviation and GAD
  // Expressed per population SD of REM slope heterogeneity.
  //
  // IMPORTANT SIGN CONVENTION:
  // The GAD model uses -beta_rem * z_slope.
  // Therefore, a NEGATIVE beta_rem corresponds to the
  // NEGATIVE association reported in the original analysis,
  // where the predictor was the negative of the REM slope.
  real beta_rem;
  // Unexplained variation in follow-up GAD after accounting
  // for baseline GAD, condition, and latent REM coupling
  real<lower=0> sigma_gad;
}

model {

// ============================================================
  // 3. PRIOR DISTRIBUTIONS
  // ============================================================

  // ---------- Nightly distress component ----------
  alpha_stress ~ normal(mean_stress, 5);
  beta_R ~ normal(0, 0.05);
  beta_E ~ normal(0, 1);

  sigma_stress ~ student_t(3, 0, 5);
  sigma_intercept ~ student_t(3, 0, 3);
  tau_R ~ normal(0, 0.02);

  // ---------- Follow-up anxiety component ----------

  alpha_gad ~ normal(mean_gad, 5);
  beta_baseline ~ normal(0, 1);
  beta_condition ~ normal(0, 3);
  beta_rem ~ normal(0, 2);
  sigma_gad ~ student_t(3, 0, 3);

  // ============================================================
  // 4. JOINT LIKELIHOOD: ONE PARTICIPANT AT A TIME
  // ============================================================
  for (j in 1:J) {
    int n = n_nights[j];
    int k = n + has_gad[j];
    // Starting row of participant j in nightly data
    int s = start[j];

    // Observed outcomes: morning distress and potentially GAD
    vector[k] y;
    // Expected values of these outcomes
    vector[k] mu;
    // Covariance matrix:
    // describes dependencies among nightly distress observations
    // and between nightly distress and follow-up GAD
    matrix[k, k] Sigma;

    // ----------------------------------------------------------
    // 4A. EXPECTED NIGHTLY DISTRESS
    // ----------------------------------------------------------
    for (a in 1:n) {
      // Row index for the current night
      int ia = s + a - 1;

      // Actual morning distress
      y[a] = stress[ia];
      
       // Expected morning distress from population-level effects
      //
      // Individual random intercepts and slopes are not
      // included explicitly here because they are integrated
      // out through the covariance matrix below
      mu[a] = alpha_stress +
              beta_R * R[ia] +
              beta_E * E[ia];
          
      // --------------------------------------------------------
      // 4B. COVARIANCE BETWEEN TWO NIGHTS
      // --------------------------------------------------------

      for (b in 1:n) {
        int ib = s + b - 1;

        // Shared random intercept:
        // generates covariance between all nights of a person.
        //
        // Shared random REM slope:
        // generates additional covariance depending on
        // REM duration on the two nights.
        //
        // Assumes random intercept and slope are independent.        
        Sigma[a, b] =
          square(sigma_intercept) +
          square(tau_R) * R[ia] * R[ib];
          
        // Add nightly residual variance only to diagonal cells
        // because independent residual errors do not create
        // covariance between different nights.
        if (a == b)
          Sigma[a, b] += square(sigma_stress);
      }
    }

// ----------------------------------------------------------
    // 4C. INCLUDE FOLLOW-UP GAD WHEN AVAILABLE
    // ----------------------------------------------------------
    if (has_gad[j] == 1) {
      // The final entry contains participant j's follow-up GAD
      y[k] = gad[j];

  // Expected follow-up GAD, adjusted for baseline and condition
      //
      // The population-average latent slope deviation is zero,
      // so no REM-slope term appears in the marginal mean.
      mu[k] = alpha_gad +
              beta_baseline * gad0[j] +
              beta_condition * condition[j];


      // --------------------------------------------------------
      // 4D. COVARIANCE BETWEEN NIGHTLY DISTRESS AND GAD
      // --------------------------------------------------------
      for (a in 1:n) {
        int ia = s + a - 1;
        // This is the key link in the entire joint model.
        //
        // If individual REM slope deviations predict GAD,
        // nightly distress and GAD become correlated in a way
        // that depends on person-centered REM duration.
        //
        // The covariance equals:
        //
        // - beta_rem * tau_R * R[ia]
        //
        // This integrates over each participant's uncertain
        // REM slope rather than using a fixed slope estimate.
        Sigma[a, k] = -beta_rem * tau_R * R[ia];
        // Covariance matrices must be symmetric

        Sigma[k, a] = Sigma[a, k];
      }

      // --------------------------------------------------------
      // 4E. TOTAL VARIANCE OF FOLLOW-UP GAD
      // --------------------------------------------------------

      // Follow-up GAD varies because of:
      //
      // 1. Differences in latent REM coupling:
      //    variance = beta_rem^2
      //
      // 2. Other unexplained anxiety variation:
      //    variance = sigma_gad^2
      //
      // Both contribute to total GAD variability.
      Sigma[k, k] =
        square(beta_rem) + square(sigma_gad);
    }
 // ==========================================================
    // 5. EVALUATE THE PARTICIPANT'S JOINT LIKELIHOOD
    // ==========================================================

    // The observed outcomes are jointly normally distributed
    // with expected values mu and covariance matrix Sigma.
    //
    // This likelihood integrates over the individual random
    // intercept and REM slope analytically.
    //
    // Cholesky decomposition is a numerically efficient and
    // stable way to evaluate the multivariate normal likelihood.
    target += multi_normal_cholesky_lpdf(
      y | mu, cholesky_decompose(Sigma)
    );
  }
}