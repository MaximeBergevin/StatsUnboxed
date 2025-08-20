#' Perform a simulation-based sensitivity analysis for a linear mixed model.
#'
#' This function iteratively searches for the Minimum Detectable Effect Size (MDES)
#' by simulating new datasets based on a fitted LMM and testing for significance.
#' It can analyze either a main effect or an interaction effect.
#'
#' @param model         A fitted object of class `lmerMod` or `lmerModLmerTest`.
#' @param effect_type   A character string, either "main" or "interaction", specifying the type of effect to test.
#' @param variables     A character vector specifying the name(s) of the variable(s) involved in the effect.
#'                      For a main effect, e.g., "condition". For an interaction, e.g., c("group", "condition").
#' @param coef_name     A character string specifying the exact name of the fixed effect coefficient from the model's summary.
#' @param start_delta   A numeric value indicating the initial effect size (in beta units) to begin the search.
#' @param step_size     A numeric value indicating the size of the step to adjust the effect size between iterations.
#' @param min_step      A numeric value indicating the smallest step size to use before stopping the search.
#' @param tolerance     A numeric value indicating the acceptable distance from the target power of 0.80.
#' @param n_sim         A numeric value indicating the number of simulations to run per iteration.
#' @param max_iter      An integer specifying the maximum number of iterations to search for the MDES.
#' @param verbose       A logical value. If TRUE, progress is printed to the console.
#'

mixed.sensitivity <- function(
    model,             # The fitted lmer/lmerTest model to use as a blueprint
    effect_type,       # The type of effect to test: "main" or "interaction"
    variables,         # The variable(s) involved in the effect, e.g., "condition" or c("group", "condition")
    coef_name,         # The exact name of the fixed effect coefficient to test
    start_delta = 0,   # The initial effect size to test (in beta units)
    step_size = 0.05,  # The step size to adjust the effect size between iterations
    min_step = 0.0005, # The smallest step size to use before stopping
    tolerance = 0.02,  # The acceptable distance from the target power (0.80)
    n_sim = 200,       # The number of simulations to run per iteration
    max_iter = 30,     # The maximum number of iterations to search for the MDES
    verbose = TRUE     # Whether to print progress during the run
) {
  
  # Extract the model's essential components once for efficiency
  outcome <- names(model.frame(model))[1]
  base_data <- model.frame(model)
  base_beta <- fixef(model)[coef_name]
  formula <- formula(model)
  
  # Initialize the variables for the iterative search
  current_delta <- start_delta
  best_fit <- NULL
  closest_gap <- 1
  direction <- NULL
  
  # Outer loop: Iterates to adjust the effect size and find the MDES
  for (iter in seq_len(max_iter)) {
    sig_count <- 0
    
    # Inner loop: Runs the specified number of simulations for the current effect size
    for (i in seq_len(n_sim)) {
      # 1. Simulate a new dataset based on the model's variance components
      sim_data <- simulate(model, nsim = 1, seed = NULL)[[1]]
      
      # 2. Create a new data frame and replace the outcome variable with the simulated data
      df <- base_data
      df[[outcome]] <- sim_data
      
      # 3. Inject the effect into the simulated data
      if (effect_type == "main") {
        # Logic to inject a main effect using sum contrasts
        lvls <- levels(df[[variables]])
        df[[outcome]] <- df[[outcome]] +
          ifelse(df[[variables]] == lvls[2], -0.5 * current_delta,
                 ifelse(df[[variables]] == lvls[1], 0.5 * current_delta, 0))
        
      } else if (effect_type == "interaction") {
        # Logic to inject an interaction effect using nested ifelse statements for sum contrasts
        lvls1 <- levels(df[[variables[1]]])
        lvls2 <- levels(df[[variables[2]]])
        df[[outcome]] <- df[[outcome]] +
          ifelse(df[[variables[1]]] == lvls1[2] & df[[variables[2]]] == lvls2[2], -0.5 * current_delta,
                 ifelse(df[[variables[1]]] == lvls1[2] & df[[variables[2]]] == lvls2[1], 0.5 * current_delta,
                        ifelse(df[[variables[1]]] == lvls1[1] & df[[variables[2]]] == lvls2[2], 0.5 * current_delta,
                               ifelse(df[[variables[1]]] == lvls1[1] & df[[variables[2]]] == lvls2[1], -0.5 * current_delta, 0))))
      }
      
      # 4. Refit the model to the new, modified dataset and handle potential errors
      fit <- suppressWarnings(
        tryCatch(
          lmerTest::lmer(formula, data = df, control = lmerControl(check.conv.grad = .makeCC("ignore", tol = 1e-3))),
          error = function(e) NULL
        )
      )
      
      # 5. If the model successfully fits, check if the p-value is significant
      if (!is.null(fit)) {
        p <- summary(fit)$coefficients[coef_name, "Pr(>|t|)"]
        if (!is.na(p) && p < 0.05) {
          sig_count <- sig_count + 1
        }
      }
    }
    
    # Calculate the power (proportion of significant simulations) and the MDES
    power <- sig_count / n_sim
    total_b <- base_beta + current_delta
    
    # Check for convergence and adjust the search direction and step size
    if (is.null(direction)) {
      direction <- if (power >= 0.80) -1 else 1
      last_direction <- direction
    }
    if (verbose) cat(sprintf("Testing effect size: %.6f | Estimated power: %.3f \n", total_b, power))
    
    gap <- abs(power - 0.80)
    if (gap < closest_gap) {
      closest_gap <- gap
      best_fit <- total_b
    }
    
    if (gap <= tolerance) {
      cat(sprintf("✅ MDES ≈ %.3f (β)\n", total_b))
      return(total_b)
    }
    
    new_direction <- if (power >= 0.80) -1 else 1
    if (!is.na(last_direction) && new_direction != last_direction) {
      step_size <- step_size / 2
    }
    current_delta <- current_delta + new_direction * step_size
    
    if (step_size < min_step) break
  }
  
  # Final output if the search does not converge
  if (!is.null(best_fit)) {
    cat(sprintf("⚠️ Closest MDES ≈ %.3f (β), but did not converge within threshold.\n", best_fit))
    return(best_fit)
  } else {
    cat("❌ MDES not found within iteration or step size limits.\n")
    return(NULL)
  }
}