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
    model,
    effect_type,
    variables,
    coef_name,
    start_delta = 0,
    step_size = 0.05,
    min_step = 0.0005,
    tolerance = 0.02,
    n_sim = 200,
    max_iter = 30,
    verbose = TRUE
) {
  # Get model components outside the loop for efficiency
  outcome <- names(model.frame(model))[1]
  base_data <- model.frame(model)
  base_beta <- fixef(model)[coef_name]
  formula <- formula(model)
  
  # Initialize iterative search variables
  current_delta <- start_delta
  best_fit <- NULL
  closest_gap <- 1
  direction <- NULL
  
  for (iter in seq_len(max_iter)) {
    sig_count <- 0
    
    for (i in seq_len(n_sim)) {
      # Simulate new data from the model
      sim_data <- simulate(model, nsim = 1, seed = NULL)[[1]]
      
      # Create a new data frame with simulated data
      df <- base_data
      df[[outcome]] <- sim_data
      
      # Inject the effect based on the effect type
      if (effect_type == "main") {
        df <- df %>%
          mutate(!!outcome := !!sym(outcome) + 
                   ifelse(!!sym(variables) == levels(!!sym(variables))[1], 0.5 * current_delta, 
                          -0.5 * current_delta))
        
      } else if (effect_type == "interaction") {
        df <- df %>%
          mutate(!!outcome := !!sym(outcome) + 
                   case_when(
                     !!sym(variables[1]) == levels(!!sym(variables[1]))[2] & !!sym(variables[2]) == levels(!!sym(variables[2]))[2] ~ -0.5 * current_delta,
                     !!sym(variables[1]) == levels(!!sym(variables[1]))[2] & !!sym(variables[2]) == levels(!!sym(variables[2]))[1] ~ 0.5 * current_delta,
                     !!sym(variables[1]) == levels(!!sym(variables[1]))[1] & !!sym(variables[2]) == levels(!!sym(variables[2]))[2] ~ 0.5 * current_delta,
                     TRUE ~ -0.5 * current_delta
                   ))
      }
      
      # Refit and check for significance
      fit <- suppressWarnings(
        tryCatch(
          lmerTest::lmer(formula, data = df),
          error = function(e) NULL
        )
      )
      
      if (!is.null(fit)) {
        p <- summary(fit)$coefficients[coef_name, "Pr(>|t|)"]
        if (!is.na(p) & p < 0.05) {
          sig_count <- sig_count + 1
        }
      }
    }
    
    # Rest of the convergence logic is the same...
    power <- sig_count / n_sim
    total_b <- base_beta + current_delta
    
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
  
  if (!is.null(best_fit)) {
    cat(sprintf("⚠️ Closest MDES ≈ %.3f (β), but did not converge within threshold.\n", best_fit))
    return(best_fit)
  } else {
    cat("❌ MDES not found within iteration or step size limits.\n")
    return(NULL)
  }
}