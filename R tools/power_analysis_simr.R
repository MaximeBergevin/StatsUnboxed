#' Perform a simulation-based power analysis for a linear mixed model using simr.
#'
#' This function calculates statistical power across a grid of participant
#' numbers and observation counts for a specific fixed effect. It uses the
#' `simr` package for stable simulation workflow.
#'
#' @param model        A fitted object of class `lmerMod` or `lmerModLmerTest`. This serves
#'                     as the template for the data structure and variance components.
#' @param cluster      Character string. The name of the random grouping variable
#'                     that represents the participants (e.g., "participant", "subject").
#'                     Must match the name in the formula.
#' @param effect_name  Character string. The name of the fixed effect to test as it
#'                     would appear in an ANOVA table (e.g., "forceLevel"). Used when `test_type` is 'anova'.
#' @param coef_name    Character string. The exact name of the coefficient in the
#'                     `fixef()` output to be modified by the SESOI (e.g., "forceLevelhighForce"). 
#'                     Used when `test_type` is 't-test'.
#' @param sesoi        Numeric. The Smallest Effect Size of Interest. This is the raw
#'                     (unstandardized) effect size you want to be able to detect.
#' @param n_range.     Numeric vector. A vector of participant sample sizes to simulate (e.g., `c(20, 30, 40, 50)`).
#' @param obs_range    Numeric vector. A vector of the number of
#'                     observations PER CONDITION to simulate for each participant.
#' @param test_type    Character string. Either 'anova' for an F-test of the overall
#'                     effect, or 't-test' for a t-test of a specific coefficient.
#' @param n_sim        Integer. The number of simulations to run per condition.
#'                     NOTE: Start with a small number (e.g., 100) for testing. Use 1000+ for final results.
#' @param alpha        Numeric. The significance level (Type I error rate).
#'
#' @return             A tibble (data frame) containing the power for each combination of
#'                     n_participants and n_observations.

power_analysis_simr <- function(
    model,
    cluster,
    effect_name,
    coef_name,
    sesoi,
    n_range,
    obs_range,
    test_type = 'anova',
    n_sim = 1000,
    alpha = 0.05
) {
  
  # Extract parameter and set up
  #' -----------------------------
  original_fixef <- lme4::fixef(model)
  
  if (!coef_name %in% names(original_fixef)) {
    stop(paste("Coefficient '", coef_name, "' not found in model's fixed effects.",
               "\nAvailable coefficients:", paste(names(original_fixef), collapse = ", ")))
  }
  
  # Create new fixed effects vector with SESOI injected
  fixed_with_sesoi <- original_fixef
  fixed_with_sesoi[coef_name] <- sesoi
  
  # Extract variance components
  variances <- lme4::VarCorr(model)
  model_sigma <- sigma(model)
  
  # Get a template of the experimental design (all predictor columns)
  outcome_var <- as.character(formula(model)[[2]])
  design_template <- model@frame %>%
    dplyr::select(-dplyr::all_of(c(outcome_var, cluster))) %>%
    dplyr::distinct()
  
  # Set up results storage
  results_list <- list()
  
  #  Main power sim loop using nested for loops
  #' -------------------------------------------
  for (n_participants in n_range) {
    for (n_obs_per_condition in obs_range) {
      
      message(paste("Simulating (simr): ", n_participants, " participants, ",
                    n_obs_per_condition, " observations per condition..."))
      
      # Create the specific data frame for this simulation run
      # This data frame has the desired N and obs count, but no outcome yet
      new_data <- tidyr::expand_grid(
        !!cluster := factor(1:n_participants),
        design_template
      ) %>%
        tidyr::uncount(n_obs_per_condition)
      
      # Create the hypothetical model for this specific design
      # makeLmer creates a model object with the desired SESOI and variance,
      # fitted to the structure of the hypothetical dataset.
      model_for_power_sim <- simr::makeLmer(
        formula = formula(model),
        fixef = fixed_with_sesoi,
        VarCorr = variances,
        sigma = model_sigma,
        data = new_data
      )
      
      # Define the statistical test to run
      if(test_type == 'anova') {
        test_to_run <- simr::fixed(effect_name, method = "anova")
      } else {
        test_to_run <- simr::fixed(coef_name, method = "t")
      }
      
      # Run the power simulation using simr
      # powerSim handles the entire "simulate > refit > test" loop internally.
      power_sim_results <- suppressMessages(
        simr::powerSim(
          model_for_power_sim,
          test = test_to_run,
          nsim = n_sim,
          alpha = alpha
          )
        )
      
      # Extract the power (mean of successful simulations)
      power <- summary(power_sim_results)$mean
      
      # Calculate the proportion of simulations with warnings
      warning_prop <- length(power_sim_results$warnings) / n_sim
      
      # Store results for this combination, including the warning proportion
      results_list[[length(results_list) + 1]] <- tibble::tibble(
        n_participants = n_participants,
        n_obs_per_condition = n_obs_per_condition,
        power = power,
        warning_prop = warning_prop # New column
      )
    }
  }
  
  # Combine all results into a single data frame
  power_results <- dplyr::bind_rows(results_list)
  
  message("✅ Simulation complete!")
  return(power_results)
}
