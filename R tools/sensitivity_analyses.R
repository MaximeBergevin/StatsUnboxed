## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## THE CENTRAL MOTOR COMMAND, BUT NO THE MUSCLE AFFERENT FEEDBACK, IS NECESSARY
## TO PERCEIVE EFFORT
## Script -- Sensitivity analyses
## 
## Author: Maxime Bergevin, MSc, PhD candidate 🧠
## School of kinesiology, Universite de Montreal, Canada
## Research center of Montreal's Geriatric Institutes, Canada
## Supervised by: 
## Benjamin Pageaux (Kinesiology, Universite de Montreal)
## Mathieu Roy (Psychology, McGill University)
## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Most packages are not explicitly loaded, save a few exceptions
# Every other package is explicit when calling their methods
library(tidyverse)      # Necessary for pipes (%>%)
library(lmerTest)       # Necessary for simr


# LOAD DATA ----
load(file = 'Data/presenceIsometric.Rdata')
load(file = 'Data/presenceDynamic.Rdata')
load(file = 'Data/magnitudeIsometric.Rdata')
load(file = 'Data/magnitudeDynamic.Rdata')



# SETITING FIXED EFFECT CONTRASTS ----
## Set R's default behavior to use treatment contrasts on lmer/rlmer
## Zero-sum contrasts are more appropriate in mixed designs
## See `An Introduction to Mixed Models for Experimental Psychology`
afex::set_sum_contrasts()



# SETTING POST-HOC CONTRASTS ----
# To be used with emmeans::emmeans to compute post-hoc analyses
# 1) Group comparisons for the same conditions (novice vs experienced)
# 2) Condition comparisons within the same group (VOL vs EMS; VOL vs EMS+VOL)
contr.int <- list(
  "Exp - Nov (VOL)"  = c(1,-1,0,0),     # VOL - Experienced vs novice
  "Exp - Nov (EMS)"  = c(0,0,1,-1),     # EMS - Experienced vs novice
  "VOL - EMS (Exp)"  = c(1,0,-1,0),     # Experienced - VOL vs EMS
  "VOL - EMS (Nov)"  = c(0,1,0,-1)      # Novice - VOL vs EMS
)

mixed.sensitivity <- function(
    model,                 # fitted lmer/lmerTest model
    effect_type,           # "main" or "interaction"
    variables,             # e.g., "condition" or c("group", "condition")
    coef_name,             # name of fixed effect in summary(model)
    start_delta = 0,
    step = 0.05,
    min_step = 0.0005,
    tolerance = 0.02,
    n_sim = 200,
    max_iter = 30,
    verbose = TRUE
){
  outcome     <- names(model.frame(model))[1]
  formula     <- formula(model)
  base_data   <- model.frame(model)
  base_beta   <- fixef(model)[coef_name]
  current_delta <- start_delta
  direction     <- NULL
  last_direction <- NA
  step_size     <- step
  best_fit      <- NULL
  closest_gap   <- 1
  
  for (iter in seq_len(max_iter)) {
    sig_count <- 0
    
    for (i in seq_len(n_sim)) {
      sim_data <- simulate(model, nsim = 1)[[1]]
      df       <- base_data
      df[[outcome]] <- sim_data
      
      # Inject delta into fixed effect
      if (effect_type == "main") {
        lvls <- levels(df[[variables]])
        df[[outcome]] <- df[[outcome]] +
          ifelse(df[[variables]] == lvls[2], -0.5 * current_delta,
                 ifelse(df[[variables]] == lvls[1],  0.5 * current_delta, 0))
        
      } else if (effect_type == "interaction") {
        lvls1 <- levels(df[[variables[1]]])
        lvls2 <- levels(df[[variables[2]]])
        df[[outcome]] <- df[[outcome]] +
          ifelse(df[[variables[1]]] == lvls1[2] & df[[variables[2]]] == lvls2[2], -0.5 * current_delta,
                 ifelse(df[[variables[1]]] == lvls1[2] & df[[variables[2]]] == lvls2[1],  0.5 * current_delta,
                        ifelse(df[[variables[1]]] == lvls1[1] & df[[variables[2]]] == lvls2[2],  0.5 * current_delta,
                               ifelse(df[[variables[1]]] == lvls1[1] & df[[variables[2]]] == lvls2[1], -0.5 * current_delta, 0))))
      }
      
      fit <- suppressWarnings(
        tryCatch(
          lmerTest::lmer(formula, data = df, control = lmerControl(check.conv.grad = .makeCC("ignore", tol = 1e-3))),
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
    
    # Set direction only on first iteration
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


# Isometric ----
magnitudeIsometric.lmer <- lmerTest::lmer(
  data = magnitudeIsometric.long,
  formula = effort ~ group * condition + (condition|id)
)
summary(magnitudeIsometric.lmer)

isometric_sensitivity <- mixed.sensitivity(
  model = magnitudeIsometric.lmer, 
  effect_type = "main",
  variables = "condition",
  coef_name = "condition1",
  n_sim = 1000             
)

isometric_sensitivity_interaction <- mixed.sensitivity(
  model = magnitudeIsometric.lmer, 
  effect_type = "interaction",
  variables = c("group", "condition"),
  coef_name = "group1:condition1",
  n_sim = 1000,
  max_iter = 50
)


# Dynamic ----
# Isometric ----
magnitudeDynamic.lmer <- lmerTest::lmer(
  data = magnitudeDynamic.long,
  formula = effort ~ group * condition + (condition|id)
)
summary(magnitudeDynamic.lmer)

dynamic_sensitivity <- mixed.sensitivity(
  model = magnitudeDynamic.lmer, 
  effect_type = "main",
  variables = "condition",
  coef_name = "condition1",
  n_sim = 1000             
)

dynamic_sensitivity_interaction <- mixed.sensitivity(
  model = magnitudeDynamic.lmer, 
  effect_type = "interaction",
  variables = c("group", "condition"),
  coef_name = "group1:condition1",
  n_sim = 1000,
  max_iter = 50
)