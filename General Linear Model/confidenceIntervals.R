library(tidyverse)
library(ggprism)

# athletes = blue; non-athletes = yellow
group.colors <- c(athletes = '#0a4c9e', nonAthletes = '#E69F00')


# DATA SIMULATION ----
  # Imagine a world where there are 1000 athletes and non-athletes
  # In this world, athletes have a pain threshold of 45 +/- 2
  # and non-athletes have a pain threshold of 42 +/- 2
  # Let's assume this is our entire population (i.e., big world)
  set.seed(123) # For reproducibility
  population.data <- tibble(
    group = rep(c('athletes', 'nonAthletes'), each = 1000),
    threshold = c(rnorm(1000, mean = 45, sd = 2),
                  rnorm(1000, mean = 42, sd = 2)
                  )
  )
  
  # The mean difference is 2.95 (so very close to 3)
  # It is not exactly 3 because of simulation variance (because we used rnorm)
  trueEffect <- population.data %>%
    group_by(group) %>%
    summarise(meanThreshold = mean(threshold)) %>%
    mutate(meanDifference = meanThreshold[1] - meanThreshold[2])

  # And if we visualize the data, we see
  population.plot <- ggplot(population.data, aes(x = threshold, fill = group)) +
    geom_density(alpha = 0.5, position = 'identity', adjust = 2) +
    scale_fill_manual(values = group.colors) +
    geom_vline(xintercept = trueEffect$meanThreshold[1], linetype = 'dotted') +
    geom_vline(xintercept = trueEffect$meanThreshold[2], linetype = 'dotted') +
    scale_y_continuous(limits = c(0, 0.20))

  ggsave(filename = 'General Linear Model/Figures/plot_ci_populationDensity.png',
         width = 6, height = 4.5, dpi = 'retina', bg = 'transparent')

  # Finally, we can observe this true effect by running a T-test...
  # The regression is -2.94733, reflecting the true difference between both groups
  summary(
    lm(data = population.data, threshold ~ group)
  )

  
# ONE SAMPLE ----
  set.seed(123)
  # Lets recruit some participants...
  # The `sample()` method here is similar to the actual recruitment process
  sample.data <- population.data %>%
    group_by(group) %>%
    sample_n(size = 30)
  # And now we perform the T-test
  sample.lm <- lm(data = sample.data, threshold ~ group)
  summary(sample.lm)
  confint(sample.lm)
  sample.emm <- data.frame(emmeans::emmeans(sample.lm, ~ group))
  
  # And now we plot :D
  sample.plot <- ggplot(sample.data , aes(group, y = threshold, color = group)) +
    # INDIVIDUAL DATA
    geom_point(alpha = 0.3, size = 2,
               position = position_dodge2(width = 0.2)) +
    # ESTIMATED MARGINAL MEANS
    geom_errorbar(aes(x = 0.8, 
                      ymin = sample.emm$lower.CL[1], ymax = sample.emm$upper.CL[1]),
                  width = 0.1, color = 'black') +
    geom_point(aes(x = 0.8, y = sample.emm$emmean[1]),
               size = 3, color = group.colors[1]) +
    geom_errorbar(aes(x = 2.2, 
                      ymin = sample.emm$lower.CL[2], ymax = sample.emm$upper.CL[2]),
                  width = 0.1, color = 'black') +
    geom_point(aes(x = 2.2, y = sample.emm$emmean[2]),
               size = 3, color = group.colors[2]) +
    # VANITY
    scale_color_manual(values = group.colors) +
    theme(
      legend.position = 'none'
    )
    
  ggsave(filename = 'General Linear Model/Figures/plot_ci_sampleDifference.png',
         width = 6, height = 4.5, dpi = 'retina', bg = 'transparent')
  

# MULTIPLE SAMPLES ----
  # Creates a function that samples numExperiment amount of time
  # and perform t.test on the resulting sample
  runExperiment <- function(population, group1, group2, sampleSize,
                            numExperiments, confidenceLevel) {
    experimentResults <- tibble()  
    
    trueEffect <- mean(population[population$group == group1,]$threshold) -
                  mean(population[population$group == group2,]$threshold)
    
    for(i in 1:numExperiments) {
      sample <- population %>%
        group_by(group) %>%
        sample_n(size = sampleSize)
      results <- t.test(threshold ~ group, data = sample, conf.level = confidenceLevel)
      ci <- results$conf.int
      
      experimentResults <- rbind(experimentResults, 
                                 tibble(experiment = i,
                                        lowerCi = ci[1],
                                        upperCi = ci[2],
                                        containsTrueEffect = (trueEffect >= ci[1] && trueEffect <= ci[2])
                                 ))
    }
    
    return(list(coverage = sum(experimentResults$containsTrueEffect) / numExperiments,
                data = experimentResults,
                trueEffect = trueEffect)) 
  }
  
  # Run the experiment (adjust parameters as needed)
  set.seed(123) # For reproducibility
  results <- runExperiment(population.data, "athletes", "nonAthletes", 500, 1000, 0.95) 
  
  # Plotting with ggplot2
  ggplot(results$data, aes(x = experiment)) +
    # TRUE EFFECT
    geom_hline(yintercept = results$trueEffect, linetype = "dashed", color = "black", linewidth = 1) +
    # INDIVIDUAL EXPERIMENT COVERAGE
    geom_errorbar(aes(ymin = lowerCi, ymax = upperCi,
                      color = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, "#E69F00", "#0a4c9e"),
                      alpha = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, 1, 0.5))) +
    geom_point(aes(y = lowerCi,
                   color = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, "#E69F00", "#0a4c9e"),
                   alpha = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, 0.8, 0.5))) +
    geom_point(aes(y = upperCi,
                   color = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, "#E69F00", "#0a4c9e"),
                   alpha = ifelse(lowerCi > results$trueEffect | upperCi < results$trueEffect, 0.8, 0.5))) +
    # VANITY
    labs(x = "Experiment Number", y = "Threshold", title = "Confidence Interval Coverage of True Effect") +
    scale_color_manual(values = c("#E69F00", "#0a4c9e")) +
    theme(legend.position = 'none') 
  
  ggsave(filename = 'General Linear Model/Figures/plot_ci_coverage99.png',
         width = 10, height = 3.5, dpi = 'retina', bg = 'transparent')


  
  