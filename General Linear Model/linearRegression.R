library(tidyverse)


# SIMPLE REGRESSION ----
## Data simulation ----
## Context 👇
## We have be provided with pain intensity on the forearm at a standard heat level,
## as well as the age and the VO2peak of 15 participants.

set.seed(42)
n <- 15 # Number of participants
# A tibble is a type of able used in R. They are similar to data frames.
simpleRegression.data <- tibble(
  # Set subject's unique ID from 1 to 15.
  id = 1:n,
  # Generate a set of 15 random numbers following a normal distribution with a
  # mean of 35 and a standard deviation of 10 (arbitrary values).
  vo2peak = rnorm(n, mean = 35, sd = 10),
  # Generate a set of 15 random numbers following a normal distribution that will
  # be negatively related to vo2peak
  pain = 7 - scale(vo2peak)[,1] + rnorm(n, mean = 0, sd = 1)
)

## Linear regression with optim() ----

# This is the function were are trying to minimize, i.e., the error sum of squares
# This actually just measures the sum of the squared residuals
sse <- function(params, x, y){
  b0 <- params[1]                      # Intercept
  b1 <- params[2]                      # Slope
  predictions <- b0 + b1*x             # Traditional regression formula
  residuals <- y - predictions         # Vector containing all residuals
  sum(residuals^2)                     # Error sum of squares
}

# Optimize the function SSE, i.e., finds params allowing sse to return the lowest value
optim(
  par = c(b0 = 10, b1 = -1),
  fn = sse,
  x = simpleRegression.data$vo2peak,
  y = simpleRegression.data$pain,
  control = list(
    reltol = 1e-16
  )
)

## Linear regression using closed-form equation ----
# We need our data as vector for further use
vo2max.vector <- simpleRegression.data$vo2peak
pain.vector <- simpleRegression.data$vo2peak

# Linear algebra
X <- cbind(1, vo2peak)   # Design matrix (X)
y <- pain.vector               # Response vector (y)

varianceCovariance <- t(X) %*% X
crossproduct <- t(X) %*% y
varianceCovarianceInverse <- solve(varianceCovariance)
betas <- varianceCovarianceInverse %*% crossproduct

# Additional challenges....
# Can you figure out why these lines give the elements of the variance-covariance matrix?
sum(X[,1] * t(X)[1,])   # First element
sum(X[,2] * t(X)[2,])   # Diagonal
sum(X[,1] * t(X)[2,])   # Off-diagonal (top-right)
sum(X[,2] * t(X)[1])    # Off-diagonal (bottom-left)
# Can you figure out why these lines give us both rows of the cross-product matrix?
sum(t(X)[1,] * y) # First row
sum(t(X)[2,] * y) # Second row
# Can you figure out why these lines give us the betas?
sum(varianceCovarianceInverse[1,] * crossproduct[,1])
sum(varianceCovarianceInverse[2,] * crossproduct[,1])


## Linear model with lm() ----
regression.lm <- lm(data = simpleRegression.data, pain ~ vo2peak)
summary(regression.lm)

## Scatterplot ----
# Average value of pain levels across all 15 subjects
regression.mean <- mean(simpleRegression.data$pain)

simpleRegressionResiduals.plot <- simpleRegression.data %>%
  mutate(predictedPain = predict(regression.lm)) %>%
  ggplot(aes(x = vo2peak, y = pain)) +
    geom_segment(aes(x = vo2max, xend = vo2max, yend = predictedPain),
                 linetype = 'dotted', color = 'red', linewidth  = 0.5) +
    geom_point(shape = 'triangle', color = '#0a4c6e', alpha = 0.9, size = 3) +
    geom_smooth(method = 'lm', se = FALSE, color = 'black') +
    labs(title = 'Pain level as a function of VO2peak',
         x = 'VO2max', y = 'Pain level') +
    theme(
      axis.title = element_text(size = 15, face = 'bold'),
      legend.position = 'none'
    )

ggsave(simpleRegressionResiduals.plot, file = 'General Linear Model/Figures/plot_regression-residuals.jpg',
       width = 4, height = 3, dpi = 'retina')


simpleRegressionMean.plot <- simpleRegression.data %>%
  ggplot(aes(x = vo2max, y = pain)) +
    geom_segment(aes(x = vo2max, xend = vo2max, yend = regression.mean),
                 linetype = 'dotted', color = 'red', linewidth = 0.5) +
    geom_point(shape = 'triangle', color = '#0a4c6e', alpha = 0.9, size = 3) +
    geom_smooth(method = 'lm', se = FALSE, color = 'black') +
    geom_hline(yintercept = regression.mean, linetype = 'dashed') +
    labs(title = 'Pain level as a function of VO2max',
         x = 'VO2max', y = 'Pain level') +
    theme(
      axis.title = element_text(size = 15, face = 'bold'),
      legend.position = 'none'
    )

ggsave(simpleRegressionMean.plot, file = 'LMM vs ANOVA/Figures/plot_regression-mean.tiff',
       width = 4, height = 3, dpi = 'retina')

# ANALYSIS OF VARIANCE ----
## Data simulation ----
## Context 👇
## Thee groups: non-athletes, non-contact athletes and contact athletes
## Pain levels measured three times on the forearm of each participants

set.seed(42)
n <- 5 # n per group
groups <- c('nonAthlete', 'nonContactAthlete', 'contactAthlete')

# Simulate thermal pain in all three groups (N = 60, n = 30)
# Create a tibble containing 3 pain scores per participant
data <- purrr::map_df(groups, ~tibble(
  id = rep(1:n, each = 3),
  group = factor(.x),
  pain = case_when(
    .x == 'nonAthlete' ~ rnorm(3 * n, mean = 5.5, sd = 1),          # Higher pain
    .x == 'nonContactAthlete' ~ rnorm(3 * n, mean = 5, sd = 1.25),   # Moderate pain
    .x == 'contactAthlete' ~ rnorm(3 * n, mean = 4.5, sd = 1.15)       # Lower pain
    )
  )) %>%
  mutate(
    id = rep(1:(n * length(groups)), each = 3),
    group = factor(group,
                      labels = c('Non-athletes', 'Non-contact athletes', 'Contact athletes'))
  )

plot.groupColors = c('Non-athletes' = 'grey',
                     'Non-contact athletes' = '#0a4c6e',
                     'Contact athletes' = '#E69F00')

plot.groupShapes = c('Non-athletes' = 'circle',
                     'Non-contact athletes' = 'square',
                     'Contact athletes' = 'triangle')



data.anova <- data %>%
group_by(id, group) %>%
summarize(pain = mean(pain)) %>%
ungroup()

anova.groupMean <- data.anova %>%
  group_by(group) %>%
  summarize(groupMean = mean(pain))

data.anova <- data.anova %>%
  left_join(anova.groupMean, by = 'group')


## ANOVA ----
anova.grandMean <- mean(data.anova$pain)                             # Grand mean
anova.sst <- sum((data.anova$pain - anova.grandMean)^2)              # Total sum of squares
anova.ssm <- sum((data.anova$groupMean - anova.grandMean)^2)         # Model sum of squares
anova.dfm <- length(anova.groupMean$group) - 1                       # Model's degree of fredoom
anova.ssr <- sum((data.anova$pain - data.anova$groupMean)^2)         # Residual sum of squares
anova.dfr <- length(data.anova$id) - length(anova.groupMean$group)   # Residual's degree of freedom
anova.msm <- anova.ssm / anova.dfm                                   # Model's mean square
anova.msr <- anova.ssr / anova.dfr                                   # Residual's mean square
anova.f   <- anova.msm / anova.msr


## Scatterplots ----
  plot.sst <- ggplot(data.anova, aes(x = id, y = pain,
                                     shape = group, color = group)) +
    # Horizontal and vertical lines
    geom_hline(yintercept = anova.grandMean, linetype = 'dashed', color = 'black') +
    geom_segment(aes(x = id, xend = id, y = pain, yend = grandMean),
                 linetype = 'dotted', color = 'red', linewidth  = 0.5) +
    # Scatterplot
    geom_point(size = 2, alpha = 0.90) +
    # Aesthethic
    scale_color_manual(values = plot.groupColors) +
    scale_shape_manual(values = plot.groupShapes) +
    facet_wrap(~group, scales = 'free_x') +
    labs(title = 'Pain level for all groups',
         x = 'Participants', y = 'Pain level') +
    theme(
      axis.title.y = element_text(size = 15, face = 'bold'),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none'
    )
  
  ggsave(plot = plot.sst, file = 'LMM vs ANOVA/Figures/plot_ANOVA-sst.tiff',
         width = 4, height = 3, dpi = 'retina')
    
  
  plot.ssm <- ggplot(data.anova, aes(x = id, y = pain,
                                     shape = group, color = group)) +
    # Horizontal and vertical lines
    geom_hline(aes(yintercept = groupMean), linetype = 'dashed') +
    geom_segment(aes(x = id, xend = id, y = groupMean, yend = grandMean),
                 linetype = "dotted", color = 'red', linewidth  = 0.5) +
    geom_hline(aes(yintercept = grandMean), linetype = 'dashed') +
    geom_point(size = 2, alpha = 0.90) +
    scale_color_manual(values = plot.groupColors) +
    scale_shape_manual(values = plot.groupShapes) +
    facet_wrap(~group, scales = 'free_x') +
    labs(title = 'Pain level for all groups',
         x = 'Participants', y = 'Pain level') +
    theme(
      axis.title.y = element_text(size = 15, face = 'bold'),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none'
    )
  
  ggsave(plot = plot.ssm, file = 'LMM vs ANOVA/Figures/plot_ANOVA-ssm.tiff',
         width = 4, height = 3, dpi = 'retina')
  
  
  plot.ssr <- ggplot(data.anova, aes(x = id, y = pain,
                                     shape = group, color = group)) +
    # Horizontal and vertical lines
    geom_hline(aes(yintercept = groupMean), linetype = 'dashed') +
    geom_segment(aes(x = id, xend = id, yend = groupMean),
                 linetype = "dotted", color = 'red', linewidth  = 0.5) +
    geom_point(size = 2, alpha = 0.90) +
    scale_color_manual(values = plot.groupColors) +
    scale_shape_manual(values = plot.groupShapes) +
    facet_wrap(~group, scales = 'free_x') +
    labs(title = 'Pain level for all groups',
         x = 'Participants', y = 'Pain level') +
    theme(
      axis.title.y = element_text(size = 15, face = 'bold'),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none'
    )
  
  ggsave(plot = plot.ssr, file = 'LMM vs ANOVA/Figures/plot_ANOVA-ssr.tiff',
         width = 4, height = 3, dpi = 'retina')
