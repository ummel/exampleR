library(gbm)
library(tidyverse)

d <- read_rds("~/Documents/Projects/CCL/Data/Results/Policy simulation results by household.rds")

# Temporary: Using total utilities expenditure at household level to proxy for monthly electricity expenditure
elec <- read_rds("~/Documents/Projects/CCL/Data/Results/Household category expenditure and emissions 2012.rds") %>%
  filter(cat == "Utilities", year == 2012, expend > 0) %>%
  mutate(elec = expend * 1400 / mean(expend),
         elec = elec / 12) %>%  # Call the variable electricity and set the mean to around 1400 per year
  select(id, elec)

# Adult pre-tax dividend amount
filter(d, np == 1, nc == 0)$div_pre[1]

# Fit GBM model predicting cost_co2, based on following user inputs:
# State of residence
# Household size
# Number of minors
# Age of householder
# Household income (slider variable)
# Rent expenditure (slider variable)

dset <- d %>%
  filter(year == 2012) %>%
  rename(state = state.name, hhsize = np, minors = nc, income = hinc) %>%
  mutate(state = factor(state)) %>%
  inner_join(elec) %>%
  filter(income > 0, income < quantile(income, 0.95), elec < quantile(elec, 0.95)) %>%  # Restrict income and electricity expenditure within reasonable bounds
  select(hw, state, hhsize, minors, age, income, elec, cost_co2) %>%
  sample_frac(0.05)

# Create vector of predictor variables
pred <- setdiff(names(dset), c("hw", "cost_co2"))

# Vector giving the var.monotone integer (-1, 0, 1) for each predictor variable
# Zero (no relationship) is the default
mono <- rep(0, length(pred))
mono[which(pred %in% c("income", "elec"))] <- 1

#----------

# Fit GBM model #0
# This model uses the full set of predictors
fitted_model0 <- gbm(formula = formula(paste0("cost_co2 ~ ", paste(pred, collapse = "+"))),
                     distribution = "gaussian",
                     data = dset,
                     weights = dset$hw,
                     var.monotone = mono,
                     n.trees = 1500,
                     interaction.depth = length(pred),
                     keep.data = FALSE)

# Save model object to disk as .rda object
save(fitted_model0, file = "data/fitted_model0.rda")

#----------

# Fit GBM model #1 - income
# This model uses ONLY the non-slider predictors  (slider predictors are assumed to be those where 'mono' == 1)
fitted_model1 <- gbm(formula = formula(paste0("income ~ ", paste(pred[mono != 1], collapse = "+"))),
                     distribution = "gaussian",
                     data = dset,
                     weights = dset$hw,
                     #var.monotone = mono,
                     n.trees = 1500,
                     interaction.depth = length(pred),
                     keep.data = FALSE)

# Save model object to disk as .rda object
save(fitted_model1, file = "data/fitted_model1.rda")

# Fit GBM model #2 - elec
# This model uses ONLY the non-slider predictors (slider predictors are assumed to be those where 'mono' == 1)
fitted_model2 <- gbm(formula = formula(paste0("elec ~ ", paste(pred[mono != 1], collapse = "+"))),
                     distribution = "gaussian",
                     data = dset,
                     weights = dset$hw,
                     #var.monotone = mono,
                     n.trees = 1500,
                     interaction.depth = length(pred),
                     keep.data = FALSE)

# Save model object to disk as .rda object
save(fitted_model2, file = "data/fitted_model2.rda")
