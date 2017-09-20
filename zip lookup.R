# Create zip lookup providing concordance between zip code and geographic predictor variables

library(sp)

d <- read_rds("~/Documents/Projects/CCL/Data/Spatial/ZCTA 2010 with state identifiers.rds")
g <- read_csv("~/Documents/Projects/CCL/Data/Spatial/Master State-Region Concordance File.csv")

zip_lookup <- d@data %>%
  rename(st = state_fips) %>%
  select(zip5, st) %>%
  filter(!duplicated(zip5)) %>%
  inner_join(select(g, st, state.name)) %>%
  rename(zip = zip5) %>%
  mutate(state = factor(state.name)) %>%  # Character predictors must be converted to factor to avoid warning in predict.gbm()
  select(zip, state)

# Save model object to disk as .rda object
save(zip_lookup, file = "data/zip_lookup.rda")
