# Copy of .onLoad script from https://github.com/rwebapps/tvscore/blob/master/R/onLoad.R
.onLoad <- function(lib, pkg){
  #automatically loads the dataset when package is loaded
  #do not use this in combination with lazydata=true

  utils::data(core_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(core_model_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(gas_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(gas_model_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(elec_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(elec_model_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(ngas_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(ngas_model_rq, package = pkg, envir = parent.env(environment()))

  utils::data(hoil_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(hoil_model_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(lpg_model_gam, package = pkg, envir = parent.env(environment()))
  utils::data(lpg_model_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(hfuel_model, package = pkg, envir = parent.env(environment()))
  
  utils::data(price_adjustment, package = pkg, envir = parent.env(environment()))
  utils::data(zip_lookup, package = pkg, envir = parent.env(environment()))
  utils::data(input_summary, package = pkg, envir = parent.env(environment()))
  
  # Added in version 2 (Aug 2019)
  utils::data(finval_model_gam, package = pkg, envir = parent.env(environment()))
  
}
