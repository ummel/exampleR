# Copy of .onLoad script from https://github.com/rwebapps/tvscore/blob/master/R/onLoad.R
.onLoad <- function(lib, pkg){
  #automatically loads the dataset when package is loaded
  #do not use this in combination with lazydata=true

  utils::data(core_model, package = pkg, envir = parent.env(environment()))
  utils::data(elec_model, package = pkg, envir = parent.env(environment()))
  utils::data(gas_model, package = pkg, envir = parent.env(environment()))
  utils::data(hoil_model, package = pkg, envir = parent.env(environment()))
  utils::data(lpg_model, package = pkg, envir = parent.env(environment()))
  utils::data(ngas_model, package = pkg, envir = parent.env(environment()))
  utils::data(zip_lookup, package = pkg, envir = parent.env(environment()))
  utils::data(input_summary, package = pkg, envir = parent.env(environment()))

}
