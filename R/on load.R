# Copy of .onLoad script from https://github.com/rwebapps/tvscore/blob/master/R/onLoad.R
.onLoad <- function(lib, pkg){
  #automatically loads the dataset when package is loaded
  #do not use this in combination with lazydata=true

  #utils::data(list = list.files("data"), package = pkg, envir = parent.env(environment()))

  utils::data(fitted_model0, package = pkg, envir = parent.env(environment()))
  utils::data(fitted_model1, package = pkg, envir = parent.env(environment()))
  utils::data(fitted_model2, package = pkg, envir = parent.env(environment()))
  utils::data(zip_lookup, package = pkg, envir = parent.env(environment()))

}
