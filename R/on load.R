# Copy of .onLoad script from https://github.com/rwebapps/tvscore/blob/master/R/onLoad.R
.onLoad <- function(lib, pkg){
  #automatically loads the dataset when package is loaded
  #do not use this in combination with lazydata=true
  utils::data(fitted_model, package = pkg, envir = parent.env(environment()))
}
