# Copy of .onLoad script from https://github.com/rwebapps/tvscore/blob/master/R/onLoad.R
.onLoad <- function(lib, pkg){
  #automatically loads the dataset when package is loaded
  #do not use this in combination with lazydata=true

  utils::data(core_cost_gam, package = pkg, envir = parent.env(environment()))
  utils::data(core_cost_rq, package = pkg, envir = parent.env(environment()))

  utils::data(div_taxrate_gam, package = pkg, envir = parent.env(environment()))
  
  utils::data(elec_cie_gam, package = pkg, envir = parent.env(environment()))
  utils::data(elec_expend_gam, package = pkg, envir = parent.env(environment()))
  utils::data(elec_expend_rq, package = pkg, envir = parent.env(environment()))

  utils::data(gas_cie_gam, package = pkg, envir = parent.env(environment()))
  utils::data(gas_expend_gam, package = pkg, envir = parent.env(environment()))
  utils::data(gas_expend_rq, package = pkg, envir = parent.env(environment()))

  utils::data(hfuel_glm, package = pkg, envir = parent.env(environment()))

  utils::data(hoil_cie_gam, package = pkg, envir = parent.env(environment()))
  utils::data(hoil_expend_gam, package = pkg, envir = parent.env(environment()))
  utils::data(hoil_expend_rq, package = pkg, envir = parent.env(environment()))
  
  utils::data(input_summary, package = pkg, envir = parent.env(environment()))

  utils::data(lpg_cie_gam, package = pkg, envir = parent.env(environment()))
  utils::data(lpg_expend_gam, package = pkg, envir = parent.env(environment()))
  utils::data(lpg_expend_rq, package = pkg, envir = parent.env(environment()))

  utils::data(ngas_cie_gam, package = pkg, envir = parent.env(environment()))
  utils::data(ngas_expend_gam, package = pkg, envir = parent.env(environment()))
  utils::data(ngas_expend_rq, package = pkg, envir = parent.env(environment()))

  utils::data(scenario_parameters, package = pkg, envir = parent.env(environment()))
  utils::data(zip_lookup, package = pkg, envir = parent.env(environment()))

}
