#' Use predict.rqion models to convert user-provided inputs into calculator return values
#'
#' @export
#' @importFrom quantreg predict.rq
#' @importFrom mgcv predict.gam
#' @importFrom utils read.csv

predictModel <- function(input) {

  # Input can either be csv file or data
  nd <- if (is.character(input) && file.exists(input)) {
    read.csv(input, stringsAsFactors = FALSE)
  } else {
    as.data.frame(input, stringsAsFactors = FALSE)
  }

  # Check that all necessary input variables are present
  nms <- names(nd)
  stopifnot("zip" %in% nms)
  stopifnot("na" %in% nms)
  stopifnot("nc" %in% nms)
  stopifnot("hinc" %in% nms)
  stopifnot("hfuel" %in% nms)
  stopifnot("veh" %in% nms)
  stopifnot("htype" %in% nms)

  #-----

  nd$id <- 1:nrow(nd)

  # Assign geographic variables to 'nd' using zip code provided
  nd <- merge(nd, zip_lookup, sort = FALSE)
  
  # Assign price adjustment factors based on assigned state
  nd <- merge(nd, price_adjustment, sort = FALSE)
  
  # Ensure that original order to maintained; probably not necessry is sort = FALSE in merge()
  nd <- nd[order(nd$id),]

  # Estimate household marginal tax rate
  # This is based on 2017 rates and brackets using household income as reported as user
  nd$mrate <- margRate(nd)
  
  # Adult pre-tax dividend amount
  # Hard-coded; based on original CCL analysis for 2008-2012 period
  # NOTE: This is NOT adjusted for inflation
  # Total revenue in analysis assumes $15 per ton at 2012 emission levels; $377 is based on those assumptions
  div.adult <- 377

  # Calculate pre-tax dividend amount, adjusted to current dollars
  nd$div_pre <- round(div.adult * nd$na + 0.5 * div.adult * pmin(2, nd$nc))

  nd$np <- nd$na + nd$na
  
  # Fuel-price-to-income ratios used in model fitting
  # Note that fuel prices are adjusted from 2012 to current price levels
  nd$elec_ratio <- nd$cents_kwh * nd$elec_adjust / (nd$hinc / 1e3)
  nd$gas_ratio <- nd$gasprice * nd$gas_adjust / (nd$hinc / 1e3)
  
  # Deflate user input household income to 2012 levels using CPI
  # This ensures that income matches currency units of original data and models
  nd$hinc <- nd$hinc / nd$cpi_adjust
  
  # Add 1 to 'nc' and 'veh' to allow log() calls in model prediction
  nd$nc <- nd$nc + 1
  nd$veh <- nd$veh + 1

  #----------------------

  # Core emissions
  # Calculate 90% CI margin of error (MOE) This is calculated by estimating the
  # standard deviation, assuming a Normal distribution with specified 25th and
  # 75th percentiles (the IQR is approximately 1.35x the standard deviation)

  # There were no 'Other or none' observations in the trimmed training dataset,
  # so this replacement assures the core_model predict() function will not fail
  # 'nd' object retains 'Other or none' so that hfuel logic further down still works

  other.id <- which(nd$hfuel == "Other or none")
  if (length(other.id) > 0) {
    nd$hfuel[other.id] <- "Natural gas"
  }
  
  # If user does not know heating fuel, make a guess
  # Fitted GLM model returns probability of Natural gas; Electricity otherwise (only two options)
  
  donotknow.id <- which(nd$hfuel == "Do not know")
  if (length(donotknow.id) > 0) {
    pred <- stats::predict.glm(hfuel_model, newdata = nd[donotknow.id,])
    # Convert from original scale to class probability (probability of Natural gas)
    # https://stats.stackexchange.com/questions/164648/output-of-logistic-regression-prediction
    prob <- exp(pred) / (1 + exp(pred))
    nd$hfuel[donotknow.id] <- ifelse(prob >= 0.5, "Natural gas", "Electricity")
  }

  # Note exp() used to convert log prediction value
  core <- as.numeric(exp(mgcv::predict.gam(core_model_gam, newdata = nd)))
  q <- exp(quantreg::predict.rq(core_model_rq, newdata = nd))
  if (nrow(nd) == 1) {
    stdev <- (q[2] - q[1]) / 1.35
  } else {
    stdev <- (q[,2] - q[,1]) / 1.35
  }
  
  #----------------------

  # Predict typical (mean) and upper-bound (97.5th percentile) expenditure values used to set the "page 2" slider preset and maximum value
  # Note that all dollar values are adjusted to reflect current price levels ("_adjust" variables)
  
  # Gasoline weekly expenditure (mean and 97.5th percentile)
  gas <- cbind(mgcv::predict.gam(gas_model_gam, newdata = nd), quantreg::predict.rq(gas_model_rq, newdata = nd))
  gas <- signif(gas * nd$gas_adjust * nd$gasprice / 52, digits = 2)
  #gas <- signif(predict(gas_model, newdata = nd) * nd$gas_adjust * nd$gasprice / 52, digits = 2)
  colnames(gas) <- c("gas", "gas_upr")
  gas[which(nd$veh == 1), "gas"] <- 0  # Set predicted gasoline expenditure to zero if Vehicles = 0 (which is actually "1" after +1 to 'veh' variable above

  # Electricity monthly expenditure (mean and 97.5th percentile)
  elec <- cbind(mgcv::predict.gam(elec_model_gam, newdata = nd), quantreg::predict.rq(elec_model_rq, newdata = nd))
  elec <- signif(gas * nd$elec_adjust * nd$cents_kwh / 12, digits = 2)
  #elec <- signif(predict(elec_model, newdata = nd) * nd$elec_adjust * nd$cents_kwh / 12, digits = 2)
  colnames(elec) <- c("elec", "elec_upr")

  # Primary heating fuel monthly expenditure (median and 95th percentile)
  # Note that expenditure values are adjusted to current price levels AND
  #  CIE values are adjusted to reflect change in fuel prices since 2012 (if price went up, CIE goes down)
  predHeatModels <- function(d) {

    if (d$hfuel[1] == "Natural gas") {
      d$heat_ratio <- d$ngasprice * d$ngas_adjust / (d$hinc / 1e3)
      out <- cbind(mgcv::predict.gam(ngas_model_gam, newdata = d), quantreg::predict.rq(ngas_model_rq, newdata = d)) * d$ngasprice * d$ngas_adjust
      out <- cbind(out, d$Natural_gas_cie / d$ngas_adjust)
    }

    if (d$hfuel[1] == "LPG/Propane") {
      d$heat_ratio <- d$lpgprice * d$lpg_adjust / (d$hinc / 1e3)
      out <- cbind(mgcv::predict.gam(lpg_model_gam, newdata = d), quantreg::predict.rq(lpg_model_rq, newdata = d)) * d$lpgprice * d$lpg_adjust
      out <- cbind(out, d$LPG_cie / d$lpg_adjust)
    }

    if (d$hfuel[1] == "Heating oil") {
      d$heat_ratio <- d$hoilprice * d$hoil_adjust / (d$hinc / 1e3)
      out <- cbind(mgcv::predict.gam(hoil_model_gam, newdata = d), quantreg::predict.rq(hoil_model_rq, newdata = d)) * d$hoilprice * d$hoil_adjust 
      out <- cbind(out, d$Heating_oil_cie / d$hoil_adjust)
    }
    
    # If heating fuel is Electricity, return zeros
    if (d$hfuel[1] %in% c("Electricity")) {
      out <- matrix(rep(0, 3 * nrow(d)), ncol = 3)
    }

    out <- cbind(d$id, out)
    colnames(out) <- c("id", "heat", "heat_upr", "heat_cie")
    out[,c("heat", "heat_upr")] <- signif(out[,c("heat", "heat_upr")] / 12, 2)
    return(out)

  }

  heat <- by(nd, nd$hfuel, predHeatModels)
  heat <- data.frame(do.call("rbind", heat))
  heat <- heat[order(heat$id), -1]

  #----------------------

  # Carbon price ($ per ton CO2; fixed)
  carbon.price <- 15

  # Estimated 90% margin of error (in dollars, annual)
  # Since emissions/tax burden from fuel use is assumed to be accurate (based on user specified values),
  #  the overall MOE is the modeled uncertainty in the "core" emissions
  nd$moe <- round(1.645 * stdev * carbon.price)

  # Annual cost equation
  nd$cost <- paste(
    round(core * carbon.price, 2),
    paste0("gas * ", signif(52 * nd$Gasoline_cie / nd$gas_adjust * carbon.price / 1e3, 4)),  # Input is weekly expenditure
    paste0("elec * ", signif(12 * nd$Electricity_cie / nd$elec_adjust * carbon.price / 1e3, 4)),  # Input is monthly expenditure
    paste0("heat * ", signif(12 * heat$heat_cie * carbon.price / 1e3, 4)), sep = " + ")  # Input is monthly expenditure

  # For users that did not know their heating fuel, replace the 'heat' component
  #  of cost equation with the default expenditure value and set 'heat' and 'heat_upr' variables to zero
  if (length(donotknow.id) > 0) {
    for (i in donotknow.id) nd$cost[i] <- sub("heat", heat$heat[i], nd$cost[i], fixed = TRUE)
    heat$heat[donotknow.id] <- 0
    heat$heat_upr[donotknow.id] <- 0
  }
  
  # For users that report heating fuel as "Other or none", replace the 'heat' component
  #  of cost equation with ZERO value and set heat-related presets to zero
  if (length(other.id) > 0) {
    for (i in other.id) nd$cost[i] <- sub("heat", 0, nd$cost[i], fixed = TRUE)
  }
  heat[other.id,] <- 0
  
  # Return results matrix
  psets <- cbind(gas, elec, heat)
  psets[psets < 0] <- 0
  result <- cbind(nd, psets)
  result <- subset(result, select = c(div_pre, mrate, cost, moe, gas, elec, heat, gas_upr, elec_upr, heat_upr))

  return(result)

}

#----------------

# INPUT FOR TESTING:
# require(quantreg, quietly = TRUE)
# require(mgcv, quietly = TRUE)
# source("R/marginal rate.R")
# for (i in list.files("data/", full.names = TRUE)) load(i)
# 
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Electricity", veh = 2, htype = "Other", stringsAsFactors = FALSE)
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Other or none", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = c("94062","80524","70032"), na = c(2, 1, 3), nc = c(2, 0, 3), hinc = c(50e3, 300e3, 100e3), hfuel = c("Do not know", "Natural gas", "Other or none"), veh = c(2, 0, 3), htype = c("Stand-alone house", "Apartment building", "Other"), stringsAsFactors = FALSE)
#
# predictModel(nd)
