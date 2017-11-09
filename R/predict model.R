#' Use prediction models to convert user-provided inputs into calculator return values
#'
#' @export
#' @importFrom quantreg predict.rq
#' @importFrom utils read.csv}

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
  stopifnot("rms" %in% nms)

  #-----

  nd$id <- 1:nrow(nd)

  # Assign geographic variables to 'nd' using zip code provided
  nd <- merge(nd, zip_lookup, sort = FALSE)
  nd <- nd[order(nd$id),]

  # Adult pre-tax dividend amount
  # Hard-coded; based on original CCL analysis for 2008-2012 period
  div.adult <- 377

  # #TO DO: Estimate household marginal tax rate
  nd$mrate <- margRate(nd)

  # Calculate pre- and post-tax dividends
  nd$div_pre <- round(div.adult * nd$na + 0.5 * div.adult * pmin(2, nd$nc))
  #nd$div_post <- round(nd$div_pre * (1 - nd$mrate))

  nd$np <- nd$na + nd$na
  nd$elec_ratio <- nd$cents_kwh / (nd$hinc / 1e3)
  nd$gas_ratio <- nd$gasprice / (nd$hinc / 1e3)

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
    #nd$hfuel = as.character(nd$hfuel)
    nd$hfuel[nd$hfuel == "Other or none"] <- "Natural gas"
  }

  q <- predict.rq(core_model, newdata = nd)
  core <- q[,2]
  stdev <- apply(q, 1, function(x) diff(x[c(1,3)])) / 1.35

  #----------------------

  # Gasoline weekly expenditure (median and 95th percentile)
  gas <- signif(predict.rq(gas_model, newdata = nd) * nd$gasprice / 52, digits = 2)
  colnames(gas) <- c("gas", "gas_upr")
  gas[which(nd$veh == "0"), "gas"] <- 0  # Set predicted gasoline expenditure to zero if Vehicles = 0 (user free to increase, if desired)

  # Electricity monthly expenditure (median and 95th percentile)
  elec <- signif(predict.rq(elec_model, newdata = nd) * nd$cents_kwh / 12, digits = 2)
  colnames(elec) <- c("elec", "elec_upr")

  # Primary heating fuel monthly expenditure (median and 95th percentile)
  predHeatModels <- function(d) {

    if (d$hfuel[1] == "Natural gas") {
      d$heat_ratio <- d$ngasprice / (d$hinc / 1e3)
      out <- predict.rq(ngas_model, newdata = d) * d$ngasprice
      out <- cbind(out, d$Natural_gas_cie)
    }

    if (d$hfuel[1] == "LPG/Propane") {
      d$heat_ratio <- d$lpgprice / (d$hinc / 1e3)
      out <- predict.rq(lpg_model, newdata = d) * d$lpgprice
      out <- cbind(out, d$LPG_cie)
    }

    if (d$hfuel[1] == "Heating oil") {
      d$heat_ratio <- d$hoilprice / (d$hinc / 1e3)
      out <- predict.rq(hoil_model, newdata = d) * d$hoilprice
      out <- cbind(out, d$Heating_oil_cie)
    }

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
  heat[other.id,] <- 0  # Set heat-related variables to zero if original hfuel='Other or none'

  #----------------------

  # Carbon price ($ per ton CO2)
  carbon.price <- 15

  # Estimated 90% margin of error (in dollars, annual)
  nd$moe <- round(1.645 * stdev * carbon.price)

  # Annual cost equation
  nd$cost <- paste(
    round(core * carbon.price, 2),
    paste0("gas * ", signif(52 * nd$Gasoline_cie * carbon.price / 1e3, 4)),
    paste0("elec * ", signif(12 * nd$Electricity_cie * carbon.price / 1e3, 4)),
    paste0("heat * ", signif(12 * heat$heat_cie * carbon.price / 1e3, 4)), sep = " + ")

  # Return results matrix
  psets <- cbind(gas, elec, heat)
  psets[psets < 0] <- 0
  result <- cbind(nd, psets)
  result <- subset(result, select = c(div_pre, mrate, cost, moe, gas, elec, heat, gas_upr, elec_upr, heat_upr))

  return(result)

}

#----------------

# INPUT FOR TESTING:
#for (i in list.files("~/Documents/Projects/exampleR/data/", full.names = TRUE)) load(i)
#require(quantreg)
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Electricity", veh = "2", rms = "7", stringsAsFactors = FALSE)
# nd <- data.frame(zip = c("94062","80524"), na = c(2, 1), nc = c(2, 0), hinc = c(50e3, 300e3), hfuel = c("Electricity", "Natural gas"), veh = c("2", "1"), rms = c("7", "5"), stringsAsFactors = FALSE)
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Other or none", veh = "2", rms = "7", stringsAsFactors = FALSE)
# predictModel(nd)
