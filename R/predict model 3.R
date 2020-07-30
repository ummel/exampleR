#' Generate results used by CCL online calculator
#' @description CCL’s
#'   \href{https://citizensclimatelobby.org/household-impact-study/}{Household
#'   Impact Study} estimates the direct financial effect of a carbon tax and
#'   dividend policy for a large, representative sample of U.S. households.
#'   Techniques, data, and assumptions are described in detail in the
#'   \href{https://11bup83sxdss1xze1i3lpol4-wpengine.netdna-ssl.com/wp-content/uploads/2016/05/Ummel-Impact-of-CCL-CFD-Policy-v1_4.pdf}{associated
#'    working paper}.\cr\cr The online calculator tool uses the study’s results
#'   to estimate a household’s additional costs under the policy (due to higher
#'   prices for goods and services), depending on a limited set of household
#'   characteristics (income, number of vehicles, etc.). It also calculates the
#'   expected dividend, which is a function of the household’s number of adults,
#'   number of minors, and expected federal marginal tax rate. The difference
#'   between the dividend and additional cost is the "net" or overall financial
#'   impact – positive if a household is likely to “come out ahead” under
#'   CF&D.\cr\cr For ease of use, a small and generally easy-to-recall set of
#'   user inputs are solicited. The calculator reports the expected average
#'   outcome for a household. The actual outcome for any specific household
#'   could vary from the average. For example, if a user household is a
#'   below-average consumer of carbon-intensive goods like air travel and meat,
#'   the calculator will understate the net impact (and vice-versa). Developing
#'   a precise estimate for every household would require many more questions
#'   and accurate recall. We have opted for simplicity over precision.\cr\cr The
#'   \code{predictModel2()} function described here translates user-provided
#'   household characteristics into results displayed on the online calculator
#'   application. The function uses a limited set of user characteristics (see
#'   'Arguments' below). The input data format is based on the OpenCPU
#'   \href{https://www.opencpu.org/posts/scoring-engine/}{'tvscore' example}.
#'   See 'Details' section below for technical details.\cr\cr The complete
#'   'exampleR' package and source code is available at:
#'   \href{https://github.com/ummel/exampleR}{https://github.com/ummel/exampleR}
#'
#' @param input A text string or .csv file passed via
#'   \href{https://www.opencpu.org/api.html}{OpenCPU API} (i.e. cURL POST) or a
#'   local R data frame (i.e. for debugging). In either case, \code{input}
#'   should contain the user-provided variables below. Calling
#'   \code{inputSummary()} or viewing the \code{input_summary} data object
#'   provided with package will provide the data types and allowable values for
#'   each of the variables. \itemize{ \item{zip: 5-digit zip code} \item{na:
#'   number of adults in household} \item{nc: number of minors in household}
#'   \item{hinc: household income} \item{hfuel: household primary heating fuel}
#'   \item{veh: number of vehicles owned by household} \item{htype: dwelling
#'   type} \item{dirfrac: fraction (0 - 1) of direct emissions tax burden that
#'   is passed through to consumer prices. If this argument is not supplied, it
#'   defaults to 0.95.} \item{indfrac: fraction (0 - 1) of indirect emissions tax
#'   burden that is passed through to consumer prices. If this argument is not
#'   supplied, it defaults to 0.495.} }
#'
#' @details Household-level results from the Household Impact Study were
#'   selected for the year 2012, resulting in a total sample of just over 1
#'   million households. Results for each household were processed to determine
#'   the expected additional financial cost (under the CF&D policy) associated
#'   with "indirect" emissions and those stemming from consumption of gasoline,
#'   electricity, and the household's primary heating fuel. A series of
#'   statistical models were fit to the household sample to determine the
#'   relationship between a limited set of household characteristics and the
#'   cost components. A wide variety of household characteristics were
#'   considered for inclusion in the models; the subset ultimatey selected are
#'   both easy for users to accurately recall and demonstrate a good ability to
#'   (collectively) predict a household's expected additional cost.\cr\cr The
#'   fitted models are capable of translating household characteristics into
#'   expected (average) additional cost (generalized additive models with
#'   smoothing terms) as well as conditional quantiles (quantile regression) for
#'   the purposes of uncertainty estimation. When estimating emissions/cost
#'   associated with a households indirect emissions component,
#'   \code{predictModels()} uses a GAM model to predict the average cost and
#'   quantile models to predict the conditional 25th and 75th percentiles. The
#'   latter are used to estimate the uncertainty around the expected value,
#'   assuming a Normal distribution.\cr\cr In the case of emissons associated
#'   with gasoline and utilities, \code{predictModel()} returns the expected
#'   (average) monthly expenditure value and a "cost formula" that can translate
#'   monthly expenditures into total annual additional cost (including cost
#'   associated with indirect emissions). This feature allows users of the
#'   calculator to adjust the "default" average expenditure values to reflect
#'   their specific situation, resulting in a more accurate overall estimate of
#'   the additonal cost for that household.\cr\cr In order to account for the
#'   fact that the data used to fit the statistical models is from 2012,
#'   \code{predictModel()} includes state-level, fuel-specific price adjustment
#'   factors to inflate or deflate (as appropriate) user-provided expenditure
#'   and income values to current price levels. This ensures that inflation and
#'   changes in fuel prices over time do not unduly affect the results. No
#'   analogous adjustment is made for changes to electricity grid
#'   carbon-intensity over time.\cr\cr A fixed carbon price of $15 per ton CO2
#'   is assumed. Per the assumptions and caveats in the Household Impact Study
#'   working paper, the household results reflect the "overnight" (i.e.
#'   short-term) direct financial impact of the CF&D proposal, ignoring dynamic
#'   economic effects and changes in employment, preferences, or
#'   technologies.\cr\cr The expected post-tax dividend for a given household is
#'   determined by the "div_pre" and "mrate" values returned by
#'   \code{predictModel()}. The former is simply a function of the number of
#'   adults and children in the household and the "full-share" dividend value.
#'   The latter is estimated via \code{\link{margRate}}.\cr\cr Validity and
#'   performance of \code{predictModel()} was tested by passing a random sample
#'   of the original Household Impact Study household-level results to the
#'   function and comparing the model-generated results to those in the original
#'   data. This quality-control test indicates that the predictive \eqn{R^2}
#'   value
#'   (\href{https://en.wikipedia.org/wiki/Coefficient_of_determination}{coefficient
#'    of determination}) is 0.55 in the event that a user relies on the
#'   model-generated gasoline and utility expenditure values. Model skill
#'   improves to a \eqn{R^2} value of 0.73 when users provide their own
#'   (accurate) expenditure values. Those interested in greater detail are asked
#'   to review the annotated public source code for the \code{predictModel2()}
#'   function \href{https://github.com/ummel/exampleR/tree/master/R}{on GitHub}.
#'
#' @return Function returns a data frame with one row of outputs for each row in
#'   \code{input}. The output variables are: \itemize{ \item{div_pre: household
#'   pre-tax annual dividend} \item{mrate: estimated marginal federal tax rate
#'   (see \code{\link{margRate}})} \item{cost: character string giving formula
#'   that (when evaluated) returns annual policy cost given annual total
#'   expenditure inputs for gasoline (gas), electricity (elec), and heating fuel
#'   (heat)} \item{gas: predicted annual gasoline expenditure for
#'   household (used as slider preset in online calculator)} \item{elec:
#'   predicted annual electricity expenditure for household (used as
#'   slider preset in online calculator)} \item{heat: predicted annual
#'   primary heating fuel expenditure for household (used as slider preset in
#'   online calculator; zero if not applicable)} \item{gas_upr: predicted
#'   maximum feasible annual gasoline expenditure for household (used as slider
#'   max value in online calculator)} \item{elec_upr: predicted maximum feasible
#'   annual electricity expenditure for household (used as slider max value in
#'   online calculator)} \item{heat_upr: predicted maximum feasible annual
#'   primary heating fuel expenditure for household (used as slider max value in
#'   online calculator; zero if not applicable)} }
#'
#' @section Example applications: The CCL calculator tool invokes
#'   \code{predictModel2()} remotely via cURL POST calls to OpenCPU. For
#'   example:\cr\cr \code{curl
#'   https://ummel.ocpu.io/exampleR/R/predictModel2/json -H "Content-Type:
#'   application/json" -d '{"input" : [ {"zip":"80524", "na":2, "nc":2,
#'   "hinc":50000, "hfuel":"Natural gas", "veh":2, "htype":"Stand-alone house"}
#'   ]}'}\cr\cr Or \code{predictModel2()} can be called locally within R:\cr\cr
#'   \code{nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel =
#'   "Electricity", veh = 2, htype = "Other")}\cr \code{predictModel2(nd)}\cr\cr
#'   A front-end developer can view input parameters details by calling:\cr\cr
#'   \code{curl https://ummel.ocpu.io/exampleR/R/inputSummary/json -H
#'   "Content-Type: application/json" -d '{}'}
#'
#' @export
#' @importFrom quantreg predict.rq
#' @importFrom mgcv predict.gam
#' @importFrom utils read.csv

predictModel3 <- function(input) {
  
  # Helper function
  enforceLimits <- function(x, limits) {
    x <- pmax(x, limits[1])
    x <- pmin(x, limits[2])
    return(x)
  }
  
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
  
  # Detect if the input is a test dataset
  # If the original training data is submitted, all 'zip' values will be NA
  test <- all(is.na(nd$zip))
  
  # NOT USED. Ensure that a single scenario is specified ('S')
  # stopifnot("scenario" %in% nms)
  # S <- nd$scenario[[1]]
  # stopifnot(all(nd$scenario == S))
  
  #-----
  
  nd$id <- 1:nrow(nd)
  
  # Assign geographic variables to 'nd' using zip code provided
  # If statement used to allow testing with the original dataset (which contains spatial variables but not zip code)
  if (!test) {
    nd <- merge(nd, zip_lookup, by = "zip", sort = FALSE, all.x = TRUE)
    if (any(is.na(nd$census_urban))) {
      miss <- which(is.na(nd$census_urban))
      j <- intersect(names(nd), names(zip_lookup))
      nd[miss, j] <- zip_lookup[match(substring(nd$zip[miss], 1, 4), substring(zip_lookup$zip, 1, 4)), j]
      if (any(is.na(nd$census_urban))) stop("Zip code not found.")  # If still missing, stop with error
    } 
  }
  
  # Ensure that original order to maintained; probably not necessary if sort = FALSE in merge()
  nd <- nd[order(nd$id),]
  
  #-----
  
  # Deflate user input household income to 2018 levels using CPI
  # This ensures that income matches currency units of original data and models
  nd$hinc <- nd$hinc * 0.98  # Hard-coded to deflate from 2020 to 2018 levels
  
  #----------------------
  
  # There are no 'Other or none' observations in the training dataset,
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
    pred <- stats::predict.glm(hfuel_glm, newdata = nd[donotknow.id,])
    # Convert from original scale to class probability (probability of Natural gas)
    # https://stats.stackexchange.com/questions/164648/output-of-logistic-regression-prediction
    prob <- exp(pred) / (1 + exp(pred))
    nd$hfuel[donotknow.id] <- ifelse(prob > 0.5, "Natural gas", "Electricity")
  }
  
  #----------------------
  
  # Pre-tax dividend amount; computed explicitly
  nd$div_pre <- round((nd$na + 0.5 * nd$nc) * scenario_parameters$dividend)
  
  # Estimate tax rate on the dividend
  nd$mrate <- as.numeric(mgcv::predict.gam(div_taxrate_gam, newdata = nd))
  nd$mrate <- enforceLimits(nd$mrate, scenario_parameters$div_taxrate)
  nd$mrate <- signif(nd$mrate, 3)
  
  #----------------------
  
  # SCENARIO 3 INDIRECT ("core") tax burden
  core <- as.numeric(mgcv::predict.gam(core_cost_gam, newdata = nd))
  core <- enforceLimits(core, scenario_parameters$bcore_3)
  
  # Associated standard deviation, assuming Normal distribution
  # Note exp() used to convert log prediction values for the 'rq' models
  q <- exp(quantreg::predict.rq(core_cost_rq, newdata = nd))
  if (nrow(nd) == 1) {
    stdev <- (q[2] - q[1]) / 1.35
  } else {
    stdev <- (q[,2] - q[,1]) / 1.35
  }
  
  #-----
  
  # SCENARIO 2 INDIRECT ("core") tax burden
  # core2 <- as.numeric(mgcv::predict.gam(core2_cost_gam, newdata = nd))
  # core2 <- enforceLimits(core2, scenario_parameters$bcore_2)
  # 
  # # Associated standard deviation, assuming Normal distribution
  # q <- exp(quantreg::predict.rq(core2_cost_rq, newdata = nd))
  # if (nrow(nd) == 1) {
  #   stdev2 <- (q[2] - q[1]) / 1.35
  # } else {
  #   stdev2 <- (q[,2] - q[,1]) / 1.35
  # }
  
  # Estimate 90% margin of error (in dollars, annual)
  # Since emissions/tax burden from fuel use is assumed to be accurate (based on user specified values),
  #  the overall MOE is the modeled uncertainty in the "core" emissions
  #nd$moe1 <- round(1.645 * stdev1)
  #nd$moe2 <- round(1.645 * stdev2)
  nd$moe <- round(1.645 * stdev)
  
  #----------------------
  
  # Predict typical (mean) and upper-bound (approximately 97.5th percentile) expenditure values used to set the "page 2" slider preset and maximum value
  
  # Gasoline annual expenditure (mean and approximately 97.5th percentile)
  
  gas <- cbind(mgcv::predict.gam(gas_expend_gam, newdata = nd),
               quantreg::predict.rq(gas_expend_rq, newdata = nd) * 1.2,
               mgcv::predict.gam(gas_cie_gam, newdata = nd))
  colnames(gas) <- c("gas", "gas_upr", "gas_cie")
  gas[, "gas"] <- enforceLimits(gas[, "gas"], scenario_parameters$gas)
  gas[, "gas_upr"] <- enforceLimits(gas[, "gas_upr"], scenario_parameters$gas)
  gas[, "gas_cie"] <- enforceLimits(gas[, "gas_cie"], scenario_parameters$gas_cie)
  gas[which(nd$veh == 0), "gas"] <- 0  # Set predicted gasoline expenditure to zero if user inputs zero vehicles
  gas <- gas * nd$seds_gas_price  # Model output is original value divided by state average price; this converts back to original units
  gas[, 2] <- pmax(gas[, 2], 2 * gas[, 1])  # The upper expenditure value is set to minimum 2x the mean prediction
  gas <- signif(gas, 3)
  
  #-----
  
  # Electricity annual expenditure (mean and approximately 97.5th percentile)
  
  elec <- cbind(mgcv::predict.gam(elec_expend_gam, newdata = nd), 
                quantreg::predict.rq(elec_expend_rq, newdata = nd) * 1.2,
                mgcv::predict.gam(elec_cie_gam, newdata = nd))
  colnames(elec) <- c("elec", "elec_upr", "elec_cie")
  elec[, "elec"] <- enforceLimits(elec[, "elec"], scenario_parameters$elec)
  elec[, "elec_upr"] <- enforceLimits(elec[, "elec_upr"], scenario_parameters$elec)
  elec[, "elec_cie"] <- enforceLimits(elec[, "elec_cie"], scenario_parameters$elec_cie)
  elec <- elec * nd$seds_elec_price  # Model output is annual expenditure divided by state average price; this converts back to original units
  elec[, 2] <- pmax(elec[, 2], 2 * elec[, 1])  # The upper expenditure value is set to minimum 2x the mean prediction
  elec <- signif(elec, 3)
  
  #-----
  
  # Primary heating fuel annual expenditure (median and approximately 95th percentile)

  predHeatModels <- function(d) {
    
    if (d$hfuel[1] == "Natural gas") {
      out <- cbind(mgcv::predict.gam(ngas_expend_gam, newdata = d), 
                   quantreg::predict.rq(ngas_expend_rq, newdata = d) * 1.2,
                   mgcv::predict.gam(ngas_cie_gam, newdata = d))
      out[, 1] <- enforceLimits(out[, 1], scenario_parameters$ngas)
      out[, 2] <- enforceLimits(out[, 2], scenario_parameters$ngas)
      out[, 3] <- enforceLimits(out[, 3], scenario_parameters$ngas_cie)
      out <- out * d$seds_ngas_price  # Model output is annual expenditure divided by state average price; this converts back to original units
    }
    
    if (d$hfuel[1] == "LPG/Propane") {
      out <- cbind(mgcv::predict.gam(lpg_expend_gam, newdata = d), 
                   quantreg::predict.rq(lpg_expend_rq, newdata = d) * 1.2,
                   mgcv::predict.gam(lpg_cie_gam, newdata = d))
      out[, 1] <- enforceLimits(out[, 1], scenario_parameters$lpg)
      out[, 2] <- enforceLimits(out[, 2], scenario_parameters$lpg)
      out[, 3] <- enforceLimits(out[, 3], scenario_parameters$lpg_cie)
      out <- out * d$seds_lpg_price  # Model output is annual expenditure divided by state average price; this converts back to original units
    }
    
    if (d$hfuel[1] == "Heating oil") {
      out <- cbind(mgcv::predict.gam(hoil_expend_gam, newdata = d), 
                   quantreg::predict.rq(hoil_expend_rq, newdata = d) * 1.2,
                   mgcv::predict.gam(hoil_cie_gam, newdata = d))
      out[, 1] <- enforceLimits(out[, 1], scenario_parameters$hoil)
      out[, 2] <- enforceLimits(out[, 2], scenario_parameters$hoil)
      out[, 3] <- enforceLimits(out[, 3], scenario_parameters$hoil_cie)
      out <- out * d$seds_fuel_price  # Model output is annual expenditure divided by state average price; this converts back to original units
    }
    
    # If heating fuel is Electricity, return zeros
    if (d$hfuel[1] %in% c("Electricity")) {
      out <- matrix(rep(0, 3 * nrow(d)), ncol = 3)
    }
    
    out <- cbind(d$id, out)
    colnames(out) <- c("id", "heat", "heat_upr", "heat_cie")
    return(out)
    
  }
  
  heat <- by(nd, nd$hfuel, predHeatModels)
  heat <- as.data.frame(do.call("rbind", heat))
  heat <- heat[order(heat$id), -1]
  heat[, 2] <- pmax(heat[, 2], 2 * heat[, 1])  # The upper expenditure value is set to minimum 2x the mean prediction
  heat <- signif(heat, 3)
  
  #-----
  
  # Annual cost equations
  
  # For each scenario, the 'cost' formula evaluates to:
  # Indirect "core" burden PLUS 
  # Household direct emissions (expend * cie)  multiplied by the average tax burden per kgCO2 direct emissions ('fct')
  
  # For each scenario, calculate the direct household tax burden ($) per kgCO2-eq (across all households)
  # This is multiplied by estimate of household direct emissions to arrive at a direct tax burden in dollars
  #fct1 <- signif(scenario_parameters$direct_cost1 / scenario_parameters$direct_kgco2, 3)  # Scenario 1
  #fct2 <- signif(scenario_parameters$direct_cost2 / scenario_parameters$direct_kgco2, 3)  # Scenario 2
  fct <- signif(scenario_parameters$direct_cost3 / scenario_parameters$direct_kgco2, 3)  # Scenario 3
  
  #---
  
  # Scenario 1 cost formula
  
  # nd$cost1 <- paste(
  #   
  #   # Burden due to sources other than direct emissions
  #   round(core1),
  #   
  #   # Burden due to direct emissions
  #   paste0(factor1, " * (gas * ", gas[, 3], " + elec * ", elec[, 3], " + heat * ", heat[, 3], ")"),
  #   
  #   sep = " + ")
  # 
  # #---
  # 
  # # Scenario 2 cost formula
  # 
  # nd$cost2 <- paste(
  #   
  #   # Burden due to sources other than direct emissions
  #   round(core2),
  #   
  #   # Burden due to direct emissions
  #   paste0(factor2, " * (gas * ", gas[,3], " + elec * ", elec[, 3], " + heat * ", heat[, 3], ")"),
  #   
  #   sep = " + ")
  
  #---
  
  # Cost formula
  
  nd$cost <- paste(
    
    # Burden due to sources other than direct emissions
    round(core),
    
    # Burden due to direct emissions
    paste0(fct, " * (gas * ", gas[,3], " + elec * ", elec[, 3], " + heat * ", heat[, 3], ")"),
    
    sep = " + ")
  
  
  #---
  
  # For users that did not know their heating fuel, replace the 'heat' component
  #  of cost equation with the default expenditure value and set 'heat' and 'heat_upr' variables to zero
  if (length(donotknow.id) > 0) {
    for (i in donotknow.id) {
      # nd$cost1[i] <- sub("heat", round(heat$heat[i]), nd$cost1[i], fixed = TRUE)
      # nd$cost2[i] <- sub("heat", round(heat$heat[i]), nd$cost2[i], fixed = TRUE)
      nd$cost[i] <- sub("heat", round(heat$heat[i]), nd$cost[i], fixed = TRUE)
    }
    heat$heat[donotknow.id] <- 0
    heat$heat_upr[donotknow.id] <- 0
  }
  
  # For users that report heating fuel as "Other or none", replace the 'heat' component
  #  of cost equation with ZERO value and set heat-related presets to zero
  if (length(other.id) > 0) {
    for (i in other.id) {
      # nd$cost1[i] <- sub("heat", 0, nd$cost1[i], fixed = TRUE)
      # nd$cost2[i] <- sub("heat", 0, nd$cost2[i], fixed = TRUE)
      nd$cost[i] <- sub("heat", 0, nd$cost[i], fixed = TRUE)
    }
  }
  heat[other.id,] <- 0
  
  #-----
  
  # Return results matrix
  psets <- cbind(gas, elec, heat)
  result <- cbind(nd, psets)
  result <- subset(result, select = c(div_pre, mrate, cost, moe, gas, elec, heat, gas_upr, elec_upr, heat_upr))
  
  return(result)
  
}

#----------------

# INPUT FOR TESTING:
# require(quantreg, quietly = TRUE)
# require(mgcv, quietly = TRUE)
# for (i in list.files("data", full.names = TRUE)) load(i)

# nd <- data.frame(zip = "80524", na = 2, nc = 2, hinc = 50e3, hfuel = "Natural gas", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Other or none", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = c("94062","80524","80521"), na = c(2, 1, 3), nc = c(2, 0, 3), hinc = c(50e3, 300e3, 100e3), hfuel = c("Do not know", "Natural gas", "Other or none"), veh = c(2, 0, 3), htype = c("Stand-alone house", "Apartment or condo", "Other"), stringsAsFactors = FALSE)
# nd <- data.frame(zip = "92626", na = 1, nc = 0, hinc = 100e3, hfuel = "Electricity", veh = 1, htype = "Apartment or condo", stringsAsFactors = FALSE)
# 
# # Manual check
# nd <- data.frame(zip = "36101", na = 2, nc = 0, hinc = 96400, hfuel = "Electricity", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = "52322", na = 3, nc = 1, hinc = 59100, hfuel = "Electricity", veh = 2, htype = "Other", stringsAsFactors = FALSE)
# 
# # Call prediction function
# test <- predictModel3(nd)
# 
# # Total monthly "cost" given by formula
# cost <- gsub("heat", "%d", gsub("elec", "%d", gsub("gas", "%d", test$cost)))
# cost <- sprintf(cost, as.integer(test$gas), as.integer(test$elec), as.integer(test$heat))
# cost <- as.numeric(sapply(cost, function(x) eval(parse(text = x)))) / 12  # Convert annual to monthly
# 
# # Monthly post-tax dividend
# div <- test$div_pre * (1 - test$mrate) / 12
# 
# # Net benefit/cost
# cost
# div
# div - cost

# To convert ?predictModel2 documentation to stand-alone HTML
#library(tools)
#Rd2HTML("man/predictModel.Rd", "Documentation.html")
