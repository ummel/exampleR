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
#'   that (when evaluated) returns annual policy cost given monthly average
#'   expenditure inputs for gasoline (gas), electricity (elec), and heating fuel
#'   (heat)} \item{gas: predicted average monthly gasoline expenditure for
#'   household (used as slider preset in online calculator)} \item{elec:
#'   predicted average monthly electricity expenditure for household (used as
#'   slider preset in online calculator)} \item{heat: predicted average monthly
#'   primary heating fuel expenditure for household (used as slider preset in
#'   online calculator; zero if not applicable)} \item{gas_upr: predicted
#'   maximum feasible monthly gasoline expenditure for household (used as slider
#'   max value in online calculator)} \item{elec_upr: predicted maximum feasible
#'   monthly electricity expenditure for household (used as slider max value in
#'   online calculator)} \item{heat_upr: predicted maximum feasible monthly
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

predictModel2 <- function(input) {
  
  # Assumed carbon price ($ per ton CO2)
  carbon.price <- 15
  
  # Assumed total tax burden
  total.tax <- 78.5 * 1e9
  
  # Indirect ("core") tax burden adjustment factor. This is the ratio of assumed tax burden
  # (total.tax) to analogous figure in original study ($93.4 billion). Since the
  # models make predictions based on original study, it is necessary to adjust
  # "core" emissions estimate to reflect fact that assumed tax base is smaller
  # than that calculated/used for original study (i.e. scale down all "core
  # emissions" model estimates)
  core.adj <- 78.5 / 93.4
  
  # Direct tax burden adjustment factor
  # This is based on calculation that total direct tax burden should be 3.5% lower than estimated using 2012 consumption levels and updated electricity CIE (see version 2 memo prepared for CCL)
  dir.adj <- 0.965
  
  # Adult pre-tax dividend amount
  # Hard-coded; originally based on original impact study for 2008-2012 period
  # Total revenue in analysis assumes $15 per ton at 2012 emission levels; $377 is based on those assumptions
  #div.adult <- 377.47  # Original study amount
  
  # REVISED June 22, 2019: Updated to reflect new assumption that total dividend amount is $73.3 billion and change in number of minors per household that are eligible for half-share
  # Revised (HR 7173) amount (see code in Projects/CCL/Scripts/Simulate household-level cost and dividend.R)
  div.adult <- 274
  
  #-----
  
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
  
  # If the direct emissions pass-through rate ('dirfrac') is not provided in 'input', set to default value of 0.95
  if ("dirfrac" %in% nms == FALSE) nd$dirfrac <- 0.950 else stopifnot(all(nd$dirfrac >= 0 & nd$dirfrac <= 1))
  
  # If the indirect emissions pass-through rate ('indfrac') is not provided in 'input', set to default value of 0.3466
  if ("indfrac" %in% nms == FALSE) nd$indfrac <- 0.495 else stopifnot(all(nd$indfrac >= 0 & nd$indfrac <= 1))
  
  # Calculate overall share of tax burden that falls on capital, assuming 45% of total emissions are direct and 55% indirect
  nd$capfrac <- 1 - (0.45 * nd$dirfrac + 0.55 * nd$indfrac)
  
  # Detect if the input is a test dataset
  # If the original training data is submitted, all 'zip' values will be NA
  test <- all(is.na(nd$zip))
  
  #-----
  
  nd$id <- 1:nrow(nd)
  
  # Assign geographic variables to 'nd' using zip code provided
  # If statement used to allow testing with the original dataset (which contains spatial variables but not zip code)
  if (!test) {
    nd <- merge(nd, zip_lookup, by = "zip", sort = FALSE, all.x = TRUE)
    if (any(is.na(nd$cents_kwh))) stop("Zip code not found.")
  }
  
  # Assign price adjustment factors based on assigned state
  # if() statement used to allow testing with the original dataset (no price adjustments necessary)
  if (!test) {
    nd <- merge(nd, price_adjustment, sort = FALSE)
  } else {
    temp <- price_adjustment
    temp[,-1] <- 1  # Sets all '_adjust' variables to 1
    nd <- merge(nd, temp, sort = FALSE)
  }
  
  # Ensure that original order to maintained; probably not necessry is sort = FALSE in merge()
  nd <- nd[order(nd$id),]
  
  # Estimate household marginal tax rate
  # This is based on 2017 rates and brackets using household income as reported as user
  nd$mrate <- margRate(nd)
  
  #-----
  
  # Calculate pre-tax dividend amount (2012 dollars)
  if (!test) {
    #nd$div_pre <- round(div.adult * nd$na + 0.5 * div.adult * pmin(2, nd$nc))  # Original study definition
    nd$div_pre <- round(div.adult * nd$na + 0.5 * div.adult * nd$nc)  # ADDED June 20, 2019: Revised (HR 7173) definition
  } else {
    # The "nc - 1" is necessary to account for fact that test dataset added 1 to both 'nc' and 'veh' to avoid errors with log() in model formulas
    # Rounding necessary to account for jitter added to test dataset numeric variables
    #nd$div_pre <- round(div.adult * round(nd$na) + 0.5 * div.adult * pmin(2, round(nd$nc - 1)))  # Original study definition
    nd$div_pre <- round(div.adult * round(nd$na) + 0.5 * div.adult * round(nd$nc - 1))  # ADDED June 20, 2019: Revised (HR 7173) definition
  }
  
  #-----
  
  # Fuel-price-to-income ratios used in model fitting
  # Note that fuel prices are adjusted from 2012 to current price levels before calculating ratios with current income
  # If statement will skip this if using test dataset (*_ratio variables already in the data frame)
  if (!test) {
    nd$elec_ratio <- nd$cents_kwh * nd$elec_adjust / (nd$hinc / 1e3)
    nd$gas_ratio <- nd$gasprice * nd$gas_adjust / (nd$hinc / 1e3)
    nd$ngas_ratio <- nd$ngasprice * nd$ngas_adjust / (nd$hinc / 1e3)
    nd$lpg_ratio <- nd$lpgprice * nd$lpg_adjust / (nd$hinc / 1e3)
    nd$hoil_ratio <- nd$hoilprice * nd$hoil_adjust / (nd$hinc / 1e3)
  }
  
  # Deflate user input household income to 2012 levels using CPI
  # This ensures that income matches currency units of original data and models
  nd$hinc <- nd$hinc / nd$cpi_adjust
  
  # Total number of people in household
  nd$np <- nd$na + nd$nc
  
  #----------------------
  
  # TEST:
  # Estimate household capital income
  nd$finval <- as.numeric(pmax(0, mgcv::predict.gam(finval_model_gam, newdata = nd)))
  
  # Calculate household tax burden associated with exposure to capital income
  # First term is total tax burden falling on capital; second terms is household-estimated share of national investment assets
  captax <- total.tax * nd$capfrac * nd$finval / (54.87 * 1e12)
  
  #----------------------
  
  # Core emissions
  # Calculate 90% CI margin of error (MOE) This is calculated by estimating the
  # standard deviation, assuming a Normal distribution with specified 25th and
  # 75th percentiles (the IQR is approximately 1.35x the standard deviation)
  
  # Add 1 to 'nc' and 'veh' to allow log() calls in model prediction
  # Not necessary when test dataset is being used
  if (!test) {
    nd$nc <- nd$nc + 1
    nd$veh <- nd$veh + 1
  }
  
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
  
  # ESTIMATE INDIRECT ("core") emissions and associated standard deviation, assuming Normal distribution
  # Note exp() used to convert log prediction value
  core <- as.numeric(exp(mgcv::predict.gam(core_model_gam, newdata = nd)))
  
  # TURNED OFF June 22, 2019 with requested revisions, as uncertainty not used by calculator and not calculated for 'captax' tax burden component
  # q <- exp(quantreg::predict.rq(core_model_rq, newdata = nd))
  # if (nrow(nd) == 1) {
  #   stdev <- (q[2] - q[1]) / 1.35
  # } else {
  #   stdev <- (q[,2] - q[,1]) / 1.35
  # }
  
  #----------------------
  
  # Predict typical (mean) and upper-bound (approximately 97.5th percentile) expenditure values used to set the "page 2" slider preset and maximum value
  # Note that all dollar values are adjusted to reflect current price levels ("_adjust" variables)
  # Note that the models output annual physical units consumed; converted to monthly expenditure at current prices
  
  # Gasoline monthly expenditure (mean and approximately 97.5th percentile)
  gas <- cbind(mgcv::predict.gam(gas_model_gam, newdata = nd), 1.1 * quantreg::predict.rq(gas_model_rq, newdata = nd))
  gas <- signif(gas * nd$gasprice * nd$gas_adjust * dir.adj / 12, digits = 2)
  colnames(gas) <- c("gas", "gas_upr")
  gas[which(nd$veh == 1), "gas"] <- 0  # Set predicted gasoline expenditure to zero if user inputs zero vehicles (which is actually "1" in nd$veh after +1 to original value (above)
  
  # Electricity monthly expenditure (mean and approximately 97.5th percentile)
  elec <- cbind(mgcv::predict.gam(elec_model_gam, newdata = nd), 1.1 * quantreg::predict.rq(elec_model_rq, newdata = nd))
  elec <- signif(elec * nd$cents_kwh * nd$elec_adjust * dir.adj / 12, digits = 2)
  colnames(elec) <- c("elec", "elec_upr")
  
  # Primary heating fuel monthly expenditure (median and approximately 95th percentile)
  # Note that expenditure values are adjusted to current price levels AND
  #  CIE values are adjusted to reflect change in fuel prices since 2012 (if price went up, CIE goes down)
  # Note that model physical consumption is converted to expenditure at current prices
  predHeatModels <- function(d) {
    
    if (d$hfuel[1] == "Natural gas") {
      if (!test) d$heat_ratio <- d$ngas_ratio
      out <- cbind(mgcv::predict.gam(ngas_model_gam, newdata = d), 1.1 * quantreg::predict.rq(ngas_model_rq, newdata = d)) * ifelse(!test, d$ngasprice, d$heatprice) * d$ngas_adjust
      out <- cbind(out, d$Natural_gas_cie / d$ngas_adjust)
    }
    
    if (d$hfuel[1] == "LPG/Propane") {
      if (!test) d$heat_ratio <- d$lpg_ratio
      out <- cbind(mgcv::predict.gam(lpg_model_gam, newdata = d), 1.1 * quantreg::predict.rq(lpg_model_rq, newdata = d)) * ifelse(!test, d$lpgprice, d$heatprice) * d$lpg_adjust
      out <- cbind(out, d$LPG_cie / d$lpg_adjust)
    }
    
    if (d$hfuel[1] == "Heating oil") {
      if (!test) d$heat_ratio <- d$hoil_ratio
      out <- cbind(mgcv::predict.gam(hoil_model_gam, newdata = d), 1.1 * quantreg::predict.rq(hoil_model_rq, newdata = d)) * ifelse(!test, d$hoilprice, d$heatprice) * d$hoil_adjust 
      out <- cbind(out, d$Heating_oil_cie / d$hoil_adjust)
    }
    
    # If heating fuel is Electricity, return zeros
    if (d$hfuel[1] %in% c("Electricity")) {
      out <- matrix(rep(0, 3 * nrow(d)), ncol = 3)
    }
    
    out <- cbind(d$id, out)
    colnames(out) <- c("id", "heat", "heat_upr", "heat_cie")
    out[,c("heat", "heat_upr")] <- signif(dir.adj * out[,c("heat", "heat_upr")] / 12, 2)  # Applies 'dir.adj' consumption adjustment to both 'heat' and 'heat_upr'
    return(out)
    
  }
  
  heat <- by(nd, nd$hfuel, predHeatModels)
  heat <- as.data.frame(do.call("rbind", heat))
  heat <- heat[order(heat$id), -1]
  
  #-----
  
  # Estimate 90% margin of error (in dollars, annual)
  # Since emissions/tax burden from fuel use is assumed to be accurate (based on user specified values),
  #  the overall MOE is the modeled uncertainty in the "core" emissions
  #nd$moe <- round(1.645 * stdev * carbon.price)
  
  #-----
  
  # Annual cost equation
  # NOTE that 'core' model estimate (indirect emissions) is adjusted by:
  #  1) 'core.adj' to reflect new assumption about total tax burden (less than in original study)
  #  2) and 'nd$indfrac', since only a portion of the tax burden is passed to consumers as higher prices for indirect emissions
  # Direct tax burden (gasoline, electricity, and heating fuel) is adjusted for passthrough by 'nd$dirfrac'
  # !!! Revised Aug 2, 2019 to reflect fact that only a portion of CIE is passed forward to consumer prices, based on values of 'dirfrac' and 'indfrac'
  # Input fuel expenditure values are MONTHLY, hence 12x applied to multiplicative term to make the resulting cost ANNUAL
  
  nd$cost <- paste(
    
    # Burden due to capital-related share of tax
    round(captax, 2),
    
    # Burden due to indirect household emissions
    round(core * core.adj * nd$indfrac * carbon.price, 2),
    
    # Burden due to gasoline consumption
    paste0("gas * ", signif(12 * nd$Gasoline_cie * nd$dirfrac / nd$gas_adjust * (carbon.price / 1e3), 4)),  # Input is monthly expenditure
    
    # Burden due to electricity consumption
    # Assumes use of electricity CIE value that is extrapolated to 2020
    # NOTE: No 'dir.adj' is applied to electricity, as the projected 2020 CIE is taken as accurate
    paste0("elec * ", signif(12 * nd$Electricity_cie * nd$dirfrac / nd$elec_adjust * (carbon.price / 1e3), 4)),  # Input is monthly expenditure
    
    # Burden due to heating fuel consumption (natural gas, LPG, etc.)
    # Note: Heating fuel CIE already divided by fuel-specific price change (above)
    paste0("heat * ", signif(12 * heat$heat_cie * nd$dirfrac * (carbon.price / 1e3), 4)), 
    
    sep = " + ")
  
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
  
  #-----
  
  # Return results matrix
  psets <- cbind(gas, elec, heat)
  psets[psets < 0] <- 0
  result <- cbind(nd, psets)
  #result <- subset(result, select = c(div_pre, mrate, cost, moe, gas, elec, heat, gas_upr, elec_upr, heat_upr))
  result <- subset(result, select = c(div_pre, mrate, cost, gas, elec, heat, gas_upr, elec_upr, heat_upr))  # 'moe' component removed June 22, 2019
  
  return(result)
  
}

#----------------

# INPUT FOR TESTING:
# require(quantreg, quietly = TRUE)
# require(mgcv, quietly = TRUE)
# source("R/marginal rate.R")
# for (i in list.files("data/", full.names = TRUE)) load(i)
# 
# nd <- data.frame(zip = "80524", na = 2, nc = 2, hinc = 50e3, hfuel = "Natural gas", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = "94062", na = 2, nc = 2, hinc = 50e3, hfuel = "Other or none", veh = 2, htype = "Stand-alone house", stringsAsFactors = FALSE)
# nd <- data.frame(zip = c("94062","80524","80521"), na = c(2, 1, 3), nc = c(2, 0, 3), hinc = c(50e3, 300e3, 100e3), hfuel = c("Do not know", "Natural gas", "Other or none"), veh = c(2, 0, 3), htype = c("Stand-alone house", "Apartment building", "Other"), stringsAsFactors = FALSE)
# 
# test <- predictModel2(nd)
# 
# # Total monthly cost given by formula
# cost <- gsub("heat", "%d", gsub("elec", "%d", gsub("gas", "%d", test$cost)))
# cost <- sprintf(cost, as.integer(test$gas), as.integer(test$elec), as.integer(test$heat))
# cost <- as.numeric(sapply(cost, function(x) eval(parse(text = x)))) / 12 
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
