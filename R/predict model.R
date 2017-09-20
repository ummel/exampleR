#' Predict a simple GBM model
#'
#' Simple model with multiple predictors
#'
#' @export
#' @importFrom gbm predict.gbm
#' @importFrom gbm plot.gbm
#' @param input data passed on as \code{newdata} to \code{\link{predict.gbm}}
#' @examples mydata <- data.frame(
#'    X1=runif(3),
#'    X2=2*runif(3)
#' )
#' predictModel(mydata)

predictModel <- function(input) {

  #input can either be csv file or data
  nd <- if (is.character(input) && file.exists(input)) {
    read.csv(input)
  } else {
    as.data.frame(input)
  }
  stopifnot("zip" %in% names(nd))
  stopifnot("hhsize" %in% names(nd))
  stopifnot("minors" %in% names(nd))
  stopifnot("age" %in% names(nd))

  #-----

  # INPUT FOR TESTING:
  #nd <- data.frame(state = "Texas", hhsize = 4, minors = 2, age = 50, income = 50e3, elec = 100)
  #nd <- data.frame(zip = "80524", hhsize = 4, minors = 2, age = 50)

  # Assign geographic variables to 'nd' using zip code provided
  nd <- merge(nd, zip_lookup)

  # Adult pre-tax dividend amount
  div.adult <- 377

  # Estimate household marginal tax rate
  nd$adults <- nd$hhsize - nd$minors
  nd$mrate <- 0.15

  # Calculate pre- and post-tax dividends
  nd$div_pre <- round(div.adult * nd$adults + 0.5 * div.adult * pmin(2, nd$minors))
  nd$div_post <- round(nd$div_pre * (1 - nd$mrate))

  #------

  # Predict the slider preset values

  for (x in ls(pattern = "fitted_model")[-1]) {
    m <- get(x)
    nd[[m$response.name]] <- round(predict.gbm(m, newdata = nd, n.trees = m$n.trees))
  }

  #------

  # How does variation in income affect the result
  slide <- plot.gbm(fitted_model0, i.var = c("income", "elec"), n.trees = fitted_model$n.trees, return.grid = TRUE, continuous.resolution = 50)

  f <- "y ~ income + elec + income:elec"
  fit0 <- lm(formula = formula(f), data = slide)

  # Create a character vector giving the equation expression to be evaluated using the slider variables as inputs
  m <- prettyNum(signif(coef(fit0), 4))
  names(m)[1] <- "1"
  eq <- gsub(" * 1", "", gsub(":", " * ", paste(sprintf("%s * %s", m, names(m)), collapse = " + ")), fixed = TRUE)

  #------

  # pred <- setdiff(names(test), "y")

  # fit0 <- lm(y ~ ., data = test)

  # fit0 <- lm(y ~ poly(income, 1, raw = T) * poly(elec, 1, raw = T), data = test)

  #fit <- lm(formula = formula(paste0("y ~ ", paste0("poly(", pred, ", degree = 3, raw = TRUE)"))), data = test)

  #fit2 <- as.stepfun(isoreg(x = test$income, y = test$y))

  #library(MonoPoly)
  #fit3 <- monpol(y ~income, data = test, degree = 5, monotone = "increasing")

  # Check fit
  # plot(test)
  # lines(x = test$income, y = predict(fit0, data.frame(income = test$income)), col = 2)
  # lines(x = test$income, y = predict(fit0, data.frame(income = test$income)), col = 3)

  #-----

  out <- subset(nd, select = c(mrate, div_pre, div_post))
  out$cost <- eq

  return(out)

}
