#' Predict a simple GBM model
#'
#' Simple model with two predictor \code{X1} and \code{X2}
#'
#' @export
#' @importFrom gbm predict.gbm
#' @param input data passed on as \code{newdata} to \code{\link{predict.gbm}}
#' @examples mydata <- data.frame(
#'    X1=runif(3),
#'    X2=2*runif(3)
#' )
#' predictModel(mydata)

predictModel <- function(input){
  #input can either be csv file or data
  newdata <- if(is.character(input) && file.exists(input)){
    read.csv(input)
  } else {
    as.data.frame(input)
  }
  stopifnot("X1" %in% names(newdata))
  stopifnot("X2" %in% names(newdata))

  #newdata$age <- as.numeric(newdata$age)

  #tv_model is included with the package
  newdata$out <- as.vector(predict.gbm(fitted_model, newdata = newdata, n.trees = fitted_model$n.trees))
  return(newdata$out)
}
