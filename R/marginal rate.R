#' Estimate household marginal tax rate
#'
#' Approximaton of household marginal federal income tax rate using tax year
#' 2017 values. Assumes gross household income is equivalent to AGI and applies
#' standard deduction in all cases. Households with >2 adults are treated as
#' married with all adults >2 counted as dependents.
#'
#' @export

margRate <- function(d) {

  # Deductions, exemptions, and tax rates/brackets based on tax year 2017
  # Source: https://taxfoundation.org/2017-tax-brackets/

  # Estimate standard deduction and exemptions
  deduct <- rep(6350, nrow(d))
  deduct[d$nc > 0] <- 9350
  deduct[d$na >= 2] <- 12700
  exempt <- 4050 * (d$na + d$nc)

  # Calculate taxable income (assumes AGI = d$hinc)
  taxable <- d$hinc - deduct - exempt

  # Determine marginal Federal income tax rate
  rates <- c(0, 0.1, 0.15, 0.25, 0.28, 0.33, 0.35, 0.396)
  mrate <- 1
  mrate[deduct == 6350] <- findInterval(taxable[deduct == 6350], c(-Inf, 0, 9325, 37950, 91900, 191650, 416700, 418400, Inf))
  mrate[deduct == 9350] <- findInterval(taxable[deduct == 9350], c(-Inf, 0, 13350, 50800, 131200, 212500, 416700, 444500, Inf))
  mrate[deduct == 12700] <- findInterval(taxable[deduct == 12700], c(-Inf, 0, 18650, 75900, 153100, 233350, 416700, 470700, Inf))
  mrate <- rates[mrate]

  return(mrate)

}
