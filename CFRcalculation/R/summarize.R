#' @title Data from simulation
#' @docType  data
#' @name summarize
#' @description Summarized survival data transformed by individual survival data from simulation.
#' @usage data(summarize)
#' @format A data frame with 151 observations on the following 8 variables.
#' \describe{
#' \item{id}{patient id number}
#' \item{date}{the cout time}
#'\item{n}{number of new admissions}
#' \item{d}{number of new deaths}
#' \item{c}{number of new discharged patients}
#' \item{N}{cumulative number of cases}
#' \item{D}{cumulative number of deaths}
#' \item{C}{cumulative number of discharged patients}
#' }
#' @examples
#' library(CFRcalculation)
#' data(summarize)
NULL
