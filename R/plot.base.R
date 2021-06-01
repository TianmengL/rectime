#' Plot the estimated baseline function
#'
#' Plot the baseline function estimated by recbase().
#'
#' @param base the output from recbase().
#' @param main title of the plot.
#' @param xlab label of the x axis.
#' @param ylab label of the y axis.

plot.base <- function(base, main=NULL, xlab="time", ylab="baseline")
{
  plot(base[,1],base[,2],type="l",main=main,xlab=xlab,ylab=ylab)
}
