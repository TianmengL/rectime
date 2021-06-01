#' Baseline estimation and plot
#'
#' Estimate baseline rates at equally spaced time points using the same method selected for
#' model estimation.
#'
#' @param fit the output from the rectime() function.
#' @param formula the same formula used in the rectime() function.
#' @param recdata the rectime data object.
#' @param tau the pre-specified time point such that the recurrent event process could be
#' potentially observed beyond tau. Usually it is the length of the study. If left unspecified,
#' the default is the maximum of the time variable.
#'
#' @param np number of time points where the baseline rate function is evaluated. The time
#' points are equally spaced between the minimum event time to tau.
#'

recbase <- function(fit, formula, recdata, tau=NULL, np=40)
{
  if (missing(fit)) stop('argument "fit" is missing, with no default value')
  if (missing(formula)) stop('argument "formula" is missing, with no default value')
  if (missing(recdata)) stop('argument "recdata" is missing, with no default value')

  d.event <- recdata[[1]]
  d.regular <- recdata[[2]]
  #data preparation
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  idall <- unique(c(d.regular$id, d.event$id))
  N <- length(idall)
  #delete the missing observations
  #in d.event, if one covariate value is missing in a row, delete the whole subject
  del.id <- c()
  for (i in 1:N){
    eventsub <- as.vector(d.event[which(d.event$id==idall[i]), covar[-1]])
    if (anyNA(eventsub)) del.id <- c(del.id, idall[i])
  }
  d.event <- d.event[which(!(d.event$id %in% del.id)),]
  #in d.regular, if one covariate value is missing in a row, just delete this row
  del.obs <- c()
  for (i in 1:nrow(d.regular)){
    obssub <- as.vector(d.regular[i, covar[-1]])
    if (anyNA(obssub)) del.obs <- c(del.obs, i)
  }
  d.regular <- d.regular[-del.obs,]

  if (is.null(tau)) {
    tau <- max(c(d.regular$cent,d.event$cent))
  }

  baseall <- c()
  tall <- seq(min(d.event[,time]), tau, length.out = np)
  if (fit$model=="prop") {
    beta <- fit$est
    if (fit$method=="kernel"){
       for (i in 1:np)
       {
         baseall <- c(baseall, base_prop(tall[i], beta=beta, formula=formula, d.event=d.event, d.regular=d.regular, method="kernel", bandwidth=as.numeric(fit$bandwidth), tau=tau))
       }
    } else {
      for (i in 1:np)
      {
        baseall <- c(baseall, base_prop(tall[i], beta=beta, formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, tau=tau))
      }
    }
  }

  if (fit$model=="add") {
    theta <- fit$est
    if (fit$method=="kernel"){
      for (i in 1:np)
      {
        baseall <- c(baseall, base_add(tall[i], theta=theta, formula=formula, d.event=d.event, d.regular=d.regular, method="kernel", bandwidth=as.numeric(fit$bandwidth), tau=tau))
      }
    } else {
      for (i in 1:np)
      {
        baseall <- c(baseall, base_add(tall[i], theta=theta, formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, tau=tau))
      }
    }
  }

  if (fit$model=="add-mul") {
    theta <- fit$est
    if (fit$method=="kernel"){
      for (i in 1:np)
      {
        baseall <- c(baseall, base_am(tall[i], theta=theta, formula=formula, d.event=d.event, d.regular=d.regular, method="kernel", bandwidth=as.numeric(fit$bandwidth), tau=tau))
      }
    } else {
      for (i in 1:np)
      {
        baseall <- c(baseall, base_am(tall[i], theta=theta, formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, tau=tau))
      }
    }
  }
  for (i in 1:length(baseall))
  {
    subbase <- baseall[1:i]
    baseall[i] <- max(subbase)
  }
  table.base <- cbind(tall, baseall)
  colnames(table.base) <- c("time", "baseline.est")
  return(invisible(table.base))
}
