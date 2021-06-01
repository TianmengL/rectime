#' Fit the marginal rates models
#'
#' Fit one of the following rates models for recurrent event data: proportional rates model, additive rates model,
#' additive-multiplicative rates model.
#'
#' @param formula a formula object. The left side of the operator '~' should be the time
#' variable in the dataset. The right side should be the terms: for the proportional rates
#' model and the additive rates model, the covariates included in the model are simply connected
#' by '+'; for the additive-multiplicative rates model, the covariates included in the
#' multiplicative part come first and are separated from the covarites included in the
#' additive part by '|'. See the 'Examples' for details.
#' @param recdata is the rectime data object generated from the recdata() function.
#' @param model the regression model to fit the data. Possible choices include: "prop" for
#' the proportional rates model; "add" for the additive rates model; "add-mul" for the
#' additive-multiplicative rates model.
#' @param method the method to deal with the intermittently observed time-dependent covariates.
#' Possible choices include: "kernel" for the kernel smoothing method; "ACCF" for the last
#' covariate carry forward method; "interp" for the linear interpolation method. The default
#' is "kernel".
#' @param SE is TRUE then the standard error estimated by bootstrap is returned. The default
#' is FALSE.
#' @param bandwidth the bandwidth used in the "kernel" method. If left unspecified, the
#' bandiwidth will be chosen automatically. See 'Details' for more information.
#' @param low the lower bound of the interval for the searching of the constant 'C' for the
#' automatically chosen bandwidth 'C*N^(-1/3)'. If left unpecified, the default is 0.
#' @param up the upper bound of the interval for the searching of the constant 'C' for the
#' automatically chosen bandwidth 'C*N^(-1/3)'. If left unpecified, the default is tau.
#' @param tau the pre-specified time point such that the recurrent event process could be
#' potentially observed beyond tau. Usually it is the length of the study. If left unspecified,
#' the default is the maximum of the time variable.
#' @param nb the resampling size of bootstrap. The default is 50.
#'
#' @details
#' When there are missing values in the covariates included in the model: for d.event, the whole
#' subject will be deleted if at least one covariate is missing in one observation; for
#' d.regular, only the observation with missing covariate values will be deleted.
#'
#'
#' @examples
#' data(simdata)
#' recd <- recdata(id="id", type="type", time="time", data=simdata)
#' #proportional rates model
#' rectime(time~z1+z2, recdata=recd, model="prop", method="kernel",bandwidth=1)
#' #additive rates model
#' rectime(time~z1+z2, recdata=recd, model="add", method="kernel",bandwidth=1)
#' #additive-multiplicative rates model
#' rectime(time~z1|z2, recdata=recd, model="add-mul", method="kernel",bandwidth=1)
#'


rectime <- function(formula, recdata, model, method="kernel", SE=FALSE, bandwidth=NULL, low=0, up=NULL, tau=NULL, nb=50)
{
  if (missing(formula)) stop('argument "formula" is missing, with no default value')
  #  if (missing(id)) stop('argument "id" is missing, with no default value')
  #  if (missing(type)) stop('argument "type" is missing, with no default value')
  if (missing(recdata)) stop('argument "recdata" is missing, with no default value')
  if (missing(model)) stop('argument "model" is missing, with no default value')

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

  if (model=="prop") {
    fit <- est_prop(formula=formula, d.event=d.event, d.regular=d.regular, method=method, bandwidth=bandwidth, low=low, up=up, tau=tau, inte.low=-10, inte.up=10)
    if (SE==TRUE) {
      if (method=="kernel"){
        se <- boot_prop(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, bandwidth=as.numeric(fit$bandwidth), tau=tau)
      } else {
        se <- boot_prop(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, tau=tau)
      }
      fit$SE <- se
      #output
      cat("Fitted proportional rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=3)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f", round(fit$est,4))
      output[,3] <- sprintf("%.4f", round(se,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.","SE")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    } else {
      #output
      cat("Fitted proportional rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=2)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f",round(fit$est,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    }
  }
  if (model=="add") {
    fit <- est_add(formula=formula, d.event=d.event, d.regular=d.regular, method=method, bandwidth=bandwidth, low=low, up=up, tau=tau)
    if (SE==TRUE) {
      if (method=="kernel"){
        se <- boot_add(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, bandwidth=as.numeric(fit$bandwidth), tau=tau)
      } else {
        se <- boot_add(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, tau=tau)
      }
      fit$SE <- se
      #output
      cat("Fitted additive rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=3)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f",round(fit$est,4))
      output[,3] <- sprintf("%.4f",round(se,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.","SE")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    } else {
      #output
      cat("Fitted additive rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=2)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f",round(fit$est,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    }
  }
  if (model=="add-mul") {
    time <- all.vars(formula[[2]])
    mulcovar <- all.names(formula[[3]][2])
    addcovar <- all.names(formula[[3]][3])
    mulcovar <- mulcovar[((length(mulcovar)-1)/2+1):length(mulcovar)]
    addcovar <- addcovar[((length(addcovar)-1)/2+1):length(addcovar)]
    p <-  length(mulcovar) + length(addcovar)

    fit <- est_am(formula=formula, d.event=d.event, d.regular=d.regular, method=method, bandwidth=bandwidth, low=low, up=up, tau=tau)
    if (SE==TRUE) {
      if (method=="kernel"){
        se <- boot_am(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, bandwidth=as.numeric(fit$bandwidth), tau=tau)
      } else {
        se <- boot_am(formula=formula, d.event=d.event, d.regular=d.regular, method=fit$method, nb=nb, tau=tau)
      }
      fit$SE <- se
      #output
      cat("Fitted additive-multiplicative rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=3)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f",round(fit$est,4))
      output[,3] <- sprintf("%.4f",round(se,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.","SE")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    } else {
      #output
      cat("Fitted additive-multiplicative rates rates model: ")
      print(formula)
      cat("Estimated coefficients:\n")
      output <- matrix(NA, nrow=p, ncol=2)
      output[,1] <- covar[-1]
      output[,2] <- sprintf("%.4f",round(fit$est,4))
      output <- as.data.frame(output)
      colnames(output) <- c("Covariate","Est.")
      print(output)
      cat("Number of subjects: ", N, "\n")
      cat("Average number of events per subject: ", nrow(d.event)/N, "\n")
      cat("Average number of regular visits per subject: ", nrow(d.regular)/N, "\n")
    }
  }
  return(invisible(fit))
}




