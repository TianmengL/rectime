#' Create the rectime data object
#'
#' Create the rectime data object used as input of other functions.
#'
#' @param id the column name for the subject id in the dataset.
#' @param type the column name for the variable which indicates the type of each observation:
#' 1: event visit; 2: regular visit; 0: censoring.
#' @param time the column name for the time variable.
#' @param data the dataset. Must include: a column for the subject id, a column for the
#' observation type, a column for the time variable and several columns for the covariates. Note
#' that we assume that time variable starts from 0 and it needs to be 'numeric'. The covariates
#' need to be 'numeric' as well.
#'
#' @details The rectime data object is a list of two datasets: d.event is for the event observations;
#' d.regular is for the observations at regular visits. The observations with a missing time variable
#' will be deleted: for d.event, the whole subject will be deleted if the time variable is
#' missing in one observation; for d.regular, only the observation with missing time variable
#' will be deleted. In addition, in d.event and d.regular, the column names for subject id,
#' type of observation and censoring time are fixed as: "id", "type", "cent".
#'
#' @return
#' A list of two datasets: d.event and d.regular.
#'
#' @examples
#' data(simdata)
#' recd <- recdata(id="id", type="type", time="time", data=simdata)


recdata <- function(id, type, time, data)
{
  if (!is.numeric(data[,time])) stop('time is not numeric')
  if (sum(!(data[,type] %in% c(0,1,2)))!=0) stop('type can only have values 0,1,2')

  idall <- unique(data[,id])
  N <- length(idall)
  data$cent <- 0
  #if a subject does not have a type=0 observation, then the censoring time is chosen as
  #the maximum of the largest event time and the largest regular visit time.
  for (i in 1:N) {
    censub <- data[which(data[,id]==idall[i] & data[,type]==0),time]
    if (length(censub)==0) {
      data[which(data[,id]==idall[i]),]$cent <- max(data[which(data[,id]==idall[i]),time])
    } else {data[which(data[,id]==idall[i]),]$cent <- censub}
  }
  d.event <- data[which(data[,type]==1),]
  d.regular <- data[which(data[,type]!=1),]
  colnames(d.event)[which(colnames(d.event)==id)] <- "id"
  colnames(d.event)[which(colnames(d.event)==type)] <- "type"
  colnames(d.regular)[which(colnames(d.regular)==id)] <- "id"
  colnames(d.regular)[which(colnames(d.regular)==type)] <- "type"

  #delete the missing observations: time variable
  #d.event
  del.id <- c()
  for (i in 1:N){
    eventsub <- as.vector(d.event[which(d.event[,id]==idall[i]), time])
    if (anyNA(eventsub)) del.id <- c(del.id, idall[i])
  }
  d.event <- d.event[which(!(d.event[,id] %in% del.id)),]
  #d.regular
  del.obs <- c()
  for (i in 1:nrow(d.regular)){
    obssub <- as.vector(d.regular[i, time])
    if (anyNA(obssub)) del.obs <- c(del.obs, i)
  }
  if (length(del.obs)!=0){
    d.regular <- d.regular[-del.obs,]
  }

  return(list(d.event,d.regular))
}
