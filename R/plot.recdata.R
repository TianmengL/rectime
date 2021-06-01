#' Plot the recurrent event data
#'
#' Plot the recurrent events and regular visits for several subjects for illustration.
#'
#' @param recdata a rectime data object generated from recdata().
#' @param time the column name for the time variable.
#' @param ids the vector of subject id that the user wants to plot. The default is NULL. If
#' ids=NULL, the first n.id subjects will be plotted.
#' @param n.id the number of subjects that the user wants to plot. The default is 5.
#' @param main the title of the plot. The default is NULL.
#' @param xlab the label of the x axis. The default is "time".
#' @param ylab the label of the y axis. The default is "id".
#'
#' @examples
#' data(simdata)
#' recd <- recdata(id="id", type="type", time="time", data=simdata)
#' plot.recdata(recd, time="time", n.id=5, xlab="time", ylab="id")


plot.recdata <- function(recdata, time, ids=NULL, n.id=5, main=NULL, xlab="time", ylab="id")
{
  d.event <- recdata[[1]]
  d.regular <- recdata[[2]]
  idall <- unique(c(d.regular$id, d.event$id))

  if (is.null(ids)) ids <- idall[1:n.id]

  subevent <- d.event[which(d.event$id %in% ids),c("id",time)]
  subregular <- d.regular[which(d.regular$id %in% ids),c("id",time,"cent")]
  for (i in 1:length(ids)) {
    if (nrow(subevent[which(subevent$id==ids[i]),])!=0) subevent[which(subevent$id==ids[i]),]$id <- i
    if (nrow(subregular[which(subregular$id==ids[i]),])!=0) subregular[which(subregular$id==ids[i]),]$id <- i
  }
  k <- length(ids)

  plot(subevent[,2], subevent[,1], pch=4, col="red", main=main, xlab=xlab, ylab=ylab,
        yaxt='n', xlim=c(0,max(subevent[,2],subregular[,2])+0.1), ylim=c(1-0.1,k+0.1))
  axis(2, at = 1:k, labels=ids,cex.axis=0.2)
  points(subregular[,2], subregular[,1],pch=20,col="black")
  for (i in 1:k)
  {
    centi <- unique(subregular[subregular$id==i,]$cent)
    start <- min(c(subregular[subregular$id==i,]$time,subevent[subevent$id==i,]$time))
    segments(0, i, centi, i)
    points(centi,i,pch=1,col="blue")
  }
}
