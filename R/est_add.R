#dyn.load("./src/rectime.so")

est_add <- function(formula, d.event, d.regular, method="kernel", bandwidth=NULL, low=NULL, up=NULL, tau=NULL)
{
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  udt <- sort(unique(d.event[,time]))
  m1 <- table(d.event[,time])

  inteup <- c()
  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)
  for (i in 1:length(idall)){
    inteup <- c(inteup, unique(c(d.regular$cent[which(d.regular$id==idall[i])],d.event$cent[which(d.event$id==idall[i])])))
  }
  inteup <- sort(inteup)

  if (method=="kernel")
  {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,covar[2],"cent")]),ncol=4)
    if (is.null(bandwidth)){
      if (is.null(up)) up <- tau
      bandwidth <- select.con_add(regularz1=regular.subz1, low=low, up=up,tau=tau,N=N)*(N^(-1/3))
    }

    part1 <- c()
    part2 <- matrix(NA, nrow=p, ncol=p)
    for (i in 1:p){
      subzi <- regular.subz1
      subzi[,3] <- regular.sub[,covar[i+1]]
      s1.i.t <- function(t){
        if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
        }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
        } else {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
        return(ker$ker[1])
      }
      s1.i.t <- Vectorize(s1.i.t)
      part1.i <- sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[!is.na(s1.i.t(udt))] * m1[!is.na(s1.i.t(udt))])
      part1 <- c(part1, part1.i)

      for (j in i:p){
        subzj <- regular.subz1
        subzj[,3] <- regular.sub[,covar[j+1]]
        s1.j.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s1.j.t <- Vectorize(s1.j.t)

        subzij <- regular.subz1
        subzij[,3] <- regular.sub[,covar[i+1]]*regular.sub[,covar[j+1]]
        s2.ij.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          u <- (s2.ij.t(t) - (s1.i.t(t))*(s1.j.t(t)))
          return(u)  }
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part2[i,j] <- part2[j,i] <- intesum.ij
      }
    }
    est <- as.vector(part1 %*% solve(part2))
  }

  if (method=="ACCF")
  {
    #prepare the data for ACCF
    pdata <- NULL
    for( k in 1:N){
      foo <- NULL
      foo1 <- d.regular[d.regular$id==idall[k],]
      foo2 <- d.event[d.event$id==idall[k],]

      if( nrow(foo1)>0 | nrow(foo2)>0 ){
        tallsub <- sort(unique(c(foo1[,time], foo2[,time])))
        centime <- unique(foo1$cent)
        tallsub <- c(tallsub, max((max(tallsub)+10^(-10)),centime))
        foo <- data.frame(cbind(id=idall[k], start=tallsub[-length(tallsub)],
                                stop=tallsub[-1]))
        for (i in 1:p)
        {
          foo[,covar[i+1]] <- sapply(foo$start,
                           function(s, t1=foo1[,time], x1=foo1[,covar[i+1]], t2=foo2[,time],
                                    x2=foo2[,covar[i+1]]){
                             if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
          }
                pdata <- rbind(pdata, foo)
      }
    }

    #estimation
    part1 <- c()
    part2 <- matrix(NA, nrow=p, ncol=p)
    for (i in 1:p){
      s1.i.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,covar[i+1]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)
      part1.i <- sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[!is.na(s1.i.t(udt))] * m1[!is.na(s1.i.t(udt))])
      part1 <- c(part1, part1.i)

      for (j in i:p){
        s1.j.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,covar[j+1]]){
          mean(x[t>=start & t<stop])}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,covar[i+1]]*pdata[,covar[j+1]]){
          mean(x[t>=start & t<stop])}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          u <- (s2.ij.t(t) - (s1.i.t(t))*(s1.j.t(t)))
          return(u)  }
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part2[i,j] <- part2[j,i] <- intesum.ij
      }
    }
    est <- as.vector(part1 %*% solve(part2))
  }

  if (method=="interp")
  {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    #estimation
    part1 <- c()
    part2 <- matrix(NA, nrow=p, ncol=p)
    for (i in 1:p){
      s1.i.t <- function(t){
        b.all.sub <- b.all[which(b.all$cent >= t),]
        b.allsub <- matrix(unlist(b.all.sub[,c("id",time,covar[i+1],"cent")]),ncol=4)
        inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
        out <- mean(inter$yinte)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)
      part1.i <- sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[!is.na(s1.i.t(udt))] * m1[!is.na(s1.i.t(udt))])
      part1 <- c(part1, part1.i)

      for (j in i:p){
        s1.j.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub <- matrix(unlist(b.all.sub[,c("id",time,covar[j+1],"cent")]),ncol=4)
          inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
          out <- mean(inter$yinte)
          return(out)}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub1 <- matrix(unlist(b.all.sub[,c("id",time,covar[i+1],"cent")]),ncol=4)
          inter1 <- .Fortran("inter",ball=as.single(b.allsub1),idall=as.single(unique(b.allsub1[,1])),dim=as.integer(nrow(b.allsub1)), length=as.integer(length(unique(b.allsub1[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub1[,1])))))
          b.allsub2 <- matrix(unlist(b.all.sub[,c("id",time,covar[j+1],"cent")]),ncol=4)
          inter2 <- .Fortran("inter",ball=as.single(b.allsub2),idall=as.single(unique(b.allsub2[,1])),dim=as.integer(nrow(b.allsub2)), length=as.integer(length(unique(b.allsub2[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub2[,1])))))
          out <- mean(inter1$yinte*inter2$yinte, na.rm=TRUE)
          return(out)}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          u <- (s2.ij.t(t) - (s1.i.t(t))*(s1.j.t(t)))
          return(u)  }
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part2[i,j] <- part2[j,i] <- intesum.ij
      }
    }
    est <- as.vector(part1 %*% solve(part2))
  }

  return(invisible(list(est=est,model="add", method=method, bandwidth=bandwidth)))
}

boot_add <- function(formula, d.event, d.regular, method, nb=50, bandwidth, tau=NULL)
{
#  if (missing(bandwidth)) stop('argument "bandwidth" is missing, with no default value')
  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)
  bootest <- NULL
  i=1
  while (i<=nb){
    #print(i)
    boot.id <- sample(idall, replace=T)
    b.event <- b.regular <- NULL
    for(k in 1:N){
      if(any(d.event$id == boot.id[k]) ){
        foo1 <- d.event[d.event$id == boot.id[k],]
        foo1$id <- k
        b.event <- rbind( b.event, foo1 )
      }
      if(any(d.regular$id == boot.id[k]) ){
        foo2 <- d.regular[d.regular$id == boot.id[k],]
        foo2$id <- k
        b.regular <- rbind( b.regular, foo2 )
      }
    }
    invisible(capture.output(est <- est_add(formula=formula, d.event=b.event, d.regular=b.regular, method=method, bandwidth=bandwidth,tau=tau)))
    bootest <- rbind(bootest, est[[1]])
    i <- i+1
  }
  bootsd <- apply(bootest, 2, sd)
  return(invisible(bootsd))
}


select.con_add <- function(regularz1, low, up,tau,N)
{
  s1.1.t.cv <- function(data,t,bandwidth){
    if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
    }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
    } else {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
    return(ker$ker[1])
  }
  mse <- function(c,regularz1)
  {
    cvsum <- 0
    ids <- unique(regularz1[,1])
    for (i in 1:length(ids))
    {
      del <- which(regularz1[,1]==ids[i])
      if (length(del)>0) {
        for (k in 1:length(del))
        {
          datasub <- matrix(unlist(regularz1[-del,]),ncol=4)
          cvsum <- cvsum + (regularz1[del[k],3]-s1.1.t.cv(datasub,regularz1[del[k],2],bandwidth=c*N^(-1/5)))^2
        }
      }
    }
    return(cvsum)
  }
  sel <- optimize(mse, interval=c(low,up),regularz1=regularz1)$minimum
  if (abs(sel-up) < 0.1) {
    up2 <- up
    while (abs(sel-up2) < 0.1)
    {
      up2 <- up+10
      sel <- optimize(mse, interval=c(up,up2),regularz1=regularz1)$minimum
      up <- up2
    }
  }
  return(sel)
}

base_add <- function(t, theta, formula, d.event, d.regular, method, bandwidth, tau=NULL)
{
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  if (t < min(d.event[,time])) {
    baseline <- 0 } else {
  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)

  #estimation
  udt <- sort(unique(d.event[,time]))
  m1 <- table(d.event[,time])

  gt <- function(t){
    centall <- c()
    for (i in 1:N)
    {
      centall <- c(centall, unique(c(d.regular[which(d.regular$id==idall[i]),]$cent, d.event[which(d.event$id==idall[i]),]$cent)))
    }
    sum(t <= centall)
  }
  gt <- Vectorize(gt)

  if (method=="kernel") {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,covar[2],"cent")]),ncol=4)

    part1 <- c()
    for (i in 1:p){
      subzi <- regular.subz1
      subzi[,3] <- regular.sub[,covar[i+1]]
      s1.i.t <- function(t){
        if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
        }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
        } else {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
        return(ker$ker[1])
      }
      s1.i.t <- Vectorize(s1.i.t)
      part1 <- c(part1, integrate(f=s1.i.t, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
    }
    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/gt(udtsub)*m1sub) - part1 %*% theta
  }

  if (method=="ACCF") {
    #prepare the data for ACCF
    pdata <- NULL
    for( k in 1:N){
      foo <- NULL
      foo1 <- d.regular[d.regular$id==idall[k],]
      foo2 <- d.event[d.event$id==idall[k],]

      if( nrow(foo1)>0 | nrow(foo2)>0 ){
        tallsub <- sort(unique(c(foo1[,time], foo2[,time])))
        centime <- unique(foo1$cent)
        tallsub <- c(tallsub, max((max(tallsub)+10^(-10)),centime))
        foo <- data.frame(cbind(id=idall[k], start=tallsub[-length(tallsub)],
                                stop=tallsub[-1]))
        for (i in 1:p)
        {
          foo[,covar[i+1]] <- sapply(foo$start,
                                     function(s, t1=foo1[,time], x1=foo1[,covar[i+1]], t2=foo2[,time],
                                              x2=foo2[,covar[i+1]]){
                                       if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
        }
        pdata <- rbind(pdata, foo)
      }
    }

    part1 <- c()
    for (i in 1:p){
      s1.i.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,covar[i+1]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)
      part1 <- c(part1, integrate(f=s1.i.t, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
    }
    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/gt(udtsub)*m1sub) - part1 %*% theta
  }

  if (method=="interp") {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    part1 <- c()
    for (i in 1:p){
      s1.i.t <- function(t){
        b.all.sub <- b.all[which(b.all$cent >= t),]
        b.allsub <- matrix(unlist(b.all.sub[,c("id",time,covar[i+1],"cent")]),ncol=4)
        inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
        out <- mean(inter$yinte)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)
      part1 <- c(part1, integrate(f=s1.i.t, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
    }
    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/gt(udtsub)*m1sub) - part1 %*% theta
  }
    }
  return(invisible(baseline))
}



