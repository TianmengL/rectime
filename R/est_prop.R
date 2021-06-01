#dyn.load("./src/rectime.so")

score_prop <- function(beta, formula, d.event, d.regular, est.method, bandwidth, low, up, tau)
{
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)
  #beta is a 'p x 1' vector
  udt <- sort(unique(d.event[,time]))
  m1 <- table(d.event[,time])

  if (est.method=="kernel")
  {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,covar[2],"cent")]),ncol=4)

    regular.subs0 <- regular.subz1
    regular.subs0[,3] <- 0
    for (i in 1:p){
      regular.subs0[,3] <- regular.subs0[,3] + beta[i]*regular.sub[,covar[i+1]]
    }
    regular.subs0[,3] <- exp(regular.subs0[,3])
    s0.t <- function(t){
      if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
      }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
      } else {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
      return(ker$ker[1])
    }
    s0.t <- Vectorize(s0.t)

    score <- 0
    for (i in 1:p){
      s1.data <- regular.subz1
      s1.data[,3] <- regular.sub[,covar[i+1]]*regular.subs0[,3]
      s1.i.t <- function(t){
        if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(s1.data),dim=as.integer(nrow(s1.data)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
        }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(s1.data),dim=as.integer(nrow(s1.data)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
        } else {ker <- .Fortran("kernel",bcarr=as.single(s1.data),dim=as.integer(nrow(s1.data)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
        return(ker$ker[1])
      }
      s1.i.t <- Vectorize(s1.i.t)
      score <- score + (sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)] * m1[which(s0.t(udt)!=0)]))^2
    }
  }

  if (est.method=="ACCF")
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
    pdatas0 <- pdata
    pdatas0$s0 <- 0
    for (i in 1:p){
      pdatas0$s0 <- pdatas0$s0 + beta[i]*pdata[,covar[i+1]]
    }
    pdatas0$s0 <- exp(pdatas0$s0)

    s0.t <- function(t, start=pdatas0$start, stop=pdatas0$stop, x=pdatas0$s0){
      mean(x[t>=start & t<stop])}
    s0.t <- Vectorize(s0.t)

    pdatas1 <- pdata
    score <- 0
    for (i in 1:p){
      pdatas1[,covar[i+1]] <- pdatas1[,covar[i+1]]*pdatas0$s0
      s1.i.t <- function(t, start=pdatas1$start, stop=pdatas1$stop, x=pdatas1[,covar[i+1]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)
      score <- score + (sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)] * m1[which(s0.t(udt)!=0)]))^2
    }
  }

  if (est.method=="interp")
  {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    #estimation
    b.alls0 <- b.all
    b.alls0$s0 <- 0
    for (i in 1:p){
      b.alls0$s0 <- b.alls0$s0 + beta[i]*b.all[,covar[i+1]]
    }
    #b.alls0$s0 <- exp(b.alls0$s0)

    s0.t <- function(t){
      b.all.sub <- b.alls0[which(b.alls0$cent >= t),]
      b.allsub <- matrix(unlist(b.all.sub[,c("id",time,"s0","cent")]),ncol=4)
      inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
      out <- mean(exp(inter$yinte))
      return(out)}
    s0.t <- Vectorize(s0.t)

    score <- 0
    for (i in 1:p){
      s1.i.t <- function(t){
        b.all.sub1 <- b.all[which(b.all$cent >= t),]
        b.allsub1 <- matrix(unlist(b.all.sub1[,c("id",time,covar[i+1],"cent")]),ncol=4)
        inter1 <- .Fortran("inter",ball=as.single(b.allsub1),idall=as.single(unique(b.allsub1[,1])),dim=as.integer(nrow(b.allsub1)), length=as.integer(length(unique(b.allsub1[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub1[,1])))))

        b.all.sub2 <- b.alls0[which(b.alls0$cent >= t),]
        b.allsub2 <- matrix(unlist(b.all.sub2[,c("id",time,"s0","cent")]),ncol=4)
        inter2 <- .Fortran("inter",ball=as.single(b.allsub2),idall=as.single(unique(b.allsub2[,1])),dim=as.integer(nrow(b.allsub2)), length=as.integer(length(unique(b.allsub2[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub2[,1])))))
        out <- mean(inter1$yinte*exp(inter2$yinte), na.rm=TRUE)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)
      score <- score + (sum(d.event[,covar[i+1]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)] * m1[which(s0.t(udt)!=0)]))^2
    }
  }

  return(score)
}

est_prop <- function(formula, d.event, d.regular, method="kernel", bandwidth=NULL, low=0, up=NULL, tau=NULL, inte.low=-10, inte.up=10)
{
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)

  if (method=="kernel")
  {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,covar[2],"cent")]),ncol=4)
    if (is.null(bandwidth)){
      if (is.null(up)) up <- tau
      bandwidth <- select.con_prop1(regularz1=regular.subz1, low=low, up=up,tau=tau,N=N)*(N^(-1/3))
    }
  }

  if (p==1) {
      est <- optimize(f=score_prop, interval=c(inte.low,inte.up),formula=formula, d.event=d.event, d.regular=d.regular, est.method=method, bandwidth=bandwidth, low=low, up=up, tau=tau)$minimum
  } else {
      est <- optim(par=rep(0,p), fn=score_prop, formula=formula, d.event=d.event, d.regular=d.regular, est.method=method, bandwidth=bandwidth, low=low, up=up, tau=tau)$par
  }
  return(invisible(list(est=est,model="prop", method=method, bandwidth=bandwidth)))
}

select.con_prop1 <- function(regularz1, low, up,tau,N=N)
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

select.con_prop2 <- function(regularz1, regular.subs0, low, up,tau,N=N)
{
  s0.t.cv <- function(t,data,bandwidth){
    if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
    }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
    } else {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
    return(ker$ker[1])
  }

  s1.1.t.cv <- function(t,data,bandwidth){
    if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
    }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
    } else {ker <- .Fortran("kernel",bcarr=as.single(data),dim=as.integer(nrow(data)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
    return(ker$ker[1])
  }

  mse <- function(c,regularz1)
  {
    cvsum <- 0
    idall <- unique(regularz1[,1])
    for (i in 1:length(idall))
    {
      del <- which(regularz1[,1]==idall[i])
      if (length(del)>0) {
        for (k in 1:length(del))
        {
          datasub1 <- matrix(unlist(regularz1[-del,]),ncol=4)
          datasub2 <- matrix(unlist(regular.subs0[-del,]),ncol=4)
          datasub3 <- datasub1
          datasub3[,3] <- datasub1[,3]*datasub2[,3]
          cvsum <- cvsum + (regularz1[del[k],3]-s1.1.t.cv(regularz1[del[k],2],datasub3,bandwidth=c*N^(-1/5))/s0.t.cv(regularz1[del[k],2],datasub2,bandwidth=c*N^(-1/5)))^2
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

boot_prop <- function(formula, d.event, d.regular, nb=50, method, bandwidth, tau=NULL)
{
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
    invisible(capture.output(est <- est_prop(formula=formula, d.event=b.event, d.regular=b.regular, method=method, bandwidth=bandwidth, tau=tau)))
    bootest <- rbind(bootest, est[[1]])
    i <- i+1
  }
  bootsd <- apply(bootest, 2, sd)
  return(invisible(bootsd))
}

base_prop <- function(t, beta, formula, d.event, d.regular, method, bandwidth, tau=NULL)
{
  covar <- all.vars(formula)
  time <- covar[1]
  p <- length(covar) - 1

  if (t < min(d.event[,time])) {
    baseline <- 0 } else {
  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)

  gt <- function(t){
    centall <- c()
    for (i in 1:N)
    {
      centall <- c(centall, unique(c(d.regular[which(d.regular$id==idall[i]),]$cent, d.event[which(d.event$id==idall[i]),]$cent)))
    }
    sum(t <= centall)
  }
  gt <- Vectorize(gt)

  udt <- sort(unique(d.event[,time]))
  m1 <- table(d.event[,time])

  if (method=="kernel") {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,covar[2],"cent")]),ncol=4)

    regular.subs0 <- regular.subz1
    regular.subs0[,3] <- 0
    for (i in 1:p){
      regular.subs0[,3] <- regular.subs0[,3] + beta[i]*regular.sub[,covar[i+1]]
    }
    regular.subs0[,3] <- exp(regular.subs0[,3])
    s0.t <- function(t){
      if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
      }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
      } else {ker <- .Fortran("kernel",bcarr=as.single(regular.subs0),dim=as.integer(nrow(regular.subs0)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
      return(ker$ker[1])
    }
    s0.t <- Vectorize(s0.t)

    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)])
  }

  if (method=="ACCF") {
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

    pdatas0 <- pdata
    pdatas0$s0 <- 0
    for (i in 1:p){
      pdatas0$s0 <- pdatas0$s0 + beta[i]*pdata[,covar[i+1]]
    }
    pdatas0$s0 <- exp(pdatas0$s0)

    s0.t <- function(t, start=pdatas0$start, stop=pdatas0$stop, x=pdatas0$s0){
      mean(x[t>=start & t<stop])}
    s0.t <- Vectorize(s0.t)

    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)])
  }

  if (method=="interp") {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    #estimation
    b.alls0 <- b.all
    b.alls0$s0 <- 0
    for (i in 1:p){
      b.alls0$s0 <- b.alls0$s0 + beta[i]*b.all[,covar[i+1]]
    }
    #b.alls0$s0 <- exp(b.alls0$s0)

    s0.t <- function(t){
      b.all.sub <- b.alls0[which(b.alls0$cent >= t),]
      b.allsub <- matrix(unlist(b.all.sub[,c("id",time,"s0","cent")]),ncol=4)
      inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
      out <- mean(exp(inter$yinte))
      return(out)}
    s0.t <- Vectorize(s0.t)

    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)])
  }
    }
  return(invisible(baseline))
}

