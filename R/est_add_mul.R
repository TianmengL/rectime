#dyn.load("./src/rectime.so")
#to use score_am, there has to be at least one covariate in the additive part and one
#covariate in the multiplicative part

score_am <- function(theta, formula, d.event, d.regular, est.method, bandwidth, low, up, tau)
{
  time <- all.vars(formula[[2]])
  mulcovar <- all.names(formula[[3]][2])
  addcovar <- all.names(formula[[3]][3])
  mulcovar <- mulcovar[((length(mulcovar)-1)/2+1):length(mulcovar)]
  addcovar <- addcovar[((length(addcovar)-1)/2+1):length(addcovar)]
  pmul <- length(mulcovar)
  padd <- length(addcovar)
  beta <- theta[1:pmul]
  gamma <- theta[(pmul+1):(pmul+padd)]

  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)

  udt <- sort(unique(d.event[,time]))
  m1 <- table(d.event[,time])

  inteup <- c()
  idall <- unique(c(d.regular$id, d.event$id))
  for (i in 1:length(idall))
  {
    inteup <- c(inteup, unique(c(d.regular$cent[which(d.regular$id==idall[i])],d.event$cent[which(d.event$id==idall[i])])))
  }
  inteup <- sort(inteup)

  if (est.method=="kernel")
  {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,addcovar[1],"cent")]),ncol=4)

    event.subsx0m <- matrix(unlist(d.event[,c("id",time,addcovar[1],"cent")]),ncol=4)

    regular.subsx0 <- regular.subz1
    regular.subsx0m <- regular.subz1
    regular.subsx0[,3] <- 0
    event.subsx0m[,3] <- 0
    for (i in 1:pmul)
    {
      regular.subsx0[,3] <- regular.subsx0[,3] + beta[i]*regular.sub[,mulcovar[i]]
      event.subsx0m[,3] <- event.subsx0m[,3] + beta[i]*d.event[,mulcovar[i]]
    }
    regular.subsx0m[,3] <- exp(-regular.subsx0[,3])
    regular.subsx0[,3] <- exp(regular.subsx0[,3])
    event.subsx0m[,3] <- exp(-event.subsx0m[,3])

    s0.t <- function(t){
      if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
      }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
      } else {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
      return(ker$ker[1])
    }
    s0.t <- Vectorize(s0.t)

    part11 <- c()
    part12 <- matrix(NA, nrow=padd, ncol=padd)
    for (i in 1:padd)
    {
      subzi <- regular.subz1
      subzi[,3] <- regular.sub[,addcovar[i]]
      s1.i.t <- function(t){
        if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
        }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
        } else {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
        return(ker$ker[1])
      }
      s1.i.t <- Vectorize(s1.i.t)

      part11.i <- sum(d.event[,addcovar[i]]*event.subsx0m[,3]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part11 <- c(part11, part11.i)

      for (j in i:padd)
      {
        subzj <- regular.subz1
        subzj[,3] <- regular.sub[,addcovar[j]]
        s1.j.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s1.j.t <- Vectorize(s1.j.t)

        subzij <- regular.subz1
        subzij[,3] <- regular.sub[,addcovar[i]]*regular.sub[,addcovar[j]]*regular.subsx0m[,3]
        s2.ij.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part12[i,j] <- part12[j,i] <- intesum.ij
      }
    }

    part21 <- c()
    part22 <- matrix(NA, nrow=pmul, ncol=padd)
    for (i in 1:pmul)
    {
      subzi <- regular.subz1
      subzi[,3] <- regular.sub[,mulcovar[i]]*regular.subsx0[,3]
      s1.i.t <- function(t){
        if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
        }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
        } else {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
        return(ker$ker[1])
      }
      s1.i.t <- Vectorize(s1.i.t)

      part21.i <- sum(d.event[,mulcovar[i]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part21 <- c(part21, part21.i)

      for (j in 1:padd)
      {
        subzj <- regular.subz1
        subzj[,3] <- regular.sub[,addcovar[j]]
        s1.j.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzj),dim=as.integer(nrow(subzj)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s1.j.t <- Vectorize(s1.j.t)

        subzij <- regular.subz1
        subzij[,3] <- regular.sub[,mulcovar[i]]*regular.sub[,addcovar[j]]
        s2.ij.t <- function(t){
          if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
          }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
          } else {ker <- .Fortran("kernel",bcarr=as.single(subzij),dim=as.integer(nrow(subzij)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
          return(ker$ker[1])
        }
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part22[i,j] <- intesum.ij
      }
    }

    part1 <- part11 - part12 %*% gamma
    part2 <- part21 - part22 %*% gamma
    res <- t(part1) %*% part1 + t(part2) %*% part2
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
        for (i in 1:pmul)
        {
          foo[,mulcovar[i]] <- sapply(foo$start,
                                     function(s, t1=foo1[,time], x1=foo1[,mulcovar[i]], t2=foo2[,time],
                                              x2=foo2[,mulcovar[i]]){
                                       if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
        }
        for (i in 1:padd)
        {
          foo[,addcovar[i]] <- sapply(foo$start,
                                      function(s, t1=foo1[,time], x1=foo1[,addcovar[i]], t2=foo2[,time],
                                               x2=foo2[,addcovar[i]]){
                                        if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
        }
        pdata <- rbind(pdata, foo)
      }
    }

    #estimation
    pdatas0 <- pdatas0m <- pdata
    event.subsx0m <- matrix(unlist(d.event[,c("id",time,addcovar[1],"cent")]),ncol=4)
    pdatas0$s0 <- event.subsx0m[,3] <- 0
    for (i in 1:pmul){
      pdatas0$s0 <- pdatas0$s0 + beta[i]*pdata[,mulcovar[i]]
      event.subsx0m[,3] <- event.subsx0m[,3] + beta[i]*d.event[,mulcovar[i]]
    }
    pdatas0m$s0 <- exp(-pdatas0$s0)
    pdatas0$s0 <- exp(pdatas0$s0)
    event.subsx0m[,3] <- exp(-event.subsx0m[,3])

    s0.t <- function(t, start=pdatas0$start, stop=pdatas0$stop, x=pdatas0$s0){
      mean(x[t>=start & t<stop])}
    s0.t <- Vectorize(s0.t)

    part11 <- c()
    part12 <- matrix(NA, nrow=padd, ncol=padd)
    for (i in 1:padd)
    {
      s1.i.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,addcovar[i]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)

      part11.i <- sum(d.event[,addcovar[i]]*event.subsx0m[,3]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part11 <- c(part11, part11.i)

      for (j in i:padd)
      {
        s1.j.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,addcovar[j]]){
          mean(x[t>=start & t<stop])}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,addcovar[i]]*pdata[,addcovar[j]]*pdatas0m$s0){
          mean(x[t>=start & t<stop])}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part12[i,j] <- part12[j,i] <- intesum.ij
      }
    }

    part21 <- c()
    part22 <- matrix(NA, nrow=pmul, ncol=padd)
    pdatas1 <- pdata
    for (i in 1:pmul)
    {
      pdatas1[,mulcovar[i]] <- pdatas1[,mulcovar[i]]*pdatas0$s0
      s1.i.t <- function(t, start=pdatas1$start, stop=pdatas1$stop, x=pdatas1[,mulcovar[i]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)

      part21.i <- sum(d.event[,mulcovar[i]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part21 <- c(part21, part21.i)

      for (j in 1:padd)
      {
        s1.j.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,addcovar[j]]){
          mean(x[t>=start & t<stop])}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,mulcovar[i]]*pdata[,addcovar[j]]){
          mean(x[t>=start & t<stop])}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part22[i,j] <- intesum.ij
      }
    }

    part1 <- part11 - part12 %*% gamma
    part2 <- part21 - part22 %*% gamma
    res <- t(part1) %*% part1 + t(part2) %*% part2
  }

  if (est.method=="interp")
  {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    #estimation
    b.alls0 <- b.alls0m <- b.all
    b.alls0$s0 <- b.alls0m$s0 <- 0
    event.subsx0m <- matrix(unlist(d.event[,c("id",time,addcovar[1],"cent")]),ncol=4)
    event.subsx0m[,3] <- 0
    for (i in 1:pmul){
      b.alls0$s0 <- b.alls0$s0 + beta[i]*b.all[,mulcovar[i]]
      event.subsx0m[,3] <- event.subsx0m[,3] + beta[i]*d.event[,mulcovar[i]]
    }
    b.alls0m$s0 <- -b.alls0$s0
    #b.alls0$s0 <- exp(b.alls0$s0)
    event.subsx0m[,3] <- exp(-event.subsx0m[,3])

    s0.t <- function(t){
      b.all.sub <- b.alls0[which(b.alls0$cent >= t),]
      b.allsub <- matrix(unlist(b.all.sub[,c("id",time,"s0","cent")]),ncol=4)
      inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
      out <- mean(exp(inter$yinte))
      return(out)}
    s0.t <- Vectorize(s0.t)

    part11 <- c()
    part12 <- matrix(NA, nrow=padd, ncol=padd)
    for (i in 1:padd)
    {
      s1.i.t <- function(t){
        b.all.sub <- b.all[which(b.all$cent >= t),]
        b.allsub <- matrix(unlist(b.all.sub[,c("id",time,addcovar[i],"cent")]),ncol=4)
        inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
        out <- mean(inter$yinte)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)

      part11.i <- sum(d.event[,addcovar[i]]*event.subsx0m[,3]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part11 <- c(part11, part11.i)

      for (j in i:padd)
      {
        s1.j.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub <- matrix(unlist(b.all.sub[,c("id",time,addcovar[j],"cent")]),ncol=4)
          inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
          out <- mean(inter$yinte)
          return(out)}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub1 <- matrix(unlist(b.all.sub[,c("id",time,addcovar[i],"cent")]),ncol=4)
          inter1 <- .Fortran("inter",ball=as.single(b.allsub1),idall=as.single(unique(b.allsub1[,1])),dim=as.integer(nrow(b.allsub1)), length=as.integer(length(unique(b.allsub1[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub1[,1])))))
          b.allsub2 <- matrix(unlist(b.all.sub[,c("id",time,addcovar[j],"cent")]),ncol=4)
          inter2 <- .Fortran("inter",ball=as.single(b.allsub2),idall=as.single(unique(b.allsub2[,1])),dim=as.integer(nrow(b.allsub2)), length=as.integer(length(unique(b.allsub2[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub2[,1])))))
          b.all.sub3 <- b.alls0m[which(b.alls0m$cent >= t),]
          b.allsub3 <- matrix(unlist(b.all.sub3[,c("id",time,"s0","cent")]),ncol=4)
          inter3 <- .Fortran("inter",ball=as.single(b.allsub3),idall=as.single(unique(b.allsub3[,1])),dim=as.integer(nrow(b.allsub3)), length=as.integer(length(unique(b.allsub3[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub3[,1])))))
          out <- mean(inter1$yinte*inter2$yinte*exp(inter3$yinte), na.rm=TRUE)
          return(out)}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part12[i,j] <- part12[j,i] <- intesum.ij
      }
    }

    part21 <- c()
    part22 <- matrix(NA, nrow=pmul, ncol=padd)
    for (i in 1:pmul)
    {
      s1.i.t <- function(t){
        b.all.sub1 <- b.all[which(b.all$cent >= t),]
        b.allsub1 <- matrix(unlist(b.all.sub1[,c("id",time,mulcovar[i],"cent")]),ncol=4)
        inter1 <- .Fortran("inter",ball=as.single(b.allsub1),idall=as.single(unique(b.allsub1[,1])),dim=as.integer(nrow(b.allsub1)), length=as.integer(length(unique(b.allsub1[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub1[,1])))))

        b.all.sub2 <- b.alls0[which(b.alls0$cent >= t),]
        b.allsub2 <- matrix(unlist(b.all.sub2[,c("id",time,"s0","cent")]),ncol=4)
        inter2 <- .Fortran("inter",ball=as.single(b.allsub2),idall=as.single(unique(b.allsub2[,1])),dim=as.integer(nrow(b.allsub2)), length=as.integer(length(unique(b.allsub2[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub2[,1])))))
        out <- mean(inter1$yinte*exp(inter2$yinte), na.rm=TRUE)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)

      part21.i <- sum(d.event[,mulcovar[i]]) - sum(s1.i.t(udt)[which(s0.t(udt)!=0)]/s0.t(udt)[which(s0.t(udt)!=0)]* m1[which(s0.t(udt)!=0)])
      part21 <- c(part21, part21.i)

      for (j in 1:padd)
      {
        s1.j.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub <- matrix(unlist(b.all.sub[,c("id",time,addcovar[j],"cent")]),ncol=4)
          inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
          out <- mean(inter$yinte)
          return(out)}
        s1.j.t <- Vectorize(s1.j.t)

        s2.ij.t <- function(t){
          b.all.sub <- b.all[which(b.all$cent >= t),]
          b.allsub1 <- matrix(unlist(b.all.sub[,c("id",time,mulcovar[i],"cent")]),ncol=4)
          inter1 <- .Fortran("inter",ball=as.single(b.allsub1),idall=as.single(unique(b.allsub1[,1])),dim=as.integer(nrow(b.allsub1)), length=as.integer(length(unique(b.allsub1[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub1[,1])))))
          b.allsub2 <- matrix(unlist(b.all.sub[,c("id",time,addcovar[j],"cent")]),ncol=4)
          inter2 <- .Fortran("inter",ball=as.single(b.allsub2),idall=as.single(unique(b.allsub2[,1])),dim=as.integer(nrow(b.allsub2)), length=as.integer(length(unique(b.allsub2[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub2[,1])))))
          out <- mean(inter1$yinte*inter2$yinte, na.rm=TRUE)
          return(out)}
        s2.ij.t <- Vectorize(s2.ij.t)
        inte.f.ij <- function(t)    {
          if (s0.t(t)!=0){
            u <- (s2.ij.t(t) - s1.i.t(t)*s1.j.t(t)/s0.t(t))
          } else {
            u <- 0
          }
          return(u)  }
        inte.f.ij <- Vectorize(inte.f.ij)
        intesum.ij <- (length(inteup))*integrate(f=inte.f.ij, 0, inteup[1], subdivisions =30000L,stop.on.error = FALSE)$value
        for (k in 1:(length(inteup)-1))  {
          intesum.ij <- intesum.ij + (length(inteup)-k)*integrate(f=inte.f.ij, inteup[k], inteup[k+1], subdivisions =30000L, stop.on.error = FALSE)$value
        }
        part22[i,j] <- intesum.ij
      }
    }

    part1 <- part11 - part12 %*% gamma
    part2 <- part21 - part22 %*% gamma
    res <- t(part1) %*% part1 + t(part2) %*% part2
  }
  return(res)
}

est_am <- function(formula, d.event, d.regular, method, bandwidth=NULL, low=NULL, up=NULL, tau=NULL)
{
  time <- all.vars(formula[[2]])
  mulcovar <- all.names(formula[[3]][2])
  addcovar <- all.names(formula[[3]][3])
  mulcovar <- mulcovar[((length(mulcovar)-1)/2+1):length(mulcovar)]
  addcovar <- addcovar[((length(addcovar)-1)/2+1):length(addcovar)]
  p <- length(mulcovar) + length(addcovar)

  idall <- unique(c(d.event$id, d.regular$id))
  N <- length(idall)

  if (method=="kernel")
  {
    regular.sub <- d.regular[which(d.regular[,time]!=0),]
    regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,addcovar[1],"cent")]),ncol=4)
    if (is.null(bandwidth)){
      if (is.null(up)) up <- tau
      bandwidth <- select.con_am(regularz1=regular.subz1, low=low, up=up,tau=tau,N=N)*(N^(-1/3))
    }
  }
  re=optim(par=rep(0,p),fn=score_am, formula=formula, d.regular=d.regular, d.event=d.event, est.method=method, bandwidth=bandwidth, low=low, up=up, tau=tau)

  est=re$par
  conver=re$convergence
  if (conver!=0) warning('the estimation does not converge')
  return(invisible(list(est=est,model="add-mul", method=method, bandwidth=bandwidth)))
}


select.con_am <- function(regularz1, low, up, tau, N)
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

boot_am <- function(formula, d.event, d.regular, nb=50, method, bandwidth, tau=NULL)
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
    invisible(capture.output(est <- est_am(formula=formula, d.event=b.event, d.regular=b.regular, method=method, bandwidth=bandwidth,tau=tau)))
    bootest <- rbind(bootest, est[[1]])
    i <- i+1
  }
  bootsd <- apply(bootest, 2, sd)
  return(invisible(bootsd))
}

base_am <- function(t, theta, formula, d.event, d.regular, method, bandwidth, tau=NULL)
{
  time <- all.vars(formula[[2]])
  mulcovar <- all.names(formula[[3]][2])
  addcovar <- all.names(formula[[3]][3])
  mulcovar <- mulcovar[((length(mulcovar)-1)/2+1):length(mulcovar)]
  addcovar <- addcovar[((length(addcovar)-1)/2+1):length(addcovar)]
  pmul <- length(mulcovar)
  padd <- length(addcovar)
  beta <- theta[1:pmul]
  gamma <- theta[(pmul+1):(pmul+padd)]

  if (t < min(d.event[,time])) {
    baseline <- 0 } else {
  idall <- unique(c(d.regular$id, d.event$id))
  N <- length(idall)

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
   regular.subz1 <- matrix(unlist(regular.sub[,c("id",time,addcovar[1],"cent")]),ncol=4)

   regular.subsx0 <- regular.subz1
   regular.subsx0[,3] <- 0
   for (i in 1:pmul)
   {
     regular.subsx0[,3] <- regular.subsx0[,3] + beta[i]*regular.sub[,mulcovar[i]]
   }
   regular.subsx0[,3] <- exp(regular.subsx0[,3])

   s0.t <- function(t){
     if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
     }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
     } else {ker <- .Fortran("kernel",bcarr=as.single(regular.subsx0),dim=as.integer(nrow(regular.subsx0)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
     return(ker$ker[1])
   }
   s0.t <- Vectorize(s0.t)

   part1 <- c()
   for (i in 1:padd)
   {
     subzi <- regular.subz1
     subzi[,3] <- regular.sub[,addcovar[i]]
     s1.i.t <- function(t){
       if (t < bandwidth) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(bandwidth),ker=as.single(0))
       }  else if (t > (tau-bandwidth)) {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(tau-bandwidth),ker=as.single(0))
       } else {ker <- .Fortran("kernel",bcarr=as.single(subzi),dim=as.integer(nrow(subzi)),bandwidth=as.single(bandwidth),tt=as.single(t),ker=as.single(0))}
       return(ker$ker[1])
     }
     s1.i.t <- Vectorize(s1.i.t)
     inte.f.i <- function(t)    {
       if (s0.t(t)!=0){
         u <- s1.i.t(t)/s0.t(t)
       } else {
         u <- 0 }
       return(u)  }
     inte.f.i <- Vectorize(inte.f.i)
     part1 <- c(part1, integrate(f=inte.f.i, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
   }
   udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
   m1sub <- table(d.event[which(d.event[,time]<= t),time])
   baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)]) - part1 %*% gamma
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
        for (i in 1:pmul)
        {
          foo[,mulcovar[i]] <- sapply(foo$start,
                                      function(s, t1=foo1[,time], x1=foo1[,mulcovar[i]], t2=foo2[,time],
                                               x2=foo2[,mulcovar[i]]){
                                        if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
        }
        for (i in 1:padd)
        {
          foo[,addcovar[i]] <- sapply(foo$start,
                                      function(s, t1=foo1[,time], x1=foo1[,addcovar[i]], t2=foo2[,time],
                                               x2=foo2[,addcovar[i]]){
                                        if( any(s==t1)){ x1[s==t1]} else { x2[s==t2]} } )
        }
        pdata <- rbind(pdata, foo)
      }
    }

    pdatas0 <- pdata
    pdatas0$s0 <- 0
    for (i in 1:pmul){
      pdatas0$s0 <- pdatas0$s0 + beta[i]*pdata[,mulcovar[i]]
    }
    pdatas0$s0 <- exp(pdatas0$s0)

    s0.t <- function(t, start=pdatas0$start, stop=pdatas0$stop, x=pdatas0$s0){
      mean(x[t>=start & t<stop])}
    s0.t <- Vectorize(s0.t)

    part1 <- c()
    for (i in 1:padd)
    {
      s1.i.t <- function(t, start=pdata$start, stop=pdata$stop, x=pdata[,addcovar[i]]){
        mean(x[t>=start & t<stop])}
      s1.i.t <- Vectorize(s1.i.t)

      s1.i.t <- Vectorize(s1.i.t)
      inte.f.i <- function(t)    {
        if (s0.t(t)!=0){
          u <- s1.i.t(t)/s0.t(t)
        } else {
          u <- 0 }
        return(u)  }
      inte.f.i <- Vectorize(inte.f.i)
      part1 <- c(part1, integrate(f=inte.f.i, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
    }
    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)]) - part1 %*% gamma
  }

  if (method=="interp")
  {
    b.all <- rbind(d.regular, d.event)
    b.all <- b.all[order(b.all[,time]),]
    b.all <- b.all[order(b.all$id),]

    #estimation
    b.alls0 <- b.all
    b.alls0$s0 <-  0
    for (i in 1:pmul){
      b.alls0$s0 <- b.alls0$s0 + beta[i]*b.all[,mulcovar[i]]
    }

    s0.t <- function(t){
      b.all.sub <- b.alls0[which(b.alls0$cent >= t),]
      b.allsub <- matrix(unlist(b.all.sub[,c("id",time,"s0","cent")]),ncol=4)
      inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
      out <- mean(exp(inter$yinte))
      return(out)}
    s0.t <- Vectorize(s0.t)

    part1 <- c()
    for (i in 1:padd)
    {
      s1.i.t <- function(t){
        b.all.sub <- b.all[which(b.all$cent >= t),]
        b.allsub <- matrix(unlist(b.all.sub[,c("id",time,addcovar[i],"cent")]),ncol=4)
        inter <- .Fortran("inter",ball=as.single(b.allsub),idall=as.single(unique(b.allsub[,1])),dim=as.integer(nrow(b.allsub)), length=as.integer(length(unique(b.allsub[,1]))),t=as.single(t),tau=as.single(tau),yinte=as.single(rep(0,length(unique(b.allsub[,1])))))
        out <- mean(inter$yinte)
        return(out)}
      s1.i.t <- Vectorize(s1.i.t)
      inte.f.i <- function(t)    {
        if (s0.t(t)!=0){
          u <- s1.i.t(t)/s0.t(t)
        } else {
          u <- 0 }
        return(u)  }
      inte.f.i <- Vectorize(inte.f.i)
      part1 <- c(part1, integrate(f=inte.f.i, 0, t, subdivisions =30000L,stop.on.error = FALSE)$value)
    }
    udtsub <- sort(unique(d.event[which(d.event[,time]<= t),time]))
    m1sub <- table(d.event[which(d.event[,time]<= t),time])
    baseline <- sum(1/s0.t(udtsub)[which(s0.t(udtsub)!=0)]/gt(udtsub)[which(s0.t(udtsub)!=0)]*m1sub[which(s0.t(udtsub)!=0)]) - part1 %*% gamma
  }
    }
  return(invisible(baseline))
}
