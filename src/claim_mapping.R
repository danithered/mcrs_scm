dirs = dir("/home/danielred/data/programs/mcrs_to_scm/OUT/", "bubble", full.names = T)

sumdata = data.frame(
  claim=rep(NA, length(dirs)),
  survived=logical(length(dirs)),
  survived_until=numeric(length(dirs)),
  max_replication=numeric(length(dirs)),
  third_alive_at=numeric(length(dirs))
                     )

for (d in dirs){
  setwd(d)
  data = read.table("output.csv", sep=";", header=T)
  #str(data)
  
  claim = strsplit(d, "_")[[1]][5]
  

  sumdata[which(is.na(sumdata$claim))[1], ] = data.frame(
    claim = as.numeric(claim),
    survived = data[nrow(data),"replicators"] > 0,
    survived_until = data$time[max(which(data$replicators > 0))],
    max_replication = max(data$percent_replicated, na.rm=T),
    third_alive_at = data[data$no_alive >= 333, "time"][1]
  )

  data = data[data$replicators>0,]
  
  par(mfrow=c(2,3))
  plot(data$time, data$replicators, type="l", main=claim)
  plot(data$time, data$no_alive, type="l", main=claim)
  plot(data$time, data$no_last_splits, type="l", main=claim)
  plot(data$time, data$mean_M, type="l", main=claim)
  plot(data$time, data$percent_replicated, type="l", main=claim)
  plot(data$time, data$percent_died, type="l", main=claim)
  par(mfrow=c(1,1))
}

sumdata = sumdata[order(sumdata$claim),]

par(mfrow=c(3, 1))
par(mai=c(0,0.82,0,0.42), oma=c(4, 0 ,1, 0))

plot(sumdata$claim, 
     ifelse(sumdata$survived, NA, sumdata$survived_until), 
     xlab=bquote(C[norep]), 
     xaxt="n",
     pch=19, type="b", las=1,
     ylab="died at [generation]",
     col=ifelse(sumdata$survived, "green", "red") )

plot(sumdata$claim, 
     sumdata$max_replication*100, 
     xlab=bquote(C[norep]), 
     xaxt="n",
     ylab="Maximum replication (% of replicators)",
     pch=19, type="b", las=1,
     col=ifelse(sumdata$survived, "green", "red")
     )

plot(sumdata$claim, 
     sumdata$third_alive_at, 
     xlab=bquote(C[norep]), 
     ylab="No vesicules alive above 1/3 of vesicules [generation]",
     las=1,
     col=ifelse(sumdata$survived, "green", "red"),
     type="b", pch=19)

mtext(bquote(C[norep]), 1, 3)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mai=c(1.02,0.82,0.82,0.42), oma=rep(0, 4))

{
par(mfrow=c(3,1))
par(mai=c(0,0.82,0,0.42), oma=c(4, 0 ,1, 0))
  
plot(data$time, data$replicators, 
     xaxt="n",
     ylab="Number of replicators",
     type="l")
abline(v=60, lty=2)

plot(data$time, data$no_alive, 
     xaxt="n",
     ylab="Number of vesicles alive",
     type="l")
abline(v=60, lty=2)

plot(data$time, data$no_last_splits, 
     # xaxt="n",
     ylab="Number of splitting events",
     type="l")
abline(v=60, lty=2)

# plot(data$time, data$mean_M, type="l")
# abline(v=60, lty=2)

mtext("Time [generations]", 1, 2)
}

claim

