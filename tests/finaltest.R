require(aa)
require(aace)
do.call(what=ceload,args=getsi("ce"))
do.call(what=ceload,c(list(da=c("2014-10-29")),getsi("ce")))
do.call(what=ceload,c(list(da=c("2014-10-22")),getsi("ce")))

da <- c("2014-10-22","2014-10-29")
putrd(rbindlist(lapply(lapply(da,ceload),data.table)),"ce test")

da <- getca()[as.Date(getca())%in%da.global]
dax <- rev(rev(da)[1:20])
ce <- rbindlist(lapply(lapply(dax,ceload),data.table))

require(snowfall)
sfInit(par=F,cpus=1)


tail(ceload(dax[1]))


if(F) {
  require(snowfall)
  sfInit(par=FALSE,cpus=1)
  dada <-c("2014-10-22","2014-10-29")
  xx <- rbindlist(lapply(sfLapplyWrap(dada,ceload),data.table))
  sfStop()
  xx[,unique(date)]
  
  sfInit(par=FALSE,cpus=2)
  dax <- rev(rev(da)[1:20])
  xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table))
  sfStop()
  xx[,unique(date)]

  sfInit(par=FALSE,cpus=5)
  dax <- rev(rev(da)[1:20])
  xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table))
  sfStop()
  xx[,unique(date)]
}
