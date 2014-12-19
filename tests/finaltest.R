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

tail(ceload(dax[1]))


if(F) {
  require(snowfall)
  sfInit(par=FALSE,cpus=1)
  sfLapplyWrap()
  dada <-c("2014-10-22","2014-10-29")
  xx<-lapply(dada,ceWrap)
  xx<-lapply(dada,loadWrapper)
  xx<-sfLapplyWrap(dada,loadWrapper)
  rm(prem.g)
  rm(vix.g)
}
