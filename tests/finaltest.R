require(aaco)
require(aa)
do.call(what=ceload,args=getsi("ce"))
do.call(what=ceload,c(list(da=c("2014-10-29")),getsi("ce")))
do.call(what=ceload,c(list(da=c("2014-10-22")),getsi("ce")))

da <- c("2014-10-22","2014-10-29")
ce <- rbindlist(lapply(lapply(da,ceload),data.table))
#putrd(rbindlist(lapply(lapply(da,ceload),data.table)),"ce test")


da.g <- getca()[as.Date(getca())%in%da.global]
dax <- rev(rev(da.g)[1:20])
ce <- rbindlist(lapply(lapply(dax,ceload),data.table))


#slower bit
if(F) {
  sfInit(par=FALSE)
  system.time(xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table)))
  sfStop()
  xx[,unique(date)]

  sfInit(par=TRUE,cpus=5)
  system.time(xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table)))
  #sfStop()
  putrd(xx,"ce 20 final dates")
}
