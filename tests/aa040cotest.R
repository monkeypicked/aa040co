require(aautil)
require(aaco)
#require(aa)
require(quadprog)
aatopselect("test")


# do.call(what=ceload,args=getsi("ce"))
# do.call(what=ceload,c(list(da=c("2014-10-29")),getsi("ce")))
# do.call(what=ceload,c(list(da=c("2014-10-22")),getsi("ce")))

da <- c("2014-10-22","2014-10-29")
bui <- buiindirs()
su.g <<- setkey(setnames(data.table(expand.grid(da,bui)),c('date','bui')),date,bui)[]

ce <- rbindlist(lapply(lapply(da,ceload,nfac=1),data.table))
expect_equal(ce[,sort(unique(date))],sort(unique(da)))

#slower bit
# if(FALSE) {
#   sfInit(par=FALSE)
#   system.time(xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table)))
#   sfStop()
#   xx[,unique(date)]
# 
#   sfInit(par=TRUE,cpus=5)
#   system.time(xx<-rbindlist(lapply(sfLapplyWrap(dax,ceload),data.table)))
#   #sfStop()
#   putrd(xx,"ce 20 final dates")
# }
