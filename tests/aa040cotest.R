require(aautil)
aatopselect("test")
require(aaco)
#require(aa)
require(quadprog)


# do.call(what=ceload,args=getsi("ce"))
# do.call(what=ceload,c(list(da=c("2014-10-29")),getsi("ce")))
# do.call(what=ceload,c(list(da=c("2014-10-22")),getsi("ce")))

da <- c("2014-10-22","2014-10-29")
bui <- buiindirs()
su.g <<- setkey(setnames(data.table(expand.grid(da,bui,stringsAsFactors=FALSE)),c('date','bui')),date,bui)[]
#su.g <<- setkey(setnames(data.table(expand.grid(da,colnames(getbdh(da=da[1])),stringsAsFactors=FALSE)),c('date','bui')),date,bui)[]   #somewhat circular

expect_false(is.factor(su.g[da[1]][,bui]))
su.g[da[1]][,bui]
buiindirs()
#expect_equal(class(getbdh(da=da[1],bui=colnames(getbdh(da=da[1])),field="prem.g",win=-230:-1)),"zoo") #prem.g does not exist yet


getbdhgl() #added this 01-11, does what is wanted ie load prem.g
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
