require(aautil)
#bui <- getrd(8)[BICS_LEVEL_CODE_ASSIGNED==17101013,bui]
aatopselect("test")
require(aaco)
require(aapa)
#require(aa)
require(quadprog)
require(testthat)


# do.call(what=ceload,args=getsi("ce"))
# do.call(what=ceload,c(list(da=c("2014-10-29")),getsi("ce")))
# do.call(what=ceload,c(list(da=c("2014-10-22")),getsi("ce")))

#da <- c("2014-10-22","2014-10-29")
#bui <- buiindirs()
#su.g <<- setkey(setnames(data.table(expand.grid(da,bui,stringsAsFactors=FALSE)),c('date','bui')),date,bui)[]
#su.g <<- setkey(setnames(data.table(expand.grid(da,colnames(getbdh(da=da[1])),stringsAsFactors=FALSE)),c('date','bui')),date,bui)[]   #somewhat circular
# .global()
# damax <- max(index(getstep()))
# su.g <<- cart(data.table(bui=nonadt()[date==max(daw.global),unique(bui)]),data.table(date=nonadt()[,unique(date)]))  #nonasu() #this is BROKEN
# bui<-su.g[,unique(bui)]
# da<-su.g[,as.Date(unique(date))]

getbdhgl()
getbdmgl()
su <- getrdatv("jo","su")
win <- -100:0
x <- getbdh(su=su,win=win)
#bui
da <-su[,max(date)]
buisu <- su[date==da,sort(unique(bui))]
expect_equal(sort(unique(colnames(x))),su[date==da,sort(unique(bui))])
#dates
expect_equal(offda(da,win),index(x))
#carry forward : test by taking the date prior to the first su date, da
da1 <- offda(da,-1)
x1 <- getbdh(su=su,da=da1,win=win)
da2 <- su[,sort(unique(date),decreasing=TRUE)[2]]
buisu <- su[date==da,sort(unique(bui))]
expect_equal(sort(unique(colnames(x1))),su[date==da2,sort(unique(bui))])
#out of range : (1) da is 2 after
x2 <- getbdh(su=su,win=win,da=offda(da,2))
expect_equal(sort(rownames(x2),decreasing=TRUE)[3],as.character(da))
expect_equal(sort(unique(colnames(x))),sort(unique(colnames(x2))))
#first date
x3 <- getbdh(su=su,win=win,da=su[,min(date)])
expect_equal(sort(unique(colnames(x))),sort(unique(colnames(x2))))
sux <- su[date==su[,min(date)]]
expect_equal(dim(x3),c(length(win),nrow(sux)))
#prior to first date
expect_warning(expect_true(is.null(getbdh(su=su,win=win,da=offda(su[,min(date)],-1)))))



expect_false(is.factor(su.g[date==da[1]][,bui]))
#su.g[da[1]][,bui]
#buiindirs()
#expect_equal(class(getbdh(da=da[1],bui=colnames(getbdh(da=da[1])),field="prem.g",win=-230:-1)),"zoo") #prem.g does not exist yet


getbdhgl() #added this 01-11, does what is wanted ie load prem.g

ce <- rbindlist(lapply(lapply(max(da),ceload,nfac=5),data.table))
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

