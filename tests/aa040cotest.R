require(aautil)
aatopselect("test")
require(aaco)
require(aapa)
require(quadprog)
require(testthat)


#getbdhgl-----
getbdhgl()
expect_true(exists("prem.g"))
expect_true(exists("vix.g"))

#getbdh-----
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

#cewrap-----
fms.df <- cewrap(getbdh(su))

#dfrce-----[internal to cewrap]
#addbench-----[internal to cewrap]
#ceload-----[not useful if exists(globals)]

#conversion back and forth
co <- getrdatv("jo","co")
co <- co[date==co[,max(date)],]
co1 <- data.table(dfrce(dtce(co))) #round trip
expect_equal(vcvce(dtce(co))$T,vcvce(dtce(co))$T)

#pruneztei
expect_true(pruneztei(nmin=3)[,.N,bcode][,all(N>=3)])
expect_true(pruneztef(nmin=3,wmin=3)[,list(.N,wgt=sum(BRPLA)),bcode][bcode!='00',all(N>=3)&all(wgt>=3)])
expect_true(pruneztef(nmin=3,wmin=3)[,length(unique(bcode))]<pruneztef(nmin=3,wmin=2)[,length(unique(bcode))])
expect_true(pruneztef(nmin=3,wmin=1)[,length(unique(bcode))]<pruneztef(nmin=2,wmin=1)[,length(unique(bcode))])
wmat <- tabtomat(data.frame(pruneztef(nmin=5,wmin=5)))
wmat[is.na(wmat)]<-0
wmat <- wmat[,order(colSums(wmat))]
wmat <- wmat[order(apply(sweep(wmat,MAR=2,STAT=1:ncol(wmat),FUN="*"),1,sum)),]
expect_equal(nrow(wmat),su[,length(unique(bui))])
myrowsums <- apply(wmat,1,sum,na.rm=TRUE)
expect_true(all(myrowsums<1.01)&all(myrowsums>0.99))
image(wmat)
barplot(sort(apply(wmat,2,sum,na.rm=T)),horiz=T)
barplot(sort(apply(wmat,1,sum,na.rm=T)),horiz=T)
barplot((wmat),col=1:100,bord=F,horiz=T,space=0,bes=F)
barplot(t(wmat)[,],col=1:100,bord=F,horiz=T,space=0,bes=F)
image((wmat[,order(apply(wmat,2,sum,na.rm=T))]))
wmat <- wmat
#pruneztef

