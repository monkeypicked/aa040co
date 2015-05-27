require(aautil)
require(aaco)
require(aapa)
require(aabd)
require(quadprog)
require(testthat)

aatopselect("test")

aabd::fixperc()

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

su <- getrdatv("jo","su",v=2)
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

#arco - autoregressive components
require(aaco)
aatopselect("t")
getbdhgl()
su <- getrdatv("jo","su",v=2)
pa <-getbdh(su=su)
phis <- 0.1
phir <- 0.1
!any(is.na(pa))


pax <- getbdh(su=su,fi='best.g')
nna <- 5000
coredata(pax)[sample(1:length(pa),nna,rep=FALSE)] <- NA
sum(is.na(pax))
image(coredata(pax))

x1 <- arco(pa=pax,phir=0,phis=0,typ='p')
x2 <- arco(pa=pax,phir=0,phis=0,typ='i')
sum(is.na(x1$fit))==0
sum(is.na(x1$act))==nna
mser <- -5:10
mses <- -5:10
for(i in seq_along(mser)) mser[i] <- arco(pa=pax,phir=mser[i]/10,phis=.6,typ='i')$MSE
for(i in seq_along(mses)) mses[i] <- arco(pa=pax,phir=1.,phis=mses[i]/10,typ='i')$MSE

par(mfrow=c(1,2))
barplot(mser)
barplot(mses)
min(mses)
min(mser)

x1 <- arco(pa=pa,phir=0,phis=0,typ='p')
x2 <- arco(pa=pa,phir=0,phis=0,typ='i')
x1$MSE==x2$MSE




x1 <- arco(pa=pa,su=su,phir=0,phis=0)






su <- getrdatv("jo","su",v=2)
pa <-getbdh(su=su)
i <- sample(1:length(pa),siz=100,rep=F)
pax <- coredata(pa)[i]
pa[i]<-NA
pa1 <- interpte(pa=pa,type='f')
pax1 <- coredata(pa1)[i]
cor(pax,pax1)
plot(pax,pax1)




arco(pa=pax,phis=1)



###Unfo this test universe is pretty illiquid so cannot really test on it... wait for prod and use liquids

#2015-05-27 test sequence prior to hatij
require(aaco)
aatopselect('t')
getbdhgl()
su0 <- getrd(140) #global ag

builiq <- getbdp()[su0[,bui]][800<CRNCY_ADJ_MKT_CAP,bui]

su <- setkey(su0,bui)[builiq]
pa <- getbdh(su)
ce <- dtce(data.table(cewrap(pa,nfac=2)))
fmp <- fmpce(ce) 
su <- setkey(su,bui)[bui] #selected no-na bui
#putrd(su,'no hols 28 ag') #156
pa <- getbdh(su)
lc <- loocvi(pa[1:230,],nfac=2)
max(abs(lc$mhatM-lc$mhatiM))
max(abs(lc$mhatS-lc$mhatiS))
max(abs(lc$mhatT-lc$mhatiT))
ilcvsumm(lc)
