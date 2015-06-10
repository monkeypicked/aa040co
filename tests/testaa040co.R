require(aautil)
require(aaco)
require(aapa)
require(aabd)
require(quadprog)
require(testthat)

aatopselect("test")

#aabd::fixperc() already done

#getbdhgl-----
getbdhgl()
expect_true(exists("prem.g"))
expect_true(exists("vix.g"))

#getpa-----
su <-getrdatv("jo","su",0)
su <- rbind(su,copy(su)[,date:=date-1000])
win <- -100:0
x <- getpa(su=su,win=win)
#bui
da <-as.Date(su[,max(date)])
buisu <- su[date==da,sort(unique(bui))]
expect_equal(sort(unique(colnames(x))),su[date==da,sort(unique(bui))])
#dates
expect_equal(offda(da,win),index(x))
#carry forward : test by taking the date prior to the first su date, da
da1 <- offda(da,-1)
x1 <- getpa(su=su,da=da1,win=win)
da2 <- su[,sort(unique(date),decreasing=TRUE)[2]]
buisu <- su[date==da,sort(unique(bui))]
expect_equal(sort(unique(colnames(x1))),su[date==da2,sort(unique(bui))])
#out of range : (1) da is 2 after
x2 <- getpa(su=su,win=win,da=offda(da,2))
expect_equal(sort(rownames(x2),decreasing=TRUE)[3],as.character(da))
expect_equal(sort(unique(colnames(x))),sort(unique(colnames(x2))))
#first date
# x3 <- getpa(su=su,win=win,da=su[,min(date)])
# expect_equal(sort(unique(colnames(x))),sort(unique(colnames(x2))))
# sux <- su[date==su[,min(date)]]
# expect_equal(dim(x3),c(length(win),nrow(sux)))
#prior to first date
expect_warning(expect_true(is.null(getpa(su=su,win=win,da=offda(su[,min(date)],-1)))))


#cewrap-----
fms.df <- cewrap(getpa(su)) #returns output from cewrap, ie dfrce data.frame
putrdatv(data.table(fms.df),'jo','co',0)

#dfrce-----[internal to cewrap]
#addbench-----[internal to cewrap]
#ceload-----[not useful if exists(globals)]

#conversion back and forth
co <- getrdatv("jo","co")
co <- co[date==co[,max(date)],]
co1 <- data.table(dfrce(dtce(co))) #round trip
expect_equal(vcvce(dtce(co))$T,vcvce(dtce(co))$T)

#getrdatv("jo","su",v=2)
#pruneztei
expect_true(pruneztei(su,nmin=3)[,.N,bcode][,all(N>=3)])
expect_true(pruneztef(su,nmin=3,wmin=3)[,list(.N,wgt=sum(BRPLA)),bcode][bcode!='00',all(N>=3)&all(wgt>=3)])
expect_true(pruneztef(su,nmin=3,wmin=3)[,length(unique(bcode))]<pruneztef(su,nmin=3,wmin=2)[,length(unique(bcode))])
expect_true(pruneztef(su,nmin=3,wmin=1)[,length(unique(bcode))]<pruneztef(su,nmin=2,wmin=1)[,length(unique(bcode))])
wmat <- tabtomat(data.frame(pruneztef(su,nmin=5,wmin=5)))
wmat[is.na(wmat)]<-0
wmat <- wmat[,order(colSums(wmat))]
wmat <- wmat[order(apply(sweep(wmat,MAR=2,STAT=1:ncol(wmat),FUN="*"),1,sum)),]
expect_equal(nrow(wmat),su[,length(unique(bui))])
myrowsums <- apply(wmat,1,sum,na.rm=TRUE)
expect_true(all(myrowsums<1.1)&all(myrowsums>0.99))
                                   
if(FALSE) {
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
pa <-getpa(su=su)
phis <- 0.1
phir <- 0.1
!any(is.na(pa))


pax <- getpa(su=su,fi='best.g')
nna <- 5000
coredata(pax)[sample(1:length(pa),nna,rep=FALSE)] <- NA
sum(is.na(pax))
image(coredata(pax))


x1 <- arco(pa=pax,co=co,phir=0,phis=0,typ='p')
#x2 <- arco(pa=pax,co=co,phir=0,phis=0,typ='i')
sum(is.na(x1$fit))==0
sum(is.na(x1$act))==nna
mser <- -5:10
mses <- -5:10
for(i in seq_along(mser)) mser[i] <- arco(pa=pax,co=co,phir=mser[i]/10,phis=.6,typ='p')$MSE
for(i in seq_along(mses)) mses[i] <- arco(pa=pax,co=co,phir=1.,phis=mses[i]/10,typ='p')$MSE

par(mfrow=c(1,2))
barplot(mser)
barplot(mses)
min(mses)
min(mser)
}