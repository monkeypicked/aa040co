require(aabd)

load(paste0(bbdir(),"/bdp/derive-000/CRNCY_ADJ_MKT_CAP.RData"))
adj <- x
load(paste0(bbdir(),"/bdp/derive-000/CUR_MKT_CAP.RData"))
raw <- x
adj[raw][bigbui]


#high market cap bui
dir(paste0(bbdir(),"/bdp/derive-000"))
load(paste0(bbdir(),"/bdp/derive-000/NAME.RData"))
namex<-x
load(paste0(bbdir(),"/bdp/derive-000/CRNCY_ADJ_MKT_CAP.RData"))
x <- x[!is.na(CRNCY_ADJ_MKT_CAP)]
bigbui <- setkey(x,CRNCY_ADJ_MKT_CAP)[nrow(x)-(seq(1,20)-1),bui]
load(paste0(bbdir(),"/bdp/derive-000/NAME.RData"))   #just to check, found a bug, now fixed
namedtab <- setkey(x[setkey(adj[raw][bigbui],'bui')],CRNCY_ADJ_MKT_CAP)[]

require(aautil)

resve("ce",896)
extract("ce","distinct date")
extract("ce","distinct bui")
ce<-data.table(gettab("ce",ret="data"))
setkey(ce,date)
#putrd(ce,desc="ice=296") #wrong, it is all dates
ce['2014-10-29',mean(as.numeric(full),na.rm=T)] #.847


require(aa)
dergu()
#############
#take a sample of stocks of different length of I

resve("ce",896)
ce<-data.table(gettab("ce",ret="data"))
bui0 <- ce[,unique(bui)]
all(bui0%in%buiindir(paste0(bbdir(),"/",bdp2dir()))) #yes

load(paste0(bbdir(),"/bdp/derive-000/NAME.RData"))
namex<-copy(x)

load(paste0(bbdir(),"/bdp/derive-000/BICS_REVENUE_%_LEVEL_ASSIGNED.RData"))
nn <- names(x)
nn[4]<-"BICS_REVENUE_PERC_LEVEL_ASSIGNED"
setnames(x,nn)
xcopy <- copy(x)

bui1 <- intersect(bui0,x[,bui]) #some bui have files but no data for segments - this needs handling
bs <-setkey(x[bui1][,seqmax:=max(seq),by=bui][,max(bui),by=seqmax],seqmax)[]
x[bs[10,V1]] #this has 12 segments and is presumably adecoagro
x[bs[2,V1]]
#for a given bui
m<-5
m<-4
bui2 <- bs[m,V1]
#select a date when it is in the universe
setkey(ce,bui)
da <- ce[bui2][,max(date)]
#check multiplicity of each of its segments in the univers
x2 <- setkey(x,bui)[bs[m,V1]]
setkey(x, BICS_LEVEL_CODE_ASSIGNED)
setkey(x2, BICS_LEVEL_CODE_ASSIGNED)
allseg <- x[x2[,BICS_LEVEL_CODE_ASSIGNED]][,.N,list(bui,BICS_LEVEL_CODE_ASSIGNED,BICS_LEVEL_NAME_ASSIGNED,BICS_REVENUE_PERC_LEVEL_ASSIGNED)]
allseg
#get the vcv
ace <- getce(da)
vcv <- vcvce(ace)
#constraint to match primary exposure
mybc <- setkey(x2,seq)[1,BICS_LEVEL_CODE_ASSIGNED]
x[mybc]
buitosolve <- x[mybc][,bui]
buiinboth <-  intersect(buitosolve,rownames(vcv$T))
A <- vcv$T[,1,drop=F]*0
LHS <- A[buiinboth,,drop=F]
RHS <- setkey(allseg[BICS_LEVEL_CODE_ASSIGNED==mybc],bui)[buiinboth]
rownames(LHS)==RHS[,bui]
A[buiinboth,] <- RHS[,BICS_REVENUE_PERC_LEVEL_ASSIGNED]
#adjust to remove target from solution to avoid trivial solution x=target
A[bui2,]
i<-match(bui2,rownames(A))
require(quadprog)
?solve.QP
Dmat<-vcv$T[-i,-i]
dvec<-rep(1,ncol(Dmat))
Amat<-A[-i,,drop=F]
bvec<-A[i,1]
sol <- solve.QP(Dmat,dvec,Amat,bvec,1)
sol$solution%*%Amat
#check solution
sol$solution%*%Amat-bvec
#generalise to multiple
######a better way to create A, buiinboth . x2[,BICS_LEVEL_CODE_ASSIGNED] ; the latter being bvec
#1 convert the irregular table of exposures to a matrix bui . BLCA; select the desired cols; remove the target's row
require(reshape2)
?cast
#get a data.table with cols bui,BLCA and then cast it into A
#raw data in x
#bui in vcv$T and in xs DONE : buix
#blca in x2[,BICS_LEVEL_CODE_ASSIGNED]
xs <- copy(xcopy) #the original
setkey(xs,bui)
buix <- intersect(rownames(vcv$T),xs[,bui]) #solution set
bui2%in%buix
xs1 <- xs[buix]
meltedx <- setnames(xs1[,list(bui,BICS_LEVEL_CODE_ASSIGNED,BICS_REVENUE_PERC_LEVEL_ASSIGNED)],c('i','j','x'))
aa<-acast(data=data.frame(meltedx),formula=i ~ j)
aa[is.na(aa)]<-0
aa[,x2[,BICS_LEVEL_CODE_ASSIGNED]]
x2[]

#check sequence of each argument 
buisol <- setdiff((intersect(rownames(vcv$T),rownames(aa))),bui2)
all(rownames(aa)%in%rownames(vcv$T))
Dmat<-vcv$T[buisol,buisol]
dvec<-matrix(rep(1,ncol(Dmat)),ncol(Dmat),1,dimnames=list(buisol,NULL))
Amat<-aa[buisol,x2[,BICS_LEVEL_CODE_ASSIGNED],drop=F]
bvec<-x2[,BICS_REVENUE_PERC_LEVEL_ASSIGNED]
sol1 <- solve.QP(Dmat,dvec,Amat,bvec,length(bvec))
sol1$solution%*%Amat
sol1$solution%*%Amat-bvec
barplot(sort(sol1$solution))

sol2<- solve.QP(Dmat,dvec/.4e5,Amat,bvec,length(bvec))
sol2$solution%*%Amat
sol2$solution%*%Amat-bvec
barplot(sort(sol2$solution))
sum(abs(sol2$solution))
sum((sol2$solution))

setkey(x,bui)[buisol[order(-sol2$solution)][1:5]]

sol2$solution[order(-sol2$solution)][1:5]

i <- order(-sol2$solution)[2]
i <-which.max(sol2$solution)
dvec<-matrix(rep(0,ncol(Dmat)),ncol(Dmat),1,dimnames=list(buisol,NULL))
dvec[i]<-140
sol3<- solve.QP(Dmat,dvec/.4e5,Amat,bvec,length(bvec))
sol3$solution%*%Amat
sol3$solution%*%Amat-bvec
barplot(sort(sol3$solution))
sum(abs(sol3$solution))
sum((sol3$solution))


namex[bui2]
namex[buisol[i]]



#unfortunately tgt... fails with an unexpected condition & dk why
# constr <- list(Am=Amat,bv=bvec,meq=length(bvec))
# tgt.solve.QP(ce=cex,dvec=dvec/.5e6,constr=constr,tgt=1)
# cex <- getce(da)
# class(ce)
