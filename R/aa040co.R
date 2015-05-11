
#' top-level derive method for covariance table
#'
#' loads data to global if not exists; loops through date estimating co, convert to data.table, save to rd
#' 
#' Details (body)
#' 
#' @param su securityuniverse object
#' @param da date of class Date
#' @param verbose flag for printing date
#' @param nfac number of factors
#' @param ... passed to cewrap
#' @family top level
#' @examples 
#' \dontrun{
#' deraaco()
#' }
#' @export
deraaco <- function(su=getrdatv("jo","su"),da=su[,sort(unique(date))],verbose=TRUE,nfac=20,...) {
  if(! exists("prem.g",envir=globalenv())) getbdhgl() #gets latest premium, vix, mcap
  x <- vector("list",length(da))
  for(i in seq_along(da)) {
    if(verbose) print(da[i])
    x[[i]] <- data.table(cewrap(pa=getbdh(su=su,da=da[i]),nfac=nfac,...))
  }
  putrdatv(rbindlist(x),"jo","co")
}

#' get global zoo panels
#'
#' loads timeseries required for co estimation
#' 
#' @param field named list of derived panels to load; they are assigned in the global environment, with the name in the list
#' @family internal
#' @examples 
#' \dontrun{
#' getbdhgl()
#' exists('prem.g') #true
#' }
#' @export
getbdhgl <- function(field=list(
                              list(prem.g="x0700redoto"),
                              list(mcap.g="x0702mcap"),
                              list(dipr.g="x0702dipr"),
                              list(erpr.g="x0702erpr"),
                              list(bopr.g="x0702bopr"),
                              list(best.g="x0702best"),
                              list(vix.g="x0502vix")),
                      myclass='zoo') {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep(unlist(field),n='001',myclass=myclass),envir=globalenv())})
}


#' return windowed zoo panel from global panel
#'
#' Description (body)
#' 
#' Details (body)
#' 
#' @param su securityuniverse object
#' @param da date of class Date, forms the datum for the lags in win
#' @param bui identifiers for inclusion
#' @param field name of global zoo panel from which to extract
#' @param win window normally negative (prior) lags referred to da
#' @family utility
#' @examples 
#' \dontrun{
#' getbdh()
#' }
#' @export
getbdh <- function(su=su.g,da=su[,max(date)],bui=su[date==su[date<=da,max(date)],bui],field="prem.g",win=-(1:230)) {
  if(0==length(bui)*length(da)) {
    return(NULL) 
  } else { 
    get(field,envir=globalenv())[as.Date(offda(da,win)),bui]
  }
}

#' derive ce
#'
#' high-level covariance estimation
#' 
#' takes panel of returns, applies normalisations, estimates ce, adds benchmark data, returns as dataframe
#' 
#' @param pa zoo panel (normally returns)
#' @param normalise flag to redistribute - see aautil::zoonorm()
#' @param nfac # factors in PCA
#' @param applyvix flag to scale historical returns to 'today's vol'
#' @param mixvix bayes adjustment to vix rescaling, strength of shrinkage to mean
#' @param bench type of benchmark weights
#' @param ... passed to fms2
#' @examples 
#' \dontrun{
#' cewrap()
#' }
#' @export
cewrap <- function (pa=getbdh(), normalise="NONE", nfac=20, applyvix=TRUE, mixvix=.5, bench="equal",...) {
  da <- max(index(pa))
  if (applyvix) {
     vix <- na.locf(vix.g )
     vixinverse <- 1/(mixvix * vix + (1 - mixvix) * 
                        rollapply(vix,6 * 52, mean, align = "right", na.rm = T, fill = 'extend'))
     pa <- sweep(pa, MARGIN = 1, FUN = "*", 
                 STAT = vixinverse[index(pa),]/as.numeric(vixinverse[as.Date(da), ]))
  }
  if(normalise!="NONE") { pa <- zoonorm(pa,dimension=normalise) }
  fm <- fms2(
          pa,
          range.factors = c(nfac,nfac), 
          ...)
  fm1 <- addbench(fm,mcap = mcap.g[match(as.Date(da), index(mcap.g)),,drop = FALSE], bench = bench)
  dfrce(fm1,da=da)
}


#' extract object of class 'ce' from a co table
#'
#'
#'
#' largely for historical reasons ce objects are converted to a dataframe - this function converts the dataframe back to a 'ce'
#' 
#' @param cetab a data.frame or data.table
#' @param dat the date to extract
#' @family conversion
#' @examples 
#' \dontrun{
#' require(aace)
#' vcvce(dtce(getrdatv("jo","co"))) #covariance matrices
#' }
#' @export
dtce <- function (cetab, dat = cetab[,max(date)]) 
{
  stopifnot(length(dat) == 1)
  #stopifnot(valda(dat) && length(dat) == 1)
  cetab <- data.frame(cetab[date==dat])
  jbui <- grep("bui",colnames(cetab))
  jloadings <- grep("loadings",colnames(cetab))
  jfmp <- grep("fmp",colnames(cetab))
  jhpl <- grep("hpl",colnames(cetab))
  jmethod <- grep("method",colnames(cetab))
  jfull <- grep("full",colnames(cetab))
  juniqueness <- grep("uniqueness",colnames(cetab))
  jsdev <- grep("sdev",colnames(cetab))
  jqua <- grep("qua",colnames(cetab))
  jmvp <- match("mvp",colnames(cetab))
  jmcp <- match("mcp",colnames(cetab))
  jevp <- match("evp",colnames(cetab))
  jbeta <- grep("beta",colnames(cetab))
  jpoco <- grep("poco",colnames(cetab))
  bui <- cetab[,jbui,drop=TRUE]
  attributes(bui) <- NULL
  ifull <- which(cetab[,jfull]=="1")
  factors <- 1:ncol(cetab[,jloadings,drop=FALSE])
  ans <- list(
    loadings=as.matrix(cetab[,jloadings,drop=FALSE]), 
    fmp=as.matrix(cetab[,jfmp,drop=FALSE]),        
    hpl=as.matrix(cetab[,jhpl,drop=FALSE]),        
    method=as.matrix(cetab[,jmethod,drop=FALSE]),                   
    full=as.matrix(cetab[,jfull,drop=FALSE]),     
    uniqueness=as.matrix(cetab[,juniqueness,drop=FALSE]),
    sdev=as.matrix(cetab[,jsdev,drop=FALSE]),
    qua=as.matrix(cetab[,jqua,drop=FALSE]),
    evp=as.matrix(cetab[,jevp,drop=FALSE]),
    mcp=as.matrix(cetab[,jmcp,drop=FALSE]),
    mvp=as.matrix(cetab[,jmvp,drop=FALSE]),
    beta=as.matrix(cetab[,jbeta,drop=FALSE]),
    poco=as.matrix(cetab[,jpoco,drop=FALSE]),
    weight=NULL,
    call="restored"
  )
  #change mode for logical
  mode(ans$full) <- "numeric"
  mode(ans$full) <- "logical"
  #label
  dimnames(ans$loadings) <- list(bui,psz("loadings",factors))
  dimnames(ans$fmp) <- list(bui,psz("fmp",factors))
  dimnames(ans$hpl) <- list(bui,psz("hpl",factors))
  dimnames(ans$method) <- list(bui,psz("method",factors))
  dimnames(ans$full) <- list(bui,"full")
  dimnames(ans$uniqueness) <- list(bui,"uniqueness")
  dimnames(ans$sdev) <- list(bui,"sdev")
  dimnames(ans$evp) <- list(bui,"evp")
  dimnames(ans$mcp) <- list(bui,"mcp")
  dimnames(ans$mvp) <- list(bui,"mvp")
  dimnames(ans$beta) <- list(bui,c("evp","mcp","mvp"))
  dimnames(ans$poco) <- list(bui,"poco")
  class(ans) <- "ce"
  #stopifnot(valce(ans))
  ans
}

#dfrce - convert ce to a dataframe for storage as co
#'
#'
#'
#' largely for historical reasons ce objects are converted to a dataframe co
#' 
#' @param x an object of class ce; see package aace
#' @param dat the date to assign it
#' @family conversion
#' @examples 
#' \dontrun{
#' require(aace)
#' dfrce(dtce(getrdatv("jo","co"))) #a round trip: get a stored co object, convert to ce for one date, convert back to co
#' }
#' @export
dfrce <- function(x,dat=max(index(x))) 
{
  colnames(x$loadings) <- 1:ncol(x$loadings)  #fix colnames because data.frame prepends object name
  colnames(x$fmp) <- 1:ncol(x$fmp)
  colnames(x$hpl) <- 1:ncol(x$fmp)
  colnames(x$method) <- 1:ncol(x$method)
  res <- data.frame(
    date=dat,
    bui=as.matrix(rownames(x$loadings)),
    sdev=x$sdev,
    uniqueness=x$uniqueness,
    full=as.character(as.numeric(x$full)),
    loadings=x$loadings,
    fmp=x$fmp,
    hpl=x$hpl,
    method=metce(x),
    qua=x$qua,
    evp=x$evp,
    mcp=x$mcp,
    mvp=x$mvp,
    poco=x$poco,
    betaevp=x$beta[,"evp"],
    betamcp=x$beta[,"mcp"],
    betamvp=x$beta[,"mvp"]
  )
  names(res) <- gsub("\\.","",names(res))
  res
}


#' @export
interpce <- function(co=getrdatv("jo","co"),pa=getbdh(su=getrdatv("jo","su")),comp=c('T','M','S')) {
  ce <- dtce(co)
  stopifnot(all(sort(colnames(pa))==sort(buice(ce))))
  i <- is.na(coredata(pa))
  paz <- pa
  coredata(paz)[i] <- 0 #this does not affect the result if pa has same na as estimation window
  coredata(pa)[i] <- coredata(mscecomp(x=ce,ret=paz)[[comp]])[i]
  pa
}

#internals----------------------------------------
#' @export
addbench <- function(
  ce,
  mcap=t(ce$loadings[,1,drop=FALSE]*NA), 
  bench=c("equal","minvar","mcap")
)
{
  bench <- match.arg(bench)
  stopifnot(all(buice(ce)%in%colnames(mcap)))
  bui <- buice(ce)
  n <- length(bui)
  mcap <- coredata(mcap)[,bui,drop=TRUE] #nb no lagging is done to mcap, should already be lagged 1 period (ie opening mcap)
  if(mean(!is.na(mcap))<0.1) mcap <- mcap*NA      #insufficient data for cap-weighted
  if(mean(!is.na(mcap))>0) mcap[is.na(mcap)] <- min(mcap,na.rm=TRUE) #na assigned minimum weight
  mcap[mcap==0] <- min(mcap[mcap!=0],na.rm=TRUE) #0 assigned minimum weight
  vv <- vcvce(ce)$T
  stopifnot(all(colnames(vv)==colnames(mcap)))
  mvp <- solve(vv,rep(1,ncol(vv)))
  evp <- rep(1/ncol(vv),ncol(vv))
  mcp <- as.numeric(mcap/sum(mcap,na.rm=TRUE))
  mvp <- mvp/sum(mvp)
  evp <- evp/sum(abs(evp))
  ce$mvp <- matrix(mvp,n,1,dimn=list(bui,NULL))
  ce$evp <- matrix(evp,n,1,dimn=list(bui,NULL))
  ce$mcp <- matrix(mcp,n,1,dimn=list(bui,NULL))
  bpf <- switch(bench,
                equal=ce$evp,
                minvar=ce$mvp,
                mcap=ce$mcp
  )
  ce$beta <- cbind(
    vv%*%evp*1./as.numeric(((t(evp)%*%vv)%*%evp)),
    vv%*%mcp*1./as.numeric(((t(mcp)%*%vv)%*%mcp)),
    vv%*%mvp*1./as.numeric(((t(mvp)%*%vv)%*%mvp))
  )
  dimnames(ce$beta) <- list(bui,c("evp","mcp","mvp"))
  pol <- as.numeric(sign(t(ldgce(ce))%*%bpf))
  pol[which(pol==0)] <- 1
  ce$bench <- bench
  ce$fmp <- sweep(ce$fmp,STAT=pol,MAR=2,FUN="*")
  ce$hpl <- sweep(ce$hpl,STAT=pol,MAR=2,FUN="*")
  ce$loadings <- sweep(ce$loadings,STAT=pol,MAR=2,FUN="*")
  ce$poco <- matrix(0,n,1,dimn=list(bui,NULL))
  ce
}


#--------na interpolation section

#' getzco - covariance-based (PCA) interpolator 
#' 
#' @export
getzco <- function(co=getrd(100)) {
  ce <- dtce(co)
  list(
    M=fmpce(ce)[,1,drop=FALSE]%*%t(ldgce(ce)[,1,drop=FALSE]),
    S=fmpce(ce)[,-1,drop=FALSE]%*%t(ldgce(ce)[,-1,drop=FALSE]),
    T=fmpce(ce)%*%t(ldgce(ce))
  )    
}


# getzte - industry-based interpolator
# 
# the te object is read from rd, and has normally been read back from directories
# this works but is cumbersome/not needed so has been replaced by getzte that uses autoprune, see below
# getzte <- function(te=getrd(103),su=getrdatv("jo","su",2),da=su[,max(date)],loocv=FALSE) {
#   buix <- su[date==da,unique(bui)]
#   te <- setkey(te[buix,list(bui,bcode,BICS_REVENUE_PERC_LEVEL_ASSIGNED)],bui)
#   setkey(te[,BRPLA:=sum(BICS_REVENUE_PERC_LEVEL_ASSIGNED),list(bui,bcode)],bui,bcode)[,BICS_REVENUE_PERC_LEVEL_ASSIGNED:=NULL]
#   te <- unique(te)
#   x <- tabtomat(data.frame(te))
#   x[is.na(x)] <- 0
#   stopifnot(all(99<apply(x,1,sum)) & all(apply(x,1,sum)<101))
#   sol <- list(T=x%*%solve(crossprod(x))%*%t(x))
#   if(loocv) {
#     for( i in 1:nrow(x) ) {
#       sol[[i+1]] <- x[-i,]%*%solve(crossprod(x[-i,]))%*%t(x)
#     }
#     names(sol)[2:(1+nrow(x))] <- paste0('x',1:nrow(x))
#   }
#   sol
# }


#' integer industries pruned to nmin/node
#' 
#' @export
pruneztei <- function(su=getrdatv("jo","su",2),da=su[,max(date)],nmin=3) {
  buix <- su[date==da,unique(bui)]
  te <- getbdp()[buix][,list(bui,bcode=BICS_LEVEL_CODE_ASSIGNED,BRPLA=1)]
  tesums <- setcolorder(setkey(te[,list(agg=sum(BRPLA)),bcode],bcode)[,startcode:=bcode][,bagg:=agg],c('startcode','agg','bcode','bagg'))
  clen <- tesums[,max(nchar(bcode))]
  while(0<clen) {
    tesums[,acode:=bcode][((bagg<nmin)&(nchar(bcode)==clen)),acode:=substr(bcode,1,nchar(bcode)-2)]
    tesums[,bcode:=acode][,acode:=NULL][,bagg:=NULL][,bagg:=sum(agg),bcode]
    clen <- clen-2
  }
  te <- setkey(te,bcode)[setkey(tesums,startcode)][,bcode:=i.bcode][,i.bcode:=NULL][,agg:=NULL][,bagg:=NULL][bcode=="",bcode:="00"]
  te
}

#' fractional industries pruned to nmin/node
#' 
#' prunes the tree, leaving nodes satisfying criteria for aggregate weight and number of segments
#' @export
pruneztef <- function(su=getrdatv("jo","su",2),da=su[,max(date)],wmin=2,nmin=4) {
  buix <- su[date==da,unique(bui)]
  te <-  getbdp(mnem=data.table(data.frame(field = c("BICS_REVENUE_PERC_LEVEL_ASSIGNED"))))[buix][,list(bui,bcode=BICS_LEVEL_CODE_ASSIGNED,BRPLA=BICS_REVENUE_PERC_LEVEL_ASSIGNED/100)][!is.na(BRPLA)]
  tesums <- setcolorder(setkey(te[,list(agg=sum(BRPLA),count=.N),bcode],bcode)[,startcode:=bcode][,bagg:=agg][,bcount:=count],c('startcode','agg','count','bcode','bagg','bcount'))
  clen <- tesums[,max(nchar(bcode))]
  while(0<clen) {
    tesums[,acode:=bcode][(  ((bagg<wmin)|(bcount<nmin))  &(nchar(bcode)==clen)),acode:=substr(bcode,1,nchar(bcode)-2)]
    tesums[,bcode:=acode][,acode:=NULL][,bagg:=NULL][,bagg:=sum(agg),bcode][,bcount:=sum(count),bcode]
    clen <- clen-2
  }
  te <- setkey(te,bcode)[setkey(tesums,startcode)][,bcode:=i.bcode][,i.bcode:=NULL][,agg:=NULL][,bagg:=NULL][,count:=NULL][,bcount:=NULL][bcode=="",bcode:="00"]
  setkey(te[,list(BRPLA=sum(BRPLA)),list(bui,bcode)],bui,bcode)
}


#' z for integer industries pruned to nmin/node
#' 
#' @export
getztei <- function(su=getrdatv("jo","su",2),da=su[,max(date)],loocv=FALSE,nmin=3,wmin=3,type=c('i','f')) {
  type <- match.arg(type)
  if(type=='i') {
    te <- pruneztei(su=su,da=da,nmin=nmin)
  } else {
    te <- pruneztef(su=su,da=da,nmin=nmin,wmin=wmin)    
  }
  x <- tabtomat(data.frame(te))
  x[is.na(x)] <- 0
  sol <- list(T=x%*%solve(crossprod(x))%*%t(x))
  if(loocv) {
    for( i in 1:nrow(x) ) {
      sol[[i+1]] <- x[-i,]%*%solve(crossprod(x[-i,]))%*%t(x)
    }
    names(sol)[2:(1+nrow(x))] <- paste0('x',1:nrow(x))
  }
  sol
}

#'loocvi - drops each row in turn and estimates new co; fits the row - this only makes sense for ce and is the only version that makes sense for ce
#' @export
loocvi <- function(pa=getbdh(su),...) {
  mhatT <- mhatS <- mhatM <- m <- pa
  comp <- 
  for(i in 1:nrow(m)) {
    mi <- m[-i,]
    ce <- dtce(data.table(cewrap(pa=mi,...)))
    mscec <- mscecomp(ce,m[i,,drop=FALSE])
    mhatM[i,] <- mscec$M
    mhatS[i,] <- mscec$S
    mhatT[i,] <- mscec$T
  }
  mhatlist <- lapply(list(M=mhatM,S=mhatS,T=mhatT),as.numeric)
  mfitlist <- lapply(mscecomp(dtce(data.table(cewrap(pa=m,...))),m),as.numeric)
  names(mhatlist) <- paste0('mhat',names(mscec))
  names(mfitlist) <- paste0('mfit',names(mscec))
  c(list(act=as.numeric(m)),mhatlist,mfitlist)
}

#'loocvj - drops each column in turn of pa, ie each identifier bui - note there is NO iteration, no substitutio of the 'left out', it is just fitted which is more correct
#' @export
loocvj <- function(pa=getbdh(su),z=getzte(te=getrd(105),loocv=TRUE),...) {
  stopifnot(all(unlist(lapply(lapply(lapply(data.frame(t(data.frame(lapply(z,colnames)))),duplicated),"[",-1),all))))
  stopifnot(all.equal(colnames(pa),colnames(z[[1]])))
  mhat <- m <- coredata(pa)
  for(j in 1:ncol(m)) { #drop each column (identifier bui) in turn and postmultiply by the (x'x)-1.X
    mj <- m[,-j]
    xxj <- z[[paste0('x',j)]]
    mhat[,j] <- mj%*%(xxj[,j,drop=FALSE])
  }
  mfit <- m%*%z[[1]]
  list(act=as.numeric(m),mhat=as.numeric(mhat),mfit=as.numeric(mfit))
}

#'iterloocv - drops each observation in turn and replaces it with an iterated
#' @export
iterloocv <- function(pa=getbdh(su),z=getzco(),niter=2,myfun=c("as.numeric","mattotab","as.matrix")) {
  myfun <- get(match.arg(myfun))
  
  stopifnot(all(unlist(lapply(lapply(lapply(data.frame(t(data.frame(lapply(z,colnames)))),duplicated),"[",-1),all))))
  stopifnot(all.equal(colnames(pa),colnames(z[[1]])))
  mhat <- m <- coredata(pa)
  comp <- names(z)
  mfitlist <- mhatlist <- NULL
  for(i in 1:nrow(m)) {
    mi <- mean(m[i,])
    for(j in 1:ncol(m)) {
      mrow <- m[i,,drop=FALSE]
      mrow[,j] <- mi
      for(l in 1:niter) {mrow[1,j] <- (mrow%*%z[['T']])[1,j]}
      mhat[i,j] <- mrow[,j]
    }
  }
  for(k in comp) {
    mhatlist[[k]] <- mhat%*%z[[k]]
    mfitlist[[k]] <- m%*%z[[k]]
  }
  names(mhatlist) <- paste0('mhat',comp)
  names(mfitlist) <- paste0('mfit',comp)
  c(list(act=lapply(list(m),myfun)[[1]]),lapply(mhatlist,myfun),lapply(mfitlist,myfun))
}

# #' @export
# ilcvsumm <- function(x=iterloocv(...),...) {
#   s1 <- summary(lm(act ~ fit,data=data.frame(x)))
#   s2 <- summary(lm(act ~ loocv,data=data.frame(x)))
#   tab <- matrix(NA,5,2,dimnames=list(c('fit','loocv','diff','loocvM','loocvS'),c('Rsq','coef')))
#   tab[1,1] <- s1$r.squared
#   tab[1,2] <- s1$coefficients[2,1]
#   tab[2,1] <- s2$r.squared
#   tab[2,2] <- s2$coefficients[2,1]
#   tab[3,] <- tab[1,]-tab[2,]
#   if('loocvM'%in%names(x)) {
#     s3 <- summary(lm(act ~ loocvM,data=data.frame(x)))
#     s4 <- summary(lm(act ~ loocvS,data=data.frame(x)))
#     tab[4,1] <- s3$r.squared
#     tab[4,2] <- s3$coefficients[2,1]
#     tab[5,1] <- s4$r.squared
#     tab[5,2] <- s4$coefficients[2,1]
#   }
#   tab
# }

#' @export
ilcvsumm <- function(x=iterloocv(...),...) {
  tab <- NULL
  for(i in 2:length(x)) {
    mo <- summary(lm((paste0("act ~ ",names(x)[i])),data=data.frame(x)))
    tab <- data.table(rbind(tab,data.frame(rsq=mo$r.squared,coef=mo$coefficients[2,1])))
  }
  setkey(tab[,component:=names(x)[-1]],component)[]
}

#'summarises correlation of 'start' with final iteration
#' @export
iterate0 <- function(n=10,niter=5,FUN=mean,...) {
  suppressWarnings(do.call(FUN,list(unlist(lapply(lapply(lapply(1:n,FUN=iterate1,niter=niter,...),FUN=cor),FUN=`[`,1,niter)))))
}

#'wrapper to iterate2, applies 'drop' ie sets NA a specified fraction of data
#' @export
iterate1 <- function(seed=1,pa=getbdh(su),z=getzco()$T,idropfraction=0,initial=mean(pa,na.rm=TRUE),niter=5) {
  m <- coredata(pa)
  if(0<idropfraction) {
    mask <- m-m+1
    idrop <- round(length(mask)*idropfraction)
    mask[1:idrop] <- NA
    for(i in 1:nrow(mask)) mask[i,] <- mask[i,sample(1:ncol(m),size=ncol(m))]
    m <- m*mask
  }
  ina <- which(is.na(m))
  if(length(ina)==0) {
    return()
  } else {
    x <- iterate2(m,z,initial,niter)
  }
  if(0<idropfraction) x[1:idrop,'start'] <- coredata(pa)[ina]
  x
}

#'iterate2 - inner function applying z 
#' @export
iterate2 <- function(m,z,initial,niter=5) {
  ina <- which(is.na(m))
  res <- matrix(0*NA,length(ina),niter+2,dimnames=list(rownames(m)[ina],c('start',c(as.character(0:niter)))))
  res[,"0"] <- initial
  m[ina] <- initial
  for(i in 1:niter) {
    m[ina] <- (m%*%z)[ina]
    res[,as.character(i)] <- m[ina]
  }
  res
}

#'runco self-documenting run of 4 flavours of testing
#' @export
runco <- function(nmin=seq(from=4,to=16,by=2)) {
  require(aapa)
  require(aate)
  aatopselect('test')
  getbdhgl()
  su <- getrd(99)
  co <- getrd(100)
  pa <- getbdh(su)
  x <- vector("list")
  
  x[[length(x)+1]] <- structure(ilcvsumm(iterloocv(pa,z=getzco()))[,nmin:=nmin[1]],desc="PCA iterated")
  x[[length(x)+1]] <- structure(ilcvsumm(iterloocv(pa,z=getzte()))[,nmin:=nmin[1]],desc="BICSF iterated")
  x[[length(x)+1]] <- structure(ilcvsumm(loocvi(pa))[,nmin:=nmin[1]],desc="PCA non-iterated")
  x[[length(x)+1]] <- structure(ilcvsumm(loocvj(pa))[,nmin:=nmin[1]],desc="BICSF non-iterated")
  for(i in seq_along(nmin)) {
    x[[length(x)+1]] <- structure(ilcvsumm(iterloocv(pa,z=getztei(nmin=nmin[i])))[,nmin:=nmin[i]],desc="BICSI iterated")
    x[[length(x)+1]] <- structure(ilcvsumm(loocvj(pa,z=getztei(loocv=T,nmin=nmin[i])))[,nmin:=nmin[i]],desc="BICSI non-iterated")
  }
  x1 <- rbindlist(x)
  putrdatv(body(runco),app='ilcv',type=length(x))
  putrdatv(x1,app='ilcv',type='nmin')
  #lapply(x,putrd,usedesc=TRUE)
}


#fitplot - in a full panel that contains NA, drop size points from the non-NA, and crossplot fitted and actual
#' @export
fitplot <- function(size=min(1000,0.1*length(iok)),vbl=list(
  book.price=zoonorm(log(1+getbdh(su=su,fi='bopr.g')^-1)),
  upside=zoonorm(log(getbdh(su=su,fi='best.g')))),
  comp="T",
  ...) {
  set.seed(123)
  par(mfrow=c(1,1))
  for(i in seq_along(vbl)) {
    print(i)
    pa <- vbl[[i]]
    iok <- which(!is.na(coredata(pa)))
    izap <- sample(iok,siz=size,rep=F)
    act <- coredata(pa)[izap]
    coredata(pa)[izap]<-NA
    xint <- interpce(pa=pa,comp=comp)
    int <- coredata(xint)[izap]
    stopifnot(all(abs(act - int)>1e-10))
    lmx <- summary(lm(act~int))
    sub <- paste0('R-squared ',round(lmx$r.squared,2),' coeff: ',round(lmx$coefficients[2,1],2),' tstat: ',round(lmx$coefficients[2,3],1))
    plot(int,act,pch=20,col='grey',ylab='actual',xlab='fitted',main=paste0(names(vbl)[i]," (",comp,")"),sub=sub,...)
    grid()
    abline(0,1)
  }
}


colm <- function(pa,z,tau=-2:2) {
  nona <- !is.na(coredata(pa))
  x <- iterloocv(pa=pa,z=z,myfun='as.m')
  x1 <- x2 <- pa
  for(i in seq_along(tau)) {
    x1 <- x1+lag(x2,k=tau[i],na.pad=TRUE)
  }
  x3 <- data.table(mattotab(x1))[!is.na(field)]
  
#get valid M,S,R,bui,t and same for t+tau
#regress 
#T on lagged MS, delta MS, lagged R
#MS on lagged MS and delta MS
#R on lagged R
}