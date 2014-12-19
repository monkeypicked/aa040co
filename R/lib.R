require(aautil)
require(aapa)
require(aa) #meantri

#getbdhgl - get global zoo panels
getbdhgl <- function(field=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap"),list(vix.g="x0502vix"))) {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep(unlist(field),n='001'),envir=globalenv())})
}

#getbdmgl - get global macro zoo panels
getbdmgl <- function(field=list(list(vix.g="VIX"))) {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep("VIX",n="000",ty='m'),envir=globalenv())})
}

#getbdh - return windowed zoo panel from global panel
getbdh <- function(da="2014-10-29",field="prem.g",win=-(1:230)) {
  get(field,envir=globalenv())[as.Date(offda(da,win)),]
}

#ceload - outer wrapper
ceload <- function(...,fieldrd=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap"))) {
  getbdhgl(field=fieldrd) 
  getbdmgl() 
  cewrap(...)
}

#dfrce from aa
dfrce <- function(x,dat=getda("no")) {
  stopifnot(valda(dat) && length(dat)==1)
  stopifnot(class(x) %in% c("ce"))
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

#addbench - totally unmodified from ce
`addbench` <- function(
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


cewrap <- function (da="2014-10-29", win=-(1:230), normalise="NONE", nfac=20,
                       applyvix=FALSE, mixvix=.5, bench="equal",...) {
  pa <- getbdh(da=da,field="prem.g",win=win)[,1:100]
  if (applyvix) {
     vix <- vix.g #should be loaded from 0502VIX
     vixinverse <- 1/(mixvix * vix + (1 - mixvix) * 
                        rollapply(vix,6 * 52, mean, align = "right", na.rm = T, fill = NA))
     pa <- sweep(pa, MARGIN = 1, FUN = "*", 
                 STAT = vixinverse[index(pa),]/as.numeric(vixinverse[as.Date(da), ]))
  }
  if(normalise!="NONE") { pa <- zoonorm(pa,dimension=normalise) }
  fm <- fms2(
          pa,
          range.factors = c(nfac,nfac), 
          ...)
  dfrce(addbench(
          fm, 
          mcap = mcap.g[match(as.Date(da), index(mcap.g)),,drop = FALSE], bench = bench
          ),da=da)
}

