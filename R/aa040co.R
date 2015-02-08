
#' get global zoo panels
#'
#' @export
getbdhgl <- function(field=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap"),list(vix.g="x0502vix"))) {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep(unlist(field),n='001'),envir=globalenv())})
}

#' get global macro zoo panels
#'
#' @export
getbdmgl <- function(field=list(list(vix.g="VIX"))) {
  x <- lapply(field,function(field) {assign(x=names(field),value=getstep("VIX",n="000",ty='m'),envir=globalenv())})
}

#' return windowed zoo panel from global panel
#'
#' @export
getbdh <- function(da=su.g[,max(date)],bui=colnames(x),field="prem.g",win=-(1:230)) {
  x <- get(field,envir=globalenv())
  x[as.Date(offda(da,win)),bui]
}

#' outer wrapper, loads absent global objects, calls cewrap
#'
#' @export
ceload <- function(...,fieldrd=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap")),fieldm="vix.g",isu=greprd()) {
  if(any(!unlist(lapply(names(unlist(fieldrd)),exists)))) {getbdhgl(field=fieldrd) }
  if(any(!exists(fieldm))) {getbdmgl() }
  if(!exists('su.g')) {su.g<<-getrd(isu)}
  cewrap(...)
}



#' derive ce
#'
#' @export
cewrap <- function (da=su.g[,max(date)], win=-(1:230), normalise="NONE", nfac=20,
                       applyvix=TRUE, mixvix=.5, bench="equal",...) {
  pa <- getbdh(da=da,bui=su.g[date==da][,bui],field="prem.g",win=win)
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


#' extract object of class 'ce' from a ce table
#'
#' @export
# dtce <- function (cetab, dat = cetab[,max(date)]) 
# {
#   stopifnot(length(dat) == 1)
#   #stopifnot(valda(dat) && length(dat) == 1)
#   cetab <- data.frame(cetab[date==dat])
#   jbui <- grep("bui",colnames(cetab))
#   jloadings <- grep("loadings",colnames(cetab))
#   jfmp <- grep("fmp",colnames(cetab))
#   jhpl <- grep("hpl",colnames(cetab))
#   jmethod <- grep("method",colnames(cetab))
#   jfull <- grep("full",colnames(cetab))
#   juniqueness <- grep("uniqueness",colnames(cetab))
#   jsdev <- grep("sdev",colnames(cetab))
#   jqua <- grep("qua",colnames(cetab))
#   jmvp <- match("mvp",colnames(cetab))
#   jmcp <- match("mcp",colnames(cetab))
#   jevp <- match("evp",colnames(cetab))
#   jbeta <- grep("beta",colnames(cetab))
#   jpoco <- grep("poco",colnames(cetab))
#   bui <- cetab[,jbui,drop=TRUE]
#   attributes(bui) <- NULL
#   ifull <- which(cetab[,jfull]=="1")
#   factors <- 1:ncol(cetab[,jloadings,drop=FALSE])
#   ans <- list(
#     loadings=as.matrix(cetab[,jloadings,drop=FALSE]), 
#     fmp=as.matrix(cetab[,jfmp,drop=FALSE]),        
#     hpl=as.matrix(cetab[,jhpl,drop=FALSE]),        
#     method=as.matrix(cetab[,jmethod,drop=FALSE]),                   
#     full=as.matrix(cetab[,jfull,drop=FALSE]),     
#     uniqueness=as.matrix(cetab[,juniqueness,drop=FALSE]),
#     sdev=as.matrix(cetab[,jsdev,drop=FALSE]),
#     qua=as.matrix(cetab[,jqua,drop=FALSE]),
#     evp=as.matrix(cetab[,jevp,drop=FALSE]),
#     mcp=as.matrix(cetab[,jmcp,drop=FALSE]),
#     mvp=as.matrix(cetab[,jmvp,drop=FALSE]),
#     beta=as.matrix(cetab[,jbeta,drop=FALSE]),
#     poco=as.matrix(cetab[,jpoco,drop=FALSE]),
#     weight=NULL,
#     call="restored"
#   )
#   #change mode for logical
#   mode(ans$full) <- "numeric"
#   mode(ans$full) <- "logical"
#   #label
#   dimnames(ans$loadings) <- list(bui,psz("loadings",factors))
#   dimnames(ans$fmp) <- list(bui,psz("fmp",factors))
#   dimnames(ans$hpl) <- list(bui,psz("hpl",factors))
#   dimnames(ans$method) <- list(bui,psz("method",factors))
#   dimnames(ans$full) <- list(bui,"full")
#   dimnames(ans$uniqueness) <- list(bui,"uniqueness")
#   dimnames(ans$sdev) <- list(bui,"sdev")
#   dimnames(ans$evp) <- list(bui,"evp")
#   dimnames(ans$mcp) <- list(bui,"mcp")
#   dimnames(ans$mvp) <- list(bui,"mvp")
#   dimnames(ans$beta) <- list(bui,c("evp","mcp","mvp"))
#   dimnames(ans$poco) <- list(bui,"poco")
#   class(ans) <- "ce"
#   #stopifnot(valce(ans))
#   ans
# }

#dfrce - returns ce as suitably labelled dataframe

#' @export
`dfrce` <- function(x,
                    dat=getda("no")) 
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
