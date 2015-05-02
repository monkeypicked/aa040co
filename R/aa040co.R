
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
getbdhgl <- function(field=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap"),list(vix.g="x0502vix")),myclass='zoo') {
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

interpce <- function(co=getrdatv("jo","co"),pa=getbdh(su=getrdatv("jo","su"))) {
  ce <- dtce(co)
  stopifnot(all(sort(colnames(pa))==sort(buice(ce))))
  i <- is.na(coredata(pa))
  paz <- pa
  coredata(paz)[i] <- 0 #this does not affect the result if pa has same na as estimation window
  coredata(pa)[i] <- coredata(msce(x=ce,ret=paz))[i]
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

#presumed junk----------------------------------------

#getpi is overcomplex; slow, probably buggy - replaced with getbdh which just uses zoo - simples
# #' Get panel for universe defined by su
# #'
# #' gets the panel corresponding to a securityuniverse, either applying locf or focb to the
# #' securityuniverse, then optionally extracting one date and applying focb to that.  The default
# #' is locf-focb, which means that the latest ex-ante universe identifies the cross-section, whose
# #' history is then extracted.  This is the correct usage for 'rolling ex-ante estimation' on a
# #' dynamic univers.
# #' @param su security universe, a data.table with columns date, bui
# #' @param mnem filename without extension
# #' @param mydir directory
# #' @param myclass 'zoo' or 'dt' is returned
# #' @param da reference date or datum for universe (this in incremented each time in a rolling estimation)
# #' @param la maximum lag to access as a positive number, so data returned starts at (da-lag)
# #' @param roll flag taking values 'locf' or 'focb', the former being correct for ex-ante estimation
# #' @param fixed logical flag indicating whether the universe is fixed as of da, or dynamic, defaults TRUE
# #' @param ... passed to dern to construct mydir
# #' @examples
# #' \dontrun{
# #' su <- getrdatv('jo','su')
# #' pa1 <- getpi(su)
# #' pa2 <- getpi(su,fix=F)
# #' pa3 <- getpi(su,da='2011-11-30')
# #' pa4 <- getpi(su,fix=F,da='2011-11-30')
# #' mean(is.na(pa2))-mean(is.na(pa1))#has more na because history screened out
# #' mean(is.na(pa4))-mean(is.na(pa3))#has more na because history screened out
# #' }
# #' @export
# #' @family accessor
# getpi <- function(su=getrdatv('jo','su'),
#                   x=prem.g, 
#                   mydir = dern(...), 
#                   myclass=c("zoo","dt"),
#                   da=max(x[,date]),
#                   la=200,
#                   roll=TRUE,
#                   fixed=TRUE,
#                   ...) {
#   myclass <- match.arg(myclass)
#   stopifnot(mode(da)=="character" && da%in%x[,unique(date)])
#   stopifnot(
#     is(su,'data.table') && 
#       all(c('bui','date')%in%colnames(su)) && 
#       all(su[,date]%in%derca()[,date]) && 
#       length(setdiff(su[,unique(bui)],x[,unique(bui)]))==0)
#   datevec <- x[,list(date=unique(date))][date<=da,list(date=sort(date,decreasing=TRUE))][1:la][,rev(date)]
#   if(su[,mode(date)!='character']) su[,date:=as.character(date)]
#   #interpolate su onto the basis of x, note the ta==1 which is then deleted
#   su1 <- unique(setkey(su[,ta:=1],bui,date))[CJ(su[,unique(bui)],x[,unique(date)]),roll=roll][,list(bui,date,ta)][ta==1][,ta:=NULL]
#   if(fixed) {
#     su2 <- setnames(CJ(setkey(su1,date)[da,bui],datevec),c('bui','date'))
#   } else {
#     su2 <- setkey(su1,date)[datevec]
#   }
#   x <- setkey(x,bui,date)[setkey(su2,bui,date)]
#   if(myclass=="zoo") {
#     mz(tabtomat(data.frame(x)))
#   } else { x }
# }

# getbdmgl <- function(field=list(list(vix.g="VIX")),myclass='zoo') {
#   x <- lapply(field,function(field) {assign(x=names(field),value=getstep("VIX",n="000",ty='m',myclass=myclass),envir=globalenv())})
# }

# ceload <- function(...,fieldrd=list(list(prem.g="x0700redoto"),list(mcap.g="x0702mcap")),fieldm="vix.g",isu=greprd()) {
#   if(any(!unlist(lapply(names(unlist(fieldrd)),exists)))) {getbdhgl(field=fieldrd) }
#   if(any(!exists(fieldm))) {getbdmgl() }
#   if(!exists('su.g')) {su.g<<-getrd(isu)}
#   cewrap(...)
# }

#--------na interpolation section

#' getzco - covariance-based (PCA) interpolator 
#' 
#' @export
#' 
getzco <- function(co=getrd(101)) {
  ce <- dtce(co)
  fmpce(ce)%*%t(ldgce(ce))
}


#' getzte - industry-based interpolator
#' 
#' @export
#' 
getzte <- function(te=getrd(103),su=getrdatv("jo","su",2),da=su[,max(date)]) {
  buix <- su[date==da,unique(bui)]
  te <- setkey(te[buix,list(bui,bcode,BICS_REVENUE_PERC_LEVEL_ASSIGNED)],bui)
  setkey(te[,BRPLA:=sum(BICS_REVENUE_PERC_LEVEL_ASSIGNED),list(bui,bcode)],bui,bcode)[,BICS_REVENUE_PERC_LEVEL_ASSIGNED:=NULL]
  te <- unique(te)
  x <- tabtomat(data.frame(te))
  x[is.na(x)] <- 0
  stopifnot(all(99<apply(x,1,sum)) & all(apply(x,1,sum)<101))
  x%*%solve(crossprod(x))%*%t(x)
}



#this function has been split into iterate1,2
# iterate <- function(pa=getbdh(su),z=getz(),idrop=NULL,initial=mean(pa,na.rm=TRUE),niter=5) {
#   m <- coredata(pa)
#   m[idrop] <- NA
#   ina <- which(is.na(m))
#   nna <- length(ina)
#   if(nna==0) return()
#   if(length(ina)>0) {
#     res <- matrix(0*NA,nna,niter+2,dimnames=list(rownames(m)[ina],c('start',c(as.character(0:niter)))))
#     res[,"0"] <- initial
#     m[ina] <- initial
#     for(i in 1:niter) {
#       m[ina] <- (m%*%z)[ina]
#       res[,as.character(i)] <- m[ina]
#     }
#   }
#   if(length(idrop)==nna) res[,'start'] <- coredata(pa)[ina]
#   res
# }


#summarises correlation of 'start' with final iteration
iterate0 <- function(n=10,niter=5,FUN=mean,...) {
  suppressWarnings(do.call(FUN,list(unlist(lapply(lapply(lapply(1:n,FUN=iterate1,niter=niter,...),FUN=cor),FUN=`[`,1,niter)))))
}

#iterloocvi - drops each row in turn and estimates new co; fits the row
iterloocvi <- function(pa=getbdh(su),...) {
  mhat <- m <- pa
  for(i in 1:nrow(m)) {
    mi <- m[-i,]
    ce <- dtce(data.table(cewrap(pa=mi,...)))
    mhat[i,] <- msce(ce,m[i,,drop=FALSE])
  }
  mfit <- msce(x=dtce(data.table(cewrap(pa=mi,...))),ret=pa)
  list(loocv=as.numeric(coredata(mhat)),fit=as.numeric(coredata(mfit)),act=as.numeric(coredata(m)))
}

#iterloocv - drops each observation in turn and replaces it with an iterated
iterloocv <- function(pa=getbdh(su),z=getzco(),niter=2) {
  mhat <- m <- coredata(pa)
  for(i in 1:nrow(m)) {
    mi <- mean(m[i,])
    for(j in 1:ncol(m)) {
      mrow <- m[i,,drop=FALSE]
      mrow[,j] <- mi
      for(k in 1:niter) {mrow[1,j] <- (mrow%*%z)[1,j]}
      mhat[i,j] <- mrow[,j]
    }
  }
  mfit <- m %*% z
  list(loocv=as.numeric(mhat),fit=as.numeric(mfit),act=as.numeric(m))
}

ilcvsumm <- function(x=iterloocv(...),...) {
  s1 <- summary(lm(act ~ fit,data=data.frame(x)))
  s2 <- summary(lm(act ~ loocv,data=data.frame(x)))
  tab <- matrix(NA,3,2,dimnames=list(c('fit','loocv','diff'),c('Rsq','coef')))
  tab[1,1] <- s1$r.squared
  tab[1,2] <- s1$coefficients[2,1]
  tab[2,1] <- s2$r.squared
  tab[2,2] <- s2$coefficients[2,1]
  tab[3,] <- tab[1,]-tab[2,]
  tab
}

#wrapper to iterate2, applies 'drop' ie sets NA a specified fraction of data
iterate1 <- function(seed=1,pa=getbdh(su),z=getzco(),idropfraction=0,initial=mean(pa,na.rm=TRUE),niter=5) {
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
  if(idrop==length(ina)) x[,'start'] <- coredata(pa)[ina]
  x
}

#iterate2 - inner function applying z 
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
