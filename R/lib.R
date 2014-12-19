require(aautil)
require(aapa)
require(aa0) #meantri, dfrce, addbench
require(aace)

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
getbdh <- function(da="2014-10-29",bui=colnames(x),field="prem.g",win=-(1:230)) {
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
cewrap <- function (da="2014-10-29", win=-(1:230), normalise="NONE", nfac=20,
                       applyvix=TRUE, mixvix=.5, bench="equal",...) {
  pa <- getbdh(da=da,bui=su.g[da][,bui],field="prem.g",win=win)
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
dtce <- function (cetab, dat = cetab[,max(date)]) 
{
  stopifnot(valda(dat) && length(dat) == 1)
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
  stopifnot(valce(ans))
  ans
}
