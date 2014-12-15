require(aautil)


#some data prep
#su <- setkey(extract("su",c("date","bui"),res=data.table),date,bui)
#putrd(su,'su isu=312')

#some functions to replace bits inside ceWrapper
#preparevix
#preparemcap



ceWrapper <- function (dat, field, win, normalise, center, weight, nfac, lambda, 
                       mcap, bench, applyvix, mixvix, vixsmooth, shrinkb, ...) 

if (applyvix) {
  vix <- rollxts(x = getvix(mylag = -1), 
                 what = "meantri", 
                 n = vixsmooth
                 )
  vixinverse <- 1/(mixvix * vix + (1 - mixvix) * 
                     rollapply(vix,6 * 52, mean, align = "right", na.rm = T, fill = NA))
  pa <- sweep(pa, MARGIN = 1, FUN = "*", 
              STAT = vixinverse[index(pa),]/as.numeric(vixinverse[as.Date(dat), ]))
}

fm <- fms2(
        pa, 
        center = center, 
        weight = weight, 
        range.factors = c(nfac,nfac), 
        lambda = lambda, 
        shrinkb = shrinkb, 
        ...)
fm <- addbench(
        fm, 
        mcap = mcap[match(as.Date(dat), index(mcap)),,drop = FALSE], bench = bench
        )
putce(x = fm, tabname = "ce", dat = dat)