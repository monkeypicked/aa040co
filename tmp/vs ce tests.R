require(aapa)
resve("ce",896)
da <- as.character(extract("ce","max(date)"))
ce<-getce(da)
all(buice(ce) %in% getstep("NAME",n='000',typ='BDP')[,bui])

prem <- getpi("prem",da)
bui <- colnames(prem)
prem1 <- getstep("x0700rrdoto",n='001',typ='BDH')
all(bui%in%colnames(prem1))
prem2 <- prem1[index(prem),bui]
identical(colnames(prem),colnames(prem2))
identical(rownames(prem),rownames(prem2))
imgzoo(prem)
imgzoo(prem2)
prem-prem2
#these are now v similar, difference may be due to a full week locf in the old version? not checked


#dates where there are a lot of hols
apply(is.na(prem2),1,sum)
da <- index(prem2)[which(apply(is.na(prem2),1,sum)>100)]
bui1 <- colnames(prem2)[which(is.na(prem2[da[4]]))]
nn <- getstep("NAME",n='000',typ='BDP')
ee <- getstep("EQY_PRIM_EXCH",n='000',typ='BDP')
tk <- getstep("TICKER_AND_EXCH_CODE",n='000',typ='BDP')
edit(tk[nn][ee][bui1]) #JP & HK



require(aapa)
buix <- 'EQ0736411200001000' #2590 jp
tk[nn][ee][buix]
getstep("x0603PL_TTU",n='001',typ='BDH')[,buix,drop=F]
getstep("x0501PL_TFU",n='001',typ='BDH')[,buix,drop=F]
getstep("x0200PL_TFU",n='001',typ='BDH')[,buix,drop=F]

rr <- getstep("x0700rrdoto",n='001',typ='BDH')[,buix,drop=F]
pp <- getstep("x0700prdo",n='001',typ='BDH')[,buix,drop=F]


x1 <- getstep("x0200PL_TTU")[,buix,drop=F]
x2 <- f0603wed(x1)

f0603wed(getstep("x0200PL_TTU"))

xx <- x1[2860:2880,drop=F]
dtlocf(xx,wd=3,roll=-5)
xx


x0601PL_TFU

f0600wed <- function(x){  #locf
  dtlocf(x,wd=3,roll=4)
}
