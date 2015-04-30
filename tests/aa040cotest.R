require(aautil)
aatopselect("test")
require(aaco)
require(aapa)
require(quadprog)
require(testthat)


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
