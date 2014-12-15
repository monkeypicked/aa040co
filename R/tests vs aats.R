require(aapa)
x <- setkey(getstep("CRNCY_ADJ_MKT_CAP",n='000',typ='BDP'),CRNCY_ADJ_MKT_CAP)
require(aa0)
prdo <- getpa("prdo")

bigbui <- x[nrow(x)-(seq(1,1000)-1),bui]

nn <- getstep("NAME",n='000',typ='BDP')
nn[bigbui]

#aats and rd match in price on wed
plt1 <- getstep("PX_LAST_TTU",n='000',typ='BDH')
plt2 <- getstep("x0603PL_TTU",n='001',typ='BDH')
bui <- intersect(intersect(colnames(plt2),colnames(prdo)),bigbui)[1]
test3 <- cbind(plt2[,bui,drop=F],prdo[,bui,drop=F]) #perfick

#return is as expected for rd  #NO LONGER WHEN APPLY AS.XTS!!!!!!!!!!!!!!!!!!!!!
x07 <- getstep("x0700redoto",n='001',typ='BDH')
test0 <-cbind(as.xts(x07[,bui]),retxts(as.xts(plt2[,bui,drop=F])))


#check return for aats
redoto1 <- getpa("redoto")
test <- retxts(as.xts(prdo[,bui,drop=F]))
test1 <- cbind(as.zoo(test),redoto1[,bui])
               

#returns match
redoto1 <- getpa("redoto")
redoto2 <- getstep("x0700redoto",n='001',typ='BDH')
xx <- cbind(redoto2[,bui],redoto1[,bui])
mean(xx[,2]-xx[,1],na.rm=T)*52


#premia
rrdoto1 <- getpa("rrdoto")
rrdoto2 <- getstep("x0700rrdoto",n='001',typ='BDH')
xx <- cbind(rrdoto2[,bui],rrdoto1[,bui])
mean(xx[,2]-xx[,1],na.rm=T)*52



#so plt2 matches prdo; plt2 -> test2; prdo -> test1; and test1 does not match test2
cbind(retxts(as.xts(test3)),as.xts(redoto1[,bui,drop=F]),as.xts(x07[,bui,drop=F]))



bui<-bigbui[1]
plt1 <- getstep("PX_LAST_TTU",n='000',typ='BDH')
plt2 <- getstep("x0603PL_TTU",n='001',typ='BDH')
cbind(plt2[,bui],prdo[,bui]) #perfick




plot(cbind(plt1[,bui],plt2[,bui]),col=1:2)



plt <- plt2
#x07 <- getstep("x0700rrdoto",n='001',typ='BDH')













require(aa0)
prdo <- getpa("prdo")
rrdoto <- getpa("rrdoto")


bui <- intersect(intersect(colnames(plt),colnames(prdo)),bigbui)[1]
plot(cbind(plt1[,bui],plt2[,bui]),scr=1,col=1:2)




plt[,bui]/prdo[,bui]

zoo(is.na(plt[,bui]),index(plt))/zoo(is.na(prdo[,bui]),index(prdo))


rrdoto[,bui]/x07[,bui]

buix<-bui[1]
xx<-cbind(rrdoto[,buix],x07[,buix])
coredata(xx)[is.na(coredata(xx))] <- 0
plot(cumsum(xx),scr=1,col=1:2)


zoo(is.na(plt[,bui]),index(plt))/zoo(is.na(prdo[,bui]),index(prdo))

#second item is much lower return
