# Analysis: Age mixing with transmission clusters


# 
# > dim(df)
# [1]   20 1877
# > dim(DF)
# [1]   48 1876

df <- read.csv("/home/niyukuri/Dropbox/2.5.2018.Simpact/study2/master.model.age.mixing.pattern.study2.csv")

DF <- read.csv("/home/niyukuri/Dropbox/2.5.2018.Simpact/study2/master.model.age.mixing.pattern.study21.csv")


stat.fun.men.women.ratio <- function(df=df){
  
  tot.pairs <- df[,7] + df[,8] + df[,9] + df[,10] + df[,11] + df[,12] + df[,13] + df[,14] + df[,15]
  
  men.15.25.w.15.25 <- df[,7]/tot.pairs
  men.15.25.w.25.40 <- df[,8]/tot.pairs
  men.15.25.w.40.50 <- df[,9]/tot.pairs
  men.25.40.w.15.25 <- df[,10]/tot.pairs
  men.25.40.w.25.40 <- df[,11]/tot.pairs
  men.25.40.w.40.50 <- df[,12]/tot.pairs
  men.40.50.w.15.25 <- df[,13]/tot.pairs
  men.40.50.w.25.40 <- df[,14]/tot.pairs
  men.40.50.w.40.50 <- df[,15]/tot.pairs
  
  return(cbind(men.15.25.w.15.25, men.15.25.w.25.40, men.15.25.w.40.50,
               men.25.40.w.15.25, men.25.40.w.25.40, men.25.40.w.40.50,
               men.40.50.w.15.25, men.40.50.w.25.40, men.40.50.w.40.50))
  
}


DF.r <- DF[,-1]

pop.stat <- DF.r[,1:3]

## I. Transmissions
####################

# MCAR

transm.CAR.cov.35 <- DF.r[,4:21] # +17

transm.CAR.cov.35.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.35)

transm.CAR.cov.40 <- DF.r[,22:39]

transm.CAR.cov.40.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.40)

transm.CAR.cov.45 <- DF.r[,40:57]

transm.CAR.cov.45.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.45)

transm.CAR.cov.50 <- DF.r[,58:75]

transm.CAR.cov.50.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.50)

transm.CAR.cov.55 <- DF.r[,76:93]

transm.CAR.cov.55.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.55)

transm.CAR.cov.60 <- DF.r[,94:111]

transm.CAR.cov.60.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.60)

transm.CAR.cov.65 <- DF.r[,112:129]

transm.CAR.cov.65.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.65)

transm.CAR.cov.70 <- DF.r[,130:147]

transm.CAR.cov.70.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.70)

transm.CAR.cov.75 <- DF.r[,148:165]

transm.CAR.cov.75.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.75)

transm.CAR.cov.80 <- DF.r[,166:182]

transm.CAR.cov.80.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.80)

transm.CAR.cov.85 <- DF.r[,184:199]

transm.CAR.cov.85.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.85)

transm.CAR.cov.90 <- DF.r[,202:219]

transm.CAR.cov.90.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.90)

transm.CAR.cov.95 <- DF.r[,220:237]

transm.CAR.cov.95.tab <- stat.fun.men.women.ratio(df=transm.CAR.cov.95)


# MAR

# a.

transm.AR.a.cov.35 <- DF.r[,238:255] # +17

transm.AR.a.cov.35.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.35)

transm.AR.a.cov.40 <- DF.r[,256:273]

transm.AR.a.cov.40.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.40)

transm.AR.a.cov.45 <- DF.r[,274:291]

transm.AR.a.cov.45.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.45)

transm.AR.a.cov.50 <- DF.r[,292:309]

transm.AR.a.cov.50.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.50)

transm.AR.a.cov.55 <- DF.r[,310:327]

transm.AR.a.cov.55.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.55)

transm.AR.a.cov.60 <- DF.r[,328:345]

transm.AR.a.cov.60.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.60)

transm.AR.a.cov.65 <- DF.r[,346:363]

transm.AR.a.cov.65.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.65)

transm.AR.a.cov.70 <- DF.r[,364:381]

transm.AR.a.cov.70.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.70)

transm.AR.a.cov.75 <- DF.r[,382:399]

transm.AR.a.cov.75.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.75)

transm.AR.a.cov.80 <- DF.r[,400:417]

transm.AR.a.cov.80.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.80)

transm.AR.a.cov.85 <- DF.r[,418:435]

transm.AR.a.cov.85.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.85)

transm.AR.a.cov.90 <- DF.r[,436:453]

transm.AR.a.cov.90.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.90)

transm.AR.a.cov.95 <- DF.r[,454:471]

transm.AR.a.cov.95.tab <- stat.fun.men.women.ratio(df=transm.AR.a.cov.95)


# b

transm.AR.b.cov.35 <- DF.r[,472:489] # +17

transm.AR.b.cov.35.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.35)

transm.AR.b.cov.40 <- DF.r[,490:507]

transm.AR.b.cov.40.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.40)

transm.AR.b.cov.45 <- DF.r[,508:525]

transm.AR.b.cov.45.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.45)

transm.AR.b.cov.50 <- DF.r[,526:543]

transm.AR.b.cov.50.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.50)

transm.AR.b.cov.55 <- DF.r[,544:561]

transm.AR.b.cov.55.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.55)

transm.AR.b.cov.60 <- DF.r[,562:579]

transm.AR.b.cov.60.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.60)

transm.AR.b.cov.65 <- DF.r[,580:597]

transm.AR.b.cov.65.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.65)

transm.AR.b.cov.70 <- DF.r[,598:615]

transm.AR.b.cov.70.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.70)

transm.AR.b.cov.75 <- DF.r[,616:633]

transm.AR.b.cov.75.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.75)

transm.AR.b.cov.80 <- DF.r[,634:651]

transm.AR.b.cov.80.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.80)

transm.AR.b.cov.85 <- DF.r[,652:669]

transm.AR.b.cov.85.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.85)

transm.AR.b.cov.90 <- DF.r[,670:687]

transm.AR.b.cov.90.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.90)

transm.AR.b.cov.95 <- DF.r[,688:705]

transm.AR.b.cov.95.tab <- stat.fun.men.women.ratio(df=transm.AR.b.cov.95)

# c

transm.AR.c.cov.35 <- DF.r[,706:723] # +17

transm.AR.c.cov.35.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.35)

transm.AR.c.cov.40 <- DF.r[,724:741]

transm.AR.c.cov.40.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.40)

transm.AR.c.cov.45 <- DF.r[,742:759]

transm.AR.c.cov.45.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.45)

transm.AR.c.cov.50 <- DF.r[,760:777]

transm.AR.c.cov.50.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.50)

transm.AR.c.cov.55 <- DF.r[,778:795]

transm.AR.c.cov.55.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.55)

transm.AR.c.cov.60 <- DF.r[,796:813]

transm.AR.c.cov.60.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.60)

transm.AR.c.cov.65 <- DF.r[,814:831]

transm.AR.c.cov.65.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.65)

transm.AR.c.cov.70 <- DF.r[,832:849]

transm.AR.c.cov.70.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.70)

transm.AR.c.cov.75 <- DF.r[,850:867]

transm.AR.c.cov.75.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.75)

transm.AR.c.cov.80 <- DF.r[,868:885]

transm.AR.c.cov.80.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.80)

transm.AR.c.cov.85 <- DF.r[,886:903]

transm.AR.c.cov.85.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.85)

transm.AR.c.cov.90 <- DF.r[,904:921]

transm.AR.c.cov.90.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.90)

transm.AR.c.cov.95 <- DF.r[,922:939]

transm.AR.c.cov.95.tab <- stat.fun.men.women.ratio(df=transm.AR.c.cov.95)


## II. Clusters
################

# MCAR

clust.CAR.cov.35 <- DF.r[,940:957] # +17

clust.CAR.cov.35.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.35)

clust.CAR.cov.40 <- DF.r[,958:975]

clust.CAR.cov.40.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.40)

clust.CAR.cov.45 <- DF.r[,976:993]

clust.CAR.cov.45.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.45)

clust.CAR.cov.50 <- DF.r[,994:1011]

clust.CAR.cov.50.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.50)

clust.CAR.cov.55 <- DF.r[,1012:1029]

clust.CAR.cov.55.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.55)

clust.CAR.cov.60 <- DF.r[,1030:1047]

clust.CAR.cov.60.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.60)

clust.CAR.cov.65 <- DF.r[,1048:1065]

clust.CAR.cov.65.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.65)

clust.CAR.cov.70 <- DF.r[,1066:1083]

clust.CAR.cov.70.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.70)

clust.CAR.cov.75 <- DF.r[,1084:1101]

clust.CAR.cov.75.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.75)

clust.CAR.cov.80 <- DF.r[,1102:1119]

clust.CAR.cov.80.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.80)

clust.CAR.cov.85 <- DF.r[,1120:1137]

clust.CAR.cov.85.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.85)

clust.CAR.cov.90 <- DF.r[,1138:1155]

clust.CAR.cov.90.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.90)

clust.CAR.cov.95 <- DF.r[,1156:1173]

clust.CAR.cov.95.tab <- stat.fun.men.women.ratio(df=clust.CAR.cov.95)


# MAR

# a.

clust.AR.a.cov.35 <- DF.r[,1174:1191] # +17

clust.AR.a.cov.35.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.35)

clust.AR.a.cov.40 <- DF.r[,1192:1209]

clust.AR.a.cov.40.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.40)

clust.AR.a.cov.45 <- DF.r[,1210:1227]

clust.AR.a.cov.45.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.45)

clust.AR.a.cov.50 <- DF.r[,1228:1245]

clust.AR.a.cov.50.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.50)

clust.AR.a.cov.55 <- DF.r[,1246:1263]

clust.AR.a.cov.55.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.55)

clust.AR.a.cov.60 <- DF.r[,1264:1281]

clust.AR.a.cov.60.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.60)

clust.AR.a.cov.65 <- DF.r[,1282:1299]

clust.AR.a.cov.65.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.65)

clust.AR.a.cov.70 <- DF.r[,1300:1317]

clust.AR.a.cov.70.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.70)

clust.AR.a.cov.75 <- DF.r[,1318:1335]

clust.AR.a.cov.75.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.75)

clust.AR.a.cov.80 <- DF.r[,1336:1353]

clust.AR.a.cov.80.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.80)

clust.AR.a.cov.85 <- DF.r[,1354:1371]

clust.AR.a.cov.85.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.85)

clust.AR.a.cov.90 <- DF.r[,1372:1389]

clust.AR.a.cov.90.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.90)

clust.AR.a.cov.95 <- DF.r[,1390:1407]

clust.AR.a.cov.95.tab <- stat.fun.men.women.ratio(df=clust.AR.a.cov.95)


# b

clust.AR.b.cov.35 <- DF.r[,1408:1425] # +17

clust.AR.b.cov.35.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.35)

clust.AR.b.cov.40 <- DF.r[,1426:1443]

clust.AR.b.cov.40.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.40)

clust.AR.b.cov.45 <- DF.r[,1444:1461]

clust.AR.b.cov.45.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.45)

clust.AR.b.cov.50 <- DF.r[,1462:1479]

clust.AR.b.cov.50.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.50)

clust.AR.b.cov.55 <- DF.r[,1480:1497]

clust.AR.b.cov.55.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.55)

clust.AR.b.cov.60 <- DF.r[,1498:1515]

clust.AR.b.cov.60.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.60)

clust.AR.b.cov.65 <- DF.r[,1516:1533]

clust.AR.b.cov.65.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.65)

clust.AR.b.cov.70 <- DF.r[,1534:1551]

clust.AR.b.cov.70.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.70)

clust.AR.b.cov.75 <- DF.r[,1552:1569]

clust.AR.b.cov.75.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.75)

clust.AR.b.cov.80 <- DF.r[,1570:1587]

clust.AR.b.cov.80.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.80)

clust.AR.b.cov.85 <- DF.r[,1588:1605]

clust.AR.b.cov.85.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.85)

clust.AR.b.cov.90 <- DF.r[,1606:1623]

clust.AR.b.cov.90.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.90)

clust.AR.b.cov.95 <- DF.r[,1624:1641]

clust.AR.b.cov.95.tab <- stat.fun.men.women.ratio(df=clust.AR.b.cov.95)


# c

clust.AR.c.cov.35 <- DF.r[,1642:1659] # +17

clust.AR.c.cov.35.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.35)

clust.AR.c.cov.40 <- DF.r[,1660:1677]

clust.AR.c.cov.40.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.40)

clust.AR.c.cov.45 <- DF.r[,1678:1695]

clust.AR.c.cov.45.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.45)

clust.AR.c.cov.50 <- DF.r[,1696:1713]

clust.AR.c.cov.50.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.50)

clust.AR.c.cov.55 <- DF.r[,1714:1731]

clust.AR.c.cov.55.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.55)

clust.AR.c.cov.60 <- DF.r[,1732:1749]

clust.AR.c.cov.60.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.60)

clust.AR.c.cov.65 <- DF.r[,1750:1767]

clust.AR.c.cov.65.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.65)

clust.AR.c.cov.70 <- DF.r[,1768:1785]

clust.AR.c.cov.70.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.70)

clust.AR.c.cov.75 <- DF.r[,1786:1803]

clust.AR.c.cov.75.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.75)

clust.AR.c.cov.80 <- DF.r[,1804:1821]

clust.AR.c.cov.80.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.80)

clust.AR.c.cov.85 <- DF.r[,1822:1839]

clust.AR.c.cov.85.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.85)

clust.AR.c.cov.90 <- DF.r[,1840:1857]

clust.AR.c.cov.90.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.90)

clust.AR.c.cov.95 <- DF.r[,1858:1875]

clust.AR.c.cov.95.tab <- stat.fun.men.women.ratio(df=clust.AR.c.cov.95)



# Plotting the input parameters of the calibrated model

# > names(as.data.frame(clust.AR.c.cov.95.tab))
# [1] "men.15.25.w.15.25" "men.15.25.w.25.40" "men.15.25.w.40.50" "men.25.40.w.15.25" "men.25.40.w.25.40" "men.25.40.w.40.50" 
# "men.40.50.w.15.25" "men.40.50.w.25.40" "men.40.50.w.40.50"

# Compare two-by-two paiirngs of true transmission and clusters
# here: [4] men.25.40.w.15.25 & [5] men.25.40.w.25.40

# Ex.1
plot(transm.AR.a.cov.95.tab[,4],
     transm.AR.a.cov.95.tab[,5],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1),
     ylim = c(0, 1))
points(clust.AR.b.cov.95.tab[,4],
       clust.AR.b.cov.95.tab[,5],
       pch = 16,
       col = "blue2")
points(clust.AR.c.cov.95.tab[,4],
       clust.AR.c.cov.95.tab[,5],
       pch = 16,
       col = "orange")

# Ex.2

plot(transm.CAR.cov.95.tab[,4],
     transm.CAR.cov.95.tab[,5],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1),
     ylim = c(0, 1))
points(clust.CAR.cov.95.tab[,4],
       clust.CAR.cov.95.tab[,5],
       pch = 16,
       col = "blue2")
points(clust.AR.a.cov.95.tab[,4],
       clust.AR.a.cov.95.tab[,5],
       pch = 16,
       col = "red")
points(clust.AR.b.cov.95.tab[,4],
       clust.AR.b.cov.95.tab[,5],
       pch = 16,
       col = "orange")
points(clust.AR.c.cov.95.tab[,4],
       clust.AR.c.cov.95.tab[,5],
       pch = 16,
       col = "cyan")





# Fix true transmission of one pairing of age-group and compare with different clusters
# here: [4] men.25.40.w.15.25

plot(transm.AR.a.cov.95.tab[,4],
     clust.AR.a.cov.95.tab[,4],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1),
     ylim = c(0, 1))
points(transm.AR.a.cov.95.tab[,4],
       clust.AR.b.cov.95.tab[,4],
       pch = 16,
       col = "blue2")
points(transm.AR.a.cov.95.tab[,4],
       clust.AR.c.cov.95.tab[,4],
       pch = 16,
       col = "orange")




################################

library(data.table)
library(lattice)
library(mice)
library(varhandle)
library(ggplot2)
library(caret)
library(nnet)
library(neuralnet)
library(MASS)
library(survival)
library(fitdistrplus)
library(ggcorrplot)

# dataset.raw.num.complete <- as.data.frame(cbind(transm.CAR.cov.95.tab[,4],
#                                                 transm.CAR.cov.95.tab[,5],
#                                                 clust.CAR.cov.95.tab[,4],
#                                                 clust.CAR.cov.95.tab[,5],
#                                                 clust.AR.a.cov.95.tab[,4],
#                                                 clust.AR.a.cov.95.tab[,5],
#                                                 clust.AR.b.cov.95.tab[,4],
#                                                 clust.AR.b.cov.95.tab[,5],
#                                                 clust.AR.c.cov.95.tab[,4],
#                                                 clust.AR.c.cov.95.tab[,5]))

dataset.raw.num.complete <- as.data.frame(cbind(transm.CAR.cov.95.tab[,4],
                                                clust.CAR.cov.95.tab[,4],
                                                clust.AR.a.cov.95.tab[,4],
                                                clust.AR.b.cov.95.tab[,4],
                                                clust.AR.c.cov.95.tab[,4]))

# names.v <- c("transm.CAR.cov1",
#              "transm.CAR.cov2",
#              "clust.CAR.cov1",
#              "clust.CAR.cov2",
#              "clust.AR.a.cov1",
#              "clust.AR.a.cov2",
#              "clust.AR.b.cov1",
#              "clust.AR.b.cov2",
#              "clust.AR.c.cov1",
#              "clust.AR.c.cov2")

names.v <- c("transm.CAR.cov1",
             "clust.CAR.cov1",
             "clust.AR.a.cov1",
             "clust.AR.b.cov1",
             "clust.AR.c.cov1")


names(dataset.raw.num.complete) <- names.v


dataset.raw.num.complete <- na.omit(dataset.raw.num.complete)

# Compute a correlation matrix
corr <- round(cor(dataset.raw.num.complete), 6)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(dataset.raw.num.complete)


# Visualize the correlation matrix
# --------------------------------
# method = "square" (default)
ggcorrplot(corr)



# Add correlation coefficients
# --------------------------------
# argument lab = TRUE
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)


# Add correlation significance level
# --------------------------------
# Argument p.mat
# Barring the no significant coefficient
ggcorrplot(corr, hc.order = TRUE,
           type = "lower", p.mat = p.mat)


# Reordering the correlation matrix
# --------------------------------
# using hierarchical clustering
ggcorrplot(corr, hc.order = TRUE, outline.col = "white")

