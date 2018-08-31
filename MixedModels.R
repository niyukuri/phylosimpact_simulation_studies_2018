

# Linear models    and    linear    mixed    effects    models    in    R    with    linguistic    applications    

# Pitch and sex

pitch = c(233,204,242,130,112,142)
sex = c(rep("female",3),rep("male",3))
my.df = data.frame(sex,pitch)

xmdl = lm(pitch ~ sex, my.df)


summary(xmdl)


age = c(14,23,35,48,52,67)
pitch = c(252,244,240,233,212,204)
my.df = data.frame(age,pitch)
xmdl = lm(pitch ~ age, my.df)
summary(xmdl)


my.df$age.c = my.df$age - mean(my.df$age)
xmdl = lm(pitch ~ age.c, my.df)
summary(xmdl)

# 2nd



library(nlme)

data(Rail)

summary(Rail)

with(Rail, tapply(travel, Rail, mean))


with(Rail, plot(travel, Rail, xlab="Travel time"))


# Rail as a fixed effect
########################

# y_ij = β_i + e_ij

r1.lm <- lm(travel ~ Rail - 1, data=Rail)

# β 2 = 31.67 , . . . , β 4 = 96.00

summary(r1.lm)

library(lattice)

with(Rail, bwplot(Rail ~ residuals(r1.lm)))


# Rail as a random effect
#########################

# y_ij = β + b_i + e_ij

library(lme4)

r2.lme <- lmer(travel ~ 1 + (1 | Rail),
               REML=FALSE, data=Rail)


summary(r2.lme)


fixef(r2.lme)+ranef(r2.lme)$Rail


# Multilevel models
###################

library(faraway)
data(jsp)

summary(jsp)


# Center raven since we would otherwise be comparing to zero

jsp.y2$ctrraven<-jsp.y2$raven-mean(jsp.y2$raven)



data(weight)


# Tuesday, 14th August, 2018 #
##############################

# URL: https://ourcodingclub.github.io/2017/03/15/mixed-models.html


######################################
#                                    #
#   Mixed effects modeling in R     #
#                                    #
######################################

## authors: Gabriela K Hajduk, based on workshop developed by Liam Bailey
## contact details: gkhajduk.github.io; email: gkhajduk@gmail.com
## date: 2017-03-09
##

###---- Explore the data -----###

## load the data and have a look at it

load("~/CC-Linear-mixed-models/dragons.RData")

head(dragons)

## Let's say we want to know how the body length affects test scores.

## Have a look at the data distribution:

hist(dragons$testScore)  # seems close to normal distribution - good!

## It is good practice to  standardise your explanatory variables before proceeding - you can use scale() to do that:

dragons$bodyLength2 <- scale(dragons$bodyLength)

## Back to our question: is test score affected by body length?

###---- Fit all data in one analysis -----###

## One way to analyse this data would be to try fitting a linear model to all our data, ignoring the sites and the mountain ranges for now.

library(lme4)

basic.lm <- lm(testScore ~ bodyLength2, data = dragons)

summary(basic.lm)

## Let's plot the data with ggplot2

library(ggplot2)

ggplot(dragons, aes(x = bodyLength, y = testScore)) +
  geom_point()+
  geom_smooth(method = "lm")


### Assumptions?

## Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(basic.lm, which = 1)  # not perfect, but look alright

## Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(basic.lm, which = 2)  # a bit off at the extremes, but that's often the case; again doesn't look too bad


## However, what about observation independence? Are our data independent?
## We collected multiple samples from eight mountain ranges
## It's perfectly plausible that the data from within each mountain range are more similar to each other than the data from different mountain ranges - they are correlated. Pseudoreplication isn't our friend.

## Have a look at the data to see if above is true
boxplot(testScore ~ mountainRange, data = dragons)  # certainly looks like something is going on here

## We could also plot it colouring points by mountain range
ggplot(dragons, aes(x = bodyLength, y = testScore, colour = mountainRange))+
  geom_point(size = 2)+
  theme_classic()+
  theme(legend.position = "none")

## From the above plots it looks like our mountain ranges vary both in the dragon body length and in their test scores. This confirms that our observations from within each of the ranges aren't independent. We can't ignore that.

## So what do we do?

###----- Run multiple analyses -----###


## We could run many separate analyses and fit a regression for each of the mountain ranges.

## Lets have a quick look at the data split by mountain range
## We use the facet_wrap to do that

ggplot(aes(bodyLength, testScore), data = dragons) + geom_point() +
  facet_wrap(~ mountainRange) +
  xlab("length") + ylab("test score")



##----- Modify the model -----###

## We want to use all the data, but account for the data coming from different mountain ranges

## let's add mountain range as a fixed effect to our basic.lm

mountain.lm <- lm(testScore ~ bodyLength2 + mountainRange, data = dragons)
summary(mountain.lm)

## now body length is not significant


###----- Mixed effects models -----###

library(lme4)

##----- First mixed model -----##

### model

### plots

### summary

### variance accounted for by mountain ranges



##-- implicit vs explicit nesting --##

head(dragons)  # we have site and mountainRange
str(dragons)  # we took samples from three sites per mountain range and eight mountain ranges in total

### create new "sample" variable


##----- Second mixed model -----##

### model

### summary

### plot



##----- Model selection for the keen -----##

### full model

### reduced model

### comparison




# Dummy data for age-mixing in transmission clusters

library(data.table)
library(nlme)
library(ggplot2)


female.age <- c(17:24, 32, 23, 21, 29, 19, 21, 23, 24, 18, 19, 19, 21, 17:24, 21:27, 19)
male.age <- c(19, 21, 23, 24, 18, 19, 19, 21, 17:24, 21:27)


female.df <- data.frame(age=female.age)
male.df <- data.frame(age=male.age)

female.df$gender <- "female"
male.df$gender <- "male"

ageTransmClust <- rbind(female.df, male.df)

ggplot(ageTransmClust, aes(age, fill = gender)) + geom_density(alpha = 0.2)



# 1
female.age1 <- c(17, 32, 23, 21, 29, 19, 21, 23, 24, 18, 19, 19, 21, 17:24, 21:27, 19)
male.age1 <- c(42, 40, 35, 49, 50, 37, 45, 49, 39, 41)


female.df1 <- data.frame(age=female.age1)
male.df1 <- data.frame(age=male.age1)

female.df1$gender <- "female"
male.df1$gender <- "male"

ageTransmClust1 <- rbind(female.df1, male.df1)

ggplot(ageTransmClust1, aes(age, fill = gender)) + geom_density(alpha = 0.2)

# 2

female.age2 <- c(17, 32, 23, 21, 22, 19, 21, 23, 24, 18, 19, 19, 21, 17:24,19)
male.age2 <- c(42, 21, 41, 35, 39, 50, 19, 27, 45, 39, 39, 21)


female.df2 <- data.frame(age=female.age2)
male.df2 <- data.frame(age=male.age2)

female.df2$gender <- "female"
male.df2$gender <- "male"

ageTransmClust2 <- rbind(female.df2, male.df2)

ggplot(ageTransmClust2, aes(age, fill = gender)) + geom_density(alpha = 0.2)

# 3

female.age3 <- c(17:34, 32, 23, 21, 22, 19, 21:27, 23, 24, 18, 19, 19, 21, 17:24,19)
male.age3 <- c(42, 21:42, 41, 35, 39, 50, 17:29, 27, 45, 39, 39, 21)


female.df3 <- data.frame(age=female.age3)
male.df3 <- data.frame(age=male.age3)

female.df3$gender <- "female"
male.df3$gender <- "male"

ageTransmClust3 <- rbind(female.df3, male.df3)


ggplot(ageTransmClust3, aes(age, fill = gender)) + geom_density(alpha = 0.2)


# 4

female.age4 <- c(17:34, 32, 23, 21, 22, 23, 24, 18, 19, 19, 21, 17:24,19)
male.age4 <- c(42, 21:42, 41, 35, 39, 50, 27, 45, 39, 39, 21)


female.df4 <- data.frame(age=female.age4)
male.df4 <- data.frame(age=male.age4)

female.df4$gender <- "female"
male.df4$gender <- "male"

ageTransmClust4 <- rbind(female.df4, male.df4)


ggplot(ageTransmClust4, aes(age, fill = gender)) + geom_density(alpha = 0.2)


# 5

female.age5 <- c(17:34, 32, 23, 21, 22, 23:40, 24, 18, 19, 19, 21, 17:24,19)
male.age5 <- c(42, 21:42, 41, 35, 39, 19, 50, 27, 45:50, 39, 39:41, 21)


female.df5 <- data.frame(age=female.age5)
male.df5 <- data.frame(age=male.age5)

female.df5$gender <- "female"
male.df5$gender <- "male"

ageTransmClust5 <- rbind(female.df5, male.df5)


ggplot(ageTransmClust5, aes(age, fill = gender)) + geom_density(alpha = 0.2)


# 6

female.age6 <- c(17:34, 32, 23, 21, 22,  24, 18, 19, 19, 21, 27:44,19)
male.age6 <- c(42, 41, 35, 39, 19, 50, 27, 45:50, 39, 39:41, 21)


female.df6 <- data.frame(age=female.age6)
male.df6 <- data.frame(age=male.age6)

female.df6$gender <- "female"
male.df6$gender <- "male"

ageTransmClust6 <- rbind(female.df6, male.df6)


ggplot(ageTransmClust6, aes(age, fill = gender)) + geom_density(alpha = 0.2)


# 7

female.age7 <- c(32, 23, 21, 22, 40:49,  24, 18)
male.age7 <- c(42, 41, 35, 39, 19, 50, 27, 45:50, 39, 39:41)


female.df7 <- data.frame(age=female.age7)
male.df7 <- data.frame(age=male.age7)

female.df7$gender <- "female"
male.df7$gender <- "male"

ageTransmClust7 <- rbind(female.df7, male.df7)

ggplot(ageTransmClust7, aes(age, fill = gender)) + geom_density(alpha = 0.2)


ageTransmClusts <- list(ageTransmClust1, ageTransmClust2, ageTransmClust3,
                        ageTransmClust4, ageTransmClust5, ageTransmClust6,
                        ageTransmClust7)



library(nlme)

ageTransmClust1$id <- 1
ageTransmClust2$id <- 2
ageTransmClust3$id <- 3
ageTransmClust4$id <- 4
ageTransmClust5$id <- 5
ageTransmClust6$id <- 6
ageTransmClust7$id <- 7

ageTransmClusts.Clean <- list(ageTransmClust1, ageTransmClust2, ageTransmClust3,
                              ageTransmClust4, ageTransmClust5, ageTransmClust6,
                              ageTransmClust7)

AGEMIX <- do.call(rbind, ageTransmClusts.Clean)

AGEMIX$id <- as.factor(AGEMIX$id)

traclust <- lme(age ~ gender, data = AGEMIX, random = ~ 1|id)

library(lme4)

traclust2 <- lmer(age ~ gender +  (1|id), data = AGEMIX)

# Error in match.arg(name) : 
#   'arg' should be one of “X”, “Z”, “Zt”, “Ztlist”, “mmList”, “y”, “mu”, “u”, “b”, 
# “Gp”, “Tp”, “L”, “Lambda”, “Lambdat”, “Lind”, “Tlist”, “A”, “RX”, “RZX”, “sigma”, 
# “flist”, “fixef”, “beta”, “theta”, “ST”, “REML”, “is_REML”, “n_rtrms”, “n_rfacs”, 
# “N”, “n”, “p”, “q”, “p_i”, “l_i”, “q_i”, “k”, “m_i”, “m”, “cnms”, “devcomp”, “offset”, 
# “lower”, “devfun”, “glmer.nb.theta”

random.effects(traclust)

fitted(traclust, level = 0:1)

coef(summary(traclust))[2] # 

ranef(traclust)[[1]][[1]]

a <- coef(summary(traclust))[1] # intercept - alpha

beta <- coef(summary(traclust))[2] # or fixed.effects(summary(traclust))[2] gender effect - gamma
# getME(traclust2,"beta")[2] < this is bridge.width 

# b <- getME(traclust2,"sigma") # between cluster var

b1 <- as.numeric(VarCorr(traclust)[3]) # between cluster variation - sigma_beta

b2 <- as.numeric(VarCorr(traclust)[4]) # within cluster variation - residual - sigma



############################## Clusters


transmission.clust.list <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster


# Check how many indiv in a cluster

# for (i in 1:length(d)) {
#   
#   clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
#   
#   t.v <- c(t.v, nrow(clus.read))
#   
# }

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  #size <- c(size, nrow(clus.read))
  
  transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
  
  transm.df.cl.dat$age <- transm.df.cl$age.samp.Rec
  transm.df.cl.dat$gender <- transm.df.cl$GenderRec
  transm.df.cl.dat$clust.id <- as.factor(rep(i, nrow(transm.df.cl)))
  
  transmission.clust.list[[i]] <- as.data.frame(transm.df.cl.dat)
  
}


clust.table.df <- as.data.frame(do.call(rbind, transmission.clust.list)) # data.table & data.frame


fit.lme.transm.clust <- lme(age ~ gender, data = clust.table.df, random = ~ 1|clust.id)



a <- coef(summary(fit.lme.transm.clust))[1] # average age in transmission clusters

beta <- coef(summary(fit.lme.transm.clust))[2] # average age difference in transmission clusters: 
# seen as bridge width which shows potential cross-generation transmission


b1 <- as.numeric(VarCorr(fit.lme.transm.clust)[3]) # between cluster variation

b2 <- as.numeric(VarCorr(fit.lme.transm.clust)[4]) # within cluster variation



toy.parallel <- function(input){
  set.seed(input[1])
  med.toy <- median(rnorm(n=input[2], mean = input[3], sd = 1))
  val.toy <- med.toy + 6
  val.ALL <- c(med.toy, val.toy)
  names(val.ALL) <- c("med.toy", "val.toy")
  return(val.ALL)
}

input=c(777, 100,12)

v <- toy.parallel(input = input)

input <- c(100,12)
reps <- 20
inputmatrix <- matrix(rep(input, reps), byrow = TRUE, nrow = reps)

library(parallel)

seed_count = 0
n_cluster = 8

cl <- makeCluster(getOption("cl.cores", n_cluster))
tab_simul_summarystat = NULL
list_param <- list(NULL)
tab_param <- NULL
paramtemp <- NULL
simultemp <- NULL

nb_simul <- nrow(inputmatrix)

for (i in 1:nb_simul) {
  l <- ncol(inputmatrix)
  param <- c((seed_count + i), inputmatrix[i, ])
  list_param[[i]] <- param
  tab_param <- rbind(tab_param, param[2:(l + 1)])
  paramtemp <- rbind(paramtemp, param[2:(l + 1)])
}
list_simul_summarystat = parLapplyLB(cl, list_param,
                                     toy.parallel)
tab_simul_summarystat <- do.call(rbind, list_simul_summarystat)
stopCluster(cl)

results <- cbind(tab_simul_summarystat, seed_count + 1:nb_simul)

write.csv(results, file = "toy.results.csv")


## another try: foreach

library(foreach)
library(doParallel)

cl<-makeCluster(n_cores)
registerDoParallel(cl)


n <- n_cores

res <- foreach(icount(n), .combine=rbind) %dopar% {
  
  v <- list_param[[i]]
  
  toy.parallel(input = v)
  
}

stopImplicitCluster()

d = as.data.frame(res)


res.z <- foreach(icount(n), .combine=cbind) %dopar% {
  
  
  fun.raba(10)
  
  
}



### Out together simulated data for age-mix in transmission

library(data.table)
library(dplyr)

# run1 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run1.csv")
# run2 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run2.csv")
# run3 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run3.csv")
# run4 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run4.csv")
# run6 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run6.csv")
# 
# 
# run10 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run10.csv")
# run12 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run12.csv")
# run15 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run15.csv")
# 
# 
# run16 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run16.csv")
# run17 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run17.csv")
# run18 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run18.csv")
# 
# run19 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run19.csv")
# run20 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run20.csv")


run3 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run3.csv")
run4 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run4.csv")
run6 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run6.csv")


run10 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run10.csv")
run12 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run12.csv")
run15 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run15.csv")


run16 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run16.csv")
run17 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run17.csv")
run18 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run18.csv")

run19 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run19.csv")
run20 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.age.mixing.run20.csv")


run.ALL <- rbind(run3, run4, run6,
                 run10, run12, run15, run16,
                 run17, run18, run19, run20)

run.ALL <- dplyr::filter(run.ALL, run.ALL$flag.lme=="1") # select simulations with age mixing

# run.ALL <- abs(run.ALL) # to make sure all values are positive - age gap

# True: population level mixing characteristics

Pop.mean.AD.v <- run.ALL$Pop.mean.AD
Pop.med.AD.v <- run.ALL$Pop.med.AD
Pop.sd.AD.v <- run.ALL$Pop.sd.AD
av.age.male.v <- run.ALL$av.age.male
av.age.diff.v <- run.ALL$av.age.diff
between.transm.var.v <- run.ALL$between.transm.var
within.transm.var.v <- run.ALL$within.transm.var

# "clust.AR.c.95.av.age.male"      
# [217] "clust.AR.c.95.av.age.diff"       "clust.AR.c.95.between.clust.var" "clust.AR.c.95.within.clust.var
# 

# MCAR

# Cov 35
MCAR.av.age.male.cov.35 <- run.ALL[,12]
MCAR.av.age.diff.cov.35 <- run.ALL[,13]
MCAR.between.transm.var.cov.35 <- run.ALL[,14]
MCAR.within.transm.var.cov.35 <- run.ALL[,15]


# Cov 40
MCAR.av.age.male.cov.40 <- run.ALL[,16]
MCAR.av.age.diff.cov.40 <- run.ALL[,17]
MCAR.between.transm.var.cov.40 <- run.ALL[,18]
MCAR.within.transm.var.cov.40 <- run.ALL[,19]


# Cov 45
MCAR.av.age.male.cov.45 <- run.ALL[,20]
MCAR.av.age.diff.cov.45 <- run.ALL[,21]
MCAR.between.transm.var.cov.45 <- run.ALL[,22]
MCAR.within.transm.var.cov.45 <- run.ALL[,23]


# Cov 50
MCAR.av.age.male.cov.50 <- run.ALL[,24]
MCAR.av.age.diff.cov.50 <- run.ALL[,25]
MCAR.between.transm.var.cov.50 <- run.ALL[,26]
MCAR.within.transm.var.cov.50 <- run.ALL[,27]


# Cov 55
MCAR.av.age.male.cov.55 <- run.ALL[,28]
MCAR.av.age.diff.cov.55 <- run.ALL[,29]
MCAR.between.transm.var.cov.55 <- run.ALL[,30]
MCAR.within.transm.var.cov.55 <- run.ALL[,31]

# Cov 60
MCAR.av.age.male.cov.60 <- run.ALL[,32]
MCAR.av.age.diff.cov.60 <- run.ALL[,33]
MCAR.between.transm.var.cov.60 <- run.ALL[,34]
MCAR.within.transm.var.cov.60 <- run.ALL[,35]


# Cov 65
MCAR.av.age.male.cov.65 <- run.ALL[,36]
MCAR.av.age.diff.cov.65 <- run.ALL[,37]
MCAR.between.transm.var.cov.65 <- run.ALL[,38]
MCAR.within.transm.var.cov.65 <- run.ALL[,39]


# Cov 70
MCAR.av.age.male.cov.70 <- run.ALL[,40]
MCAR.av.age.diff.cov.70 <- run.ALL[,41]
MCAR.between.transm.var.cov.70 <- run.ALL[,42]
MCAR.within.transm.var.cov.70 <- run.ALL[,43]


# Cov 75
MCAR.av.age.male.cov.75 <- run.ALL[,44]
MCAR.av.age.diff.cov.75 <- run.ALL[,45]
MCAR.between.transm.var.cov.75 <- run.ALL[,46]
MCAR.within.transm.var.cov.75 <- run.ALL[,47]

# Cov 80
MCAR.av.age.male.cov.80 <- run.ALL[,48]
MCAR.av.age.diff.cov.80 <- run.ALL[,49]
MCAR.between.transm.var.cov.80 <- run.ALL[,50]
MCAR.within.transm.var.cov.80 <- run.ALL[,51]


# Cov 85
MCAR.av.age.male.cov.85 <- run.ALL[,52]
MCAR.av.age.diff.cov.85 <- run.ALL[,53]
MCAR.between.transm.var.cov.85 <- run.ALL[,54]
MCAR.within.transm.var.cov.85 <- run.ALL[,55]

# Cov 90
MCAR.av.age.male.cov.90 <- run.ALL[,56]
MCAR.av.age.diff.cov.90 <- run.ALL[,57]
MCAR.between.transm.var.cov.90 <- run.ALL[,58]
MCAR.within.transm.var.cov.90 <- run.ALL[,59]


# Cov 95
MCAR.av.age.male.cov.95 <- run.ALL[,60]
MCAR.av.age.diff.cov.95 <- run.ALL[,61]
MCAR.between.transm.var.cov.95 <- run.ALL[,62]
MCAR.within.transm.var.cov.95 <- run.ALL[,63]


# MAR

# a.

# Cov 35
MAR.a.av.age.male.cov.35 <- run.ALL[,64]
MAR.a.av.age.diff.cov.35 <- run.ALL[,65]
MAR.a.between.transm.var.cov.35 <- run.ALL[,66]
MAR.a.within.transm.var.cov.35 <- run.ALL[,67]


# Cov 40
MAR.a.av.age.male.cov.40 <- run.ALL[,68]
MAR.a.av.age.diff.cov.40 <- run.ALL[,69]
MAR.a.between.transm.var.cov.40 <- run.ALL[,70]
MAR.a.within.transm.var.cov.40 <- run.ALL[,71]


# Cov 45
MAR.a.av.age.male.cov.45 <- run.ALL[,72]
MAR.a.av.age.diff.cov.45 <- run.ALL[,73]
MAR.a.between.transm.var.cov.45 <- run.ALL[,74]
MAR.a.within.transm.var.cov.45 <- run.ALL[,75]


# Cov 50
MAR.a.av.age.male.cov.50 <- run.ALL[,76]
MAR.a.av.age.diff.cov.50 <- run.ALL[,77]
MAR.a.between.transm.var.cov.50 <- run.ALL[,78]
MAR.a.within.transm.var.cov.50 <- run.ALL[,79]


# Cov 55
MAR.a.av.age.male.cov.55 <- run.ALL[,80]
MAR.a.av.age.diff.cov.55 <- run.ALL[,81]
MAR.a.between.transm.var.cov.55 <- run.ALL[,82]
MAR.a.within.transm.var.cov.55 <- run.ALL[,83]

# Cov 60
MAR.a.av.age.male.cov.60 <- run.ALL[,84]
MAR.a.av.age.diff.cov.60 <- run.ALL[,85]
MAR.a.between.transm.var.cov.60 <- run.ALL[,86]
MAR.a.within.transm.var.cov.60 <- run.ALL[,87]


# Cov 65
MAR.a.av.age.male.cov.65 <- run.ALL[,88]
MAR.a.av.age.diff.cov.65 <- run.ALL[,89]
MAR.a.between.transm.var.cov.65 <- run.ALL[,90]
MAR.a.within.transm.var.cov.65 <- run.ALL[,91]


# Cov 70
MAR.a.av.age.male.cov.70 <- run.ALL[,92]
MAR.a.av.age.diff.cov.70 <- run.ALL[,93]
MAR.a.between.transm.var.cov.70 <- run.ALL[,94]
MAR.a.within.transm.var.cov.70 <- run.ALL[,95]


# Cov 75
MAR.a.av.age.male.cov.75 <- run.ALL[,96]
MAR.a.av.age.diff.cov.75 <- run.ALL[,97]
MAR.a.between.transm.var.cov.75 <- run.ALL[,98]
MAR.a.within.transm.var.cov.75 <- run.ALL[,99]

# Cov 80
MAR.a.av.age.male.cov.80 <- run.ALL[,100]
MAR.a.av.age.diff.cov.80 <- run.ALL[,101]
MAR.a.between.transm.var.cov.80 <- run.ALL[,102]
MAR.a.within.transm.var.cov.80 <- run.ALL[,103]


# Cov 85
MAR.a.av.age.male.cov.85 <- run.ALL[,104]
MAR.a.av.age.diff.cov.85 <- run.ALL[,105]
MAR.a.between.transm.var.cov.85 <- run.ALL[,106]
MAR.a.within.transm.var.cov.85 <- run.ALL[,107]

# Cov 90
MAR.a.av.age.male.cov.90 <- run.ALL[,108]
MAR.a.av.age.diff.cov.90 <- run.ALL[,109]
MAR.a.between.transm.var.cov.90 <- run.ALL[,110]
MAR.a.within.transm.var.cov.90 <- run.ALL[,111]


# Cov 95
MAR.a.av.age.male.cov.95 <- run.ALL[,112]
MAR.a.av.age.diff.cov.95 <- run.ALL[,113]
MAR.a.between.transm.var.cov.95 <- run.ALL[,114]
MAR.a.within.transm.var.cov.95 <- run.ALL[,115]


# b

# Cov 35
MAR.b.av.age.male.cov.35 <- run.ALL[,116]
MAR.b.av.age.diff.cov.35 <- run.ALL[,117]
MAR.b.between.transm.var.cov.35 <- run.ALL[,118]
MAR.b.within.transm.var.cov.35 <- run.ALL[,119]


# Cov 40
MAR.b.av.age.male.cov.40 <- run.ALL[,120]
MAR.b.av.age.diff.cov.40 <- run.ALL[,121]
MAR.b.between.transm.var.cov.40 <- run.ALL[,122]
MAR.b.within.transm.var.cov.40 <- run.ALL[,123]


# Cov 45
MAR.b.av.age.male.cov.45 <- run.ALL[,124]
MAR.b.av.age.diff.cov.45 <- run.ALL[,125]
MAR.b.between.transm.var.cov.45 <- run.ALL[,126]
MAR.b.within.transm.var.cov.45 <- run.ALL[,127]


# Cov 50
MAR.b.av.age.male.cov.50 <- run.ALL[,128]
MAR.b.av.age.diff.cov.50 <- run.ALL[,129]
MAR.b.between.transm.var.cov.50 <- run.ALL[,130]
MAR.b.within.transm.var.cov.50 <- run.ALL[,131]


# Cov 55
MAR.b.av.age.male.cov.55 <- run.ALL[,132]
MAR.b.av.age.diff.cov.55 <- run.ALL[,133]
MAR.b.between.transm.var.cov.55 <- run.ALL[,134]
MAR.b.within.transm.var.cov.55 <- run.ALL[,135]

# Cov 60
MAR.b.av.age.male.cov.60 <- run.ALL[,136]
MAR.b.av.age.diff.cov.60 <- run.ALL[,137]
MAR.b.between.transm.var.cov.60 <- run.ALL[,138]
MAR.b.within.transm.var.cov.60 <- run.ALL[,139]


# Cov 65
MAR.b.av.age.male.cov.65 <- run.ALL[,140]
MAR.b.av.age.diff.cov.65 <- run.ALL[,141]
MAR.b.between.transm.var.cov.65 <- run.ALL[,142]
MAR.b.within.transm.var.cov.65 <- run.ALL[,143]


# Cov 70
MAR.b.av.age.male.cov.70 <- run.ALL[,144]
MAR.b.av.age.diff.cov.70 <- run.ALL[,145]
MAR.b.between.transm.var.cov.70 <- run.ALL[,146]
MAR.b.within.transm.var.cov.70 <- run.ALL[,147]


# Cov 75
MAR.b.av.age.male.cov.75 <- run.ALL[,148]
MAR.b.av.age.diff.cov.75 <- run.ALL[,149]
MAR.b.between.transm.var.cov.75 <- run.ALL[,150]
MAR.b.within.transm.var.cov.75 <- run.ALL[,151]

# Cov 80
MAR.b.av.age.male.cov.80 <- run.ALL[,152]
MAR.b.av.age.diff.cov.80 <- run.ALL[,153]
MAR.b.between.transm.var.cov.80 <- run.ALL[,154]
MAR.b.within.transm.var.cov.80 <- run.ALL[,155]


# Cov 85
MAR.b.av.age.male.cov.85 <- run.ALL[,156]
MAR.b.av.age.diff.cov.85 <- run.ALL[,157]
MAR.b.between.transm.var.cov.85 <- run.ALL[,158]
MAR.b.within.transm.var.cov.85 <- run.ALL[,159]

# Cov 90
MAR.b.av.age.male.cov.90 <- run.ALL[,160]
MAR.b.av.age.diff.cov.90 <- run.ALL[,161]
MAR.b.between.transm.var.cov.90 <- run.ALL[,162]
MAR.b.within.transm.var.cov.90 <- run.ALL[,163]


# Cov 95
MAR.b.av.age.male.cov.95 <- run.ALL[,164]
MAR.b.av.age.diff.cov.95 <- run.ALL[,165]
MAR.b.between.transm.var.cov.95 <- run.ALL[,166]
MAR.b.within.transm.var.cov.95 <- run.ALL[,167]


# c

# Cov 35
MAR.c.av.age.male.cov.35 <- run.ALL[,168]
MAR.c.av.age.diff.cov.35 <- run.ALL[,169]
MAR.c.between.transm.var.cov.35 <- run.ALL[,170]
MAR.c.within.transm.var.cov.35 <- run.ALL[,171]


# Cov 40
MAR.c.av.age.male.cov.40 <- run.ALL[,172]
MAR.c.av.age.diff.cov.40 <- run.ALL[,173]
MAR.c.between.transm.var.cov.40 <- run.ALL[,174]
MAR.c.within.transm.var.cov.40 <- run.ALL[,175]


# Cov 45
MAR.c.av.age.male.cov.45 <- run.ALL[,176]
MAR.c.av.age.diff.cov.45 <- run.ALL[,177]
MAR.c.between.transm.var.cov.45 <- run.ALL[,178]
MAR.c.within.transm.var.cov.45 <- run.ALL[,179]


# Cov 50
MAR.c.av.age.male.cov.50 <- run.ALL[,180]
MAR.c.av.age.diff.cov.50 <- run.ALL[,181]
MAR.c.between.transm.var.cov.50 <- run.ALL[,182]
MAR.c.within.transm.var.cov.50 <- run.ALL[,183]


# Cov 55
MAR.c.av.age.male.cov.55 <- run.ALL[,184]
MAR.c.av.age.diff.cov.55 <- run.ALL[,185]
MAR.c.between.transm.var.cov.55 <- run.ALL[,186]
MAR.c.within.transm.var.cov.55 <- run.ALL[,187]

# Cov 60
MAR.c.av.age.male.cov.60 <- run.ALL[,188]
MAR.c.av.age.diff.cov.60 <- run.ALL[,189]
MAR.c.between.transm.var.cov.60 <- run.ALL[,190]
MAR.c.within.transm.var.cov.60 <- run.ALL[,191]


# Cov 65
MAR.c.av.age.male.cov.65 <- run.ALL[,192]
MAR.c.av.age.diff.cov.65 <- run.ALL[,193]
MAR.c.between.transm.var.cov.65 <- run.ALL[,194]
MAR.c.within.transm.var.cov.65 <- run.ALL[,195]


# Cov 70
MAR.c.av.age.male.cov.70 <- run.ALL[,196]
MAR.c.av.age.diff.cov.70 <- run.ALL[,197]
MAR.c.between.transm.var.cov.70 <- run.ALL[,198]
MAR.c.within.transm.var.cov.70 <- run.ALL[,199]


# Cov 75
MAR.c.av.age.male.cov.75 <- run.ALL[,200]
MAR.c.av.age.diff.cov.75 <- run.ALL[,201]
MAR.c.between.transm.var.cov.75 <- run.ALL[,202]
MAR.c.within.transm.var.cov.75 <- run.ALL[,203]

# Cov 80
MAR.c.av.age.male.cov.80 <- run.ALL[,204]
MAR.c.av.age.diff.cov.80 <- run.ALL[,205]
MAR.c.between.transm.var.cov.80 <- run.ALL[,206]
MAR.c.within.transm.var.cov.80 <- run.ALL[,207]


# Cov 85
MAR.c.av.age.male.cov.85 <- run.ALL[,208]
MAR.c.av.age.diff.cov.85 <- run.ALL[,209]
MAR.c.between.transm.var.cov.85 <- run.ALL[,210]
MAR.c.within.transm.var.cov.85 <- run.ALL[,211]

# Cov 90
MAR.c.av.age.male.cov.90 <- run.ALL[,212]
MAR.c.av.age.diff.cov.90 <- run.ALL[,213]
MAR.c.between.transm.var.cov.90 <- run.ALL[,214]
MAR.c.within.transm.var.cov.90 <- run.ALL[,215]


# Cov 95
MAR.c.av.age.male.cov.95 <- run.ALL[,216]
MAR.c.av.age.diff.cov.95 <- run.ALL[,217]
MAR.c.between.transm.var.cov.95 <- run.ALL[,218]
MAR.c.within.transm.var.cov.95 <- run.ALL[,219]


# Comparisons between population-level stats and sequence coverage scenarios

# Population level
av.age.male.pop <- av.age.male.v
av.age.diff.pop <- av.age.diff.v
between.transm.var.pop <- between.transm.var.v
within.transm.var.pop <- within.transm.var.v

my.df.pop <- data.frame(av.age.male.pop, av.age.diff.pop, between.transm.var.pop, within.transm.var.pop)


# MCAR

MCAR.df.35 <- data.frame(MCAR.av.age.male.cov.35, MCAR.av.age.diff.cov.35, MCAR.between.transm.var.cov.35, MCAR.within.transm.var.cov.35)
MCAR.df.40 <- data.frame(MCAR.av.age.male.cov.40, MCAR.av.age.diff.cov.40, MCAR.between.transm.var.cov.40, MCAR.within.transm.var.cov.40)
MCAR.df.45 <- data.frame(MCAR.av.age.male.cov.45, MCAR.av.age.diff.cov.45, MCAR.between.transm.var.cov.45, MCAR.within.transm.var.cov.45)
MCAR.df.50 <- data.frame(MCAR.av.age.male.cov.50, MCAR.av.age.diff.cov.50, MCAR.between.transm.var.cov.50, MCAR.within.transm.var.cov.50)
MCAR.df.55 <- data.frame(MCAR.av.age.male.cov.55, MCAR.av.age.diff.cov.55, MCAR.between.transm.var.cov.55, MCAR.within.transm.var.cov.55)
MCAR.df.60 <- data.frame(MCAR.av.age.male.cov.60, MCAR.av.age.diff.cov.60, MCAR.between.transm.var.cov.60, MCAR.within.transm.var.cov.60)
MCAR.df.65 <- data.frame(MCAR.av.age.male.cov.65, MCAR.av.age.diff.cov.65, MCAR.between.transm.var.cov.65, MCAR.within.transm.var.cov.65)
MCAR.df.70 <- data.frame(MCAR.av.age.male.cov.70, MCAR.av.age.diff.cov.70, MCAR.between.transm.var.cov.70, MCAR.within.transm.var.cov.70)
MCAR.df.75 <- data.frame(MCAR.av.age.male.cov.75, MCAR.av.age.diff.cov.75, MCAR.between.transm.var.cov.75, MCAR.within.transm.var.cov.75)
MCAR.df.80 <- data.frame(MCAR.av.age.male.cov.80, MCAR.av.age.diff.cov.80, MCAR.between.transm.var.cov.80, MCAR.within.transm.var.cov.80)
MCAR.df.85 <- data.frame(MCAR.av.age.male.cov.85, MCAR.av.age.diff.cov.85, MCAR.between.transm.var.cov.85, MCAR.within.transm.var.cov.85)
MCAR.df.90 <- data.frame(MCAR.av.age.male.cov.90, MCAR.av.age.diff.cov.90, MCAR.between.transm.var.cov.90, MCAR.within.transm.var.cov.90)
MCAR.df.95 <- data.frame(MCAR.av.age.male.cov.95, MCAR.av.age.diff.cov.95, MCAR.between.transm.var.cov.95, MCAR.within.transm.var.cov.95)

# MAR

# a

MAR.a.df.35 <- data.frame(MAR.a.av.age.male.cov.35, MAR.a.av.age.diff.cov.35, MAR.a.between.transm.var.cov.35, MAR.a.within.transm.var.cov.35)
MAR.a.df.40 <- data.frame(MAR.a.av.age.male.cov.40, MAR.a.av.age.diff.cov.40, MAR.a.between.transm.var.cov.40, MAR.a.within.transm.var.cov.40)
MAR.a.df.45 <- data.frame(MAR.a.av.age.male.cov.45, MAR.a.av.age.diff.cov.45, MAR.a.between.transm.var.cov.45, MAR.a.within.transm.var.cov.45)
MAR.a.df.50 <- data.frame(MAR.a.av.age.male.cov.50, MAR.a.av.age.diff.cov.50, MAR.a.between.transm.var.cov.50, MAR.a.within.transm.var.cov.50)
MAR.a.df.55 <- data.frame(MAR.a.av.age.male.cov.55, MAR.a.av.age.diff.cov.55, MAR.a.between.transm.var.cov.55, MAR.a.within.transm.var.cov.55)
MAR.a.df.60 <- data.frame(MAR.a.av.age.male.cov.60, MAR.a.av.age.diff.cov.60, MAR.a.between.transm.var.cov.60, MAR.a.within.transm.var.cov.60)
MAR.a.df.65 <- data.frame(MAR.a.av.age.male.cov.65, MAR.a.av.age.diff.cov.65, MAR.a.between.transm.var.cov.65, MAR.a.within.transm.var.cov.65)
MAR.a.df.70 <- data.frame(MAR.a.av.age.male.cov.70, MAR.a.av.age.diff.cov.70, MAR.a.between.transm.var.cov.70, MAR.a.within.transm.var.cov.70)
MAR.a.df.75 <- data.frame(MAR.a.av.age.male.cov.75, MAR.a.av.age.diff.cov.75, MAR.a.between.transm.var.cov.75, MAR.a.within.transm.var.cov.75)
MAR.a.df.80 <- data.frame(MAR.a.av.age.male.cov.80, MAR.a.av.age.diff.cov.80, MAR.a.between.transm.var.cov.80, MAR.a.within.transm.var.cov.80)
MAR.a.df.85 <- data.frame(MAR.a.av.age.male.cov.85, MAR.a.av.age.diff.cov.85, MAR.a.between.transm.var.cov.85, MAR.a.within.transm.var.cov.85)
MAR.a.df.90 <- data.frame(MAR.a.av.age.male.cov.90, MAR.a.av.age.diff.cov.90, MAR.a.between.transm.var.cov.90, MAR.a.within.transm.var.cov.90)
MAR.a.df.95 <- data.frame(MAR.a.av.age.male.cov.95, MAR.a.av.age.diff.cov.95, MAR.a.between.transm.var.cov.95, MAR.a.within.transm.var.cov.95)

# b

MAR.b.df.35 <- data.frame(MAR.b.av.age.male.cov.35, MAR.b.av.age.diff.cov.35, MAR.b.between.transm.var.cov.35, MAR.b.within.transm.var.cov.35)
MAR.b.df.40 <- data.frame(MAR.b.av.age.male.cov.40, MAR.b.av.age.diff.cov.40, MAR.b.between.transm.var.cov.40, MAR.b.within.transm.var.cov.40)
MAR.b.df.45 <- data.frame(MAR.b.av.age.male.cov.45, MAR.b.av.age.diff.cov.45, MAR.b.between.transm.var.cov.45, MAR.b.within.transm.var.cov.45)
MAR.b.df.50 <- data.frame(MAR.b.av.age.male.cov.50, MAR.b.av.age.diff.cov.50, MAR.b.between.transm.var.cov.50, MAR.b.within.transm.var.cov.50)
MAR.b.df.55 <- data.frame(MAR.b.av.age.male.cov.55, MAR.b.av.age.diff.cov.55, MAR.b.between.transm.var.cov.55, MAR.b.within.transm.var.cov.55)
MAR.b.df.60 <- data.frame(MAR.b.av.age.male.cov.60, MAR.b.av.age.diff.cov.60, MAR.b.between.transm.var.cov.60, MAR.b.within.transm.var.cov.60)
MAR.b.df.65 <- data.frame(MAR.b.av.age.male.cov.65, MAR.b.av.age.diff.cov.65, MAR.b.between.transm.var.cov.65, MAR.b.within.transm.var.cov.65)
MAR.b.df.70 <- data.frame(MAR.b.av.age.male.cov.70, MAR.b.av.age.diff.cov.70, MAR.b.between.transm.var.cov.70, MAR.b.within.transm.var.cov.70)
MAR.b.df.75 <- data.frame(MAR.b.av.age.male.cov.75, MAR.b.av.age.diff.cov.75, MAR.b.between.transm.var.cov.75, MAR.b.within.transm.var.cov.75)
MAR.b.df.80 <- data.frame(MAR.b.av.age.male.cov.80, MAR.b.av.age.diff.cov.80, MAR.b.between.transm.var.cov.80, MAR.b.within.transm.var.cov.80)
MAR.b.df.85 <- data.frame(MAR.b.av.age.male.cov.85, MAR.b.av.age.diff.cov.85, MAR.b.between.transm.var.cov.85, MAR.b.within.transm.var.cov.85)
MAR.b.df.90 <- data.frame(MAR.b.av.age.male.cov.90, MAR.b.av.age.diff.cov.90, MAR.b.between.transm.var.cov.90, MAR.b.within.transm.var.cov.90)
MAR.b.df.95 <- data.frame(MAR.b.av.age.male.cov.95, MAR.b.av.age.diff.cov.95, MAR.b.between.transm.var.cov.95, MAR.b.within.transm.var.cov.95)

# c
# c

MAR.c.df.35 <- data.frame(MAR.c.av.age.male.cov.35, MAR.c.av.age.diff.cov.35, MAR.c.between.transm.var.cov.35, MAR.c.within.transm.var.cov.35)
MAR.c.df.40 <- data.frame(MAR.c.av.age.male.cov.40, MAR.c.av.age.diff.cov.40, MAR.c.between.transm.var.cov.40, MAR.c.within.transm.var.cov.40)
MAR.c.df.45 <- data.frame(MAR.c.av.age.male.cov.45, MAR.c.av.age.diff.cov.45, MAR.c.between.transm.var.cov.45, MAR.c.within.transm.var.cov.45)
MAR.c.df.50 <- data.frame(MAR.c.av.age.male.cov.50, MAR.c.av.age.diff.cov.50, MAR.c.between.transm.var.cov.50, MAR.c.within.transm.var.cov.50)
MAR.c.df.55 <- data.frame(MAR.c.av.age.male.cov.55, MAR.c.av.age.diff.cov.55, MAR.c.between.transm.var.cov.55, MAR.c.within.transm.var.cov.55)
MAR.c.df.60 <- data.frame(MAR.c.av.age.male.cov.60, MAR.c.av.age.diff.cov.60, MAR.c.between.transm.var.cov.60, MAR.c.within.transm.var.cov.60)
MAR.c.df.65 <- data.frame(MAR.c.av.age.male.cov.65, MAR.c.av.age.diff.cov.65, MAR.c.between.transm.var.cov.65, MAR.c.within.transm.var.cov.65)
MAR.c.df.70 <- data.frame(MAR.c.av.age.male.cov.70, MAR.c.av.age.diff.cov.70, MAR.c.between.transm.var.cov.70, MAR.c.within.transm.var.cov.70)
MAR.c.df.75 <- data.frame(MAR.c.av.age.male.cov.75, MAR.c.av.age.diff.cov.75, MAR.c.between.transm.var.cov.75, MAR.c.within.transm.var.cov.75)
MAR.c.df.80 <- data.frame(MAR.c.av.age.male.cov.80, MAR.c.av.age.diff.cov.80, MAR.c.between.transm.var.cov.80, MAR.c.within.transm.var.cov.80)
MAR.c.df.85 <- data.frame(MAR.c.av.age.male.cov.85, MAR.c.av.age.diff.cov.85, MAR.c.between.transm.var.cov.85, MAR.c.within.transm.var.cov.85)
MAR.c.df.90 <- data.frame(MAR.c.av.age.male.cov.90, MAR.c.av.age.diff.cov.90, MAR.c.between.transm.var.cov.90, MAR.c.within.transm.var.cov.90)
MAR.c.df.95 <- data.frame(MAR.c.av.age.male.cov.95, MAR.c.av.age.diff.cov.95, MAR.c.between.transm.var.cov.95, MAR.c.within.transm.var.cov.95)



# x = rnorm(10)
# y = rnorm(10)
# t.test(x,y)
# 
# confidence_interval <- function(vector, interval) {
#   # Standard deviation of sample
#   vec_sd <- sd(vector)
#   # Sample size
#   n <- length(vector)
#   # Mean of sample
#   vec_mean <- mean(vector)
#   # Error according to t distribution
#   error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
#   # Confidence interval as a vector
#   result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
#   return(result)
# }
# 
# 
# vector <- c(12, 17, 24, 35, 23, 34, 56)
# confidence_interval(vector, 0.90)
# 
# 
# confidence_interval(av.age.male.pop, 0.95)
# 
# confidence_interval(MCAR.df.95$MCAR.av.age.male.cov.95, 0.95)

# Or use gmodels

# Estimates averages values and confidence intervals

library(gmodels)

ci(av.age.male.pop, confidence=0.95)

# Population level
Est.age.male.pop <- ci(my.df.pop$av.age.male.pop, confidence=0.95)
Est.age.diff.pop <- ci(my.df.pop$av.age.diff.pop, confidence = 0.95)

Est.between.transm.pop <- ci(my.df.pop$between.transm.var.pop, confidence=0.95)
Est.within.transm.pop <- ci(my.df.pop$within.transm.var.pop, confidence = 0.95)

# MCAR
Est.MCAR.age.male.35 <- ci(MCAR.df.35$MCAR.av.age.male.cov.35, confidence=0.95)
Est.MCAR.age.diff.35 <- ci(MCAR.df.35$MCAR.av.age.diff.cov.35, confidence = 0.95)

Est.MCAR.between.transm.35 <- ci(MCAR.df.35$MCAR.between.transm.var.cov.35, confidence=0.95)
Est.MCAR.within.transm.35 <- ci(MCAR.df.35$MCAR.within.transm.var.cov.35, confidence = 0.95)

Est.MCAR.age.male.40 <- ci(MCAR.df.40$MCAR.av.age.male.cov.40, confidence=0.95)
Est.MCAR.age.diff.40 <- ci(MCAR.df.40$MCAR.av.age.diff.cov.40, confidence = 0.95)

Est.MCAR.between.transm.40 <- ci(MCAR.df.40$MCAR.between.transm.var.cov.40, confidence=0.95)
Est.MCAR.within.transm.40 <- ci(MCAR.df.40$MCAR.within.transm.var.cov.40, confidence = 0.95)

Est.MCAR.age.male.45 <- ci(MCAR.df.45$MCAR.av.age.male.cov.45, confidence=0.95)
Est.MCAR.age.diff.45 <- ci(MCAR.df.45$MCAR.av.age.diff.cov.45, confidence = 0.95)

Est.MCAR.between.transm.45 <- ci(MCAR.df.45$MCAR.between.transm.var.cov.45, confidence=0.95)
Est.MCAR.within.transm.45 <- ci(MCAR.df.45$MCAR.within.transm.var.cov.45, confidence = 0.95)

Est.MCAR.age.male.50 <- ci(MCAR.df.50$MCAR.av.age.male.cov.50, confidence=0.95)
Est.MCAR.age.diff.50 <- ci(MCAR.df.50$MCAR.av.age.diff.cov.50, confidence = 0.95)

Est.MCAR.between.transm.50 <- ci(MCAR.df.50$MCAR.between.transm.var.cov.50, confidence=0.95)
Est.MCAR.within.transm.50 <- ci(MCAR.df.50$MCAR.within.transm.var.cov.50, confidence = 0.95)

Est.MCAR.age.male.55 <- ci(MCAR.df.55$MCAR.av.age.male.cov.55, confidence=0.95)
Est.MCAR.age.diff.55 <- ci(MCAR.df.55$MCAR.av.age.diff.cov.55, confidence = 0.95)

Est.MCAR.between.transm.55 <- ci(MCAR.df.55$MCAR.between.transm.var.cov.55, confidence=0.95)
Est.MCAR.within.transm.55 <- ci(MCAR.df.55$MCAR.within.transm.var.cov.55, confidence = 0.95)

Est.MCAR.age.male.60 <- ci(MCAR.df.60$MCAR.av.age.male.cov.60, confidence=0.95)
Est.MCAR.age.diff.60 <- ci(MCAR.df.60$MCAR.av.age.diff.cov.60, confidence = 0.95)

Est.MCAR.between.transm.60 <- ci(MCAR.df.60$MCAR.between.transm.var.cov.60, confidence=0.95)
Est.MCAR.within.transm.60 <- ci(MCAR.df.60$MCAR.within.transm.var.cov.60, confidence = 0.95)

Est.MCAR.age.male.65 <- ci(MCAR.df.65$MCAR.av.age.male.cov.65, confidence=0.95)
Est.MCAR.age.diff.65 <- ci(MCAR.df.65$MCAR.av.age.diff.cov.65, confidence = 0.95)

Est.MCAR.between.transm.65 <- ci(MCAR.df.65$MCAR.between.transm.var.cov.65, confidence=0.95)
Est.MCAR.within.transm.65 <- ci(MCAR.df.65$MCAR.within.transm.var.cov.65, confidence = 0.95)

Est.MCAR.age.male.70 <- ci(MCAR.df.70$MCAR.av.age.male.cov.70, confidence=0.95)
Est.MCAR.age.diff.70 <- ci(MCAR.df.70$MCAR.av.age.diff.cov.70, confidence = 0.95)

Est.MCAR.between.transm.70 <- ci(MCAR.df.70$MCAR.between.transm.var.cov.70, confidence=0.95)
Est.MCAR.within.transm.70 <- ci(MCAR.df.70$MCAR.within.transm.var.cov.70, confidence = 0.95)

Est.MCAR.age.male.75 <- ci(MCAR.df.75$MCAR.av.age.male.cov.75, confidence=0.95)
Est.MCAR.age.diff.75 <- ci(MCAR.df.75$MCAR.av.age.diff.cov.75, confidence = 0.95)

Est.MCAR.between.transm.75 <- ci(MCAR.df.75$MCAR.between.transm.var.cov.75, confidence=0.95)
Est.MCAR.within.transm.75 <- ci(MCAR.df.75$MCAR.within.transm.var.cov.75, confidence = 0.95)

Est.MCAR.age.male.80 <- ci(MCAR.df.80$MCAR.av.age.male.cov.80, confidence=0.95)
Est.MCAR.age.diff.80 <- ci(MCAR.df.80$MCAR.av.age.diff.cov.80, confidence = 0.95)

Est.MCAR.between.transm.80 <- ci(MCAR.df.80$MCAR.between.transm.var.cov.80, confidence=0.95)
Est.MCAR.within.transm.80 <- ci(MCAR.df.80$MCAR.within.transm.var.cov.80, confidence = 0.95)

Est.MCAR.age.male.85 <- ci(MCAR.df.85$MCAR.av.age.male.cov.85, confidence=0.95)
Est.MCAR.age.diff.85 <- ci(MCAR.df.85$MCAR.av.age.diff.cov.85, confidence = 0.95)

Est.MCAR.between.transm.85 <- ci(MCAR.df.85$MCAR.between.transm.var.cov.85, confidence=0.95)
Est.MCAR.within.transm.85 <- ci(MCAR.df.85$MCAR.within.transm.var.cov.85, confidence = 0.95)

Est.MCAR.age.male.90 <- ci(MCAR.df.90$MCAR.av.age.male.cov.90, confidence=0.95)
Est.MCAR.age.diff.90 <- ci(MCAR.df.90$MCAR.av.age.diff.cov.90, confidence = 0.95)

Est.MCAR.between.transm.90 <- ci(MCAR.df.90$MCAR.between.transm.var.cov.90, confidence=0.95)
Est.MCAR.within.transm.90 <- ci(MCAR.df.90$MCAR.within.transm.var.cov.90, confidence = 0.95)

Est.MCAR.age.male.95 <- ci(MCAR.df.95$MCAR.av.age.male.cov.95, confidence=0.95)
Est.MCAR.age.diff.95 <- ci(MCAR.df.95$MCAR.av.age.diff.cov.95, confidence = 0.95)

Est.MCAR.between.transm.95 <- ci(MCAR.df.95$MCAR.between.transm.var.cov.95, confidence=0.95)
Est.MCAR.within.transm.95 <- ci(MCAR.df.95$MCAR.within.transm.var.cov.95, confidence = 0.95)

# MAR

# a
Est.MAR.a.age.male.35 <- ci(MAR.a.df.35$MAR.a.av.age.male.cov.35, confidence=0.95)
Est.MAR.a.age.diff.35 <- ci(MAR.a.df.35$MAR.a.av.age.diff.cov.35, confidence = 0.95)

Est.MAR.a.between.transm.35 <- ci(MAR.a.df.35$MAR.a.between.transm.var.cov.35, confidence=0.95)
Est.MAR.a.within.transm.35 <- ci(MAR.a.df.35$MAR.a.within.transm.var.cov.35, confidence = 0.95)

Est.MAR.a.age.male.40 <- ci(MAR.a.df.40$MAR.a.av.age.male.cov.40, confidence=0.95)
Est.MAR.a.age.diff.40 <- ci(MAR.a.df.40$MAR.a.av.age.diff.cov.40, confidence = 0.95)

Est.MAR.a.between.transm.40 <- ci(MAR.a.df.40$MAR.a.between.transm.var.cov.40, confidence=0.95)
Est.MAR.a.within.transm.40 <- ci(MAR.a.df.40$MAR.a.within.transm.var.cov.40, confidence = 0.95)

Est.MAR.a.age.male.45 <- ci(MAR.a.df.45$MAR.a.av.age.male.cov.45, confidence=0.95)
Est.MAR.a.age.diff.45 <- ci(MAR.a.df.45$MAR.a.av.age.diff.cov.45, confidence = 0.95)

Est.MAR.a.between.transm.45 <- ci(MAR.a.df.45$MAR.a.between.transm.var.cov.45, confidence=0.95)
Est.MAR.a.within.transm.45 <- ci(MAR.a.df.45$MAR.a.within.transm.var.cov.45, confidence = 0.95)

Est.MAR.a.age.male.50 <- ci(MAR.a.df.50$MAR.a.av.age.male.cov.50, confidence=0.95)
Est.MAR.a.age.diff.50 <- ci(MAR.a.df.50$MAR.a.av.age.diff.cov.50, confidence = 0.95)

Est.MAR.a.between.transm.50 <- ci(MAR.a.df.50$MAR.a.between.transm.var.cov.50, confidence=0.95)
Est.MAR.a.within.transm.50 <- ci(MAR.a.df.50$MAR.a.within.transm.var.cov.50, confidence = 0.95)

Est.MAR.a.age.male.55 <- ci(MAR.a.df.55$MAR.a.av.age.male.cov.55, confidence=0.95)
Est.MAR.a.age.diff.55 <- ci(MAR.a.df.55$MAR.a.av.age.diff.cov.55, confidence = 0.95)

Est.MAR.a.between.transm.55 <- ci(MAR.a.df.55$MAR.a.between.transm.var.cov.55, confidence=0.95)
Est.MAR.a.within.transm.55 <- ci(MAR.a.df.55$MAR.a.within.transm.var.cov.55, confidence = 0.95)

Est.MAR.a.age.male.60 <- ci(MAR.a.df.60$MAR.a.av.age.male.cov.60, confidence=0.95)
Est.MAR.a.age.diff.60 <- ci(MAR.a.df.60$MAR.a.av.age.diff.cov.60, confidence = 0.95)

Est.MAR.a.between.transm.60 <- ci(MAR.a.df.60$MAR.a.between.transm.var.cov.60, confidence=0.95)
Est.MAR.a.within.transm.60 <- ci(MAR.a.df.60$MAR.a.within.transm.var.cov.60, confidence = 0.95)

Est.MAR.a.age.male.65 <- ci(MAR.a.df.65$MAR.a.av.age.male.cov.65, confidence=0.95)
Est.MAR.a.age.diff.65 <- ci(MAR.a.df.65$MAR.a.av.age.diff.cov.65, confidence = 0.95)

Est.MAR.a.between.transm.65 <- ci(MAR.a.df.65$MAR.a.between.transm.var.cov.65, confidence=0.95)
Est.MAR.a.within.transm.65 <- ci(MAR.a.df.65$MAR.a.within.transm.var.cov.65, confidence = 0.95)

Est.MAR.a.age.male.70 <- ci(MAR.a.df.70$MAR.a.av.age.male.cov.70, confidence=0.95)
Est.MAR.a.age.diff.70 <- ci(MAR.a.df.70$MAR.a.av.age.diff.cov.70, confidence = 0.95)

Est.MAR.a.between.transm.70 <- ci(MAR.a.df.70$MAR.a.between.transm.var.cov.70, confidence=0.95)
Est.MAR.a.within.transm.70 <- ci(MAR.a.df.70$MAR.a.within.transm.var.cov.70, confidence = 0.95)

Est.MAR.a.age.male.75 <- ci(MAR.a.df.75$MAR.a.av.age.male.cov.75, confidence=0.95)
Est.MAR.a.age.diff.75 <- ci(MAR.a.df.75$MAR.a.av.age.diff.cov.75, confidence = 0.95)

Est.MAR.a.between.transm.75 <- ci(MAR.a.df.75$MAR.a.between.transm.var.cov.75, confidence=0.95)
Est.MAR.a.within.transm.75 <- ci(MAR.a.df.75$MAR.a.within.transm.var.cov.75, confidence = 0.95)

Est.MAR.a.age.male.80 <- ci(MAR.a.df.80$MAR.a.av.age.male.cov.80, confidence=0.95)
Est.MAR.a.age.diff.80 <- ci(MAR.a.df.80$MAR.a.av.age.diff.cov.80, confidence = 0.95)

Est.MAR.a.between.transm.80 <- ci(MAR.a.df.80$MAR.a.between.transm.var.cov.80, confidence=0.95)
Est.MAR.a.within.transm.80 <- ci(MAR.a.df.80$MAR.a.within.transm.var.cov.80, confidence = 0.95)

Est.MAR.a.age.male.85 <- ci(MAR.a.df.85$MAR.a.av.age.male.cov.85, confidence=0.95)
Est.MAR.a.age.diff.85 <- ci(MAR.a.df.85$MAR.a.av.age.diff.cov.85, confidence = 0.95)

Est.MAR.a.between.transm.85 <- ci(MAR.a.df.85$MAR.a.between.transm.var.cov.85, confidence=0.95)
Est.MAR.a.within.transm.85 <- ci(MAR.a.df.85$MAR.a.within.transm.var.cov.85, confidence = 0.95)

Est.MAR.a.age.male.90 <- ci(MAR.a.df.90$MAR.a.av.age.male.cov.90, confidence=0.95)
Est.MAR.a.age.diff.90 <- ci(MAR.a.df.90$MAR.a.av.age.diff.cov.90, confidence = 0.95)

Est.MAR.a.between.transm.90 <- ci(MAR.a.df.90$MAR.a.between.transm.var.cov.90, confidence=0.95)
Est.MAR.a.within.transm.90 <- ci(MAR.a.df.90$MAR.a.within.transm.var.cov.90, confidence = 0.95)

Est.MAR.a.age.male.95 <- ci(MAR.a.df.95$MAR.a.av.age.male.cov.95, confidence=0.95)
Est.MAR.a.age.diff.95 <- ci(MAR.a.df.95$MAR.a.av.age.diff.cov.95, confidence = 0.95)

Est.MAR.a.between.transm.95 <- ci(MAR.a.df.95$MAR.a.between.transm.var.cov.95, confidence=0.95)
Est.MAR.a.within.transm.95 <- ci(MAR.a.df.95$MAR.a.within.transm.var.cov.95, confidence = 0.95)


# b
Est.MAR.b.age.male.35 <- ci(MAR.b.df.35$MAR.b.av.age.male.cov.35, confidence=0.95)
Est.MAR.b.age.diff.35 <- ci(MAR.b.df.35$MAR.b.av.age.diff.cov.35, confidence = 0.95)

Est.MAR.b.between.transm.35 <- ci(MAR.b.df.35$MAR.b.between.transm.var.cov.35, confidence=0.95)
Est.MAR.b.within.transm.35 <- ci(MAR.b.df.35$MAR.b.within.transm.var.cov.35, confidence = 0.95)

Est.MAR.b.age.male.40 <- ci(MAR.b.df.40$MAR.b.av.age.male.cov.40, confidence=0.95)
Est.MAR.b.age.diff.40 <- ci(MAR.b.df.40$MAR.b.av.age.diff.cov.40, confidence = 0.95)

Est.MAR.b.between.transm.40 <- ci(MAR.b.df.40$MAR.b.between.transm.var.cov.40, confidence=0.95)
Est.MAR.b.within.transm.40 <- ci(MAR.b.df.40$MAR.b.within.transm.var.cov.40, confidence = 0.95)

Est.MAR.b.age.male.45 <- ci(MAR.b.df.45$MAR.b.av.age.male.cov.45, confidence=0.95)
Est.MAR.b.age.diff.45 <- ci(MAR.b.df.45$MAR.b.av.age.diff.cov.45, confidence = 0.95)

Est.MAR.b.between.transm.45 <- ci(MAR.b.df.45$MAR.b.between.transm.var.cov.45, confidence=0.95)
Est.MAR.b.within.transm.45 <- ci(MAR.b.df.45$MAR.b.within.transm.var.cov.45, confidence = 0.95)

Est.MAR.b.age.male.50 <- ci(MAR.b.df.50$MAR.b.av.age.male.cov.50, confidence=0.95)
Est.MAR.b.age.diff.50 <- ci(MAR.b.df.50$MAR.b.av.age.diff.cov.50, confidence = 0.95)

Est.MAR.b.between.transm.50 <- ci(MAR.b.df.50$MAR.b.between.transm.var.cov.50, confidence=0.95)
Est.MAR.b.within.transm.50 <- ci(MAR.b.df.50$MAR.b.within.transm.var.cov.50, confidence = 0.95)

Est.MAR.b.age.male.55 <- ci(MAR.b.df.55$MAR.b.av.age.male.cov.55, confidence=0.95)
Est.MAR.b.age.diff.55 <- ci(MAR.b.df.55$MAR.b.av.age.diff.cov.55, confidence = 0.95)

Est.MAR.b.between.transm.55 <- ci(MAR.b.df.55$MAR.b.between.transm.var.cov.55, confidence=0.95)
Est.MAR.b.within.transm.55 <- ci(MAR.b.df.55$MAR.b.within.transm.var.cov.55, confidence = 0.95)

Est.MAR.b.age.male.60 <- ci(MAR.b.df.60$MAR.b.av.age.male.cov.60, confidence=0.95)
Est.MAR.b.age.diff.60 <- ci(MAR.b.df.60$MAR.b.av.age.diff.cov.60, confidence = 0.95)

Est.MAR.b.between.transm.60 <- ci(MAR.b.df.60$MAR.b.between.transm.var.cov.60, confidence=0.95)
Est.MAR.b.within.transm.60 <- ci(MAR.b.df.60$MAR.b.within.transm.var.cov.60, confidence = 0.95)

Est.MAR.b.age.male.65 <- ci(MAR.b.df.65$MAR.b.av.age.male.cov.65, confidence=0.95)
Est.MAR.b.age.diff.65 <- ci(MAR.b.df.65$MAR.b.av.age.diff.cov.65, confidence = 0.95)

Est.MAR.b.between.transm.65 <- ci(MAR.b.df.65$MAR.b.between.transm.var.cov.65, confidence=0.95)
Est.MAR.b.within.transm.65 <- ci(MAR.b.df.65$MAR.b.within.transm.var.cov.65, confidence = 0.95)

Est.MAR.b.age.male.70 <- ci(MAR.b.df.70$MAR.b.av.age.male.cov.70, confidence=0.95)
Est.MAR.b.age.diff.70 <- ci(MAR.b.df.70$MAR.b.av.age.diff.cov.70, confidence = 0.95)

Est.MAR.b.between.transm.70 <- ci(MAR.b.df.70$MAR.b.between.transm.var.cov.70, confidence=0.95)
Est.MAR.b.within.transm.70 <- ci(MAR.b.df.70$MAR.b.within.transm.var.cov.70, confidence = 0.95)

Est.MAR.b.age.male.75 <- ci(MAR.b.df.75$MAR.b.av.age.male.cov.75, confidence=0.95)
Est.MAR.b.age.diff.75 <- ci(MAR.b.df.75$MAR.b.av.age.diff.cov.75, confidence = 0.95)

Est.MAR.b.between.transm.75 <- ci(MAR.b.df.75$MAR.b.between.transm.var.cov.75, confidence=0.95)
Est.MAR.b.within.transm.75 <- ci(MAR.b.df.75$MAR.b.within.transm.var.cov.75, confidence = 0.95)

Est.MAR.b.age.male.80 <- ci(MAR.b.df.80$MAR.b.av.age.male.cov.80, confidence=0.95)
Est.MAR.b.age.diff.80 <- ci(MAR.b.df.80$MAR.b.av.age.diff.cov.80, confidence = 0.95)

Est.MAR.b.between.transm.80 <- ci(MAR.b.df.80$MAR.b.between.transm.var.cov.80, confidence=0.95)
Est.MAR.b.within.transm.80 <- ci(MAR.b.df.80$MAR.b.within.transm.var.cov.80, confidence = 0.95)

Est.MAR.b.age.male.85 <- ci(MAR.b.df.85$MAR.b.av.age.male.cov.85, confidence=0.95)
Est.MAR.b.age.diff.85 <- ci(MAR.b.df.85$MAR.b.av.age.diff.cov.85, confidence = 0.95)

Est.MAR.b.between.transm.85 <- ci(MAR.b.df.85$MAR.b.between.transm.var.cov.85, confidence=0.95)
Est.MAR.b.within.transm.85 <- ci(MAR.b.df.85$MAR.b.within.transm.var.cov.85, confidence = 0.95)

Est.MAR.b.age.male.90 <- ci(MAR.b.df.90$MAR.b.av.age.male.cov.90, confidence=0.95)
Est.MAR.b.age.diff.90 <- ci(MAR.b.df.90$MAR.b.av.age.diff.cov.90, confidence = 0.95)

Est.MAR.b.between.transm.90 <- ci(MAR.b.df.90$MAR.b.between.transm.var.cov.90, confidence=0.95)
Est.MAR.b.within.transm.90 <- ci(MAR.b.df.90$MAR.b.within.transm.var.cov.90, confidence = 0.95)

Est.MAR.b.age.male.95 <- ci(MAR.b.df.95$MAR.b.av.age.male.cov.95, confidence=0.95)
Est.MAR.b.age.diff.95 <- ci(MAR.b.df.95$MAR.b.av.age.diff.cov.95, confidence = 0.95)

Est.MAR.b.between.transm.95 <- ci(MAR.b.df.95$MAR.b.between.transm.var.cov.95, confidence=0.95)
Est.MAR.b.within.transm.95 <- ci(MAR.b.df.95$MAR.b.within.transm.var.cov.95, confidence = 0.95)

# c
Est.MAR.c.age.male.35 <- ci(MAR.c.df.35$MAR.c.av.age.male.cov.35, confidence=0.95)
Est.MAR.c.age.diff.35 <- ci(MAR.c.df.35$MAR.c.av.age.diff.cov.35, confidence = 0.95)

Est.MAR.c.between.transm.35 <- ci(MAR.c.df.35$MAR.c.between.transm.var.cov.35, confidence=0.95)
Est.MAR.c.within.transm.35 <- ci(MAR.c.df.35$MAR.c.within.transm.var.cov.35, confidence = 0.95)

Est.MAR.c.age.male.40 <- ci(MAR.c.df.40$MAR.c.av.age.male.cov.40, confidence=0.95)
Est.MAR.c.age.diff.40 <- ci(MAR.c.df.40$MAR.c.av.age.diff.cov.40, confidence = 0.95)

Est.MAR.c.between.transm.40 <- ci(MAR.c.df.40$MAR.c.between.transm.var.cov.40, confidence=0.95)
Est.MAR.c.within.transm.40 <- ci(MAR.c.df.40$MAR.c.within.transm.var.cov.40, confidence = 0.95)

Est.MAR.c.age.male.45 <- ci(MAR.c.df.45$MAR.c.av.age.male.cov.45, confidence=0.95)
Est.MAR.c.age.diff.45 <- ci(MAR.c.df.45$MAR.c.av.age.diff.cov.45, confidence = 0.95)

Est.MAR.c.between.transm.45 <- ci(MAR.c.df.45$MAR.c.between.transm.var.cov.45, confidence=0.95)
Est.MAR.c.within.transm.45 <- ci(MAR.c.df.45$MAR.c.within.transm.var.cov.45, confidence = 0.95)

Est.MAR.c.age.male.50 <- ci(MAR.c.df.50$MAR.c.av.age.male.cov.50, confidence=0.95)
Est.MAR.c.age.diff.50 <- ci(MAR.c.df.50$MAR.c.av.age.diff.cov.50, confidence = 0.95)

Est.MAR.c.between.transm.50 <- ci(MAR.c.df.50$MAR.c.between.transm.var.cov.50, confidence=0.95)
Est.MAR.c.within.transm.50 <- ci(MAR.c.df.50$MAR.c.within.transm.var.cov.50, confidence = 0.95)

Est.MAR.c.age.male.55 <- ci(MAR.c.df.55$MAR.c.av.age.male.cov.55, confidence=0.95)
Est.MAR.c.age.diff.55 <- ci(MAR.c.df.55$MAR.c.av.age.diff.cov.55, confidence = 0.95)

Est.MAR.c.between.transm.55 <- ci(MAR.c.df.55$MAR.c.between.transm.var.cov.55, confidence=0.95)
Est.MAR.c.within.transm.55 <- ci(MAR.c.df.55$MAR.c.within.transm.var.cov.55, confidence = 0.95)

Est.MAR.c.age.male.60 <- ci(MAR.c.df.60$MAR.c.av.age.male.cov.60, confidence=0.95)
Est.MAR.c.age.diff.60 <- ci(MAR.c.df.60$MAR.c.av.age.diff.cov.60, confidence = 0.95)

Est.MAR.c.between.transm.60 <- ci(MAR.c.df.60$MAR.c.between.transm.var.cov.60, confidence=0.95)
Est.MAR.c.within.transm.60 <- ci(MAR.c.df.60$MAR.c.within.transm.var.cov.60, confidence = 0.95)

Est.MAR.c.age.male.65 <- ci(MAR.c.df.65$MAR.c.av.age.male.cov.65, confidence=0.95)
Est.MAR.c.age.diff.65 <- ci(MAR.c.df.65$MAR.c.av.age.diff.cov.65, confidence = 0.95)

Est.MAR.c.between.transm.65 <- ci(MAR.c.df.65$MAR.c.between.transm.var.cov.65, confidence=0.95)
Est.MAR.c.within.transm.65 <- ci(MAR.c.df.65$MAR.c.within.transm.var.cov.65, confidence = 0.95)

Est.MAR.c.age.male.70 <- ci(MAR.c.df.70$MAR.c.av.age.male.cov.70, confidence=0.95)
Est.MAR.c.age.diff.70 <- ci(MAR.c.df.70$MAR.c.av.age.diff.cov.70, confidence = 0.95)

Est.MAR.c.between.transm.70 <- ci(MAR.c.df.70$MAR.c.between.transm.var.cov.70, confidence=0.95)
Est.MAR.c.within.transm.70 <- ci(MAR.c.df.70$MAR.c.within.transm.var.cov.70, confidence = 0.95)

Est.MAR.c.age.male.75 <- ci(MAR.c.df.75$MAR.c.av.age.male.cov.75, confidence=0.95)
Est.MAR.c.age.diff.75 <- ci(MAR.c.df.75$MAR.c.av.age.diff.cov.75, confidence = 0.95)

Est.MAR.c.between.transm.75 <- ci(MAR.c.df.75$MAR.c.between.transm.var.cov.75, confidence=0.95)
Est.MAR.c.within.transm.75 <- ci(MAR.c.df.75$MAR.c.within.transm.var.cov.75, confidence = 0.95)

Est.MAR.c.age.male.80 <- ci(MAR.c.df.80$MAR.c.av.age.male.cov.80, confidence=0.95)
Est.MAR.c.age.diff.80 <- ci(MAR.c.df.80$MAR.c.av.age.diff.cov.80, confidence = 0.95)

Est.MAR.c.between.transm.80 <- ci(MAR.c.df.80$MAR.c.between.transm.var.cov.80, confidence=0.95)
Est.MAR.c.within.transm.80 <- ci(MAR.c.df.80$MAR.c.within.transm.var.cov.80, confidence = 0.95)

Est.MAR.c.age.male.85 <- ci(MAR.c.df.85$MAR.c.av.age.male.cov.85, confidence=0.95)
Est.MAR.c.age.diff.85 <- ci(MAR.c.df.85$MAR.c.av.age.diff.cov.85, confidence = 0.95)

Est.MAR.c.between.transm.85 <- ci(MAR.c.df.85$MAR.c.between.transm.var.cov.85, confidence=0.95)
Est.MAR.c.within.transm.85 <- ci(MAR.c.df.85$MAR.c.within.transm.var.cov.85, confidence = 0.95)

Est.MAR.c.age.male.90 <- ci(MAR.c.df.90$MAR.c.av.age.male.cov.90, confidence=0.95)
Est.MAR.c.age.diff.90 <- ci(MAR.c.df.90$MAR.c.av.age.diff.cov.90, confidence = 0.95)

Est.MAR.c.between.transm.90 <- ci(MAR.c.df.90$MAR.c.between.transm.var.cov.90, confidence=0.95)
Est.MAR.c.within.transm.90 <- ci(MAR.c.df.90$MAR.c.within.transm.var.cov.90, confidence = 0.95)

Est.MAR.c.age.male.95 <- ci(MAR.c.df.95$MAR.c.av.age.male.cov.95, confidence=0.95)
Est.MAR.c.age.diff.95 <- ci(MAR.c.df.95$MAR.c.av.age.diff.cov.95, confidence = 0.95)

Est.MAR.c.between.transm.95 <- ci(MAR.c.df.95$MAR.c.between.transm.var.cov.95, confidence=0.95)
Est.MAR.c.within.transm.95 <- ci(MAR.c.df.95$MAR.c.within.transm.var.cov.95, confidence = 0.95)

require(ggplot2)

# Age male - MCAR
df.age.male <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                          F = c(Est.age.male.pop[[1]], Est.MCAR.age.male.35[[1]], Est.MCAR.age.male.40[[1]], Est.MCAR.age.male.45[[1]],
                                Est.MCAR.age.male.50[[1]], Est.MCAR.age.male.55[[1]], Est.MCAR.age.male.60[[1]], Est.MCAR.age.male.65[[1]],
                                Est.MCAR.age.male.70[[1]], Est.MCAR.age.male.75[[1]], Est.MCAR.age.male.80[[1]], Est.MCAR.age.male.85[[1]],
                                Est.MCAR.age.male.90[[1]], Est.MCAR.age.male.95[[1]]),
                          L = c(Est.age.male.pop[[2]], Est.MCAR.age.male.35[[2]], Est.MCAR.age.male.40[[2]], Est.MCAR.age.male.45[[2]],
                                Est.MCAR.age.male.50[[2]], Est.MCAR.age.male.55[[2]], Est.MCAR.age.male.60[[2]], Est.MCAR.age.male.65[[2]],
                                Est.MCAR.age.male.70[[2]], Est.MCAR.age.male.75[[2]], Est.MCAR.age.male.80[[2]], Est.MCAR.age.male.85[[2]],
                                Est.MCAR.age.male.90[[2]], Est.MCAR.age.male.95[[2]]),
                          U = c(Est.age.male.pop[[3]], Est.MCAR.age.male.35[[3]], Est.MCAR.age.male.40[[3]], Est.MCAR.age.male.45[[3]],
                                Est.MCAR.age.male.50[[3]], Est.MCAR.age.male.55[[3]], Est.MCAR.age.male.60[[3]], Est.MCAR.age.male.65[[3]],
                                Est.MCAR.age.male.70[[3]], Est.MCAR.age.male.75[[3]], Est.MCAR.age.male.80[[3]], Est.MCAR.age.male.85[[3]],
                                Est.MCAR.age.male.90[[3]], Est.MCAR.age.male.95[[3]]))

# Age diff - MCAR
df.age.diff <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                          F = c(Est.age.diff.pop[[1]], Est.MCAR.age.diff.35[[1]], Est.MCAR.age.diff.40[[1]], Est.MCAR.age.diff.45[[1]],
                                Est.MCAR.age.diff.50[[1]], Est.MCAR.age.diff.55[[1]], Est.MCAR.age.diff.60[[1]], Est.MCAR.age.diff.65[[1]],
                                Est.MCAR.age.diff.70[[1]], Est.MCAR.age.diff.75[[1]], Est.MCAR.age.diff.80[[1]], Est.MCAR.age.diff.85[[1]],
                                Est.MCAR.age.diff.90[[1]], Est.MCAR.age.diff.95[[1]]),
                          L = c(Est.age.diff.pop[[2]], Est.MCAR.age.diff.35[[2]], Est.MCAR.age.diff.40[[2]], Est.MCAR.age.diff.45[[2]],
                                Est.MCAR.age.diff.50[[2]], Est.MCAR.age.diff.55[[2]], Est.MCAR.age.diff.60[[2]], Est.MCAR.age.diff.65[[2]],
                                Est.MCAR.age.diff.70[[2]], Est.MCAR.age.diff.75[[2]], Est.MCAR.age.diff.80[[2]], Est.MCAR.age.diff.85[[2]],
                                Est.MCAR.age.diff.90[[2]], Est.MCAR.age.diff.95[[2]]),
                          U = c(Est.age.diff.pop[[3]], Est.MCAR.age.diff.35[[3]], Est.MCAR.age.diff.40[[3]], Est.MCAR.age.diff.45[[3]],
                                Est.MCAR.age.diff.50[[3]], Est.MCAR.age.diff.55[[3]], Est.MCAR.age.diff.60[[3]], Est.MCAR.age.diff.65[[3]],
                                Est.MCAR.age.diff.70[[3]], Est.MCAR.age.diff.75[[3]], Est.MCAR.age.diff.80[[3]], Est.MCAR.age.diff.85[[3]],
                                Est.MCAR.age.diff.90[[3]], Est.MCAR.age.diff.95[[3]]))


# Age between groups - MCAR
df.age.between <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                             F = c(Est.between.transm.pop[[1]], Est.MCAR.between.transm.35[[1]], Est.MCAR.between.transm.40[[1]], Est.MCAR.between.transm.45[[1]],
                                   Est.MCAR.between.transm.50[[1]], Est.MCAR.between.transm.55[[1]], Est.MCAR.between.transm.60[[1]], Est.MCAR.between.transm.65[[1]],
                                   Est.MCAR.between.transm.70[[1]], Est.MCAR.between.transm.75[[1]], Est.MCAR.between.transm.80[[1]], Est.MCAR.between.transm.85[[1]],
                                   Est.MCAR.between.transm.90[[1]], Est.MCAR.between.transm.95[[1]]),
                             L = c(Est.between.transm.pop[[2]], Est.MCAR.between.transm.35[[2]], Est.MCAR.between.transm.40[[2]], Est.MCAR.between.transm.45[[2]],
                                   Est.MCAR.between.transm.50[[2]], Est.MCAR.between.transm.55[[2]], Est.MCAR.between.transm.60[[2]], Est.MCAR.between.transm.65[[2]],
                                   Est.MCAR.between.transm.70[[2]], Est.MCAR.between.transm.75[[2]], Est.MCAR.between.transm.80[[2]], Est.MCAR.between.transm.85[[2]],
                                   Est.MCAR.between.transm.90[[2]], Est.MCAR.between.transm.95[[2]]),
                             U = c(Est.between.transm.pop[[3]], Est.MCAR.between.transm.35[[3]], Est.MCAR.between.transm.40[[3]], Est.MCAR.between.transm.45[[3]],
                                   Est.MCAR.between.transm.50[[3]], Est.MCAR.between.transm.55[[3]], Est.MCAR.between.transm.60[[3]], Est.MCAR.between.transm.65[[3]],
                                   Est.MCAR.between.transm.70[[3]], Est.MCAR.between.transm.75[[3]], Est.MCAR.between.transm.80[[3]], Est.MCAR.between.transm.85[[3]],
                                   Est.MCAR.between.transm.90[[3]], Est.MCAR.between.transm.95[[3]]))

# Age within groups - MCAR
df.age.within <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                            F = c(Est.within.transm.pop[[1]], Est.MCAR.within.transm.35[[1]], Est.MCAR.within.transm.40[[1]], Est.MCAR.within.transm.45[[1]],
                                  Est.MCAR.within.transm.50[[1]], Est.MCAR.within.transm.55[[1]], Est.MCAR.within.transm.60[[1]], Est.MCAR.within.transm.65[[1]],
                                  Est.MCAR.within.transm.70[[1]], Est.MCAR.within.transm.75[[1]], Est.MCAR.within.transm.80[[1]], Est.MCAR.within.transm.85[[1]],
                                  Est.MCAR.within.transm.90[[1]], Est.MCAR.within.transm.95[[1]]),
                            L = c(Est.within.transm.pop[[2]], Est.MCAR.within.transm.35[[2]], Est.MCAR.within.transm.40[[2]], Est.MCAR.within.transm.45[[2]],
                                  Est.MCAR.within.transm.50[[2]], Est.MCAR.within.transm.55[[2]], Est.MCAR.within.transm.60[[2]], Est.MCAR.within.transm.65[[2]],
                                  Est.MCAR.within.transm.70[[2]], Est.MCAR.within.transm.75[[2]], Est.MCAR.within.transm.80[[2]], Est.MCAR.within.transm.85[[2]],
                                  Est.MCAR.within.transm.90[[2]], Est.MCAR.within.transm.95[[2]]),
                            U = c(Est.within.transm.pop[[3]], Est.MCAR.within.transm.35[[3]], Est.MCAR.within.transm.40[[3]], Est.MCAR.within.transm.45[[3]],
                                  Est.MCAR.within.transm.50[[3]], Est.MCAR.within.transm.55[[3]], Est.MCAR.within.transm.60[[3]], Est.MCAR.within.transm.65[[3]],
                                  Est.MCAR.within.transm.70[[3]], Est.MCAR.within.transm.75[[3]], Est.MCAR.within.transm.80[[3]], Est.MCAR.within.transm.85[[3]],
                                  Est.MCAR.within.transm.90[[3]], Est.MCAR.within.transm.95[[3]]))


ggplot(df.age.male, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) + 
  ggtitle("Population level average man age & Sequence coverage - MCAR") +
  xlab("Sequence coverage scenarios") + ylab("Average man age")

ggplot(df.age.diff, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Gender effect & Sequence coverage - MCAR") +
  xlab("Sequence coverage scenarios") + ylab("Age difference")


ggplot(df.age.between, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Between clusters age variation & Sequence coverage - MCAR") +
  xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")


ggplot(df.age.within, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Within clusters age variation & Sequence coverage - MCAR") +
  xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")





# Age male - MAR - a
df.age.male.MAR.a <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.male.pop[[1]], Est.MAR.a.age.male.35[[1]], Est.MAR.a.age.male.40[[1]], Est.MAR.a.age.male.45[[1]],
                                      Est.MAR.a.age.male.50[[1]], Est.MAR.a.age.male.55[[1]], Est.MAR.a.age.male.60[[1]], Est.MAR.a.age.male.65[[1]],
                                      Est.MAR.a.age.male.70[[1]], Est.MAR.a.age.male.75[[1]], Est.MAR.a.age.male.80[[1]], Est.MAR.a.age.male.85[[1]],
                                      Est.MAR.a.age.male.90[[1]], Est.MAR.a.age.male.95[[1]]),
                                L = c(Est.age.male.pop[[2]], Est.MAR.a.age.male.35[[2]], Est.MAR.a.age.male.40[[2]], Est.MAR.a.age.male.45[[2]],
                                      Est.MAR.a.age.male.50[[2]], Est.MAR.a.age.male.55[[2]], Est.MAR.a.age.male.60[[2]], Est.MAR.a.age.male.65[[2]],
                                      Est.MAR.a.age.male.70[[2]], Est.MAR.a.age.male.75[[2]], Est.MAR.a.age.male.80[[2]], Est.MAR.a.age.male.85[[2]],
                                      Est.MAR.a.age.male.90[[2]], Est.MAR.a.age.male.95[[2]]),
                                U = c(Est.age.male.pop[[3]], Est.MAR.a.age.male.35[[3]], Est.MAR.a.age.male.40[[3]], Est.MAR.a.age.male.45[[3]],
                                      Est.MAR.a.age.male.50[[3]], Est.MAR.a.age.male.55[[3]], Est.MAR.a.age.male.60[[3]], Est.MAR.a.age.male.65[[3]],
                                      Est.MAR.a.age.male.70[[3]], Est.MAR.a.age.male.75[[3]], Est.MAR.a.age.male.80[[3]], Est.MAR.a.age.male.85[[3]],
                                      Est.MAR.a.age.male.90[[3]], Est.MAR.a.age.male.95[[3]]))

# Age diff - MAR - a
df.age.diff.MAR.a <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.diff.pop[[1]], Est.MAR.a.age.diff.35[[1]], Est.MAR.a.age.diff.40[[1]], Est.MAR.a.age.diff.45[[1]],
                                      Est.MAR.a.age.diff.50[[1]], Est.MAR.a.age.diff.55[[1]], Est.MAR.a.age.diff.60[[1]], Est.MAR.a.age.diff.65[[1]],
                                      Est.MAR.a.age.diff.70[[1]], Est.MAR.a.age.diff.75[[1]], Est.MAR.a.age.diff.80[[1]], Est.MAR.a.age.diff.85[[1]],
                                      Est.MAR.a.age.diff.90[[1]], Est.MAR.a.age.diff.95[[1]]),
                                L = c(Est.age.diff.pop[[2]], Est.MAR.a.age.diff.35[[2]], Est.MAR.a.age.diff.40[[2]], Est.MAR.a.age.diff.45[[2]],
                                      Est.MAR.a.age.diff.50[[2]], Est.MAR.a.age.diff.55[[2]], Est.MAR.a.age.diff.60[[2]], Est.MAR.a.age.diff.65[[2]],
                                      Est.MAR.a.age.diff.70[[2]], Est.MAR.a.age.diff.75[[2]], Est.MAR.a.age.diff.80[[2]], Est.MAR.a.age.diff.85[[2]],
                                      Est.MAR.a.age.diff.90[[2]], Est.MAR.a.age.diff.95[[2]]),
                                U = c(Est.age.diff.pop[[3]], Est.MAR.a.age.diff.35[[3]], Est.MAR.a.age.diff.40[[3]], Est.MAR.a.age.diff.45[[3]],
                                      Est.MAR.a.age.diff.50[[3]], Est.MAR.a.age.diff.55[[3]], Est.MAR.a.age.diff.60[[3]], Est.MAR.a.age.diff.65[[3]],
                                      Est.MAR.a.age.diff.70[[3]], Est.MAR.a.age.diff.75[[3]], Est.MAR.a.age.diff.80[[3]], Est.MAR.a.age.diff.85[[3]],
                                      Est.MAR.a.age.diff.90[[3]], Est.MAR.a.age.diff.95[[3]]))


# Age between groups - MAR - a
df.age.between.MAR.a <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                   F = c(Est.between.transm.pop[[1]], Est.MAR.a.between.transm.35[[1]], Est.MAR.a.between.transm.40[[1]], Est.MAR.a.between.transm.45[[1]],
                                         Est.MAR.a.between.transm.50[[1]], Est.MAR.a.between.transm.55[[1]], Est.MAR.a.between.transm.60[[1]], Est.MAR.a.between.transm.65[[1]],
                                         Est.MAR.a.between.transm.70[[1]], Est.MAR.a.between.transm.75[[1]], Est.MAR.a.between.transm.80[[1]], Est.MAR.a.between.transm.85[[1]],
                                         Est.MAR.a.between.transm.90[[1]], Est.MAR.a.between.transm.95[[1]]),
                                   L = c(Est.between.transm.pop[[2]], Est.MAR.a.between.transm.35[[2]], Est.MAR.a.between.transm.40[[2]], Est.MAR.a.between.transm.45[[2]],
                                         Est.MAR.a.between.transm.50[[2]], Est.MAR.a.between.transm.55[[2]], Est.MAR.a.between.transm.60[[2]], Est.MAR.a.between.transm.65[[2]],
                                         Est.MAR.a.between.transm.70[[2]], Est.MAR.a.between.transm.75[[2]], Est.MAR.a.between.transm.80[[2]], Est.MAR.a.between.transm.85[[2]],
                                         Est.MAR.a.between.transm.90[[2]], Est.MAR.a.between.transm.95[[2]]),
                                   U = c(Est.between.transm.pop[[3]], Est.MAR.a.between.transm.35[[3]], Est.MAR.a.between.transm.40[[3]], Est.MAR.a.between.transm.45[[3]],
                                         Est.MAR.a.between.transm.50[[3]], Est.MAR.a.between.transm.55[[3]], Est.MAR.a.between.transm.60[[3]], Est.MAR.a.between.transm.65[[3]],
                                         Est.MAR.a.between.transm.70[[3]], Est.MAR.a.between.transm.75[[3]], Est.MAR.a.between.transm.80[[3]], Est.MAR.a.between.transm.85[[3]],
                                         Est.MAR.a.between.transm.90[[3]], Est.MAR.a.between.transm.95[[3]]))

# Age within groups - MAR - a
df.age.within.MAR.a <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                  F = c(Est.within.transm.pop[[1]], Est.MAR.a.within.transm.35[[1]], Est.MAR.a.within.transm.40[[1]], Est.MAR.a.within.transm.45[[1]],
                                        Est.MAR.a.within.transm.50[[1]], Est.MAR.a.within.transm.55[[1]], Est.MAR.a.within.transm.60[[1]], Est.MAR.a.within.transm.65[[1]],
                                        Est.MAR.a.within.transm.70[[1]], Est.MAR.a.within.transm.75[[1]], Est.MAR.a.within.transm.80[[1]], Est.MAR.a.within.transm.85[[1]],
                                        Est.MAR.a.within.transm.90[[1]], Est.MAR.a.within.transm.95[[1]]),
                                  L = c(Est.within.transm.pop[[2]], Est.MAR.a.within.transm.35[[2]], Est.MAR.a.within.transm.40[[2]], Est.MAR.a.within.transm.45[[2]],
                                        Est.MAR.a.within.transm.50[[2]], Est.MAR.a.within.transm.55[[2]], Est.MAR.a.within.transm.60[[2]], Est.MAR.a.within.transm.65[[2]],
                                        Est.MAR.a.within.transm.70[[2]], Est.MAR.a.within.transm.75[[2]], Est.MAR.a.within.transm.80[[2]], Est.MAR.a.within.transm.85[[2]],
                                        Est.MAR.a.within.transm.90[[2]], Est.MAR.a.within.transm.95[[2]]),
                                  U = c(Est.within.transm.pop[[3]], Est.MAR.a.within.transm.35[[3]], Est.MAR.a.within.transm.40[[3]], Est.MAR.a.within.transm.45[[3]],
                                        Est.MAR.a.within.transm.50[[3]], Est.MAR.a.within.transm.55[[3]], Est.MAR.a.within.transm.60[[3]], Est.MAR.a.within.transm.65[[3]],
                                        Est.MAR.a.within.transm.70[[3]], Est.MAR.a.within.transm.75[[3]], Est.MAR.a.within.transm.80[[3]], Est.MAR.a.within.transm.85[[3]],
                                        Est.MAR.a.within.transm.90[[3]], Est.MAR.a.within.transm.95[[3]]))




ggplot(df.age.male.MAR.a, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) + 
  ggtitle("Population level average man age & Sequence coverage - MAR (ratio = 0.7)") +
  xlab("Sequence coverage scenarios") + ylab("Average man age")

ggplot(df.age.diff.MAR.a, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Gender effect & Sequence coverage - MAR (ratio = 0.7)") +
  xlab("Sequence coverage scenarios") + ylab("Age difference")


ggplot(df.age.between.MAR.a, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Between clusters age variation & Sequence coverage - MAR (ratio = 0.7)") +
  xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")


ggplot(df.age.within.MAR.a, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Within clusters age variation & Sequence coverage - MAR (ratio = 0.7)") +
  xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")



# Age male - MAR - b
df.age.male.MAR.b <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.male.pop[[1]], Est.MAR.b.age.male.35[[1]], Est.MAR.b.age.male.40[[1]], Est.MAR.b.age.male.45[[1]],
                                      Est.MAR.b.age.male.50[[1]], Est.MAR.b.age.male.55[[1]], Est.MAR.b.age.male.60[[1]], Est.MAR.b.age.male.65[[1]],
                                      Est.MAR.b.age.male.70[[1]], Est.MAR.b.age.male.75[[1]], Est.MAR.b.age.male.80[[1]], Est.MAR.b.age.male.85[[1]],
                                      Est.MAR.b.age.male.90[[1]], Est.MAR.b.age.male.95[[1]]),
                                L = c(Est.age.male.pop[[2]], Est.MAR.b.age.male.35[[2]], Est.MAR.b.age.male.40[[2]], Est.MAR.b.age.male.45[[2]],
                                      Est.MAR.b.age.male.50[[2]], Est.MAR.b.age.male.55[[2]], Est.MAR.b.age.male.60[[2]], Est.MAR.b.age.male.65[[2]],
                                      Est.MAR.b.age.male.70[[2]], Est.MAR.b.age.male.75[[2]], Est.MAR.b.age.male.80[[2]], Est.MAR.b.age.male.85[[2]],
                                      Est.MAR.b.age.male.90[[2]], Est.MAR.b.age.male.95[[2]]),
                                U = c(Est.age.male.pop[[3]], Est.MAR.b.age.male.35[[3]], Est.MAR.b.age.male.40[[3]], Est.MAR.b.age.male.45[[3]],
                                      Est.MAR.b.age.male.50[[3]], Est.MAR.b.age.male.55[[3]], Est.MAR.b.age.male.60[[3]], Est.MAR.b.age.male.65[[3]],
                                      Est.MAR.b.age.male.70[[3]], Est.MAR.b.age.male.75[[3]], Est.MAR.b.age.male.80[[3]], Est.MAR.b.age.male.85[[3]],
                                      Est.MAR.b.age.male.90[[3]], Est.MAR.b.age.male.95[[3]]))

# Age diff - MAR - b
df.age.diff.MAR.b <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.diff.pop[[1]], Est.MAR.b.age.diff.35[[1]], Est.MAR.b.age.diff.40[[1]], Est.MAR.b.age.diff.45[[1]],
                                      Est.MAR.b.age.diff.50[[1]], Est.MAR.b.age.diff.55[[1]], Est.MAR.b.age.diff.60[[1]], Est.MAR.b.age.diff.65[[1]],
                                      Est.MAR.b.age.diff.70[[1]], Est.MAR.b.age.diff.75[[1]], Est.MAR.b.age.diff.80[[1]], Est.MAR.b.age.diff.85[[1]],
                                      Est.MAR.b.age.diff.90[[1]], Est.MAR.b.age.diff.95[[1]]),
                                L = c(Est.age.diff.pop[[2]], Est.MAR.b.age.diff.35[[2]], Est.MAR.b.age.diff.40[[2]], Est.MAR.b.age.diff.45[[2]],
                                      Est.MAR.b.age.diff.50[[2]], Est.MAR.b.age.diff.55[[2]], Est.MAR.b.age.diff.60[[2]], Est.MAR.b.age.diff.65[[2]],
                                      Est.MAR.b.age.diff.70[[2]], Est.MAR.b.age.diff.75[[2]], Est.MAR.b.age.diff.80[[2]], Est.MAR.b.age.diff.85[[2]],
                                      Est.MAR.b.age.diff.90[[2]], Est.MAR.b.age.diff.95[[2]]),
                                U = c(Est.age.diff.pop[[3]], Est.MAR.b.age.diff.35[[3]], Est.MAR.b.age.diff.40[[3]], Est.MAR.b.age.diff.45[[3]],
                                      Est.MAR.b.age.diff.50[[3]], Est.MAR.b.age.diff.55[[3]], Est.MAR.b.age.diff.60[[3]], Est.MAR.b.age.diff.65[[3]],
                                      Est.MAR.b.age.diff.70[[3]], Est.MAR.b.age.diff.75[[3]], Est.MAR.b.age.diff.80[[3]], Est.MAR.b.age.diff.85[[3]],
                                      Est.MAR.b.age.diff.90[[3]], Est.MAR.b.age.diff.95[[3]]))


# Age between groups - MAR - b
df.age.between.MAR.b <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                   F = c(Est.between.transm.pop[[1]], Est.MAR.b.between.transm.35[[1]], Est.MAR.b.between.transm.40[[1]], Est.MAR.b.between.transm.45[[1]],
                                         Est.MAR.b.between.transm.50[[1]], Est.MAR.b.between.transm.55[[1]], Est.MAR.b.between.transm.60[[1]], Est.MAR.b.between.transm.65[[1]],
                                         Est.MAR.b.between.transm.70[[1]], Est.MAR.b.between.transm.75[[1]], Est.MAR.b.between.transm.80[[1]], Est.MAR.b.between.transm.85[[1]],
                                         Est.MAR.b.between.transm.90[[1]], Est.MAR.b.between.transm.95[[1]]),
                                   L = c(Est.between.transm.pop[[2]], Est.MAR.b.between.transm.35[[2]], Est.MAR.b.between.transm.40[[2]], Est.MAR.b.between.transm.45[[2]],
                                         Est.MAR.b.between.transm.50[[2]], Est.MAR.b.between.transm.55[[2]], Est.MAR.b.between.transm.60[[2]], Est.MAR.b.between.transm.65[[2]],
                                         Est.MAR.b.between.transm.70[[2]], Est.MAR.b.between.transm.75[[2]], Est.MAR.b.between.transm.80[[2]], Est.MAR.b.between.transm.85[[2]],
                                         Est.MAR.b.between.transm.90[[2]], Est.MAR.b.between.transm.95[[2]]),
                                   U = c(Est.between.transm.pop[[3]], Est.MAR.b.between.transm.35[[3]], Est.MAR.b.between.transm.40[[3]], Est.MAR.b.between.transm.45[[3]],
                                         Est.MAR.b.between.transm.50[[3]], Est.MAR.b.between.transm.55[[3]], Est.MAR.b.between.transm.60[[3]], Est.MAR.b.between.transm.65[[3]],
                                         Est.MAR.b.between.transm.70[[3]], Est.MAR.b.between.transm.75[[3]], Est.MAR.b.between.transm.80[[3]], Est.MAR.b.between.transm.85[[3]],
                                         Est.MAR.b.between.transm.90[[3]], Est.MAR.b.between.transm.95[[3]]))

# Age within groups - MAR - b
df.age.within.MAR.b <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                  F = c(Est.within.transm.pop[[1]], Est.MAR.b.within.transm.35[[1]], Est.MAR.b.within.transm.40[[1]], Est.MAR.b.within.transm.45[[1]],
                                        Est.MAR.b.within.transm.50[[1]], Est.MAR.b.within.transm.55[[1]], Est.MAR.b.within.transm.60[[1]], Est.MAR.b.within.transm.65[[1]],
                                        Est.MAR.b.within.transm.70[[1]], Est.MAR.b.within.transm.75[[1]], Est.MAR.b.within.transm.80[[1]], Est.MAR.b.within.transm.85[[1]],
                                        Est.MAR.b.within.transm.90[[1]], Est.MAR.b.within.transm.95[[1]]),
                                  L = c(Est.within.transm.pop[[2]], Est.MAR.b.within.transm.35[[2]], Est.MAR.b.within.transm.40[[2]], Est.MAR.b.within.transm.45[[2]],
                                        Est.MAR.b.within.transm.50[[2]], Est.MAR.b.within.transm.55[[2]], Est.MAR.b.within.transm.60[[2]], Est.MAR.b.within.transm.65[[2]],
                                        Est.MAR.b.within.transm.70[[2]], Est.MAR.b.within.transm.75[[2]], Est.MAR.b.within.transm.80[[2]], Est.MAR.b.within.transm.85[[2]],
                                        Est.MAR.b.within.transm.90[[2]], Est.MAR.b.within.transm.95[[2]]),
                                  U = c(Est.within.transm.pop[[3]], Est.MAR.b.within.transm.35[[3]], Est.MAR.b.within.transm.40[[3]], Est.MAR.b.within.transm.45[[3]],
                                        Est.MAR.b.within.transm.50[[3]], Est.MAR.b.within.transm.55[[3]], Est.MAR.b.within.transm.60[[3]], Est.MAR.b.within.transm.65[[3]],
                                        Est.MAR.b.within.transm.70[[3]], Est.MAR.b.within.transm.75[[3]], Est.MAR.b.within.transm.80[[3]], Est.MAR.b.within.transm.85[[3]],
                                        Est.MAR.b.within.transm.90[[3]], Est.MAR.b.within.transm.95[[3]]))




ggplot(df.age.male.MAR.b, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) + 
  ggtitle("Population level average man age & Sequence coverage - MAR (ratio = 0.3)") +
  xlab("Sequence coverage scenarios") + ylab("Average man age")

ggplot(df.age.diff.MAR.b, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Gender effect & Sequence coverage - MAR (ratio = 0.3)") +
  xlab("Sequence coverage scenarios") + ylab("Age difference")


ggplot(df.age.between.MAR.b, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Between clusters age variation & Sequence coverage - MAR (ratio = 0.3)") +
  xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")


ggplot(df.age.within.MAR.b, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Within clusters age variation & Sequence coverage - MAR (ratio = 0.3)") +
  xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")



# Age male - MAR - c
df.age.male.MAR.c <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.male.pop[[1]], Est.MAR.c.age.male.35[[1]], Est.MAR.c.age.male.40[[1]], Est.MAR.c.age.male.45[[1]],
                                      Est.MAR.c.age.male.50[[1]], Est.MAR.c.age.male.55[[1]], Est.MAR.c.age.male.60[[1]], Est.MAR.c.age.male.65[[1]],
                                      Est.MAR.c.age.male.70[[1]], Est.MAR.c.age.male.75[[1]], Est.MAR.c.age.male.80[[1]], Est.MAR.c.age.male.85[[1]],
                                      Est.MAR.c.age.male.90[[1]], Est.MAR.c.age.male.95[[1]]),
                                L = c(Est.age.male.pop[[2]], Est.MAR.c.age.male.35[[2]], Est.MAR.c.age.male.40[[2]], Est.MAR.c.age.male.45[[2]],
                                      Est.MAR.c.age.male.50[[2]], Est.MAR.c.age.male.55[[2]], Est.MAR.c.age.male.60[[2]], Est.MAR.c.age.male.65[[2]],
                                      Est.MAR.c.age.male.70[[2]], Est.MAR.c.age.male.75[[2]], Est.MAR.c.age.male.80[[2]], Est.MAR.c.age.male.85[[2]],
                                      Est.MAR.c.age.male.90[[2]], Est.MAR.c.age.male.95[[2]]),
                                U = c(Est.age.male.pop[[3]], Est.MAR.c.age.male.35[[3]], Est.MAR.c.age.male.40[[3]], Est.MAR.c.age.male.45[[3]],
                                      Est.MAR.c.age.male.50[[3]], Est.MAR.c.age.male.55[[3]], Est.MAR.c.age.male.60[[3]], Est.MAR.c.age.male.65[[3]],
                                      Est.MAR.c.age.male.70[[3]], Est.MAR.c.age.male.75[[3]], Est.MAR.c.age.male.80[[3]], Est.MAR.c.age.male.85[[3]],
                                      Est.MAR.c.age.male.90[[3]], Est.MAR.c.age.male.95[[3]]))

# Age diff - MAR - c
df.age.diff.MAR.c <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                F = c(Est.age.diff.pop[[1]], Est.MAR.c.age.diff.35[[1]], Est.MAR.c.age.diff.40[[1]], Est.MAR.c.age.diff.45[[1]],
                                      Est.MAR.c.age.diff.50[[1]], Est.MAR.c.age.diff.55[[1]], Est.MAR.c.age.diff.60[[1]], Est.MAR.c.age.diff.65[[1]],
                                      Est.MAR.c.age.diff.70[[1]], Est.MAR.c.age.diff.75[[1]], Est.MAR.c.age.diff.80[[1]], Est.MAR.c.age.diff.85[[1]],
                                      Est.MAR.c.age.diff.90[[1]], Est.MAR.c.age.diff.95[[1]]),
                                L = c(Est.age.diff.pop[[2]], Est.MAR.c.age.diff.35[[2]], Est.MAR.c.age.diff.40[[2]], Est.MAR.c.age.diff.45[[2]],
                                      Est.MAR.c.age.diff.50[[2]], Est.MAR.c.age.diff.55[[2]], Est.MAR.c.age.diff.60[[2]], Est.MAR.c.age.diff.65[[2]],
                                      Est.MAR.c.age.diff.70[[2]], Est.MAR.c.age.diff.75[[2]], Est.MAR.c.age.diff.80[[2]], Est.MAR.c.age.diff.85[[2]],
                                      Est.MAR.c.age.diff.90[[2]], Est.MAR.c.age.diff.95[[2]]),
                                U = c(Est.age.diff.pop[[3]], Est.MAR.c.age.diff.35[[3]], Est.MAR.c.age.diff.40[[3]], Est.MAR.c.age.diff.45[[3]],
                                      Est.MAR.c.age.diff.50[[3]], Est.MAR.c.age.diff.55[[3]], Est.MAR.c.age.diff.60[[3]], Est.MAR.c.age.diff.65[[3]],
                                      Est.MAR.c.age.diff.70[[3]], Est.MAR.c.age.diff.75[[3]], Est.MAR.c.age.diff.80[[3]], Est.MAR.c.age.diff.85[[3]],
                                      Est.MAR.c.age.diff.90[[3]], Est.MAR.c.age.diff.95[[3]]))


# Age between groups - MAR - c
df.age.between.MAR.c <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                   F = c(Est.between.transm.pop[[1]], Est.MAR.c.between.transm.35[[1]], Est.MAR.c.between.transm.40[[1]], Est.MAR.c.between.transm.45[[1]],
                                         Est.MAR.c.between.transm.50[[1]], Est.MAR.c.between.transm.55[[1]], Est.MAR.c.between.transm.60[[1]], Est.MAR.c.between.transm.65[[1]],
                                         Est.MAR.c.between.transm.70[[1]], Est.MAR.c.between.transm.75[[1]], Est.MAR.c.between.transm.80[[1]], Est.MAR.c.between.transm.85[[1]],
                                         Est.MAR.c.between.transm.90[[1]], Est.MAR.c.between.transm.95[[1]]),
                                   L = c(Est.between.transm.pop[[2]], Est.MAR.c.between.transm.35[[2]], Est.MAR.c.between.transm.40[[2]], Est.MAR.c.between.transm.45[[2]],
                                         Est.MAR.c.between.transm.50[[2]], Est.MAR.c.between.transm.55[[2]], Est.MAR.c.between.transm.60[[2]], Est.MAR.c.between.transm.65[[2]],
                                         Est.MAR.c.between.transm.70[[2]], Est.MAR.c.between.transm.75[[2]], Est.MAR.c.between.transm.80[[2]], Est.MAR.c.between.transm.85[[2]],
                                         Est.MAR.c.between.transm.90[[2]], Est.MAR.c.between.transm.95[[2]]),
                                   U = c(Est.between.transm.pop[[3]], Est.MAR.c.between.transm.35[[3]], Est.MAR.c.between.transm.40[[3]], Est.MAR.c.between.transm.45[[3]],
                                         Est.MAR.c.between.transm.50[[3]], Est.MAR.c.between.transm.55[[3]], Est.MAR.c.between.transm.60[[3]], Est.MAR.c.between.transm.65[[3]],
                                         Est.MAR.c.between.transm.70[[3]], Est.MAR.c.between.transm.75[[3]], Est.MAR.c.between.transm.80[[3]], Est.MAR.c.between.transm.85[[3]],
                                         Est.MAR.c.between.transm.90[[3]], Est.MAR.c.between.transm.95[[3]]))

# Age within groups - MAR - c
df.age.within.MAR.c <- data.frame(x=c("pop", "cov.35", "cov.40", "cov.45", "cov.50", "cov.55", "cov.60", "cov.65", "cov.70", "cov.75", "cov.80", "cov.85", "cov.90", "cov.95"),
                                  F = c(Est.within.transm.pop[[1]], Est.MAR.c.within.transm.35[[1]], Est.MAR.c.within.transm.40[[1]], Est.MAR.c.within.transm.45[[1]],
                                        Est.MAR.c.within.transm.50[[1]], Est.MAR.c.within.transm.55[[1]], Est.MAR.c.within.transm.60[[1]], Est.MAR.c.within.transm.65[[1]],
                                        Est.MAR.c.within.transm.70[[1]], Est.MAR.c.within.transm.75[[1]], Est.MAR.c.within.transm.80[[1]], Est.MAR.c.within.transm.85[[1]],
                                        Est.MAR.c.within.transm.90[[1]], Est.MAR.c.within.transm.95[[1]]),
                                  L = c(Est.within.transm.pop[[2]], Est.MAR.c.within.transm.35[[2]], Est.MAR.c.within.transm.40[[2]], Est.MAR.c.within.transm.45[[2]],
                                        Est.MAR.c.within.transm.50[[2]], Est.MAR.c.within.transm.55[[2]], Est.MAR.c.within.transm.60[[2]], Est.MAR.c.within.transm.65[[2]],
                                        Est.MAR.c.within.transm.70[[2]], Est.MAR.c.within.transm.75[[2]], Est.MAR.c.within.transm.80[[2]], Est.MAR.c.within.transm.85[[2]],
                                        Est.MAR.c.within.transm.90[[2]], Est.MAR.c.within.transm.95[[2]]),
                                  U = c(Est.within.transm.pop[[3]], Est.MAR.c.within.transm.35[[3]], Est.MAR.c.within.transm.40[[3]], Est.MAR.c.within.transm.45[[3]],
                                        Est.MAR.c.within.transm.50[[3]], Est.MAR.c.within.transm.55[[3]], Est.MAR.c.within.transm.60[[3]], Est.MAR.c.within.transm.65[[3]],
                                        Est.MAR.c.within.transm.70[[3]], Est.MAR.c.within.transm.75[[3]], Est.MAR.c.within.transm.80[[3]], Est.MAR.c.within.transm.85[[3]],
                                        Est.MAR.c.within.transm.90[[3]], Est.MAR.c.within.transm.95[[3]]))




ggplot(df.age.male.MAR.c, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) + 
  ggtitle("Population level average man age & Sequence coverage - MAR (ratio = 0.5)") +
  xlab("Sequence coverage scenarios") + ylab("Average man age")

ggplot(df.age.diff.MAR.c, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Gender effect & Sequence coverage - MAR (ratio = 0.5)") +
  xlab("Sequence coverage scenarios") + ylab("Age difference")


ggplot(df.age.between.MAR.c, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Between clusters age variation & Sequence coverage - MAR (ratio = 0.5)") +
  xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")


ggplot(df.age.within.MAR.c, aes(x = x, y = F)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  ggtitle("Within clusters age variation & Sequence coverage - MAR (ratio = 0.5)") +
  xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")



# Plotting the input parameters of the calibrated model

# av.age.male.v <- run.ALL$av.age.male
# av.age.diff.v <- run.ALL$av.age.diff
# between.transm.var.v <- run.ALL$between.transm.var
# within.transm.var.v <- run.ALL$within.transm.var

plot(av.age.male.v,
     av.age.diff.v,
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,60),
     ylim = c(0, 20))
points(MCAR.av.age.male.cov.95,
       abs(MCAR.av.age.diff.cov.95),
       pch = 16,
       col = "blue2")
points(MCAR.av.age.male.cov.35,
       abs(MCAR.av.age.diff.cov.35),
       pch = 16,
       col = "orange")

##

plot(between.transm.var.v,
     within.transm.var.v,
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,10),
     ylim = c(0, 10))
points(MCAR.between.transm.var.cov.35,
       abs(MCAR.within.transm.var.cov.35),
       pch = 16,
       col = "blue2")
points(MCAR.between.transm.var.cov.35,
       abs(MCAR.within.transm.var.cov.35),
       pch = 16,
       col = "orange")



## Variables correlation
av.age.male.pop <- av.age.male.v
av.age.diff.pop <- av.age.diff.v
between.transm.var.pop <- between.transm.var.v
within.transm.var.pop <- within.transm.var.v

my.df.pop <- data.frame(av.age.male.pop, av.age.diff.pop, between.transm.var.pop, within.transm.var.pop)


# Correlation - MCAR 
age.male.df.MCAR <- data.frame(av.age.male.pop, 
                               MCAR.av.age.male.cov.35, MCAR.av.age.male.cov.40, MCAR.av.age.male.cov.45,
                               MCAR.av.age.male.cov.50, MCAR.av.age.male.cov.55, MCAR.av.age.male.cov.60,
                               MCAR.av.age.male.cov.65, MCAR.av.age.male.cov.70, MCAR.av.age.male.cov.75,
                               MCAR.av.age.male.cov.80, MCAR.av.age.male.cov.85, MCAR.av.age.male.cov.90,
                               MCAR.av.age.male.cov.95)

diff.df.MCAR <- data.frame(av.age.diff.pop,
                           MCAR.av.age.diff.cov.35, MCAR.av.age.diff.cov.40, MCAR.av.age.diff.cov.45,
                           MCAR.av.age.diff.cov.50, MCAR.av.age.diff.cov.55, MCAR.av.age.diff.cov.60,
                           MCAR.av.age.diff.cov.65, MCAR.av.age.diff.cov.70, MCAR.av.age.diff.cov.75,
                           MCAR.av.age.diff.cov.80, MCAR.av.age.diff.cov.85, MCAR.av.age.diff.cov.90,
                           MCAR.av.age.diff.cov.95)

between.df.MCAR <- data.frame(between.transm.var.pop,
                              MCAR.between.transm.var.cov.35, MCAR.between.transm.var.cov.40, MCAR.between.transm.var.cov.45,
                              MCAR.between.transm.var.cov.50, MCAR.between.transm.var.cov.55, MCAR.between.transm.var.cov.60,
                              MCAR.between.transm.var.cov.65, MCAR.between.transm.var.cov.70, MCAR.between.transm.var.cov.75,
                              MCAR.between.transm.var.cov.80, MCAR.between.transm.var.cov.85, MCAR.between.transm.var.cov.90,
                              MCAR.between.transm.var.cov.95)

within.df.MCAR <- data.frame(within.transm.var.pop,
                             MCAR.within.transm.var.cov.35, MCAR.within.transm.var.cov.40, MCAR.within.transm.var.cov.45,
                             MCAR.within.transm.var.cov.50, MCAR.within.transm.var.cov.55, MCAR.within.transm.var.cov.60,
                             MCAR.within.transm.var.cov.65, MCAR.within.transm.var.cov.70, MCAR.within.transm.var.cov.75,
                             MCAR.within.transm.var.cov.80, MCAR.within.transm.var.cov.85, MCAR.within.transm.var.cov.90,
                             MCAR.within.transm.var.cov.95)


# Correlation - MAR.a 
age.male.df.MAR.a <- data.frame(av.age.male.pop, 
                                MAR.a.av.age.male.cov.35, MAR.a.av.age.male.cov.40, MAR.a.av.age.male.cov.45,
                                MAR.a.av.age.male.cov.50, MAR.a.av.age.male.cov.55, MAR.a.av.age.male.cov.60,
                                MAR.a.av.age.male.cov.65, MAR.a.av.age.male.cov.70, MAR.a.av.age.male.cov.75,
                                MAR.a.av.age.male.cov.80, MAR.a.av.age.male.cov.85, MAR.a.av.age.male.cov.90,
                                MAR.a.av.age.male.cov.95)

diff.df.MAR.a <- data.frame(av.age.diff.pop,
                            MAR.a.av.age.diff.cov.35, MAR.a.av.age.diff.cov.40, MAR.a.av.age.diff.cov.45,
                            MAR.a.av.age.diff.cov.50, MAR.a.av.age.diff.cov.55, MAR.a.av.age.diff.cov.60,
                            MAR.a.av.age.diff.cov.65, MAR.a.av.age.diff.cov.70, MAR.a.av.age.diff.cov.75,
                            MAR.a.av.age.diff.cov.80, MAR.a.av.age.diff.cov.85, MAR.a.av.age.diff.cov.90,
                            MAR.a.av.age.diff.cov.95)

between.df.MAR.a <- data.frame(between.transm.var.pop,
                               MAR.a.between.transm.var.cov.35, MAR.a.between.transm.var.cov.40, MAR.a.between.transm.var.cov.45,
                               MAR.a.between.transm.var.cov.50, MAR.a.between.transm.var.cov.55, MAR.a.between.transm.var.cov.60,
                               MAR.a.between.transm.var.cov.65, MAR.a.between.transm.var.cov.70, MAR.a.between.transm.var.cov.75,
                               MAR.a.between.transm.var.cov.80, MAR.a.between.transm.var.cov.85, MAR.a.between.transm.var.cov.90,
                               MAR.a.between.transm.var.cov.95)

within.df.MAR.a <- data.frame(within.transm.var.pop,
                              MAR.a.within.transm.var.cov.35, MAR.a.within.transm.var.cov.40, MAR.a.within.transm.var.cov.45,
                              MAR.a.within.transm.var.cov.50, MAR.a.within.transm.var.cov.55, MAR.a.within.transm.var.cov.60,
                              MAR.a.within.transm.var.cov.65, MAR.a.within.transm.var.cov.70, MAR.a.within.transm.var.cov.75,
                              MAR.a.within.transm.var.cov.80, MAR.a.within.transm.var.cov.85, MAR.a.within.transm.var.cov.90,
                              MAR.a.within.transm.var.cov.95)



# Correlation - MAR.b 
age.male.df.MAR.b <- data.frame(av.age.male.pop, 
                                MAR.b.av.age.male.cov.35, MAR.b.av.age.male.cov.40, MAR.b.av.age.male.cov.45,
                                MAR.b.av.age.male.cov.50, MAR.b.av.age.male.cov.55, MAR.b.av.age.male.cov.60,
                                MAR.b.av.age.male.cov.65, MAR.b.av.age.male.cov.70, MAR.b.av.age.male.cov.75,
                                MAR.b.av.age.male.cov.80, MAR.b.av.age.male.cov.85, MAR.b.av.age.male.cov.90,
                                MAR.b.av.age.male.cov.95)

diff.df.MAR.b <- data.frame(av.age.diff.pop,
                            MAR.b.av.age.diff.cov.35, MAR.b.av.age.diff.cov.40, MAR.b.av.age.diff.cov.45,
                            MAR.b.av.age.diff.cov.50, MAR.b.av.age.diff.cov.55, MAR.b.av.age.diff.cov.60,
                            MAR.b.av.age.diff.cov.65, MAR.b.av.age.diff.cov.70, MAR.b.av.age.diff.cov.75,
                            MAR.b.av.age.diff.cov.80, MAR.b.av.age.diff.cov.85, MAR.b.av.age.diff.cov.90,
                            MAR.b.av.age.diff.cov.95)

between.df.MAR.b <- data.frame(between.transm.var.pop,
                               MAR.b.between.transm.var.cov.35, MAR.b.between.transm.var.cov.40, MAR.b.between.transm.var.cov.45,
                               MAR.b.between.transm.var.cov.50, MAR.b.between.transm.var.cov.55, MAR.b.between.transm.var.cov.60,
                               MAR.b.between.transm.var.cov.65, MAR.b.between.transm.var.cov.70, MAR.b.between.transm.var.cov.75,
                               MAR.b.between.transm.var.cov.80, MAR.b.between.transm.var.cov.85, MAR.b.between.transm.var.cov.90,
                               MAR.b.between.transm.var.cov.95)

within.df.MAR.b <- data.frame(within.transm.var.pop,
                              MAR.b.within.transm.var.cov.35, MAR.b.within.transm.var.cov.40, MAR.b.within.transm.var.cov.45,
                              MAR.b.within.transm.var.cov.50, MAR.b.within.transm.var.cov.55, MAR.b.within.transm.var.cov.60,
                              MAR.b.within.transm.var.cov.65, MAR.b.within.transm.var.cov.70, MAR.b.within.transm.var.cov.75,
                              MAR.b.within.transm.var.cov.80, MAR.b.within.transm.var.cov.85, MAR.b.within.transm.var.cov.90,
                              MAR.b.within.transm.var.cov.95)




# Correlation - MAR.c 
age.male.df.MAR.c <- data.frame(av.age.male.pop, 
                                MAR.c.av.age.male.cov.35, MAR.c.av.age.male.cov.40, MAR.c.av.age.male.cov.45,
                                MAR.c.av.age.male.cov.50, MAR.c.av.age.male.cov.55, MAR.c.av.age.male.cov.60,
                                MAR.c.av.age.male.cov.65, MAR.c.av.age.male.cov.70, MAR.c.av.age.male.cov.75,
                                MAR.c.av.age.male.cov.80, MAR.c.av.age.male.cov.85, MAR.c.av.age.male.cov.90,
                                MAR.c.av.age.male.cov.95)

diff.df.MAR.c <- data.frame(av.age.diff.pop,
                            MAR.c.av.age.diff.cov.35, MAR.c.av.age.diff.cov.40, MAR.c.av.age.diff.cov.45,
                            MAR.c.av.age.diff.cov.50, MAR.c.av.age.diff.cov.55, MAR.c.av.age.diff.cov.60,
                            MAR.c.av.age.diff.cov.65, MAR.c.av.age.diff.cov.70, MAR.c.av.age.diff.cov.75,
                            MAR.c.av.age.diff.cov.80, MAR.c.av.age.diff.cov.85, MAR.c.av.age.diff.cov.90,
                            MAR.c.av.age.diff.cov.95)

between.df.MAR.c <- data.frame(between.transm.var.pop,
                               MAR.c.between.transm.var.cov.35, MAR.c.between.transm.var.cov.40, MAR.c.between.transm.var.cov.45,
                               MAR.c.between.transm.var.cov.50, MAR.c.between.transm.var.cov.55, MAR.c.between.transm.var.cov.60,
                               MAR.c.between.transm.var.cov.65, MAR.c.between.transm.var.cov.70, MAR.c.between.transm.var.cov.75,
                               MAR.c.between.transm.var.cov.80, MAR.c.between.transm.var.cov.85, MAR.c.between.transm.var.cov.90,
                               MAR.c.between.transm.var.cov.95)

within.df.MAR.c <- data.frame(within.transm.var.pop,
                              MAR.c.within.transm.var.cov.35, MAR.c.within.transm.var.cov.40, MAR.c.within.transm.var.cov.45,
                              MAR.c.within.transm.var.cov.50, MAR.c.within.transm.var.cov.55, MAR.c.within.transm.var.cov.60,
                              MAR.c.within.transm.var.cov.65, MAR.c.within.transm.var.cov.70, MAR.c.within.transm.var.cov.75,
                              MAR.c.within.transm.var.cov.80, MAR.c.within.transm.var.cov.85, MAR.c.within.transm.var.cov.90,
                              MAR.c.within.transm.var.cov.95)


# Compute a correlation matrix

# c(age.male.df.MCAR, diff.df.MCAR, between.df.MCAR, within.df.MCAR,
#   age.male.df.MAR.a, diff.df.MAR.a, between.df.MAR.a, within.df.MAR.a,
#   age.male.df.MAR.b, diff.df.MAR.b, between.df.MAR.b, within.df.MAR.b,
#   age.male.df.MAR.c, diff.df.MAR.c, between.df.MAR.c, within.df.MAR.c)

library(ggcorrplot)

dataset.raw.num.complete <- within.df.MCAR 


corr <- round(cor(dataset.raw.num.complete), 6)

ggcorrplot(corr)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)


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

