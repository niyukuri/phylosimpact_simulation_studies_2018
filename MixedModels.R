# Linear models	 and	 linear	mixed	effects	models	in	R	with	linguistic	applications	

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

female.age <- c(17:24, 32, 23, 21, 29, 19, 21, 23, 24, 18, 19, 19, 21, 17:24, 21:27, 19)
male.age <- c(19, 21, 23, 24, 18, 19, 19, 21, 17:24, 21:27)


female.df <- data.frame(age=female.age)
male.df <- data.frame(age=male.age)

female.df$gender <- "female"
male.df$gender <- "male"

ageTransmClust <- rbind(female.df, male.df)

ggplot(ageTransmClust, aes(age, fill = gender)) + geom_density(alpha = 0.2)


library(ggplot2)

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

traclust <- lme(age ~ val, data = AGEMIX, random = ~ 1|id)

library(lme4)

traclust2 <- lmer(age ~ val +  (1|id), data = AGEMIX)

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

a <- coef(summary(traclust))[1] # 

beta <- coef(summary(traclust))[2] # or fixed.effects(summary(traclust))[2] 
# getME(traclust2,"beta")[2] < this is bridge.width 

# b <- getME(traclust2,"sigma") # between cluster var

b1 <- as.numeric(VarCorr(traclust)[3]) # between cluster variation

b2 <- as.numeric(VarCorr(traclust)[4]) # within cluster variation



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
