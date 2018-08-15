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

