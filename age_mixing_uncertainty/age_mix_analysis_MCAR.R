

### Out together simulated data for age-mix in transmission

library(data.table)
library(dplyr)
library(gmodels)
require(ggplot2)


run1 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.NonAge.mix.run1.csv")
run2 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.NonAge.mix.run2.csv")
run5 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.NonAge.mix.run5.csv")
run6 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.NonAge.mix.run6.csv")

NonAgeMix.run.ALL <- rbind(run1, run2, run5, run6)

NonAgeMix.run.ALL <- dplyr::filter(NonAgeMix.run.ALL, NonAgeMix.run.ALL$flag.lme=="0") # select simulations with Non-Age mixing

NonAgeMix.run.ALL <- na.omit(NonAgeMix.run.ALL)

# NonAgeMix.run.ALL <- abs(NonAgeMix.run.ALL) # to make sure all values are positive - age gap

# True: population level mixing characteristics

Pop.mean.AD.v <- NonAgeMix.run.ALL$Pop.mean.AD
Pop.med.AD.v <- NonAgeMix.run.ALL$Pop.med.AD
Pop.sd.AD.v <- NonAgeMix.run.ALL$Pop.sd.AD
av.age.male.v <- NonAgeMix.run.ALL$av.age.male
av.age.diff.v <- NonAgeMix.run.ALL$av.age.diff
between.transm.var.v <- NonAgeMix.run.ALL$between.transm.var
within.transm.var.v <- NonAgeMix.run.ALL$within.transm.var

# "clust.AR.c.95.av.age.male"      
# [217] "clust.AR.c.95.av.age.diff"       "clust.AR.c.95.between.clust.var" "clust.AR.c.95.within.clust.var
# 

# MCAR

# Cov 35
MCAR.av.age.male.cov.35 <- NonAgeMix.run.ALL[,12]
MCAR.av.age.diff.cov.35 <- NonAgeMix.run.ALL[,13]
MCAR.between.transm.var.cov.35 <- NonAgeMix.run.ALL[,14]
MCAR.within.transm.var.cov.35 <- NonAgeMix.run.ALL[,15]


# Cov 40
MCAR.av.age.male.cov.40 <- NonAgeMix.run.ALL[,16]
MCAR.av.age.diff.cov.40 <- NonAgeMix.run.ALL[,17]
MCAR.between.transm.var.cov.40 <- NonAgeMix.run.ALL[,18]
MCAR.within.transm.var.cov.40 <- NonAgeMix.run.ALL[,19]


# Cov 45
MCAR.av.age.male.cov.45 <- NonAgeMix.run.ALL[,20]
MCAR.av.age.diff.cov.45 <- NonAgeMix.run.ALL[,21]
MCAR.between.transm.var.cov.45 <- NonAgeMix.run.ALL[,22]
MCAR.within.transm.var.cov.45 <- NonAgeMix.run.ALL[,23]


# Cov 50
MCAR.av.age.male.cov.50 <- NonAgeMix.run.ALL[,24]
MCAR.av.age.diff.cov.50 <- NonAgeMix.run.ALL[,25]
MCAR.between.transm.var.cov.50 <- NonAgeMix.run.ALL[,26]
MCAR.within.transm.var.cov.50 <- NonAgeMix.run.ALL[,27]


# Cov 55
MCAR.av.age.male.cov.55 <- NonAgeMix.run.ALL[,28]
MCAR.av.age.diff.cov.55 <- NonAgeMix.run.ALL[,29]
MCAR.between.transm.var.cov.55 <- NonAgeMix.run.ALL[,30]
MCAR.within.transm.var.cov.55 <- NonAgeMix.run.ALL[,31]

# Cov 60
MCAR.av.age.male.cov.60 <- NonAgeMix.run.ALL[,32]
MCAR.av.age.diff.cov.60 <- NonAgeMix.run.ALL[,33]
MCAR.between.transm.var.cov.60 <- NonAgeMix.run.ALL[,34]
MCAR.within.transm.var.cov.60 <- NonAgeMix.run.ALL[,35]


# Cov 65
MCAR.av.age.male.cov.65 <- NonAgeMix.run.ALL[,36]
MCAR.av.age.diff.cov.65 <- NonAgeMix.run.ALL[,37]
MCAR.between.transm.var.cov.65 <- NonAgeMix.run.ALL[,38]
MCAR.within.transm.var.cov.65 <- NonAgeMix.run.ALL[,39]


# Cov 70
MCAR.av.age.male.cov.70 <- NonAgeMix.run.ALL[,40]
MCAR.av.age.diff.cov.70 <- NonAgeMix.run.ALL[,41]
MCAR.between.transm.var.cov.70 <- NonAgeMix.run.ALL[,42]
MCAR.within.transm.var.cov.70 <- NonAgeMix.run.ALL[,43]


# Cov 75
MCAR.av.age.male.cov.75 <- NonAgeMix.run.ALL[,44]
MCAR.av.age.diff.cov.75 <- NonAgeMix.run.ALL[,45]
MCAR.between.transm.var.cov.75 <- NonAgeMix.run.ALL[,46]
MCAR.within.transm.var.cov.75 <- NonAgeMix.run.ALL[,47]

# Cov 80
MCAR.av.age.male.cov.80 <- NonAgeMix.run.ALL[,48]
MCAR.av.age.diff.cov.80 <- NonAgeMix.run.ALL[,49]
MCAR.between.transm.var.cov.80 <- NonAgeMix.run.ALL[,50]
MCAR.within.transm.var.cov.80 <- NonAgeMix.run.ALL[,51]


# Cov 85
MCAR.av.age.male.cov.85 <- NonAgeMix.run.ALL[,52]
MCAR.av.age.diff.cov.85 <- NonAgeMix.run.ALL[,53]
MCAR.between.transm.var.cov.85 <- NonAgeMix.run.ALL[,54]
MCAR.within.transm.var.cov.85 <- NonAgeMix.run.ALL[,55]

# Cov 90
MCAR.av.age.male.cov.90 <- NonAgeMix.run.ALL[,56]
MCAR.av.age.diff.cov.90 <- NonAgeMix.run.ALL[,57]
MCAR.between.transm.var.cov.90 <- NonAgeMix.run.ALL[,58]
MCAR.within.transm.var.cov.90 <- NonAgeMix.run.ALL[,59]


# Cov 95
MCAR.av.age.male.cov.95 <- NonAgeMix.run.ALL[,60]
MCAR.av.age.diff.cov.95 <- NonAgeMix.run.ALL[,61]
MCAR.between.transm.var.cov.95 <- NonAgeMix.run.ALL[,62]
MCAR.within.transm.var.cov.95 <- NonAgeMix.run.ALL[,63]



# Cov 100
MCAR.av.age.male.cov.100 <- NonAgeMix.run.ALL[,220]
MCAR.av.age.diff.cov.100 <- NonAgeMix.run.ALL[,221]
MCAR.between.transm.var.cov.100 <- NonAgeMix.run.ALL[,222]
MCAR.within.transm.var.cov.100 <- NonAgeMix.run.ALL[,223]




# Comparisons between population-level stats and sequence coverage scenarios

# Population level
av.age.male.pop <- MCAR.av.age.male.cov.100 # av.age.male.v
av.age.diff.pop <- MCAR.av.age.diff.cov.100 # av.age.diff.v
between.transm.var.pop <- MCAR.between.transm.var.cov.100 # between.transm.var.v
within.transm.var.pop <- MCAR.within.transm.var.cov.100 # within.transm.var.v

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




# Age male - MCAR
df.age.male.NonAgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
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
df.age.diff.NonAgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
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
df.age.between.NonAgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
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
df.age.within.NonAgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
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


# ggplot(df.age.male.NonAgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) + 
#   ggtitle("Population level average man age & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Average man age")
# 
# ggplot(df.age.diff.NonAgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Gender effect & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Age difference")
# 
# 
# ggplot(df.age.between.NonAgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Between clusters age variation & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")
# 
# 
# ggplot(df.age.within.NonAgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Within clusters age variation & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")
# 




##

runAM2 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.Age.mix.run2.csv")
runAM5 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.Age.mix.run5.csv")
runAM8 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.Age.mix.run8.csv")
runAM9 <- read.csv("~/Dropbox/2.5.2018.Simpact/age_mixing_uncertainty/LMEM.master.Age.mix.run9.csv")

AgeMix.run.ALL <- rbind(runAM2, runAM5, runAM8, runAM9)

AgeMix.run.ALL <- dplyr::filter(AgeMix.run.ALL, AgeMix.run.ALL$flag.lme=="1") # select simulations with Age mixing

AgeMix.run.ALL <- na.omit(AgeMix.run.ALL)

# AgeMix.run.ALL <- abs(AgeMix.run.ALL) # to make sure all values are positive - age gap

# True: population level mixing characteristics

Pop.mean.AD.v.AM <- AgeMix.run.ALL$Pop.mean.AD
Pop.med.AD.v.AM <- AgeMix.run.ALL$Pop.med.AD
Pop.sd.AD.v.AM <- AgeMix.run.ALL$Pop.sd.AD
av.age.male.v.AM <- AgeMix.run.ALL$av.age.male
av.age.diff.v.AM <- AgeMix.run.ALL$av.age.diff
between.transm.var.v.AM <- AgeMix.run.ALL$between.transm.var
within.transm.var.v.AM <- AgeMix.run.ALL$within.transm.var

# "clust.AR.c.95.av.age.male"      
# [217] "clust.AR.c.95.av.age.diff"       "clust.AR.c.95.between.clust.var" "clust.AR.c.95.within.clust.var
# 

# MCAR

# Cov 35
AM.MCAR.av.age.male.cov.35 <- AgeMix.run.ALL[,12]
AM.MCAR.av.age.diff.cov.35 <- AgeMix.run.ALL[,13]
AM.MCAR.between.transm.var.cov.35 <- AgeMix.run.ALL[,14]
AM.MCAR.within.transm.var.cov.35 <- AgeMix.run.ALL[,15]


# Cov 40
AM.MCAR.av.age.male.cov.40 <- AgeMix.run.ALL[,16]
AM.MCAR.av.age.diff.cov.40 <- AgeMix.run.ALL[,17]
AM.MCAR.between.transm.var.cov.40 <- AgeMix.run.ALL[,18]
AM.MCAR.within.transm.var.cov.40 <- AgeMix.run.ALL[,19]


# Cov 45
AM.MCAR.av.age.male.cov.45 <- AgeMix.run.ALL[,20]
AM.MCAR.av.age.diff.cov.45 <- AgeMix.run.ALL[,21]
AM.MCAR.between.transm.var.cov.45 <- AgeMix.run.ALL[,22]
AM.MCAR.within.transm.var.cov.45 <- AgeMix.run.ALL[,23]


# Cov 50
AM.MCAR.av.age.male.cov.50 <- AgeMix.run.ALL[,24]
AM.MCAR.av.age.diff.cov.50 <- AgeMix.run.ALL[,25]
AM.MCAR.between.transm.var.cov.50 <- AgeMix.run.ALL[,26]
AM.MCAR.within.transm.var.cov.50 <- AgeMix.run.ALL[,27]


# Cov 55
AM.MCAR.av.age.male.cov.55 <- AgeMix.run.ALL[,28]
AM.MCAR.av.age.diff.cov.55 <- AgeMix.run.ALL[,29]
AM.MCAR.between.transm.var.cov.55 <- AgeMix.run.ALL[,30]
AM.MCAR.within.transm.var.cov.55 <- AgeMix.run.ALL[,31]

# Cov 60
AM.MCAR.av.age.male.cov.60 <- AgeMix.run.ALL[,32]
AM.MCAR.av.age.diff.cov.60 <- AgeMix.run.ALL[,33]
AM.MCAR.between.transm.var.cov.60 <- AgeMix.run.ALL[,34]
AM.MCAR.within.transm.var.cov.60 <- AgeMix.run.ALL[,35]

# Cov 65
AM.MCAR.av.age.male.cov.65 <- AgeMix.run.ALL[,36]
AM.MCAR.av.age.diff.cov.65 <- AgeMix.run.ALL[,37]
AM.MCAR.between.transm.var.cov.65 <- AgeMix.run.ALL[,38]
AM.MCAR.within.transm.var.cov.65 <- AgeMix.run.ALL[,39]


# Cov 70
AM.MCAR.av.age.male.cov.70 <- AgeMix.run.ALL[,40]
AM.MCAR.av.age.diff.cov.70 <- AgeMix.run.ALL[,41]
AM.MCAR.between.transm.var.cov.70 <- AgeMix.run.ALL[,42]
AM.MCAR.within.transm.var.cov.70 <- AgeMix.run.ALL[,43]


# Cov 75
AM.MCAR.av.age.male.cov.75 <- AgeMix.run.ALL[,44]
AM.MCAR.av.age.diff.cov.75 <- AgeMix.run.ALL[,45]
AM.MCAR.between.transm.var.cov.75 <- AgeMix.run.ALL[,46]
AM.MCAR.within.transm.var.cov.75 <- AgeMix.run.ALL[,47]

# Cov 80
AM.MCAR.av.age.male.cov.80 <- AgeMix.run.ALL[,48]
AM.MCAR.av.age.diff.cov.80 <- AgeMix.run.ALL[,49]
AM.MCAR.between.transm.var.cov.80 <- AgeMix.run.ALL[,50]
AM.MCAR.within.transm.var.cov.80 <- AgeMix.run.ALL[,51]


# Cov 85
AM.MCAR.av.age.male.cov.85 <- AgeMix.run.ALL[,52]
AM.MCAR.av.age.diff.cov.85 <- AgeMix.run.ALL[,53]
AM.MCAR.between.transm.var.cov.85 <- AgeMix.run.ALL[,54]
AM.MCAR.within.transm.var.cov.85 <- AgeMix.run.ALL[,55]

# Cov 90
AM.MCAR.av.age.male.cov.90 <- AgeMix.run.ALL[,56]
AM.MCAR.av.age.diff.cov.90 <- AgeMix.run.ALL[,57]
AM.MCAR.between.transm.var.cov.90 <- AgeMix.run.ALL[,58]
AM.MCAR.within.transm.var.cov.90 <- AgeMix.run.ALL[,59]


# Cov 95
AM.MCAR.av.age.male.cov.95 <- AgeMix.run.ALL[,60]
AM.MCAR.av.age.diff.cov.95 <- AgeMix.run.ALL[,61]
AM.MCAR.between.transm.var.cov.95 <- AgeMix.run.ALL[,62]
AM.MCAR.within.transm.var.cov.95 <- AgeMix.run.ALL[,63]


# Cov 100
AM.MCAR.av.age.male.cov.100 <- AgeMix.run.ALL[,220]
AM.MCAR.av.age.diff.cov.100 <- AgeMix.run.ALL[,221]
AM.MCAR.between.transm.var.cov.100 <- AgeMix.run.ALL[,222]
AM.MCAR.within.transm.var.cov.100 <- AgeMix.run.ALL[,223]



# Comparisons between population-level stats and sequence coverage scenarios

# Population level
AM.av.age.male.pop <- AM.MCAR.av.age.male.cov.100 # av.age.male.v
AM.av.age.diff.pop <- AM.MCAR.av.age.diff.cov.100 # av.age.diff.v
AM.between.transm.var.pop <- AM.MCAR.between.transm.var.cov.100 # between.transm.var.v
AM.within.transm.var.pop <- AM.MCAR.within.transm.var.cov.100 # within.transm.var.v

AM.my.df.pop <- data.frame(AM.av.age.male.pop, AM.av.age.diff.pop, AM.between.transm.var.pop, AM.within.transm.var.pop)





# MCAR

AM.MCAR.df.35 <- data.frame(AM.MCAR.av.age.male.cov.35, AM.MCAR.av.age.diff.cov.35, AM.MCAR.between.transm.var.cov.35, AM.MCAR.within.transm.var.cov.35)
AM.MCAR.df.40 <- data.frame(AM.MCAR.av.age.male.cov.40, AM.MCAR.av.age.diff.cov.40, AM.MCAR.between.transm.var.cov.40, AM.MCAR.within.transm.var.cov.40)
AM.MCAR.df.45 <- data.frame(AM.MCAR.av.age.male.cov.45, AM.MCAR.av.age.diff.cov.45, AM.MCAR.between.transm.var.cov.45, AM.MCAR.within.transm.var.cov.45)
AM.MCAR.df.50 <- data.frame(AM.MCAR.av.age.male.cov.50, AM.MCAR.av.age.diff.cov.50, AM.MCAR.between.transm.var.cov.50, AM.MCAR.within.transm.var.cov.50)
AM.MCAR.df.55 <- data.frame(AM.MCAR.av.age.male.cov.55, AM.MCAR.av.age.diff.cov.55, AM.MCAR.between.transm.var.cov.55, AM.MCAR.within.transm.var.cov.55)
AM.MCAR.df.60 <- data.frame(AM.MCAR.av.age.male.cov.60, AM.MCAR.av.age.diff.cov.60, AM.MCAR.between.transm.var.cov.60, AM.MCAR.within.transm.var.cov.60)
AM.MCAR.df.65 <- data.frame(AM.MCAR.av.age.male.cov.65, AM.MCAR.av.age.diff.cov.65, AM.MCAR.between.transm.var.cov.65, AM.MCAR.within.transm.var.cov.65)
AM.MCAR.df.70 <- data.frame(AM.MCAR.av.age.male.cov.70, AM.MCAR.av.age.diff.cov.70, AM.MCAR.between.transm.var.cov.70, AM.MCAR.within.transm.var.cov.70)
AM.MCAR.df.75 <- data.frame(AM.MCAR.av.age.male.cov.75, AM.MCAR.av.age.diff.cov.75, AM.MCAR.between.transm.var.cov.75, AM.MCAR.within.transm.var.cov.75)
AM.MCAR.df.80 <- data.frame(AM.MCAR.av.age.male.cov.80, AM.MCAR.av.age.diff.cov.80, AM.MCAR.between.transm.var.cov.80, AM.MCAR.within.transm.var.cov.80)
AM.MCAR.df.85 <- data.frame(AM.MCAR.av.age.male.cov.85, AM.MCAR.av.age.diff.cov.85, AM.MCAR.between.transm.var.cov.85, AM.MCAR.within.transm.var.cov.85)
AM.MCAR.df.90 <- data.frame(AM.MCAR.av.age.male.cov.90, AM.MCAR.av.age.diff.cov.90, AM.MCAR.between.transm.var.cov.90, AM.MCAR.within.transm.var.cov.90)
AM.MCAR.df.95 <- data.frame(AM.MCAR.av.age.male.cov.95, AM.MCAR.av.age.diff.cov.95, AM.MCAR.between.transm.var.cov.95, AM.MCAR.within.transm.var.cov.95)

# Population level
AM.Est.age.male.pop <- ci(AM.my.df.pop$AM.av.age.male.pop, confidence=0.95)
AM.Est.age.diff.pop <- ci(AM.my.df.pop$AM.av.age.diff.pop, confidence = 0.95)

AM.Est.between.transm.pop <- ci(AM.my.df.pop$AM.between.transm.var.pop, confidence=0.95)
AM.Est.within.transm.pop <- ci(AM.my.df.pop$AM.within.transm.var.pop, confidence = 0.95)



# MCAR

AM.Est.MCAR.age.male.35 <- ci(AM.MCAR.df.35$AM.MCAR.av.age.male.cov.35, confidence=0.95)
AM.Est.MCAR.age.diff.35 <- ci(AM.MCAR.df.35$AM.MCAR.av.age.diff.cov.35, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.35 <- ci(AM.MCAR.df.35$AM.MCAR.between.transm.var.cov.35, confidence=0.95)
AM.Est.AM.MCAR.within.transm.35 <- ci(AM.MCAR.df.35$AM.MCAR.within.transm.var.cov.35, confidence = 0.95)

AM.Est.MCAR.age.male.40 <- ci(AM.MCAR.df.40$AM.MCAR.av.age.male.cov.40, confidence=0.95)
AM.Est.MCAR.age.diff.40 <- ci(AM.MCAR.df.40$AM.MCAR.av.age.diff.cov.40, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.40 <- ci(AM.MCAR.df.40$AM.MCAR.between.transm.var.cov.40, confidence=0.95)
AM.Est.AM.MCAR.within.transm.40 <- ci(AM.MCAR.df.40$AM.MCAR.within.transm.var.cov.40, confidence = 0.95)

AM.Est.MCAR.age.male.45 <- ci(AM.MCAR.df.45$AM.MCAR.av.age.male.cov.45, confidence=0.95)
AM.Est.MCAR.age.diff.45 <- ci(AM.MCAR.df.45$AM.MCAR.av.age.diff.cov.45, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.45 <- ci(AM.MCAR.df.45$AM.MCAR.between.transm.var.cov.45, confidence=0.95)
AM.Est.AM.MCAR.within.transm.45 <- ci(AM.MCAR.df.45$AM.MCAR.within.transm.var.cov.45, confidence = 0.95)

AM.Est.MCAR.age.male.50 <- ci(AM.MCAR.df.50$AM.MCAR.av.age.male.cov.50, confidence=0.95)
AM.Est.MCAR.age.diff.50 <- ci(AM.MCAR.df.50$AM.MCAR.av.age.diff.cov.50, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.50 <- ci(AM.MCAR.df.50$AM.MCAR.between.transm.var.cov.50, confidence=0.95)
AM.Est.AM.MCAR.within.transm.50 <- ci(AM.MCAR.df.50$AM.MCAR.within.transm.var.cov.50, confidence = 0.95)

AM.Est.MCAR.age.male.55 <- ci(AM.MCAR.df.55$AM.MCAR.av.age.male.cov.55, confidence=0.95)
AM.Est.MCAR.age.diff.55 <- ci(AM.MCAR.df.55$AM.MCAR.av.age.diff.cov.55, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.55 <- ci(AM.MCAR.df.55$AM.MCAR.between.transm.var.cov.55, confidence=0.95)
AM.Est.AM.MCAR.within.transm.55 <- ci(AM.MCAR.df.55$AM.MCAR.within.transm.var.cov.55, confidence = 0.95)

AM.Est.MCAR.age.male.60 <- ci(AM.MCAR.df.60$AM.MCAR.av.age.male.cov.60, confidence=0.95)
AM.Est.MCAR.age.diff.60 <- ci(AM.MCAR.df.60$AM.MCAR.av.age.diff.cov.60, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.60 <- ci(AM.MCAR.df.60$AM.MCAR.between.transm.var.cov.60, confidence=0.95)
AM.Est.AM.MCAR.within.transm.60 <- ci(AM.MCAR.df.60$AM.MCAR.within.transm.var.cov.60, confidence = 0.95)

AM.Est.MCAR.age.male.65 <- ci(AM.MCAR.df.65$AM.MCAR.av.age.male.cov.65, confidence=0.95)
AM.Est.MCAR.age.diff.65 <- ci(AM.MCAR.df.65$AM.MCAR.av.age.diff.cov.65, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.65 <- ci(AM.MCAR.df.65$AM.MCAR.between.transm.var.cov.65, confidence=0.95)
AM.Est.AM.MCAR.within.transm.65 <- ci(AM.MCAR.df.65$AM.MCAR.within.transm.var.cov.65, confidence = 0.95)

AM.Est.MCAR.age.male.70 <- ci(AM.MCAR.df.70$AM.MCAR.av.age.male.cov.70, confidence=0.95)
AM.Est.MCAR.age.diff.70 <- ci(AM.MCAR.df.70$AM.MCAR.av.age.diff.cov.70, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.70 <- ci(AM.MCAR.df.70$AM.MCAR.between.transm.var.cov.70, confidence=0.95)
AM.Est.AM.MCAR.within.transm.70 <- ci(AM.MCAR.df.70$AM.MCAR.within.transm.var.cov.70, confidence = 0.95)

AM.Est.MCAR.age.male.75 <- ci(AM.MCAR.df.75$AM.MCAR.av.age.male.cov.75, confidence=0.95)
AM.Est.MCAR.age.diff.75 <- ci(AM.MCAR.df.75$AM.MCAR.av.age.diff.cov.75, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.75 <- ci(AM.MCAR.df.75$AM.MCAR.between.transm.var.cov.75, confidence=0.95)
AM.Est.AM.MCAR.within.transm.75 <- ci(AM.MCAR.df.75$AM.MCAR.within.transm.var.cov.75, confidence = 0.95)

AM.Est.MCAR.age.male.80 <- ci(AM.MCAR.df.80$AM.MCAR.av.age.male.cov.80, confidence=0.95)
AM.Est.MCAR.age.diff.80 <- ci(AM.MCAR.df.80$AM.MCAR.av.age.diff.cov.80, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.80 <- ci(AM.MCAR.df.80$AM.MCAR.between.transm.var.cov.80, confidence=0.95)
AM.Est.AM.MCAR.within.transm.80 <- ci(AM.MCAR.df.80$AM.MCAR.within.transm.var.cov.80, confidence = 0.95)

AM.Est.MCAR.age.male.85 <- ci(AM.MCAR.df.85$AM.MCAR.av.age.male.cov.85, confidence=0.95)
AM.Est.MCAR.age.diff.85 <- ci(AM.MCAR.df.85$AM.MCAR.av.age.diff.cov.85, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.85 <- ci(AM.MCAR.df.85$AM.MCAR.between.transm.var.cov.85, confidence=0.95)
AM.Est.AM.MCAR.within.transm.85 <- ci(AM.MCAR.df.85$AM.MCAR.within.transm.var.cov.85, confidence = 0.95)

AM.Est.MCAR.age.male.90 <- ci(AM.MCAR.df.90$AM.MCAR.av.age.male.cov.90, confidence=0.95)
AM.Est.MCAR.age.diff.90 <- ci(AM.MCAR.df.90$AM.MCAR.av.age.diff.cov.90, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.90 <- ci(AM.MCAR.df.90$AM.MCAR.between.transm.var.cov.90, confidence=0.95)
AM.Est.AM.MCAR.within.transm.90 <- ci(AM.MCAR.df.90$AM.MCAR.within.transm.var.cov.90, confidence = 0.95)

AM.Est.MCAR.age.male.95 <- ci(AM.MCAR.df.95$AM.MCAR.av.age.male.cov.95, confidence=0.95)
AM.Est.MCAR.age.diff.95 <- ci(AM.MCAR.df.95$AM.MCAR.av.age.diff.cov.95, confidence = 0.95)

AM.Est.AM.MCAR.between.transm.95 <- ci(AM.MCAR.df.95$AM.MCAR.between.transm.var.cov.95, confidence=0.95)
AM.Est.AM.MCAR.within.transm.95 <- ci(AM.MCAR.df.95$AM.MCAR.within.transm.var.cov.95, confidence = 0.95)



# Age male - MCAR
df.age.male.AgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
                                 F = c(AM.Est.age.male.pop[[1]], AM.Est.MCAR.age.male.35[[1]], AM.Est.MCAR.age.male.40[[1]], AM.Est.MCAR.age.male.45[[1]],
                                       AM.Est.MCAR.age.male.50[[1]], AM.Est.MCAR.age.male.55[[1]], AM.Est.MCAR.age.male.60[[1]], AM.Est.MCAR.age.male.65[[1]],
                                       AM.Est.MCAR.age.male.70[[1]], AM.Est.MCAR.age.male.75[[1]], AM.Est.MCAR.age.male.80[[1]], AM.Est.MCAR.age.male.85[[1]],
                                       AM.Est.MCAR.age.male.90[[1]], AM.Est.MCAR.age.male.95[[1]]),
                                 L = c(AM.Est.age.male.pop[[2]], AM.Est.MCAR.age.male.35[[2]], AM.Est.MCAR.age.male.40[[2]], AM.Est.MCAR.age.male.45[[2]],
                                       AM.Est.MCAR.age.male.50[[2]], AM.Est.MCAR.age.male.55[[2]], AM.Est.MCAR.age.male.60[[2]], AM.Est.MCAR.age.male.65[[2]],
                                       AM.Est.MCAR.age.male.70[[2]], AM.Est.MCAR.age.male.75[[2]], AM.Est.MCAR.age.male.80[[2]], AM.Est.MCAR.age.male.85[[2]],
                                       AM.Est.MCAR.age.male.90[[2]], AM.Est.MCAR.age.male.95[[2]]),
                                 U = c(AM.Est.age.male.pop[[3]], AM.Est.MCAR.age.male.35[[3]], AM.Est.MCAR.age.male.40[[3]], AM.Est.MCAR.age.male.45[[3]],
                                       AM.Est.MCAR.age.male.50[[3]], AM.Est.MCAR.age.male.55[[3]], AM.Est.MCAR.age.male.60[[3]], AM.Est.MCAR.age.male.65[[3]],
                                       AM.Est.MCAR.age.male.70[[3]], AM.Est.MCAR.age.male.75[[3]], AM.Est.MCAR.age.male.80[[3]], AM.Est.MCAR.age.male.85[[3]],
                                       AM.Est.MCAR.age.male.90[[3]], AM.Est.MCAR.age.male.95[[3]]))

# Age diff - MCAR
df.age.diff.AgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
                                 F = c(AM.Est.age.diff.pop[[1]], AM.Est.MCAR.age.diff.35[[1]], AM.Est.MCAR.age.diff.40[[1]], AM.Est.MCAR.age.diff.45[[1]],
                                       AM.Est.MCAR.age.diff.50[[1]], AM.Est.MCAR.age.diff.55[[1]], AM.Est.MCAR.age.diff.60[[1]], AM.Est.MCAR.age.diff.65[[1]],
                                       AM.Est.MCAR.age.diff.70[[1]], AM.Est.MCAR.age.diff.75[[1]], AM.Est.MCAR.age.diff.80[[1]], AM.Est.MCAR.age.diff.85[[1]],
                                       AM.Est.MCAR.age.diff.90[[1]], AM.Est.MCAR.age.diff.95[[1]]),
                                 L = c(AM.Est.age.diff.pop[[2]], AM.Est.MCAR.age.diff.35[[2]], AM.Est.MCAR.age.diff.40[[2]], AM.Est.MCAR.age.diff.45[[2]],
                                       AM.Est.MCAR.age.diff.50[[2]], AM.Est.MCAR.age.diff.55[[2]], AM.Est.MCAR.age.diff.60[[2]], AM.Est.MCAR.age.diff.65[[2]],
                                       AM.Est.MCAR.age.diff.70[[2]], AM.Est.MCAR.age.diff.75[[2]], AM.Est.MCAR.age.diff.80[[2]], AM.Est.MCAR.age.diff.85[[2]],
                                       AM.Est.MCAR.age.diff.90[[2]], AM.Est.MCAR.age.diff.95[[2]]),
                                 U = c(AM.Est.age.diff.pop[[3]], AM.Est.MCAR.age.diff.35[[3]], AM.Est.MCAR.age.diff.40[[3]], AM.Est.MCAR.age.diff.45[[3]],
                                       AM.Est.MCAR.age.diff.50[[3]], AM.Est.MCAR.age.diff.55[[3]], AM.Est.MCAR.age.diff.60[[3]], AM.Est.MCAR.age.diff.65[[3]],
                                       AM.Est.MCAR.age.diff.70[[3]], AM.Est.MCAR.age.diff.75[[3]], AM.Est.MCAR.age.diff.80[[3]], AM.Est.MCAR.age.diff.85[[3]],
                                       AM.Est.MCAR.age.diff.90[[3]], AM.Est.MCAR.age.diff.95[[3]]))


# Age between groups - MCAR
df.age.between.AgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
                                    F = c(AM.Est.between.transm.pop[[1]], AM.Est.AM.MCAR.between.transm.35[[1]], AM.Est.AM.MCAR.between.transm.40[[1]], AM.Est.AM.MCAR.between.transm.45[[1]],
                                          AM.Est.AM.MCAR.between.transm.50[[1]], AM.Est.AM.MCAR.between.transm.55[[1]], AM.Est.AM.MCAR.between.transm.60[[1]], AM.Est.AM.MCAR.between.transm.65[[1]],
                                          AM.Est.AM.MCAR.between.transm.70[[1]], AM.Est.AM.MCAR.between.transm.75[[1]], AM.Est.AM.MCAR.between.transm.80[[1]], AM.Est.AM.MCAR.between.transm.85[[1]],
                                          AM.Est.AM.MCAR.between.transm.90[[1]], AM.Est.AM.MCAR.between.transm.95[[1]]),
                                    L = c(AM.Est.between.transm.pop[[2]], AM.Est.AM.MCAR.between.transm.35[[2]], AM.Est.AM.MCAR.between.transm.40[[2]], AM.Est.AM.MCAR.between.transm.45[[2]],
                                          AM.Est.AM.MCAR.between.transm.50[[2]], AM.Est.AM.MCAR.between.transm.55[[2]], AM.Est.AM.MCAR.between.transm.60[[2]], AM.Est.AM.MCAR.between.transm.65[[2]],
                                          AM.Est.AM.MCAR.between.transm.70[[2]], AM.Est.AM.MCAR.between.transm.75[[2]], AM.Est.AM.MCAR.between.transm.80[[2]], AM.Est.AM.MCAR.between.transm.85[[2]],
                                          AM.Est.AM.MCAR.between.transm.90[[2]], AM.Est.AM.MCAR.between.transm.95[[2]]),
                                    U = c(AM.Est.between.transm.pop[[3]], AM.Est.AM.MCAR.between.transm.35[[3]], AM.Est.AM.MCAR.between.transm.40[[3]], AM.Est.AM.MCAR.between.transm.45[[3]],
                                          AM.Est.AM.MCAR.between.transm.50[[3]], AM.Est.AM.MCAR.between.transm.55[[3]], AM.Est.AM.MCAR.between.transm.60[[3]], AM.Est.AM.MCAR.between.transm.65[[3]],
                                          AM.Est.AM.MCAR.between.transm.70[[3]], AM.Est.AM.MCAR.between.transm.75[[3]], AM.Est.AM.MCAR.between.transm.80[[3]], AM.Est.AM.MCAR.between.transm.85[[3]],
                                          AM.Est.AM.MCAR.between.transm.90[[3]], AM.Est.AM.MCAR.between.transm.95[[3]]))

# Age within groups - MCAR
df.age.within.AgeMix <- data.frame(x=c("pop", seq(from=35, to=95, by=5)),
                                   F = c(AM.Est.within.transm.pop[[1]], AM.Est.AM.MCAR.within.transm.35[[1]], AM.Est.AM.MCAR.within.transm.40[[1]], AM.Est.AM.MCAR.within.transm.45[[1]],
                                         AM.Est.AM.MCAR.within.transm.50[[1]], AM.Est.AM.MCAR.within.transm.55[[1]], AM.Est.AM.MCAR.within.transm.60[[1]], AM.Est.AM.MCAR.within.transm.65[[1]],
                                         AM.Est.AM.MCAR.within.transm.70[[1]], AM.Est.AM.MCAR.within.transm.75[[1]], AM.Est.AM.MCAR.within.transm.80[[1]], AM.Est.AM.MCAR.within.transm.85[[1]],
                                         AM.Est.AM.MCAR.within.transm.90[[1]], AM.Est.AM.MCAR.within.transm.95[[1]]),
                                   L = c(AM.Est.within.transm.pop[[2]], AM.Est.AM.MCAR.within.transm.35[[2]], AM.Est.AM.MCAR.within.transm.40[[2]], AM.Est.AM.MCAR.within.transm.45[[2]],
                                         AM.Est.AM.MCAR.within.transm.50[[2]], AM.Est.AM.MCAR.within.transm.55[[2]], AM.Est.AM.MCAR.within.transm.60[[2]], AM.Est.AM.MCAR.within.transm.65[[2]],
                                         AM.Est.AM.MCAR.within.transm.70[[2]], AM.Est.AM.MCAR.within.transm.75[[2]], AM.Est.AM.MCAR.within.transm.80[[2]], AM.Est.AM.MCAR.within.transm.85[[2]],
                                         AM.Est.AM.MCAR.within.transm.90[[2]], AM.Est.AM.MCAR.within.transm.95[[2]]),
                                   U = c(AM.Est.within.transm.pop[[3]], AM.Est.AM.MCAR.within.transm.35[[3]], AM.Est.AM.MCAR.within.transm.40[[3]], AM.Est.AM.MCAR.within.transm.45[[3]],
                                         AM.Est.AM.MCAR.within.transm.50[[3]], AM.Est.AM.MCAR.within.transm.55[[3]], AM.Est.AM.MCAR.within.transm.60[[3]], AM.Est.AM.MCAR.within.transm.65[[3]],
                                         AM.Est.AM.MCAR.within.transm.70[[3]], AM.Est.AM.MCAR.within.transm.75[[3]], AM.Est.AM.MCAR.within.transm.80[[3]], AM.Est.AM.MCAR.within.transm.85[[3]],
                                         AM.Est.AM.MCAR.within.transm.90[[3]], AM.Est.AM.MCAR.within.transm.95[[3]]))


# ggplot(df.age.male.AgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) + 
#   ggtitle("Population level average man age & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Average man age")
# 
# ggplot(df.age.diff.AgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Gender effect & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Age difference")
# 
# 
# ggplot(df.age.between.AgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Between clusters age variation & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Between clusters age variation")
# 
# 
# ggplot(df.age.within.AgeMix, aes(x = x, y = F)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L)) +
#   ggtitle("Within clusters age variation & Sequence coverage - MCAR") +
#   xlab("Sequence coverage scenarios") + ylab("Within clusters age variation")


## Comparing plots

# Average of age of male

ggplot(data = df.age.male.NonAgeMix, aes(x = x, y = F)) + 
  geom_point(colour = 'blue', size = 4) + 
  geom_errorbar(data = df.age.male.NonAgeMix, aes(ymax = U, ymin = L)) +
  geom_point(data = df.age.male.AgeMix, colour = 'red' , size = 4) +   
  geom_errorbar(data = df.age.male.AgeMix, aes(ymax = U, ymin = L)) +
  # ggtitle("Average of age of male - MCAR") +
  xlab("Sequence coverage scenarios") + ylab("Average age of men") + 
  theme(axis.text = element_text(face = 'bold', size = 18), 
        axis.title.x = element_text(face = 'bold', size = 28),
        axis.title.y = element_text(face = 'bold', size = 28),
        title = element_text(face = 'bold', size = 28))


# Gender effect

ggplot(data = df.age.diff.NonAgeMix, aes(x = x, y = F)) + 
  geom_point(colour = 'blue', size = 4) + 
  geom_errorbar(data = df.age.diff.NonAgeMix, aes(ymax = U, ymin = L)) +
  geom_point(data = df.age.diff.AgeMix, colour = 'red' , size = 4) +   
  geom_errorbar(data = df.age.diff.AgeMix, aes(ymax = U, ymin = L)) +
  # ggtitle("Age difference - MCAR") +
  xlab("Sequence coverage") + ylab("Gender effect") + 
  theme(axis.text = element_text(face = 'bold', size = 18), 
        axis.title.x = element_text(face = 'bold', size = 28),
        axis.title.y = element_text(face = 'bold', size = 28),
        title = element_text(face = 'bold', size = 28))


# Between cluster variation

ggplot(data = df.age.between.NonAgeMix, aes(x = x, y = F)) + 
  geom_point(colour = 'blue', size = 4) + 
  geom_errorbar(data = df.age.between.NonAgeMix, aes(ymax = U, ymin = L)) +
  geom_point(data = df.age.between.AgeMix, colour = 'red' , size = 4) +   
  geom_errorbar(data = df.age.between.AgeMix, aes(ymax = U, ymin = L)) +
  # ggtitle("Gender effect & Sequence coverage - MCAR") +
  xlab("Sequence coverage") + ylab("Age variation") + 
  theme(axis.text = element_text(face = 'bold', size = 18), 
        axis.title.x = element_text(face = 'bold', size = 28),
        axis.title.y = element_text(face = 'bold', size = 28),
        title = element_text(face = 'bold', size = 28))

# Within cluster variation

ggplot(data = df.age.within.NonAgeMix, aes(x = x, y = F)) + 
  geom_point(colour = 'blue', size = 4) + 
  geom_errorbar(data = df.age.within.NonAgeMix, aes(ymax = U, ymin = L)) +
  geom_point(data = df.age.within.AgeMix, colour = 'red' , size = 4) +   
  geom_errorbar(data = df.age.within.AgeMix, aes(ymax = U, ymin = L)) +
  # ggtitle("Gender effect & Sequence coverage - MCAR") +
  xlab("Sequence coverage") + ylab("Age variation") + 
  theme(axis.text = element_text(face = 'bold', size = 18), 
        axis.title.x = element_text(face = 'bold', size = 28),
        axis.title.y = element_text(face = 'bold', size = 28),
        title = element_text(face = 'bold', size = 28))

