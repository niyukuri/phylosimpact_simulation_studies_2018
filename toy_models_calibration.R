library(devtools)
# library(devtools)
# install_github("wdelva/RSimpactHelp")
library(RSimpactHelper)
# THIS IS ESSENTIAL: So that the imputation functions are in the global environment, and can be found.
mice.impute.norm <- mice::mice.impute.norm
mice.impute.rf <- mice::mice.impute.rf


# 2 parameters
##################

toy_model <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + rnorm(1,0,0.1) ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + rnorm(1,0,0.1) )
}

sum_stat_obs <- c(1.5, 0.5)
lls = c(0, 0)
uls = c(1.5, 2)

# To get more info about the arguments of the MaC function:
# help(MaC)

# MaC.toy <- MaC(targets.empirical = sum_stat_obs,
#                RMSD.tol.max = 2,
#                min.givetomice = 200,
#                n.experiments = 2000,
#                lls = lls,
#                uls = uls,
#                model = toy_model,
#                strict.positive.params = 0,
#                probability.params = 0,
#                method = "norm",
#                predictorMatrix = "complete",
#                maxit = 20,
#                maxwaves = 6,
#                n_cluster = 8)

# debugonce(MaC.weighted) # To go into debug mode

MaC.toy.weighted <- MaC.weighted(targets.empirical = sum_stat_obs,
                                 RMSD.tol.max = 2,
                                 min.givetomice = 200,
                                 n.experiments = 2000,
                                 lls = lls,
                                 uls = uls,
                                 model = toy_model,
                                 strict.positive.params = 0,
                                 probability.params = 0,
                                 inside_prior = TRUE,
                                 method = "norm",
                                 predictorMatrix = "complete",
                                 maxit = 20,
                                 maxwaves = 6,
                                 n_cluster = 8)



# 3 parameters
##############

toy_model.three <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1))
}

sum_stat_obs <- c(1.5, 0.5)

lls.three = c(0, 1, 1)
uls.three = c(1.5, 2, 3)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.three <- MaC.weighted(targets.empirical = sum_stat_obs,
                              RMSD.tol.max = 2,
                              min.givetomice = 200,
                              n.experiments = 2000,
                              lls = lls.three,
                              uls = uls.three,
                              model = toy_model.three,
                              strict.positive.params = 0,
                              probability.params = 0,
                              inside_prior = TRUE,
                              method = "norm",
                              predictorMatrix = "complete",
                              maxit = 20,
                              maxwaves = 6,
                              n_cluster = 8)


# 4 parameters
##############

toy_model.four <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3] + x[4] + rnorm(1,0,0.1) + x[5],   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5])
}

sum_stat_obs <- c(1.5, 0.5)

lls.four = c(0, 1, 1, 0)
uls.four = c(1.5, 2, 3, 1)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.four.weighted <- MaC.weighted(targets.empirical = sum_stat_obs,
                                      RMSD.tol.max = 2,
                                      min.givetomice = 200,
                                      n.experiments = 2000,
                                      lls = lls.four,
                                      uls = uls.four,
                                      model = toy_model.four,
                                      strict.positive.params = 0,
                                      probability.params = 0,
                                      inside_prior = TRUE,
                                      method = "norm",
                                      predictorMatrix = "complete",
                                      maxit = 20,
                                      maxwaves = 6,
                                      n_cluster = 8)




# 5 parameters
##############

toy_model.five <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5]+ x[6])
}

sum_stat_obs <- c(1.5, 0.5)

lls.five = c(0, 1, 1, 0, 3)
uls.five = c(1.5, 2, 3, 1, 5)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.five <- MaC.weighted(targets.empirical = sum_stat_obs,
                             RMSD.tol.max = 2,
                             min.givetomice = 200,
                             n.experiments = 2000,
                             lls = lls.five,
                             uls = uls.five,
                             model = toy_model.five,
                             strict.positive.params = 0,
                             probability.params = 0,
                             inside_prior = TRUE,
                             method = "norm",
                             predictorMatrix = "complete",
                             maxit = 20,
                             maxwaves = 6,
                             n_cluster = 8)




# 6 parameters
##############

toy_model.six <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] +x[7],   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5]+ x[6]*x[7])
}

sum_stat_obs <- c(1.5, 0.5)

lls.six = c(0, 1, 1, 0, 3,1)
uls.six = c(1.5, 2, 3, 1, 5,3)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.six <- MaC.weighted(targets.empirical = sum_stat_obs,
                            RMSD.tol.max = 2,
                            min.givetomice = 200,
                            n.experiments = 2000,
                            lls = lls.six,
                            uls = uls.six,
                            model = toy_model.six,
                            strict.positive.params = 0,
                            probability.params = 0,
                            inside_prior = TRUE,
                            method = "norm",
                            predictorMatrix = "complete",
                            maxit = 20,
                            maxwaves = 6,
                            n_cluster = 8)



# 7 parameters
##############

toy_model.seven <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] + x[7] + x[8] ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5] + x[7] + x[6] + x[8])
}

sum_stat_obs <- c(1.5, 0.5)

lls.seven = c(0, 1, 1, 0, 3, 5 , 0)
uls.seven = c(1.5, 2, 3, 1, 5, 10, 1)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.seven <- MaC.weighted(targets.empirical = sum_stat_obs,
                              RMSD.tol.max = 2,
                              min.givetomice = 200,
                              n.experiments = 2000,
                              lls = lls.seven,
                              uls = uls.seven,
                              model = toy_model.seven,
                              strict.positive.params = 0,
                              probability.params = 0,
                              inside_prior = TRUE,
                              method = "norm",
                              predictorMatrix = "complete",
                              maxit = 20,
                              maxwaves = 6,
                              n_cluster = 8)



##########################################################################################################
##########################################################################################################


# Let's compare this to accept-reject ABC
library(EasyABC)

# help("ABC_rejection")

toy_prior <- list(c("unif", lls[1], uls[1]),
                  c("unif", lls[2], uls[2]))

Rej.toy <- ABC_rejection(model = toy_model,
                         prior = toy_prior,
                         summary_stat_target = sum_stat_obs,
                         nb_simul = 12000,
                         use_seed = TRUE,
                         seed_count = 1,
                         n_cluster = 8,
                         tol = 200/12000)

# Let's compare this to sequential ABC

# help("ABC_sequential")

Seq.toy <- ABC_sequential(model = toy_model,
                          method = "Lenormand",
                          prior = toy_prior,
                          summary_stat_target = sum_stat_obs,
                          nb_simul = 2000,
                          alpha = 0.1,
                          p_acc_min = 0.03,
                          use_seed = TRUE,
                          seed_count = 1,
                          n_cluster = 8,
                          inside_prior = FALSE)
# To see how many waves were done:
1 + (Seq.toy$nsim - 2000) / 1800

# Plotting the input parameters of the calibrated model
plot(Rej.toy$param[, 1],
     Rej.toy$param[, 2],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1.5),
     ylim = c(0, 2))
points(Seq.toy$param[, 1],
       Seq.toy$param[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy$selected.experiments[[6]][, 1],
       MaC.toy$selected.experiments[[6]][, 2],
       pch = 16,
       col = "orange")
points(MaC.toy.weighted$selected.experiments[[1]][, 1],
       MaC.toy.weighted$selected.experiments[[1]][, 2],
       pch = 16,
       col = "red")

# Plotting the summary statistics of the calibrated model
plot(Rej.toy$stats[, 1],
     Rej.toy$stats[, 2],
     pch = 16,
     col = "black",
     xlab = "summary statistic 1",
     ylab = "summary statistic 2",
     xlim = c(1.4, 1.6),
     ylim = c(0.3, 0.7))
points(Seq.toy$stats[, 1],
       Seq.toy$stats[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy$selected.experiments[[6]][, 3],
       MaC.toy$selected.experiments[[6]][, 4],
       pch = 16,
       col = "orange")
points(MaC.toy.weighted$selected.experiments[[6]][, 3],
       MaC.toy.weighted$selected.experiments[[6]][, 4],
       pch = 16,
       col = "red")

points(sum_stat_obs[1],
       sum_stat_obs[2],
       pch = "+",
       cex = 3,
       col = "red3")

# Comparing simulation time:
MaC.toy$secondspassed[3]
Rej.toy$computime
Seq.toy$computime


# 3 parameters
##############

toy_model.three <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1))
}

sum_stat_obs <- c(1.5, 0.5)

lls.three = c(0, 1, 1)
uls.three = c(1.5, 2, 3)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.three <- MaC.weighted(targets.empirical = sum_stat_obs,
                              RMSD.tol.max = 2,
                              min.givetomice = 200,
                              n.experiments = 2000,
                              lls = lls.three,
                              uls = uls.three,
                              model = toy_model.three,
                              strict.positive.params = 0,
                              probability.params = 0,
                              inside_prior = TRUE,
                              method = "norm",
                              predictorMatrix = "complete",
                              maxit = 20,
                              maxwaves = 6,
                              n_cluster = 8)


# 4 parameters
##############

toy_model.four <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3] + x[4] + rnorm(1,0,0.1) + x[5],   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5])
}

sum_stat_obs <- c(1.5, 0.5)

lls.four = c(0, 1, 1, 0)
uls.four = c(1.5, 2, 3, 1)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.four.weighted <- MaC.weighted(targets.empirical = sum_stat_obs,
                                      RMSD.tol.max = 2,
                                      min.givetomice = 200,
                                      n.experiments = 2000,
                                      lls = lls.four,
                                      uls = uls.four,
                                      model = toy_model.four,
                                      strict.positive.params = 0,
                                      probability.params = 0,
                                      inside_prior = TRUE,
                                      method = "norm",
                                      predictorMatrix = "complete",
                                      maxit = 20,
                                      maxwaves = 6,
                                      n_cluster = 8)
# 
# Error in MaC(targets.empirical = sum_stat_obs, RMSD.tol.max = 2, min.givetomice = 200,  : 
#                unused arguments (strict.positive.fourams = 0, probability.fourams = 0)


# Let's compare this to accept-reject ABC
library(EasyABC)

# help("ABC_rejection")

toy_prior.four <- list(c("unif", lls.four[1], uls.four[1]),
                       c("unif", lls.four[2], uls.four[2]),
                       c("unif", lls.four[3], uls.four[3]),
                       c("unif", lls.four[4], uls.four[4]))


Rej.toy.four <- ABC_rejection(model = toy_model.four,
                              prior = toy_prior.four,
                              summary_stat_target = sum_stat_obs,
                              nb_simul = 12000,
                              use_seed = TRUE,
                              seed_count = 1,
                              n_cluster = 8,
                              tol = 200/12000)



# Let's compare this to sequential ABC

# help("ABC_sequential")

Seq.toy.four <- ABC_sequential(model = toy_model.four,
                               method = "Lenormand",
                               prior = toy_prior.four,
                               summary_stat_target = sum_stat_obs,
                               nb_simul = 2000,
                               alpha = 0.1,
                               p_acc_min = 0.03,
                               use_seed = TRUE,
                               seed_count = 1,
                               n_cluster = 8,
                               inside_prior = FALSE)
# To see how many waves were done:
1 + (Seq.toy.four$nsim - 2000) / 1800

# Plotting the input parameters of the calibrated model
plot(Rej.toy.four$param[, 1],
     Rej.toy.four$param[, 2],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1.5),
     ylim = c(0, 2))
points(Seq.toy.four$param[, 1],
       Seq.toy.four$param[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy.four$selected.experiments[[6]][, 1],
       MaC.toy.four$selected.experiments[[6]][, 2],
       pch = 16,
       col = "orange")


# Plotting the summary statistics of the calibrated model
plot(Rej.toy.four$stats[, 1],
     Rej.toy.four$stats[, 2],
     pch = 16,
     col = "black",
     xlab = "summary statistic 1",
     ylab = "summary statistic 2",
     xlim = c(1.4, 1.6),
     ylim = c(0.3, 0.7))
points(Seq.toy.four$stats[, 1],
       Seq.toy.four$stats[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy.four$selected.experiments[[6]][, 3],
       MaC.toy.four$selected.experiments[[6]][, 4],
       pch = 16,
       col = "orange")

points(sum_stat_obs[1],
       sum_stat_obs[2],
       pch = "+",
       cex = 3,
       col = "red3")

# Comparing simulation time:
MaC.toy.four$secondspassed[3]
Rej.toy.four$computime
Seq.toy.four$computime



# 5 parameters
##############

toy_model.five <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5]+ x[6])
}

sum_stat_obs <- c(1.5, 0.5)

lls.five = c(0, 1, 1, 0, 3)
uls.five = c(1.5, 2, 3, 1, 5)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.five <- MaC.weighted(targets.empirical = sum_stat_obs,
                             RMSD.tol.max = 2,
                             min.givetomice = 200,
                             n.experiments = 2000,
                             lls = lls.five,
                             uls = uls.five,
                             model = toy_model.five,
                             strict.positive.params = 0,
                             probability.params = 0,
                             inside_prior = TRUE,
                             method = "norm",
                             predictorMatrix = "complete",
                             maxit = 20,
                             maxwaves = 6,
                             n_cluster = 8)




# 6 parameters
##############

toy_model.six <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] +x[7],   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5]+ x[6]*x[7])
}

sum_stat_obs <- c(1.5, 0.5)

lls.six = c(0, 1, 1, 0, 3,1)
uls.six = c(1.5, 2, 3, 1, 5,3)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.six <- MaC.weighted(targets.empirical = sum_stat_obs,
                            RMSD.tol.max = 2,
                            min.givetomice = 200,
                            n.experiments = 2000,
                            lls = lls.six,
                            uls = uls.six,
                            model = toy_model.six,
                            strict.positive.params = 0,
                            probability.params = 0,
                            inside_prior = TRUE,
                            method = "norm",
                            predictorMatrix = "complete",
                            maxit = 20,
                            maxwaves = 6,
                            n_cluster = 8)



# 7 parameters
##############

toy_model.seven <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + x[3]*x[4] + rnorm(1,0,0.1) + x[5] + x[6] + x[7] + x[8] ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + x[4] + rnorm(1,0,0.1)  + x[5] + x[7] + x[6] + x[8])
}

sum_stat_obs <- c(1.5, 0.5)

lls.seven = c(0, 1, 1, 0, 3, 5 , 0)
uls.seven = c(1.5, 2, 3, 1, 5, 10, 1)

# To get more info about the arguments of the MaC function:
# help(MaC)

MaC.toy.seven <- MaC.weighted(targets.empirical = sum_stat_obs,
                              RMSD.tol.max = 2,
                              min.givetomice = 200,
                              n.experiments = 2000,
                              lls = lls.seven,
                              uls = uls.seven,
                              model = toy_model.seven,
                              strict.positive.params = 0,
                              probability.params = 0,
                              inside_prior = TRUE,
                              method = "norm",
                              predictorMatrix = "complete",
                              maxit = 20,
                              maxwaves = 6,
                              n_cluster = 8)


# Let's compare this to accept-reject ABC
library(EasyABC)

# help("ABC_rejection")

toy_prior.seven <- list(c("unif", lls.seven[1], uls.seven[1]),
                        c("unif", lls.seven[2], uls.seven[2]),
                        c("unif", lls.seven[3], uls.seven[3]),
                        c("unif", lls.seven[4], uls.seven[4]),
                        c("unif", lls.seven[5], uls.seven[5]),
                        c("unif", lls.seven[6], uls.seven[6]),
                        c("unif", lls.seven[7], uls.seven[7]))


Rej.toy.seven <- ABC_rejection(model = toy_model.seven,
                               prior = toy_prior.seven,
                               summary_stat_target = sum_stat_obs,
                               nb_simul = 12000,
                               use_seed = TRUE,
                               seed_count = 1,
                               n_cluster = 8,
                               tol = 200/12000)



# Let's compare this to sequential ABC

# help("ABC_sequential")

Seq.toy.seven <- ABC_sequential(model = toy_model.seven,
                                method = "Lenormand",
                                prior = toy_prior.seven,
                                summary_stat_target = sum_stat_obs,
                                nb_simul = 2000,
                                alpha = 0.1,
                                p_acc_min = 0.03,
                                use_seed = TRUE,
                                seed_count = 1,
                                n_cluster = 8,
                                inside_prior = FALSE)
# To see how many waves were done:
1 + (Seq.toy.seven$nsim - 2000) / 1800

# Plotting the input parameters of the calibrated model
plot(Rej.toy.seven$param[, 1],
     Rej.toy.seven$param[, 2],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1.5),
     ylim = c(0, 2))
points(Seq.toy.seven$param[, 1],
       Seq.toy.seven$param[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy.seven$selected.experiments[[6]][, 1],
       MaC.toy.seven$selected.experiments[[6]][, 2],
       pch = 16,
       col = "orange")


# Plotting the summary statistics of the calibrated model
plot(Rej.toy.seven$stats[, 1],
     Rej.toy.seven$stats[, 2],
     pch = 16,
     col = "black",
     xlab = "summary statistic 1",
     ylab = "summary statistic 2",
     xlim = c(1.4, 1.6),
     ylim = c(0.3, 0.7))
points(Seq.toy.seven$stats[, 1],
       Seq.toy.seven$stats[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy.seven$selected.experiments[[6]][, 3],
       MaC.toy.seven$selected.experiments[[6]][, 4],
       pch = 16,
       col = "orange")

points(sum_stat_obs[1],
       sum_stat_obs[2],
       pch = "+",
       cex = 3,
       col = "red3")

# Comparing simulation time:
MaC.toy.seven$secondspassed[3]
Rej.toy.seven$computime
Seq.toy.seven$computime


