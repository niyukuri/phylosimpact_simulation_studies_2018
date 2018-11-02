#######   CALIBRATION    WITH COMBINED CLASSIC STATISTICS  ##########

library(EasyABC)
library(RSimpactCyan)
library(RSimpactHelper)
library(Rcpp)
library(ape)
library(expoTree)
library(data.table)
library(readr)
# library(phangorn)
# library(lme4)
# library(nlme)
library(minque) # with lmer
library(dplyr)
library(adephylo)
library(treedater)
library(geiger)
library(picante)
library(igraph)
library(phyloTop)
library(phytools)
library(Rsamtools) # select IDs sequences in a file
library(robustbase) # colMedians
library(lme4)

library(lhs)
library(abc)


work.dir <- "/home/niyukuri/Desktop/calibration" # on laptop

# work.dir <- "/home/niyukuri/Desktop/calibration" # on PC


setwd(paste0(work.dir))


simpact4ABC.classic <- function(inputvector){
  
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
  setwd(paste0(work.dir))
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.classic.R")
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(readr)
  # library(phangorn)
  # library(lme4)
  # library(nlme)
  library(minque) # with lmer
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  library(Rsamtools) # select IDs sequences in a file
  library(robustbase) # colMedians
  library(lme4)
  
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 4000, 
                                   population.numwomen = 4000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 20, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  ABC_DestDir.classic <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Error function when computing summary statistics
  
  err.functionGEN <- function(e){
    return(chunk.summary.stats = rep(NA,28))
    stop(e)
  }
  
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = ABC_DestDir.classic,
                                  agedist = age.distr,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 28) # 37 + 82 + 33 + 1 + 1 + 53 = 207
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 28)
    }else{
      
      
      datalist <- readthedata(results)
      
      simpact.trans.net <-  transmission.network.builder(datalist = datalist, endpoint = 40)
      
      
      summary.stat.classic <- tryCatch(compute.summary.statistics.classic(datalist = datalist,
                                                                          simpact.trans.net = simpact.trans.net,
                                                                          work.dir = work.dir,
                                                                          sub.dir.rename = ABC_DestDir.classic,
                                                                          dirfasttree = work.dir,
                                                                          limitTransmEvents = 7,
                                                                          seq.cov = 100,
                                                                          age.group.15.25 = c(15,25),
                                                                          age.group.25.40 = c(25,40),
                                                                          age.group.40.50 = c(40,50),
                                                                          endpoint = 40,
                                                                          timewindow = c(30,40)),
                                       error=function(e) return(rep(NA, 28)))
      #                                   error = err.functionGEN)
      
      #  error=function(e) return(rep(NA, 36))
      
      # 
      # summary.stat.phylo <- compute.summary.statistics.phylo(simpact.trans.net = simpact.trans.net,
      #                                                        work.dir = work.dir,
      #                                                        sub.dir.rename = ABC_DestDir.classic,
      #                                                        dirfasttree = work.dir,
      #                                                        limitTransmEvents = 7,
      #                                                        seq.cov = 100,
      #                                                        age.group.15.25 = c(15,25),
      #                                                        age.group.25.40 = c(25,40),
      #                                                        age.group.40.50 = c(40,50),
      #                                                        endpoint = 40,
      #                                                        timewindow = c(30,40))
      
      
      # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
      # agegapsd <- sd(datalist$rtable$AgeGap)
      # outputvector <- c(relsperpersonperyear, agegapsd)
      
      
      outputvector <- summary.stat.classic
      
    }
    
  }
  
  unlink(paste0(ABC_DestDir.classic), recursive = TRUE)
  
  return(outputvector)
}





wrapper.simpact4ABC.classic <- function(inputvector = inputvector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/simpact4ABC.classic.R")
  
  res.simpact4ABC.classic <- tryCatch(simpact4ABC.classic(inputvector = inputvector),
                                      error=function(e) return(rep(NA, 28)))
  
  
  return(res.simpact4ABC.classic)
  
}


source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/calibration.ABC.R")


cal.par <- calibration.ABC(model.sim = wrapper.simpact4ABC.classic,
                           sum_stat_obs = sum_stat_obs,
                           simpact_prior = simpact_prior, 
                           design.points = 200,
                           seed.val = 1,
                           n_cores = 8)

save(cal.par, file = "cal.par.RData")

cal.par <- get(load("cal.par.RData"))

sum.neuralnet.cal.par <- summary(cal.par, intvl = .9)

par.run.models <- cal.par$adj.values

source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metrics.R")

after.cal.res <- simpact.parallel(model = complete.master.epic.metrics,
                                  actual.input.matrix = par.run.models,
                                  seed_count = 124,
                                  n_cluster = 8)


par.run.models=as.matrix(par.run.models[1:8,])

v = complete.master.epic.metrics(inputvector = inputvector)

# Read results from the master model and compute values for observed summary statistics

# These results were obtained by inputvector

inputvector <- c(123,-0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) 


df <- read.csv(file = "~/Desktop/mastermodeltest/epi.Metrics.features.Save.csv")


# place holder sum_stat_obs

sum_stat_obs <- c(-0.217246,  0.51472275,  0.08759751,  2.18014706,  0.7844332,  1.1012360,  0.0017159763,  0.99035413,  1.42268384, 9.05340867,  7.04394701,
                  0.24145866,  6.94741730,  5.46243176, -1.50937223,  1.06993007, 8.39706516,  7.32921946,  2.37464164,  7.92930670,  4.55751, -9.78684982,
                  45.41647034, -1.06226184,  9.80972900,  7.48176900,  12.48176900,  7.74260912)


# Compute sum_stat_obs according to different scenarios


# Calibration with ABC approach

# (i). Having plausible ranges for the parameters, sample the parameter spaces by latin hypercube several times
# (ii). Run the default model and compute the summary statistics
# (iii). Use different ABC-based methods to fit the model, this will give parameters estimates and associated summary statistics


# (i)

design.points <- 100 # number of simulations repeats


simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 2, 7), c("unif", 5, 10), c("unif", 2, 4), c("unif", 0, 1),
                      c("unif", -1, 0), c("unif", -0.9, 0), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), c("unif", 0, 0.5),
                      c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9), c("unif", 10, 14), c("unif", -3.5, -1.7))

min.v <- vector()
max.v <- vector()

for( i in 1:length(simpact_prior)){
  
  min.v <- c(min.v, as.numeric(simpact_prior[[i]][[2]]))
  max.v <- c(max.v, as.numeric(simpact_prior[[i]][[3]]))  
  
}


variables <- length(inputvector)

set.seed(1)

rlhs <- randomLHS(design.points, variables)


lhs.df <- rlhs

for (j in 1:ncol(lhs.df)){
  
  min.var <- min.v[j]
  max.var <- max.v[j]
  
  for(k in 1:nrow(lhs.df)){
    
    lhs.df[k,j] <- qunif(lhs.df[k,j], min = min.var, max = max.var)
    
  }
  
}



par.sim <- inputmatrix <- lhs.df # results of (i): parameter matrix



# (ii)

stat.sim <- simpact.parallel(model = simpact4ABC.classic,
                             actual.input.matrix = par.sim,
                             seed_count = 124,
                             n_cluster = 8)


stat.sim <- read.csv("stat.sim.csv")


stat.sim <- stat.sim[,2:29] # results of (ii): summary statistics matrix obtained from simulations done with parameter matrix

stat.obs <- sum_stat_obs


# (iii)

# Condition: dim(par.sim) == dim(stat.sim)

rej <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, method = "rejection") 


neuralnet <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, method = "neuralnet") 

param.vect <- as.data.frame(neuralnet$adj.values)

# Summary
sum.neuralnet <- summary(neuralnet, intvl = .9)

## posterior histograms
hist(neuralnet)















################    Examples with ABC

library(abc)

require(abc.data)
data(musigma2)


## The rejection algorithm
##
rej <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, method =
             "rejection") 

## ABC with local linear regression correction without/with correction
## for heteroscedasticity 
##
lin <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, hcorr =
             FALSE, method = "loclinear", transf=c("none","log"))

linhc <- abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=.1, method =
               "loclinear", transf=c("none","log")) 

## posterior summaries
##
linsum <- summary(linhc, intvl = .9)
linsum
## compare with the rejection sampling
summary(linhc, unadj = TRUE, intvl = .9)

## posterior histograms
##
hist(linhc, breaks=30, caption=c(expression(mu),
                                 expression(sigma^2))) 

## or send histograms to a pdf file
hist(linhc, file="linhc", breaks=30, caption=c(expression(mu),
                                               expression(sigma^2)))

## diagnostic plots: compare the 2 'abc' objects: "loclinear",
## "loclinear" with correction for heteroscedasticity
##
plot(lin, param=par.sim)
plot(linhc, param=par.sim)

## example illustrates how to add "true" parameter values to a plot
##
postmod <- c(post.mu[match(max(post.mu[,2]), post.mu[,2]),1],
             post.sigma2[match(max(post.sigma2[,2]), post.sigma2[,2]),1])
plot(linhc, param=par.sim, true=postmod)


## artificial example to show how to use the logit tranformations
##
myp <- data.frame(par1=runif(1000,-1,1),par2=rnorm(1000),par3=runif(1000,0,2))
mys <- myp+rnorm(1000,sd=.1)
myt <- c(0,0,1.5)
lin2 <- abc(target=myt, param=myp, sumstat=mys, tol=.1, method =
              "loclinear", transf=c("logit","none","logit"),logit.bounds = rbind(c(-1,
                                                                                   1), c(NA, NA), c(0, 2)))
summary(lin2)

