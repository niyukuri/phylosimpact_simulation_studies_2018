
#######   CALIBRATION

library(EasyABC)
library(RSimpactCyan)
library(RSimpactHelper)
library(Rcpp)
library(ape)
library(expoTree)
library(data.table)
library(readr)
library(phangorn)
library(lme4)
library(nlme)
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


work.dir <- "/home/david/Desktop/calibration" # on laptop


setwd(paste0(work.dir))

# ABC_DestDir <- "/home/david/Desktop/mastermodeltest/calibration"

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

ABC_DestDir <- paste0(work.dir,"/temp/",generate.filename(10))



simpact4ABC <- function(inputvector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/gender.men50.women50.compute.summary.statistics.combined.65.R")
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 400, 
                                   population.numwomen = 400,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 20, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Sexual behaviour
  # ###################
  #
  cfg.list["dissolution.alpha_0"] <- inputvector[1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector[2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[8] # -0.1 c("unif", -0.9, 0)
  cfg.list["population.eyecap.fraction"] <- inputvector[9] # 0.2 c("unif", 0, 0.5)
  #
  # # HIV transmission
  # ###################
  #
  
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # -0.1386294 c("unif", -0.5, 0)
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # 12 c("unif", 10, 14)
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # -2.7 c("unif", -3.5, -1.7)
  
  #
  #
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  #
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  #
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
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
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.4
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
  
  #
  # Let's introduce ART,
  art.intro <- list()
  art.intro["time"] <- 25 #25
  art.intro["diagnosis.baseline"] <- 100
  art.intro["monitoring.cd4.threshold"] <- 100
  
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 + 5 #25 + 5 = 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 25 + 8 #25 + 8 = 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 25 + 11 #25 + 11 = 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 25 + 13 #25 + 13 = 38
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  cfg <- cfg.list
  
  results <- simpact.run(cfg, ABC_DestDir) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  summary.stat <- gender.men50.women50.compute.summary.statistics.combined.65(datalist = datalist,
                                                                              work.dir = work.dir,
                                                                              sub.dir.rename = ABC_DestDir)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat)
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


#We specify the prior distributions for the input parameters
# The relationship formation rate must roughly double: hF' = 2hF
# hF = exp(a0 + a1X1 + a2X2 + .... anXn)
# 2hF = exp(a0 + a1X1 + a2X2 + .... anXn) * 2
# 2hF = exp(a0 + a1X1 + a2X2 + .... anXn) * exp(log(2))
# 2hF = exp(a0 + log(2) + a1X1 + a2X2 + .... anXn)
# So we would naively expect that the baseline parameter (0.1) should be increased by log(2) ~ 0.7 to 0.8
# However, we are also adjusting the "gap factors" and making relationships with large age gaps less
# likely will result in an overall decrease in the number of relationships formed per time unit.

simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 1, 3), c("unif", -0.5, 0.5),
                      c("unif", 2, 4), c("unif", 0, 1), c("unif", -1, 0), c("unif", -0.9, 0), 
                      c("unif", 0, 0.5), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), 
                      c("unif", 0, 0.5), c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9),
                      c("unif", 10, 14),  c("unif", -3.5, -1.7))
# Lastly, we specify the target summary statistic
# sum_stat_obs <- c(0.003916345 , 0.000000000 , 0.000000000 , 0.000000000 , 0.000000000 , 0.083333333 , 0.000000000 , 2.033333333, 22.050000000,
#                   0.066818182 , 5.006713645 , 4.920962084 , 3.942774744 , 0.613515378 , 2.449657529 , 1.033541863 , -1.386704665,  0.062130178,
#                   6.003394756 , 0.285714286 , 0.588235294 , 4.666666667 , 3.000000000 , 0.428571429 , 1.142857143 , 1.000000000)


sum_stat_obs <- c(0.003916345 , 0.000000000 , 0.000000000 , 0.000000000 , 
                  0.000000000 , 0.083333333 , 0.000000000 , 
                  # 2.033333333, 
                  # 22.050000000, 
                  # 
                  # 0.066818182 , 5.006713645 , 4.920962084 , 3.942774744 , 
                  # 0.613515378 , 2.449657529 , 
                  1.033541863 , -1.386704665,  0.062130178, 6.003394756 , 
                  0.285714286 , 0.588235294 , 4.666666667 , 3.000000000 , 
                  0.428571429 , 1.142857143 , 1.000000000)


# sum_stat_obs <- c(0.003916345)
# Now we try to run a sequential ABC scheme, according to the method proposed by Lenormand et al. 2013
# Maxime Lenormand, Franck Jabot and Guillaume Deffuant. Adaptive approximate Bayesian computation for complex models. Comput Stat (2013) 28:2777–2796 DOI 10.1007/s00180-013-0428-3


# Initial number of simulations
# n_init <- 4 #40
# alpha <- 0.75 #0.5 # This is the proportion of particles kept at each step
# pacc <- 0.9 #0.5 # This is the stopping criterion of the algorithm: a small number ensures a better convergence of the algorithm, but at a cost in computing time. Must be 0 < p_acc_min < 1. The smaller, the more strict the criterion.
# 
# ABC_LenormandResult0 <- ABC_sequential(method="Lenormand",
#                                        model=simpact4ABC,
#                                        prior=simpact_prior,
#                                        nb_simul=n_init,
#                                        summary_stat_target=sum_stat_obs,
#                                        alpha=alpha,
#                                        p_acc_min=pacc,
#                                        verbose=FALSE)
# 
# # Time to get a coffee and a biscuit, this will take a while.
# 
# ABC_LenormandResult0

n=3
p=0.7
ABC_rej<-ABC_rejection(model=simpact4ABC, prior=simpact_prior, nb_simul=n,
                       summary_stat_target=sum_stat_obs, tol=p)

# Error: grouping factors must have > 1 sampled level
# Called from: checkNlevels(reTrms$flist, n = n, control)

# Error: Invalid grouping factor specification, ID
# Called from: FUN(X[[i]], ...)

# Error in names(summary.df) <- features.names : 
#   'names' attribute [26] must be the same length as the vector [24]
# In addition: Warning message:
#   In min(smallest.branches, na.rm = TRUE) :
#   no non-missing arguments to min; returning Inf
# Called from: compute.summary.statistics.combined.95(datalist = datalist, work.dir = work.dir, 
#                                                     sub.dir.rename = ABC_DestDir)
# n=3
# ABC_Marjoram_original<-ABC_mcmc(method="Wegmann", model=simpact4ABC,
#                                 prior=simpact_prior, summary_stat_target=sum_stat_obs, n_rec=n)
