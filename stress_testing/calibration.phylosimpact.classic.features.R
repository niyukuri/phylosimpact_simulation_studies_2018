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


work.dir <- "/home/niyukuri/Desktop/calibration" # on laptop

# work.dir <- "/home/niyukuri/Desktop/calibration" # on PC


setwd(paste0(work.dir))


simpact4ABC.classic <- function(inputvector){
  
  
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
  
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/compute.summary.statistics.classic.R")
  
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 50, 
                                   population.nummen = 2000, 
                                   population.numwomen = 2000,
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
  
  
  cfg.list["dissolution.alpha_0"] <- inputvector[1] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [2] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[3] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[6] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[6] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[7] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[7] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[8] # [8] # -0.1 c("unif", -0.9, 0)
  
  # cfg.list["population.eyecap.fraction"] <- inputvector[10] # [9] # 0.2 c("unif", 0, 0.5) # REMOVED in params 
  #
  # # HIV transmission
  # ###################
  #
  
  
  cfg.list["hivtransmission.param.a"] <- inputvector[9] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[10] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[11] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[12] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[13] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[14] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[15] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[16] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[17] # [18] # -2.7 c("unif", -3.5, -1.7)
  
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

simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 1, 3), c("unif", -0.5, 0.5), c("unif", 2, 4), c("unif", 0, 1),
                      c("unif", -1, 0), c("unif", -0.9, 0), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), c("unif", 0, 0.5),
                      c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9), c("unif", 10, 14), c("unif", -3.5, -1.7))

# inputvector <- c(-0.52, -0.05, 2, 0, 2, 0.25, -0.3, -0.1,
#                  # 0.2,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18

sum_stat_obs <- c(-0.02181803 , 0.000001 , 0.000001,  0.00000002,  0.00000002,  0.00000002,  
                  0.00000002 ,  0.05265748,  3.65380066 , 2.49587156,  4.48834260,  0.10748877,  
                  2.13390420,  0.00000000,  1.03904202,  0.01818182, 2.49587156 , 4.48834260 , 
                  0.10748877, 2.13390420, 0.00000000 , 1.03904202 , 35.28732255 ,-2.64059723 , 
                  6.97765300 , 4.00243500, 4.00243500,  4.36124880)





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


# # Time to get a coffee and a biscuit, this will take a while.
# 
# ABC_LenormandResult0

n <- 4
p=0.7

# Can't work when using clusters

# ABC_rej.classic <- ABC_rejection(model = simpact4ABC.classic,
#                                  prior = simpact_prior, 
#                                  nb_simul=n,
#                                  summary_stat_target=sum_stat_obs,
#                                  tol=p,
#                                  use_seed = TRUE,
#                                  seed_count = 0,
#                                  n_cluster = 8)

ABC_rej.classic  <- ABC_rejection(model = simpact4ABC.classic,
                                  prior = simpact_prior, 
                                  nb_simul=n,
                                  summary_stat_target=sum_stat_obs, 
                                  tol=p,
                                  verbose=TRUE)




n_init <- 4 #40
alpha <- 0.75 #0.5 # This is the proportion of particles kept at each step
pacc <- 0.9 #0.5 # This is the stopping criterion of the algorithm: a small number ensures a better convergence of the algorithm, but at a cost in computing time. Must be 0 < p_acc_min < 1. The smaller, the more strict the criterion.



ABC_LenormandResult0 <- ABC_sequential(model=simpact4ABC.classic,
                                       method="Lenormand",
                                       prior=simpact_prior,
                                       nb_simul=n_init,
                                       summary_stat_target=sum_stat_obs,
                                       alpha=alpha,
                                       p_acc_min=pacc,
                                       verbose=TRUE)


# ABC_sequential(model = toy_model,
#                method = "Lenormand",
#                prior = toy_prior,
#                summary_stat_target = sum_stat_obs,
#                nb_simul = 2000,
#                alpha = 0.1,
#                p_acc_min = 0.03,
#                use_seed = TRUE,
#                seed_count = 1,
#                n_cluster = 8,
#                inside_prior = FALSE)


output.params.classic <- as.data.table(ABC_rej.classic$param)

write.csv(output.params.classic, file = "output.params.classic.csv")


