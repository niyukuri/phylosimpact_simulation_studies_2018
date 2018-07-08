#######   CALIBRATION    WITH EPIDEMIOLOGICAL AND BEHAVIOURAL FEATURES  ##########


simpact4ABC.classic <- function(inputvector){
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  setwd(paste0(work.dir))
  
  
  set.seed(inputvector[1])
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(expoTree)
  library(data.table)
  library(readr)
  library(lme4)
  library(nlme)
  library(minque) # with lmer
  library(dplyr)
  library(igraph)
  library(robustbase)
  
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
  
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 50, 
                                   population.nummen = 600, 
                                   population.numwomen = 600,
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
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector[3] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # -0.1 c("unif", -0.9, 0)
  cfg.list["population.eyecap.fraction"] <- inputvector[10] # 0.2 c("unif", 0, 0.5)
  #
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[11] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[12] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[13] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[14] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[15] # -0.1386294 c("unif", -0.5, 0)
  cfg.list["person.vsp.toacute.x"] <- inputvector[16] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[17] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[18] # 12 c("unif", 10, 14)
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[19] # -2.7 c("unif", -3.5, -1.7)
  
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
  
  cfg <- cfg.list
  
  
  results <- simpact.run(cfg, 
                         destDir = ABC_DestDir.classic) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  summary.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic)
  
  
  outputvector <- as.numeric(summary.stat) 
  
  
  return(outputvector)
  
}


# Priors

simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 1, 3), c("unif", -0.5, 0.5),
                      c("unif", 2, 4), c("unif", 0, 1), c("unif", -1, 0), c("unif", -0.9, 0), 
                      c("unif", 0, 0.5), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), 
                      c("unif", 0, 0.5), c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9),
                      c("unif", 10, 14),  c("unif", -3.5, -1.7))


# Observed features

# library(robustbase)

source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")

targets.stat <- read.csv("~/Desktop/mastermodeltest/features.matrix.csv")

targets.stat <- targets.stat[,2:length(targets.stat)]

median.targets.stat <-  colMedians(targets.stat)  # Ok # library(robustbase)

classic.target <- colMedians(targets.stat[,11:26])


sum_stat_obs <- as.numeric(classic.target)

# ABC approach

library(EasyABC)

ABC_rej.classic <- ABC_rejection(model = simpact4ABC.classic,
                                 prior = simpact_prior,
                                 summary_stat_target = sum_stat_obs,
                                 nb_simul = 12,
                                 use_seed = TRUE,
                                 seed_count = 1,
                                 n_cluster = 4,
                                 tol = 2/12)

output.params.classic <- as.data.table(ABC_rej.classic$param)
write.csv(output.params.classic, file = "output.params.classic.csv")

output.params.classic <- as.data.frame(as.matrix(output.params.classic))

<<<<<<< HEAD
=======
median.targets.stat.classic <-  colMedians((output.params.classic))  # Ok # library(robustbase)

>>>>>>> master
# 
# Seq.simp <- ABC_sequential(model = simpact4ABC.classic,
#                           method = "Lenormand",
#                           prior = simpact_prior,
#                           summary_stat_target = sum_stat_obs,
#                           nb_simul = 4,
#                           alpha = 0.1,
#                           p_acc_min = 0.03,
#                           use_seed = TRUE,
#                           seed_count = 1,
#                           n_cluster = 4,
#                           inside_prior = FALSE)
# 
# 
# ## MaC approach
<<<<<<< HEAD
# 
# library(RSimpactHelper)
# library(mice)
# library(gsubfn)
# library(data.table)
# library(readcsvcolumns)
# library(randtoolbox)
# library(pcaPP)
# library(glmnet)
# 
# lls <- c(-1, -0.5, 1, -0.5, 2, 0, -1, -0.9, 0, -2, -100, 0, 0, -0.5, 3, 5, 10, -3.5)
# 
# uls <- c(0, 0, 3, 0.5, 4, 1, 0, 0, 0.5, 0, -80, 1, 0.5, 0, 7, 9, 14, -1.7)
# 
# obs.targets <- as.numeric(sum_stat_obs)
# 
# 
# MaC.simp <- MaC(targets.empirical = obs.targets,
#                 RMSD.tol.max = 2,
#                 min.givetomice = 2,
#                 n.experiments = 20,
#                 lls = lls,
#                 uls = uls,
#                 model = simpact4ABC.classic,
#                 strict.positive.params = 0,
#                 probability.params = 0,
#                 method = "norm",
#                 predictorMatrix = "complete",
#                 maxit = 20,
#                 maxwaves = 1,
#                 n_cluster = 8)
=======

library(RSimpactHelper)
library(mice)
library(gsubfn)
library(data.table)
library(readcsvcolumns)
library(randtoolbox)
library(pcaPP)
library(glmnet)

lls <- c(-1, -0.5, 1, -0.5, 2, 0, -1, -0.9, 0, -2, -100, 0, 0, -0.5, 3, 5, 10, -3.5)

uls <- c(0, 0, 3, 0.5, 4, 1, 0, 0, 0.5, 0, -80, 1, 0.5, 0, 7, 9, 14, -1.7)

obs.targets <- as.numeric(sum_stat_obs)


MaC.simp <- MaC(targets.empirical = obs.targets,
                RMSD.tol.max = 2,
                min.givetomice = 20,
                n.experiments = 200,
                lls = lls,
                uls = uls,
                model = simpact4ABC.classic,
                strict.positive.params = c(5, 12, 14:16),
                probability.params = 6,
                method = "norm",
                predictorMatrix = "complete",
                maxit = 20,
                maxwaves = 1,
                n_cluster = 8)
>>>>>>> master



# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.classic <- as.numeric(median.targets.stat.classic)

reps <- 4

# # Input parameters in matrix form reps times (rows).

inputmatrix.classic <- matrix(rep(inputvector.classic, reps), byrow = TRUE, nrow = reps)


epi.metric.calibrates.classic <- simpact.parallel(model = wrapper.test.study.1,
                                               actual.input.matrix = inputmatrix.classic,
                                               seed_count = 124,
                                               n_cluster = 4)


# 
# # save features in the working directory

write.csv(epi.metric.calibrates.classic, file = paste0(work.dir,"/epi.metric.calibrates.classic.csv"))


