
# Define directory

work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop

setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel)


wrapper.phylo.simpact.master.model <- function(inputvector = input.vector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop
  
  
  
  destDir <- "/home/david/Desktop/mastermodeltest/temp"
  
  
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
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  ######################### Parameter space ###################
  
  # 1.1. Parameter space:
  
  #   1.1.1. Sexual behaviour parameters
  
  # [1] #' @param dissolution.alpha_0: Baseline parameter for relationship dissolution rate. (0.1)
  # [2] #' @param dissolution.alpha_4 : Effect of increasing mean age of the couple on relationship dissolution rate (-0.05)
  # [XX] #' @param formation.hazard.type: Type of hazard function for relationship formation. Choose between "simple", "agegap" and "agegapry".
  # [3] #' @param person.agegap.man.dist.normal.mu: Mean of preferred age differences distribution for men (-4)
  # [4] #' @param person.agegap.woman.dist.normal.mu: Mean of preferred age differences distribution for women (-4)
  # [5] #' @param person.agegap.man.dist.normal.sigma: Standard deviation of preferred age differences distribution for men (3)
  # [6] #' @param person.agegap.woman.dist.normal.sigma: Standard deviation of preferred age differences distribution for women (3)
  # [7] #' @param formation.hazard.agegapry.gap_agescale_man: Effect of male age on preferred age difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
  # [8] #' @param formation.hazard.agegapry.gap_agescale_woman: Effect of male age on preferred age difference (~ 1 - slope in regression model FemaleAge ~ MaleAge) (0.3)
  # [9] #' @param formation.hazard.agegapry.numrel_man: Effect of number of ongoing relationships on the relationship formation rate for men (-0.2)
  # [10] #' @param formation.hazard.agegapry.numrel_woman: Effect of number of ongoing relationships on the relationship formation rate for women (-0.2)
  # [11] #' @param formation.hazard.agegapry.numrel_diff: Effect of absolute difference in number of ongoing relationships on the relationship formation rate (-0.1)
  # [12] #' @param population.eyecap.fraction: Allow for the indication of how many people can a person possible engage in a relation with. (0.2)
  # 
  # 1.1.2. HIV transmission and diasese progression
  
  # [13] #' @param hivtransmission.param.a: Baseline parameter for HIV transmission rate in serodiscordant couples (-1.0352239)
  # [14] #' @param hivtransmission.param.b: Parameter "b" for the linear component of the effect of viral load on the HIV transmission rate in serodiscordant couples (-89.339994)
  # [15] #' @param hivtransmission.param.c: Parameter "c" for the exponential component of the effect of viral load on the HIV transmission rate in serodiscordant couples (0.4948478)
  # [16] #' @param hivtransmission.param.f1: Effect of youngest age on HIV susceptibility (log(5) ~1.6 such that the hazard is x 5 in 15 year olds)
  # [17] #' @param hivtransmission.param.f2: Effect of female age on HIV susceptibility (log(log(2.5) / log(5)) / 5 ~-0.11 such that the hazard is x 2.5 in 20 year olds, compared to the reference (>>25 year olds)
  # [18] #' @param person.vsp.toacute.x: Effect of acute versus chronic HIV infection on infectiousness (5)  # See Bellan PLoS Medicine
  # [19] #' @param person.vsp.toaids.x: Effect of "initial" AIDS stage versus chronic HIV infection on infectiousness (7)
  # [20] #' @param person.vsp.tofinalaids.x: Effect of "final" AIDS stage versus chronic HIV infection on infectiousness (12)
  
  # 1.1.3. Demographic parameters
  
  # [21] #' @param conception.alpha_base: Baseline parameter for conception rate. (-3)
  
  # inputvector: vector of 22 parameters, 1st being seed
  
  
  # age.distr <- agedistr.creator(shape = 5, scale = 65)
  # 
  # cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
  #                                  population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
  #                                  population.nummen = 600, #3000, #600, # 3800, #2500,
  #                                  population.numwomen = 600, # 3000, #600, #4200, #2500,
  #                                  hivseed.time = 10, # 20,
  #                                  hivseed.type = "amount",
  #                                  hivseed.amount = 20, #30,
  #                                  hivseed.age.min = 20,
  #                                  hivseed.age.max = 50,
  #                                  hivtransmission.param.a = -1, # -1,
  #                                  hivtransmission.param.b = -90,
  #                                  hivtransmission.param.c = 0.5,
  #                                  hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
  #                                  hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
  #                                  formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
  #                                  formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
  #                                  formation.hazard.agegapry.meanage = -0.025,
  #                                  formation.hazard.agegapry.gap_factor_man_const = 0,
  #                                  formation.hazard.agegapry.gap_factor_woman_const = 0,
  #                                  formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
  #                                  formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
  #                                  formation.hazard.agegapry.gap_agescale_man = 0.25, #inputvector[3], # 0.25,
  #                                  formation.hazard.agegapry.gap_agescale_woman = 0.25, #inputvector[3], # 0.25,#-0.30000007,#-0.03,
  #                                  debut.debutage = 15,
  #                                  conception.alpha_base = -2.5#inputvector[14]#-2.5#,
  #                                  #person.art.accept.threshold.dist.fixed.value = 0
  # )
  # 
  # # 
  # # inputvector <- c(123, 0.1, -0.05, -4, -4, 3, 3,
  # #                  0.3, 0.3, -0.2, -0.2, -0.1, 0.2,
  # #                  -1.0352239, -89.339994, 0.4948478,
  # #                  1.6, -0.11, 5, 7, 12, -3)
  # 
  # 
  # # Seed for reproducability
  # ##########################
  # 
  # seedid <- inputvector[1]
  # 
  # 
  # # Sexual behaviour
  # ###################
  # 
  # cfg.list["dissolution.alpha_0"] <- inputvector[2]
  # cfg.list["dissolution.alpha_4"] <- inputvector[3]
  # cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  # cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5]
  # cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6]
  # cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[7]
  # cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[8]
  # cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[9]
  # cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[10]
  # cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[11]
  # cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[12]
  # cfg.list["population.eyecap.fraction"] <- inputvector[13]
  # 
  # # HIV transmission
  # ###################
  # 
  # cfg.list["hivtransmission.param.a"] <- inputvector[14]
  # cfg.list["hivtransmission.param.b"] <- inputvector[15]
  # cfg.list["hivtransmission.param.c"] <- inputvector[16]
  # cfg.list["hivtransmission.param.f1"] <- inputvector[17]
  # cfg.list["hivtransmission.param.f2"] <- inputvector[18]
  # cfg.list["person.vsp.toacute.x"] <- inputvector[19]
  # cfg.list["person.vsp.toaids.x"] <- inputvector[20]
  # cfg.list["person.vsp.tofinalaids.x"] <- inputvector[21]
  # 
  # # Demographic
  # ##############
  # 
  # cfg.list["conception.alpha_base"] <- inputvector[22]
  # 
  # 
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # 
  # cfg.list["monitoring.fraction.log_viralload"] <- 0
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  # 
  # # Assumption of nature of sexual network
  # #########################################
  # 
  # cfg.list["population.msm"] = "no"
  # 
  # 
  # ## Add-ons
  # 
  # cfg.list["formation.hazard.agegapry.baseline"] <- 2
  # cfg.list["mortality.aids.survtime.C"] <- 65
  # cfg.list["mortality.aids.survtime.k"] <- -0.2
  # cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  # cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  # cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  # cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  # cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  # 
  # cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  # #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  # cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  # #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  # 
  # 
  # # ART intervention
  # ###################
  # 
  # # ART acceptability paramter and the ART  interventions
  # 
  # cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6 
  # 
  # 
  # # Let's introduce ART, 
  # art.intro <- list()
  # art.intro["time"] <- 25 #25
  # art.intro["diagnosis.baseline"] <- 100
  # art.intro["monitoring.cd4.threshold"] <- 100 
  # 
  # # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  # 
  # art.intro2 <- list()
  # art.intro2["time"] <- 25 + 5 #25 + 5 = 30
  # art.intro2["monitoring.cd4.threshold"] <- 200
  # 
  # art.intro3 <- list()
  # art.intro3["time"] <- 25 + 8 #25 + 8 = 33
  # art.intro3["monitoring.cd4.threshold"] <- 350
  # 
  # art.intro4 <- list()
  # art.intro4["time"] <- 25 + 11 #25 + 11 = 36
  # art.intro4["monitoring.cd4.threshold"] <- 500
  # 
  # art.intro5 <- list()
  # art.intro5["time"] <- 25 + 13 #25 + 13 = 38
  # art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  # 
  # interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)
  # 
  # intervention <- interventionlist
  # 
  # 
  
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                   # population.msm = "no",
                                   population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                   population.nummen = 100, #3000, #600, # 3800, #2500,
                                   population.numwomen = 100, # 3000, #600, #4200, #2500,
                                   hivseed.time = 10, # 20,
                                   hivseed.type = "amount",
                                   hivseed.amount = 20, #30,
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   hivtransmission.param.a = -1, # -1,
                                   hivtransmission.param.b = -90,
                                   hivtransmission.param.c = 0.5,
                                   hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
                                   hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                   formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                   formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   formation.hazard.agegapry.gap_factor_man_const = 0,
                                   formation.hazard.agegapry.gap_factor_woman_const = 0,
                                   formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                   formation.hazard.agegapry.gap_agescale_man = 0.25, #inputvector[3], # 0.25,
                                   formation.hazard.agegapry.gap_agescale_woman = 0.25, #inputvector[3], # 0.25,#-0.30000007,#-0.03,
                                   debut.debutage = 15,
                                   conception.alpha_base = -2.5#inputvector[14]#-2.5#,
                                   #person.art.accept.threshold.dist.fixed.value = 0
  )
  
  
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
  
  intervention <- interventionlist # scenario(interventionlist, tasp.indicator)
  
  cfg.list["hivtransmission.param.f1"] <- log(inputvector[2])
  cfg.list["hivtransmission.param.f2"] <- log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[3]
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[3]
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
  cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
  cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]
  
  #cfg <- cfg.list
  
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  # cfg["monitoring.fraction.log_viralload"] <- 0.3
  cfg.list["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine
  
  seedid <- inputvector[1]
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
  cfg.list["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
  cfg.list["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[11]
  
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
  cfg.list["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
  cfg.list["dissolution.alpha_0"] <- inputvector[15]
  cfg.list["dissolution.alpha_4"] <- inputvector[16]
  
  
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
  
  sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Run Simpact
  ##############
  
  
  results <- simpact.run(configParams = cfg.list,
                         destDir = sub.dir.rename,
                         agedist = age.distr,
                         seed = seedid,
                         intervention = intervention)
  
  
  
  datalist.agemix <- readthedata(results)
  
  
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  
  
  ############################ METRICS: TRANSMISSION NETWORK CHARACTERISTICS #####################
  
  # 1. Incidence trend #
  ######################
  
  incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25), timewindow = c(10, 40))
  
  METRICS.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
  
  incidence.df.25.34 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 35), timewindow = c(10, 40))
  
  METRICS.incidence.df.25.34 <- incidence.df.25.34$incidence[3]
  
  incidence.df.35.44 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(35, 45), timewindow = c(10, 40))
  
  METRICS.incidence.df.35.44 <- incidence.df.35.44$incidence[3]
  
  # c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
  #   METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average)
  # 
  # 2. Age mixing in transmissions #
  ##################################
  
  agemix.df <- agemixing.trans.df(datalist = datalist.agemix, 
                                  trans.network = simpact.trans.net)
  
  agemix.fit <- fit.agemix.trans(datatable = agemix.df)
  
  coef.inter <- fixef(agemix.fit)
  
  METRICS.age.mix.trans.interc <- coef.inter[[1]]
  METRICS.age.mix.slope <- coef.inter[2]
  
  # 3. Onward transmissions #
  ###########################
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = simpact.trans.net)
  
  METRICS.transm.average <- mean(transm.count)
  
  
  
  #################################### Features #############################
  
  # 1.2. Features from sexual and transmission network
  
  
  # 
  # 1.2.1. Demographic feature:
  
  #   (i) Population growth rate (pop.growth.calculator function)
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = c(0, datalist.agemix$itable$population.simtime[1]))
  
  
  # 1.2.2. Transmission features:	
  
  #   (i) Prevalence (prevalence.calculator function)
  
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  hiv.prev.25.34.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 35),
                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.25.34.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 35),
                                              timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  hiv.prev.35.44.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(35, 45),
                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  hiv.prev.35.44.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(35, 45),
                                              timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[1]
  
  
  # (ii) Transmission 	rate (transmission.rate.calculator function)
  
  transm.rate <- transmission.rate.calculator(datalist = datalist.agemix,
                                              timewindow = c(10, 40), 
                                              int = FALSE, by=1)
  
  
  
  # (iii) ART coverage
  
  # cov.vector <- ART.coverage.vector.creator(datalist = datalist.agemix,
  #                                          agegroup = c(15, 50))
  # plot(cov.vector)
  
  
  
  # 1.2.3. Sexual behaviour features:
  
  #  (i) Relationship 	rate (relationship.rate.calculator function) - rate of new relationship formation (partner turnover rate):  	
  
  
  relas.rate <- relationship.rate.calculator(datalist = datalist.agemix,
                                             timewindow = c(10, 40), 
                                             int = FALSE, by=1)
  
  # (ii) Relationship per person per year
  
  relsperpersonperyear <- nrow(datalist.agemix$rtable) / (nrow(datalist.agemix$ptable)/2) / cfg.list$population.simtime
  
  # (iv) SD age gap between couples
  
  agegapsd <- sd(datalist.agemix$rtable$AgeGap)
  
  
  # (v) Age mixing in relationships
  
  # 
  # agemix.df <- agemix.df.maker(datalist.agemix)
  # 
  # agemix.model <- pattern.modeller(dataframe = agemix.df,
  #                                  agegroup = c(15, 50),
  #                                  timepoint = datalist.agemix$itable$population.simtime[1],
  #                                  timewindow = 3)#1)#3)
  # 
  # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  # 
  # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
  # 
  # bignumber <- NA # let's try if NA works (instead of 9999 for example)
  # 
  # AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  # SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  # #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
  # slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
  # WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
  # 
  # BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
  # 
  # intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
  # 
  
  # age.scatter.df <- agemix.model[[1]]
  
  #  (iii) Point 	prevalence of concurrency in the adult population:
  
  # Concurrency point prevalence 6 months before a survey, among men
  
  pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                     timepoint = datalist.agemix$itable$population.simtime[1] - 0.5)
  
  # c(growthrate, hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
  #   hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, transm.rate, # cov.vector
  #   relas.rate,  relsperpersonperyear, agegapsd,
  #   
  #   # AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
  #   pp.cp.6months.male
  #   )
  # 
  
  

  ###############################
  # Step 3: Sequence simulation #
  ###############################
  
  
  trans.net <- simpact.trans.net # all transmission networks
  
  
  dirseqgen <- work.dir
  
  seeds.num <- inputvector[1]
  
  # Sequence simulation is done for at least a transmission network with 6 individuals
  # This means that limitTransmEvents equal at least 7
  
  sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                 sub.dir.rename = sub.dir.rename,
                                 simpact.trans.net = simpact.trans.net,
                                 seq.gen.tool = "seq-gen",
                                 datalist = datalist,
                                 seeds.num = seeds.num,
                                 endpoint = 40,
                                 limitTransmEvents = 3, # no less than 7
                                 hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                 clust = FALSE) # hiv.seq.file lodged in work.dir
  
  
  
  #####################################################
  # Step 4: Construct time stamped phylogenetic trees #
  #####################################################
  
  
  dirfasttree <- work.dir
  
  
  if(file.exists(paste0(sub.dir.rename, "/C.Epidemic_seed.seq.bis.sim.nwk.fasta"))==TRUE){
    
    # check if the run has simulated sequence,
    # in other words: e have a transmission network with at least 6 individuals
    
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = sub.dir.rename,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(sub.dir.rename,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    ########################################
    ### FEATURES FROM PHYLOGENETIC TREE ####
    ########################################
    
    
    # 1.3. Features from phylogenetic tree:
    
    # library(phytools)
    
    tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
    
    
    # tree <- read.tree("calibrated.tree.save.nwk")
    
    # Mean height of internal nodes
    
    H <- nodeHeights(tree.cal) # similar to node.depth.edgelength(tree)
    
    # It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
    # is represented twice and only twice. Thus, if we exclude the root node (zero height), 
    # we can just take the mean of H[,1].
    
    mean.feature <- mean(sort(H[,1])[3:nrow(H)]) # important
    
    # library(phyloTop)
    
    colless.feature <- colless.phylo(tree.cal, normalise = TRUE)
    
    sackin.feature <- sackin.phylo(tree.cal, normalise = TRUE)
    
    
    Depths <- getDepths(tree.cal) # depth of tips and nodes
    
    mean.tipsDepths.feature <- mean(Depths$tipDepths)
    
    mean.nodesDepths.feature <- mean(Depths$nodeDepths)
    
    maxHeight.feature <- maxHeight(tree.cal, normalise = TRUE)
    
    # Estimating confidence intervals for rates and dates using a parametric bootstrap
    pb <- parboot.treedater(tree.calib.LTT) # Lineage Through Time
    
    # Lineages through time
    LTT <- plot.parboot.ltt.dat(pb)
    
    lb.mean.feature <- mean(LTT$lb) # mean of low values of LTT
    lb.median.feature <- median(LTT$lb) # median of low values of LTT
    
    ub.mean.feature <- mean(LTT$ub) # mean of upper values of LTT
    ub.median.feature <- median(LTT$ub) # median of upper values of LTT
    
    median.mean.feature <- mean(LTT$median) # mean of medians of values of LTT
    median.median.feature <- median(LTT$median) # median of medians of values of LTT
    
    
    summary.df <- c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
                    METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average,
                    
                    growthrate, 
                    
                    hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
                    hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, transm.rate, # cov.vector
                    relas.rate,  relsperpersonperyear, agegapsd, 
                    #AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                    pp.cp.6months.male,
                    
                    mean.feature, colless.feature, sackin.feature, mean.tipsDepths.feature, mean.nodesDepths.feature,
                    maxHeight.feature, lb.mean.feature, lb.median.feature, ub.mean.feature, ub.median.feature,
                    median.mean.feature, median.median.feature)
    
    
    features.names <- c("METRICS.incidence.df.15.24", "METRICS.incidence.df.25.34", "METRICS.incidence.df.35.44",
                        "METRICS.age.mix.trans.interc", "METRICS.age.mix.slope", "METRICS.transm.average",
                        
                        "Pop.growthrate", 
                        
                        "hiv.prev.lt25.women", "hiv.prev.lt25.men", "hiv.prev.25.34.women",
                        "hiv.prev.25.34.men", "hiv.prev.35.44.women", "hiv.prev.35.44.men", "transm.rate", # cov.vector
                        "relas.rate",  "relsperpersonperyear", "agegapsd",
                        #"AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        "pp.cp.6months.male",
                        
                        "meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", "mean.nodesDepths.feature",
                        "maxHeight.feature", "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
                        "LTT.median.mean.feature", "LTT.median.median.feature")
    
    names(summary.df) <- features.names
    
    return(summary.df)
    
    
  }else{
    
    features.names <- c("METRICS.incidence.df.15.24", "METRICS.incidence.df.25.34", "METRICS.incidence.df.35.44",
                        "METRICS.age.mix.trans.interc", "METRICS.age.mix.slope", "METRICS.transm.average",
                        
                        "Pop.growthrate", 
                        
                        "hiv.prev.lt25.women", "hiv.prev.lt25.men", "hiv.prev.25.34.women",
                        "hiv.prev.25.34.men", "hiv.prev.35.44.women", "hiv.prev.35.44.men", "transm.rate", # cov.vector
                        "relas.rate",   "relsperpersonperyear", "agegapsd",
                        #"AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        "pp.cp.6months.male",
                        
                        "meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", "mean.nodesDepths.feature",
                        "maxHeight.feature", "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
                        "LTT.median.mean.feature", "LTT.median.median.feature")
    
    summary.NA <- rep(NA,12)
    
    summary.df.classic <- c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
                            METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average,
                            
                            growthrate, 
                            
                            hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
                            hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, transm.rate, # cov.vector
                            relas.rate, 
                            # AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                            pp.cp.6months.male)
    
    summary.df <- c(summary.df.classic,  summary.NA)
    
    names(summary.df) <- features.names
    
    return(summary.df)
    
  }
  
  unlink(paste0(sub.dir.rename), recursive = FALSE) # sub.dir.rename
  
}


# 
# inputvector <- c(123, 0.1, -0.05, -4, -4, 3, 3,
#                  0.3, 0.3, -0.2, -0.2, -0.1, 0.2,
#                  -1.0352239, -89.339994, 0.4948478,
#                  1.6, -0.11, 5, 7, 12, -3)


inputvector <- c(1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)



#####   PARALLELIZATION   ######
################################


# replication number

reps <- 4


# Input parameters in matrix form reps times (rows).
inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)



sim.start.time <- proc.time()[3]

features.matrix <- phylo.simpact.parallel(model = wrapper.phylo.simpact.master.model,
                                          actual.input.matrix = inputmatrix,
                                          seed_count = 123,
                                          n_cluster = 4)

sim.end.time <- proc.time()[3] - sim.start.time

print(paste0("Simulation time: ", round(sim.end.time/60,2), " minutes"))



# save features in the working directory

write.csv(features.matrix, file = paste0(work.dir,"/features.matrix.csv"))



#######   CALIBRATION

library(EasyABC)

simpact4ABC <- function(inputvector.cal){
  cfg <- cfg.list
  cfg["formation.hazard.agegap.baseline"] <- inputvector.cal[1]
  cfg["formation.hazard.agegap.gap_factor_man"] <- inputvector.cal[2]
  cfg["formation.hazard.agegap.gap_factor_woman"] <- inputvector.cal[2]
  results <- simpact.run(cfg, ABC_DestDir) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  datalist <- readthedata(results)
  relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  agegapsd <- sd(datalist$rtable$AgeGap)
  outputvector <- c(relsperpersonperyear, agegapsd)
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

simpact_prior <- list(c("unif", 0.8, 5), c("unif", -1, 0))
# Lastly, we specify the target summary statistic
sum_stat_obs <- c(1, 2)

# Now we try to run a sequential ABC scheme, according to the method proposed by Lenormand et al. 2013
# Maxime Lenormand, Franck Jabot and Guillaume Deffuant. Adaptive approximate Bayesian computation for complex models. Comput Stat (2013) 28:2777–2796 DOI 10.1007/s00180-013-0428-3


# Initial number of simulations
n_init <- 10 #40
alpha <- 0.25 #0.5 # This is the proportion of particles kept at each step
pacc <- 0.2 #0.5 # This is the stopping criterion of the algorithm: a small number ensures a better convergence of the algorithm, but at a cost in computing time. Must be 0 < p_acc_min < 1. The smaller, the more strict the criterion.

ABC_LenormandResult0 <- ABC_sequential(method="Lenormand",
                                       model=simpact4ABC,
                                       prior=simpact_prior,
                                       nb_simul=n_init,
                                       summary_stat_target=sum_stat_obs,
                                       alpha=alpha,
                                       p_acc_min=pacc,
                                       verbose=FALSE)

# Time to get a coffee and a biscuit, this will take a while.

ABC_LenormandResult0

