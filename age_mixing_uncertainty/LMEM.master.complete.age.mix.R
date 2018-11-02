# Master model for simulation of age-mixing patterns


# Define directory

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


# inputvector <- c(777, -0.52, -0.05, 2.8, 0, 3, 0.25, -0.3, -0.1,
#                  # 0.2,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18


# 
# inputvector <- c(10, -0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
#                  # 0.2,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18


LMEM.master.model.age.mixing.toy1 <- function(inputvector = input.vector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/LMEMphylo.AR.groups.fun.agemix.R")
  source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/CAR.groups.fun.agemixBIS.R")
  source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/AR.groups.fun.agemixBIS.R")
  source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/LMEMphylo.CAR.groups.fun.agemix.R") 
  
  
  # work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  # destDir <- "/home/david/Desktop/mastermodeltest/temp" # on laptop
  
  # destDir <- "/home/niyukuri/Desktop/mastermodeltest/temp" # on PC
  
  
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
  
  library(Rsamtools)
  
  library(robustbase)
  
  
  
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
  # cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1a
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
  # art.intro2["monitoring.cd4.threshold"] <- 200a
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
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 50, 
                                   population.nummen = 1200, 
                                   population.numwomen = 1200,
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
  
  # cfg.list["population.eyecap.fraction"] <- inputvector[10] # [9] # 0.2 c("unif", 0, 0.5) # REMOVED in params 
  #
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
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  
  DataListALL <- readthedata(results)
  
  datalist.agemix <- DataListALL
  
  # datalist.agemix <- get(load("datalist.agemix.RData"))
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
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
                                 seeds.num = seeds.num,
                                 endpoint = 40,
                                 limitTransmEvents = 7, # no less than 7
                                 hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                 clust = FALSE) # hiv.seq.file lodged in work.dir
  
  # Transform the sequence format to be handled by ClusterPicker
  sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
  write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")
  
  
  
  
  # (i) Age mixing in relationships
  
  # 
  
  agemix.rels.df <- agemix.df.maker(datalist.agemix)
  
  # 
  agemix.model <- pattern.modeller(dataframe = agemix.rels.df,
                                   agegroup = c(15, 50),
                                   timepoint = 40, # datalist.agemix$itable$population.simtime[1],
                                   timewindow = 10)#1)#3)
  # 
  # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  #
  # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
  
  data = dplyr::filter(agemix.model[[1]], Gender =="male")
  
  if( nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
    
    men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                     data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                     REML = TRUE,
                     control=lmerControl(check.nobs.vs.nlev = "ignore",
                                         check.nobs.vs.rankZ = "ignore",
                                         check.nobs.vs.nRE="ignore"))
    
    bignumber <- NA # let's try if NA works (instead of 9999 for example)
    AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
    slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
    WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
    
    BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
    
    intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
    
    # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
    ## AAD: average age difference across all relationship
    ## VAD: variance of these age differences
    ## SDAD: standard deviation of age differences
    ## BSD: between-subject standard deviation of age differences
    
    mix.rels.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
  }
  
  # age.scatter.df <- agemix.model[[1]]
  
  #  (ii) Point 	prevalence of concurrency in the adult population:
  
  # Concurrency point prevalence 6 months before a survey, among men
  
  
  pp.cp.6months.male <- tryCatch(concurr.pointprev.calculator(datalist = datalist.agemix,
                                                              timepoint = 40 - 0.5), error=function(e) return(NA))
  
  
  # (iii) Prevalence
  
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = 40)$pointprevalence[2]
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = 40)$pointprevalence[1]
  
  hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 40),
                                                timepoint = 40)$pointprevalence[2]
  hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 40),
                                              timepoint = 40)$pointprevalence[1]
  
  hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(40, 50),
                                                timepoint = 40)$pointprevalence[2]
  hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(40, 50),
                                              timepoint = 40)$pointprevalence[1]
  
  
  # (iv) Incidence
  
  incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25), timewindow = c(30, 40))
  
  METRICS.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
  
  METRICS.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
  METRICS.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]
  
  
  incidence.df.25.39 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 40), timewindow = c(30, 40))
  
  METRICS.incidence.df.25.39 <- incidence.df.25.39$incidence[3]
  
  METRICS.incidence.df.25.39.men <- incidence.df.25.39$incidence[1]
  METRICS.incidence.df.25.39.women <- incidence.df.25.39$incidence[2]
  
  
  incidence.df.40.49 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 40), timewindow = c(30, 40))
  
  METRICS.incidence.df.40.49 <- incidence.df.40.49$incidence[3]
  
  METRICS.incidence.df.40.49.men <- incidence.df.40.49$incidence[1] # res
  METRICS.incidence.df.40.49.women <- incidence.df.40.49$incidence[2] # res
  
  
  
  summary.epidemic.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                           hiv.prev.25.40.women, hiv.prev.25.40.men,
                           hiv.prev.40.50.women, hiv.prev.40.50.men, 
                           mix.rels.dat,
                           pp.cp.6months.male,
                           
                           METRICS.incidence.df.15.24.men, METRICS.incidence.df.15.24.women, 
                           METRICS.incidence.df.25.39.men, METRICS.incidence.df.25.39.women,
                           METRICS.incidence.df.40.49.men, METRICS.incidence.df.40.49.women)
  
  
  
  
  ##########################################################
  # Step 3: Data and age-mixing in transmissions #
  ##########################################################
  
  
  # Make a data table from all transmissions networks with at least limitTransmEvents transmission events
  
  agemixing.df <- agemixing.trans.df(trans.network = simpact.trans.net,
                                     limitTransmEvents = 7)
  
  
  # IDs of individuals infected in the time window of the study
  
  IDs.study <- new.transmissions.dat(datalist = datalist.agemix, 
                                     time.window=c(30,40))
  
  agemixing.df.IDs <- dplyr::filter(agemixing.df, agemixing.df$RecId%in%IDs.study)
  
  
  # True age mixing in transmissions within differewnt age groups seen by MELM #
  ##############################################################################
  
  
  if( nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 ){
    
    
    # Simple LMM
    #############
    # 
    # fit.lme.agemixing <- lme(AgeInfecRec ~ GenderRec, data = agemixing.df.IDs, random = ~ 1|DonId)
    # 
    # 
    # 
    # a <- coef(summary(fit.lme.agemixing))[1] # average age in transmission clusters
    # 
    # beta <- coef(summary(fit.lme.agemixing))[2] # average age difference in transmission clusters: 
    # # seen as bridge width which shows potential cross-generation transmission
    # 
    # 
    # b1 <- as.numeric(VarCorr(fit.lme.agemixing)[3]) # between cluster variation
    # 
    # b2 <- as.numeric(VarCorr(fit.lme.agemixing)[4]) # within cluster variation
    # 
    # lme.val <- c(a, beta, b1, b2)
    # 
    # names(lme.val) <-  c("av.age.male", "gendEffect.clust", "between.transm.var", "within.transm.var")
    # 
    # 
    # 
    
    # Heteroscedasticity: model variance component of within-group errors
    ######################################################################
    
    
    het.fit.lme.agemixing <- lme(AgeInfecRec ~ GenderRec, data = agemixing.df.IDs, random = ~ 1|DonId,
                                 weights = varIdent( c("1" = 0.5), ~ 1 |GenderRec ))
    
    
    het.a <- coef(summary(het.fit.lme.agemixing))[1] # average age in transmission clusters
    
    het.beta <- coef(summary(het.fit.lme.agemixing))[2] # average age difference in transmission clusters: 
    # seen as bridge width which shows potential cross-generation transmission
    
    
    het.b1 <- as.numeric(VarCorr(het.fit.lme.agemixing)[3]) # between cluster variation
    
    het.b2 <- as.numeric(VarCorr(het.fit.lme.agemixing)[4]) # within cluster variation
    
    
    # SD for the two strata
    
    unique.val.strat <- unique(attributes(het.fit.lme.agemixing$modelStruct$varStruct)$weights)
    
    het.fit.lme.agemixing$modelStruct$varStruct
    
    # reference group: female == 1
    delta.female <- 1
    
    female.val <- unique.val.strat[1]
    male.val <- unique.val.strat[2]
    
    delta.male <- female.val/male.val # delta_ref_group / val
    
    SD.female <- as.numeric(VarCorr(het.fit.lme.agemixing)[4])
    SD.male <- delta.male * SD.female 
    
    
    het.lme.val <- c(het.a, het.beta, het.b1, het.b2, SD.female, SD.male)
    
    names(het.lme.val) <-  c("het.av.age.male", "het.gendEffect.clust", "het.between.transm.var", "het.within.transm.var", "SD.female", "SD.male")
    
    
    flag.lme <- NA
    
    if(abs(het.lme.val[[2]]) > 5){ # If average age difference is greater than 5, there is a cross-generation transmission
      flag.lme <- 1
    }else{
      flag.lme <- 0
    }
    
  }else{
    
    names.lme.val <-  c("het.av.age.male", "het.gendEffect.clust", "het.between.transm.var", "het.within.transm.var", "SD.female", "SD.male")
    
    # c("av.age.male", "gendEffect.clust", "between.transm.var", "within.transm.var")
    
    het.lme.val <- rep(NA, length(names.lme.val))
    
    names(het.lme.val) <- names.lme.val
    
    flag.lme <- NA
  }
  
  
  
  # val.names <- c("num.men.15.25", "num.women.15.25",
  #                "num.men.25.40", "num.women.25.40",
  #                "num.men.40.50", "num.women.40.50",
  #                
  #                "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
  #                "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
  #                "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
  #                
  #                "mean.AD", "median.AD", "sd.AD")
  
  CAR.100 <- CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,
                                      datalist = datalist.agemix,
                                      limitTransmEvents = 7,
                                      timewindow = c(30,40),
                                      seq.cov = 100,
                                      age.group.15.25 = c(15,25),
                                      age.group.25.40 = c(25,40),
                                      age.group.40.50 = c(40,50))
  
  flag.women.val <- CAR.100[13] # partners.men.40.50.w.15.25
  
  flag.men.val <- CAR.100[9] # partners.men.15.25.w.40.50
  
  flag.women <- NA
  flag.men <- NA
  
  if(flag.women.val >=1){ # If we have at least one transmission from older men to younger women 
    flag.women <- 1
  }else{
    flag.women <- 0
  }
  
  if(flag.men.val >=1){  # If we have at least one transmission from older women to younger men 
    flag.men <- 1
  }else{
    flag.men <- 0
  }
  
  
  
  # Age difference statistics #
  #############################
  AD <- abs(abs(agemixing.df.IDs$TOBDon) - abs(agemixing.df.IDs$TOBRec))
  mean.AD <- mean(AD)
  med.AD <- median(AD)
  sd.AD <- sd(AD)
  
  
  
  # III. Transmission clusters with MCAR
  
  dirfasttree <- work.dir
  
  
  
  transm.clust.MCAR.cov.35 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 35,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.35.val <- sapply(transm.clust.MCAR.cov.35, mean)
  
  
  
  
  
  transm.clust.MCAR.cov.40 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree,  
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 40,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.40.val <- sapply(transm.clust.MCAR.cov.40, mean)
  
  
  
  transm.clust.MCAR.cov.45 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,   dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 45,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.45.val <- sapply(transm.clust.MCAR.cov.45, mean)
  
  transm.clust.MCAR.cov.50 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 50,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.50.val <- sapply(transm.clust.MCAR.cov.50, mean)
  
  
  transm.clust.MCAR.cov.55 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree,  
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 55,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.55.val <- sapply(transm.clust.MCAR.cov.55, mean)
  
  
  
  transm.clust.MCAR.cov.60 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree,  
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 60,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.60.val <- sapply(transm.clust.MCAR.cov.60, mean)
  
  
  
  transm.clust.MCAR.cov.65 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 65,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.65.val <- sapply(transm.clust.MCAR.cov.65, mean)
  
  
  transm.clust.MCAR.cov.70 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree,  
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 70,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.70.val <- sapply(transm.clust.MCAR.cov.70, mean)
  
  transm.clust.MCAR.cov.75 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 75,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.75.val <- sapply(transm.clust.MCAR.cov.75, mean)
  
  transm.clust.MCAR.cov.80 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 80,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.80.val <- sapply(transm.clust.MCAR.cov.80, mean)
  
  transm.clust.MCAR.cov.85 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 85,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.85.val <- sapply(transm.clust.MCAR.cov.85, mean)
  
  transm.clust.MCAR.cov.90 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                       work.dir = work.dir,   
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 90,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.90.val <- sapply(transm.clust.MCAR.cov.90, mean)
  
  transm.clust.MCAR.cov.95 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                       work.dir = work.dir,  
                                                                       dirfasttree = dirfasttree, 
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       limitTransmEvents = 7,
                                                                       timewindow = c(30,40),
                                                                       seq.cov = 95,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.95.val <- sapply(transm.clust.MCAR.cov.95, mean)
  
  
  transm.clust.MCAR.cov.100 <- tryCatch(LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                        work.dir = work.dir,  
                                                                        dirfasttree = dirfasttree, 
                                                                        sub.dir.rename = sub.dir.rename,
                                                                        limitTransmEvents = 7,
                                                                        timewindow = c(30,40),
                                                                        seq.cov = 100,
                                                                        age.group.15.25 = c(15,25),
                                                                        age.group.25.40 = c(25,40),
                                                                        age.group.40.50 = c(40,50)),
                                        error=function(e) return(rep(NA, 8)))
  transm.clust.MCAR.cov.100.val <- sapply(transm.clust.MCAR.cov.100, mean)
  
  
  
  
  
  
  
  # IV. Transmission clusters with MAR
  
  # IV.a 
  
  
  transm.clust.AR.a.cov.35 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 35,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.35.val <- sapply(transm.clust.AR.a.cov.35, mean)
  
  transm.clust.AR.a.cov.40 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 40,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.40.val <- sapply(transm.clust.AR.a.cov.40, mean)
  
  transm.clust.AR.a.cov.45 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 45,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.45.val <- sapply(transm.clust.AR.a.cov.45, mean)
  
  transm.clust.AR.a.cov.50 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 50,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.50.val <- sapply(transm.clust.AR.a.cov.50, mean)
  
  transm.clust.AR.a.cov.55 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 55,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.55.val <- sapply(transm.clust.AR.a.cov.55, mean)
  
  transm.clust.AR.a.cov.60 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 60,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.60.val <- sapply(transm.clust.AR.a.cov.60, mean)
  
  
  
  transm.clust.AR.a.cov.65 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 65,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.65.val <- sapply(transm.clust.AR.a.cov.65, mean)
  
  transm.clust.AR.a.cov.70 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 70,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.70.val <- sapply(transm.clust.AR.a.cov.70, mean)
  
  transm.clust.AR.a.cov.75 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 75,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.75.val <- sapply(transm.clust.AR.a.cov.75, mean)
  
  
  
  transm.clust.AR.a.cov.80 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 80,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.80.val <- sapply(transm.clust.AR.a.cov.80, mean)
  
  
  
  transm.clust.AR.a.cov.85 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 85,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.85.val <- sapply(transm.clust.AR.a.cov.85, mean)
  
  
  transm.clust.AR.a.cov.90 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 90,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.90.val <- sapply(transm.clust.AR.a.cov.90, mean)
  
  
  
  transm.clust.AR.a.cov.95 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 95,
                                                                      seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.a.cov.95.val <- sapply(transm.clust.AR.a.cov.95, mean)
  
  
  
  
  # IV.b
  
  transm.clust.AR.b.cov.35 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 35,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.35.val <- sapply(transm.clust.AR.b.cov.35, mean)
  
  transm.clust.AR.b.cov.40 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 40,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.40.val <- sapply(transm.clust.AR.b.cov.40, mean)
  
  transm.clust.AR.b.cov.45 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 45,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.45.val <- sapply(transm.clust.AR.b.cov.45, mean)
  
  transm.clust.AR.b.cov.50 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 50,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.50.val <- sapply(transm.clust.AR.b.cov.50, mean)
  
  transm.clust.AR.b.cov.55 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 55,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.55.val <- sapply(transm.clust.AR.b.cov.55, mean)
  
  transm.clust.AR.b.cov.60 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 60,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.60.val <- sapply(transm.clust.AR.b.cov.60, mean)
  
  
  
  transm.clust.AR.b.cov.65 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 65,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.65.val <- sapply(transm.clust.AR.b.cov.65, mean)
  
  transm.clust.AR.b.cov.70 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 70,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.70.val <- sapply(transm.clust.AR.b.cov.70, mean)
  
  transm.clust.AR.b.cov.75 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 75,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.75.val <- sapply(transm.clust.AR.b.cov.75, mean)
  
  
  
  transm.clust.AR.b.cov.80 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 80,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.80.val <- sapply(transm.clust.AR.b.cov.80, mean)
  
  
  
  transm.clust.AR.b.cov.85 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 85,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.85.val <- sapply(transm.clust.AR.b.cov.85, mean)
  
  
  transm.clust.AR.b.cov.90 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 90,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.90.val <- sapply(transm.clust.AR.b.cov.90, mean)
  
  
  
  transm.clust.AR.b.cov.95 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 95,
                                                                      seq.gender.ratio = 0.3, # where women have 30% of being sampled & men have 70%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.b.cov.95.val <- sapply(transm.clust.AR.b.cov.95, mean)
  
  
  
  # IV.c
  
  transm.clust.AR.c.cov.35 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 35,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.35.val <- sapply(transm.clust.AR.c.cov.35, mean)
  
  transm.clust.AR.c.cov.40 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 40,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.40.val <- sapply(transm.clust.AR.c.cov.40, mean)
  
  transm.clust.AR.c.cov.45 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir, 
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 45,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.45.val <- sapply(transm.clust.AR.c.cov.45, mean)
  
  transm.clust.AR.c.cov.50 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 50,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.50.val <- sapply(transm.clust.AR.c.cov.50, mean)
  
  transm.clust.AR.c.cov.55 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 55,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.55.val <- sapply(transm.clust.AR.c.cov.55, mean)
  
  transm.clust.AR.c.cov.60 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 60,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.60.val <- sapply(transm.clust.AR.c.cov.60, mean)
  
  
  
  transm.clust.AR.c.cov.65 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 65,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.65.val <- sapply(transm.clust.AR.c.cov.65, mean)
  
  transm.clust.AR.c.cov.70 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 70,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.70.val <- sapply(transm.clust.AR.c.cov.70, mean)
  
  transm.clust.AR.c.cov.75 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 75,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.75.val <- sapply(transm.clust.AR.c.cov.75, mean)
  
  
  
  transm.clust.AR.c.cov.80 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,   
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 80,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.80.val <- sapply(transm.clust.AR.c.cov.80, mean)
  
  
  
  transm.clust.AR.c.cov.85 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree, 
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 85,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.85.val <- sapply(transm.clust.AR.c.cov.85, mean)
  
  
  transm.clust.AR.c.cov.90 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net, 
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 90,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.90.val <- sapply(transm.clust.AR.c.cov.90, mean)
  
  
  
  transm.clust.AR.c.cov.95 <- tryCatch(LMEMphylo.AR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
                                                                      work.dir = work.dir,  
                                                                      dirfasttree = dirfasttree,  
                                                                      sub.dir.rename = sub.dir.rename,
                                                                      limitTransmEvents = 7,
                                                                      timewindow = c(30,40),
                                                                      seq.cov = 95,
                                                                      seq.gender.ratio = 0.5, # where women have 50% of being sampled & men have 50%
                                                                      age.group.15.25 = c(15,25),
                                                                      age.group.25.40 = c(25,40),
                                                                      age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 8)))
  transm.clust.AR.c.cov.95.val <- sapply(transm.clust.AR.c.cov.95, mean)
  
  
  
  # V : True age mixing in transmission network with MCAR
  
  
  true.transm.MCAR.cov.35 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 35,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.35.val <- sapply(true.transm.MCAR.cov.35, mean)
  
  
  
  true.transm.MCAR.cov.40 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 40,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.40.val <- sapply(true.transm.MCAR.cov.40, mean)
  
  
  
  true.transm.MCAR.cov.45 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 45,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.45.val <- sapply(true.transm.MCAR.cov.45, mean)
  
  true.transm.MCAR.cov.50 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 50,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.50.val <- sapply(true.transm.MCAR.cov.50, mean)
  
  
  true.transm.MCAR.cov.55 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 55,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.55.val <- sapply(true.transm.MCAR.cov.55, mean)
  
  
  
  true.transm.MCAR.cov.60 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 60,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.60.val <- sapply(true.transm.MCAR.cov.60, mean)
  
  
  
  true.transm.MCAR.cov.65 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 65,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.65.val <- sapply(true.transm.MCAR.cov.65, mean)
  
  
  true.transm.MCAR.cov.70 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 70,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.70.val <- sapply(true.transm.MCAR.cov.70, mean)
  
  true.transm.MCAR.cov.75 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,  datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 75,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.75.val <- sapply(true.transm.MCAR.cov.75, mean)
  
  true.transm.MCAR.cov.80 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 80,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.80.val <- sapply(true.transm.MCAR.cov.80, mean)
  
  true.transm.MCAR.cov.85 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 85,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.85.val <- sapply(true.transm.MCAR.cov.85, mean)
  
  true.transm.MCAR.cov.90 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 90,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.90.val <- sapply(true.transm.MCAR.cov.90, mean)
  
  true.transm.MCAR.cov.95 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 95,
                                                               age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.95.val <- sapply(true.transm.MCAR.cov.95, mean)
  
  
  true.transm.MCAR.cov.100 <- tryCatch(CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 100,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 30)))
  true.transm.MCAR.cov.100.val <- sapply(true.transm.MCAR.cov.100, mean)
  
  
  # VI:  True age mixing in transmission network with MAR
  
  # VI.a
  
  true.transm.AR.a.cov.35 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.35.val <- sapply(true.transm.AR.a.cov.35, mean)
  
  
  
  true.transm.AR.a.cov.40 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.40.val <- sapply(true.transm.AR.a.cov.40, mean)
  
  
  
  true.transm.AR.a.cov.45 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.45.val <- sapply(true.transm.AR.a.cov.45, mean)
  
  true.transm.AR.a.cov.50 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.50.val <- sapply(true.transm.AR.a.cov.50, mean)
  
  
  true.transm.AR.a.cov.55 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.55.val <- sapply(true.transm.AR.a.cov.55, mean)
  
  
  
  true.transm.AR.a.cov.60 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.60.val <- sapply(true.transm.AR.a.cov.60, mean)
  
  
  
  true.transm.AR.a.cov.65 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.65.val <- sapply(true.transm.AR.a.cov.65, mean)
  
  
  true.transm.AR.a.cov.70 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.70.val <- sapply(true.transm.AR.a.cov.70, mean)
  
  true.transm.AR.a.cov.75 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.75.val <- sapply(true.transm.AR.a.cov.75, mean)
  
  true.transm.AR.a.cov.80 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.80.val <- sapply(true.transm.AR.a.cov.80, mean)
  
  true.transm.AR.a.cov.85 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.85.val <- sapply(true.transm.AR.a.cov.85, mean)
  
  true.transm.AR.a.cov.90 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,    datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.90.val <- sapply(true.transm.AR.a.cov.90, mean)
  
  true.transm.AR.a.cov.95 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.95.val <- sapply(true.transm.AR.a.cov.95, mean)
  
  
  true.transm.AR.a.cov.100 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 100,
                                                               seq.gender.ratio = 0.7, age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 30)))
  true.transm.AR.a.cov.100.val <- sapply(true.transm.AR.a.cov.100, mean)
  
  
  # VI.b
  
  
  true.transm.AR.b.cov.35 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.35.val <- sapply(true.transm.AR.b.cov.35, mean)
  
  
  
  true.transm.AR.b.cov.40 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.40.val <- sapply(true.transm.AR.b.cov.40, mean)
  
  
  
  true.transm.AR.b.cov.45 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.45.val <- sapply(true.transm.AR.b.cov.45, mean)
  
  true.transm.AR.b.cov.50 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.50.val <- sapply(true.transm.AR.b.cov.50, mean)
  
  
  true.transm.AR.b.cov.55 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.55.val <- sapply(true.transm.AR.b.cov.55, mean)
  
  
  
  true.transm.AR.b.cov.60 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.60.val <- sapply(true.transm.AR.b.cov.60, mean)
  
  
  
  true.transm.AR.b.cov.65 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.65.val <- sapply(true.transm.AR.b.cov.65, mean)
  
  
  true.transm.AR.b.cov.70 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.70.val <- sapply(true.transm.AR.b.cov.70, mean)
  
  true.transm.AR.b.cov.75 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.75.val <- sapply(true.transm.AR.b.cov.75, mean)
  
  true.transm.AR.b.cov.80 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.80.val <- sapply(true.transm.AR.b.cov.80, mean)
  
  true.transm.AR.b.cov.85 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.85.val <- sapply(true.transm.AR.b.cov.85, mean)
  
  true.transm.AR.b.cov.90 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.90.val <- sapply(true.transm.AR.b.cov.90, mean)
  
  true.transm.AR.b.cov.95 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.95.val <- sapply(true.transm.AR.b.cov.95, mean)
  
  
  true.transm.AR.b.cov.100 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 100,
                                                               seq.gender.ratio = 0.3, age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 30)))
  true.transm.AR.b.cov.100.val <- sapply(true.transm.AR.b.cov.100, mean)
  
  
  
  # VI.c.
  
  
  true.transm.AR.c.cov.35 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 35,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.35.val <- sapply(true.transm.AR.c.cov.35, mean)
  
  
  
  true.transm.AR.c.cov.40 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 40,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.40.val <- sapply(true.transm.AR.c.cov.40, mean)
  
  
  
  true.transm.AR.c.cov.45 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 45,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.45.val <- sapply(true.transm.AR.c.cov.45, mean)
  
  true.transm.AR.c.cov.50 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 50,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.50.val <- sapply(true.transm.AR.c.cov.50, mean)
  
  
  true.transm.AR.c.cov.55 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 55,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.55.val <- sapply(true.transm.AR.c.cov.55, mean)
  
  
  
  true.transm.AR.c.cov.60 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 60,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.60.val <- sapply(true.transm.AR.c.cov.60, mean)
  
  
  
  true.transm.AR.c.cov.65 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 65,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.65.val <- sapply(true.transm.AR.c.cov.65, mean)
  
  
  true.transm.AR.c.cov.70 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 70,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.70.val <- sapply(true.transm.AR.c.cov.70, mean)
  
  true.transm.AR.c.cov.75 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 75,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.75.val <- sapply(true.transm.AR.c.cov.75, mean)
  
  true.transm.AR.c.cov.80 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 80,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.80.val <- sapply(true.transm.AR.c.cov.80, mean)
  
  true.transm.AR.c.cov.85 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 85,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.85.val <- sapply(true.transm.AR.c.cov.85, mean)
  
  true.transm.AR.c.cov.90 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net, datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),  
                                                              seq.cov = 90,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.90.val <- sapply(true.transm.AR.c.cov.90, mean)
  
  true.transm.AR.c.cov.95 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 95,
                                                              seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50)),
                                      error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.95.val <- sapply(true.transm.AR.c.cov.95, mean)
  
  
  true.transm.AR.c.cov.100 <- tryCatch(AR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,   datalist = DataListALL,
                                                               limitTransmEvents = 7,
                                                               timewindow = c(30,40),
                                                               seq.cov = 100,
                                                               seq.gender.ratio = 0.5, age.group.15.25 = c(15,25),
                                                               age.group.25.40 = c(25,40),
                                                               age.group.40.50 = c(40,50)),
                                       error=function(e) return(rep(NA, 30)))
  true.transm.AR.c.cov.100.val <- sapply(true.transm.AR.c.cov.100, mean)
  
  
  
  # Names
  
  
  names.epidemic.df <- c("hiv.prev.lt25.women", "hiv.prev.lt25.men", 
                         "hiv.prev.25.40.women", "hiv.prev.25.40.men",
                         "hiv.prev.40.50.women", "hiv.prev.40.50.men", 
                         "Pop.AAD.male", "Pop.SDAD.male", "Pop.slope.male", "Pop.WSD.male", "Pop.BSD.male", "Pop.intercept.male",
                         
                         "pp.cp.6months.male",
                         
                         "METRICS.incidence.df.15.24.men", "METRICS.incidence.df.15.24.women", 
                         "METRICS.incidence.df.25.39.men", "METRICS.incidence.df.25.39.women",
                         "METRICS.incidence.df.40.49.men", "METRICS.incidence.df.40.49.women")
  
  
  
  
  name.lme <- names(het.lme.val)
  
  
  het.lme.val <- c(het.a, het.beta, het.b1, het.b2, SD.female, SD.male)
  
  names(het.lme.val) <-  c("het.av.age.male", "het.gendEffect.clust", "het.between.transm.var", "het.within.transm.var", "SD.female", "SD.male")
  
  
  
  # same names all and we may have NA for low coverage
  
  
  ## Names for stat in clusters
  
  name.clust.MCAR.35 <- paste0("clust.MCAR.35.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.40 <- paste0("clust.MCAR.40.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.45 <- paste0("clust.MCAR.45.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.50 <- paste0("clust.MCAR.50.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.55 <- paste0("clust.MCAR.55.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.60 <- paste0("clust.MCAR.60.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.65 <- paste0("clust.MCAR.65.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.70 <- paste0("clust.MCAR.70.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.75 <- paste0("clust.MCAR.75.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.80 <- paste0("clust.MCAR.80.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.85 <- paste0("clust.MCAR.85.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.90 <- paste0("clust.MCAR.90.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.MCAR.95 <- paste0("clust.MCAR.95.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  
  
  name.clust.MCAR.100 <- paste0("clust.MCAR.100.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  
  
  
  name.clust.MCAR.scenari <- c(name.clust.MCAR.35, name.clust.MCAR.40, name.clust.MCAR.45, name.clust.MCAR.50,
                               name.clust.MCAR.55, name.clust.MCAR.60, name.clust.MCAR.65, name.clust.MCAR.70,
                               name.clust.MCAR.75, name.clust.MCAR.80, name.clust.MCAR.85,
                               name.clust.MCAR.90, name.clust.MCAR.95)
  
  
  name.clust.AR.a.35 <- paste0("clust.AR.a.35.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.40 <- paste0("clust.AR.a.40.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.45 <- paste0("clust.AR.a.45.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.50 <- paste0("clust.AR.a.50.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.55 <- paste0("clust.AR.a.55.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.60 <- paste0("clust.AR.a.60.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.65 <- paste0("clust.AR.a.65.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.70 <- paste0("clust.AR.a.70.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.75 <- paste0("clust.AR.a.75.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.80 <- paste0("clust.AR.a.80.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.85 <- paste0("clust.AR.a.85.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.90 <- paste0("clust.AR.a.90.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.a.95 <- paste0("clust.AR.a.95.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  
  
  
  name.clust.AR.a.scenari <- c(name.clust.AR.a.35, name.clust.AR.a.40, name.clust.AR.a.45, name.clust.AR.a.50,
                               name.clust.AR.a.55, name.clust.AR.a.60, name.clust.AR.a.65, name.clust.AR.a.70,
                               name.clust.AR.a.75, name.clust.AR.a.80, name.clust.AR.a.85,
                               name.clust.AR.a.90, name.clust.AR.a.95)
  
  
  name.clust.AR.b.35 <- paste0("clust.AR.b.35.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.40 <- paste0("clust.AR.b.40.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.45 <- paste0("clust.AR.b.45.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.50 <- paste0("clust.AR.b.50.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.55 <- paste0("clust.AR.b.55.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.60 <- paste0("clust.AR.b.60.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.65 <- paste0("clust.AR.b.65.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.70 <- paste0("clust.AR.b.70.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.75 <- paste0("clust.AR.b.75.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.80 <- paste0("clust.AR.b.80.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.85 <- paste0("clust.AR.b.85.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.90 <- paste0("clust.AR.b.90.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.b.95 <- paste0("clust.AR.b.95.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  
  
  
  
  name.clust.AR.b.scenari <- c(name.clust.AR.b.35, name.clust.AR.b.40, name.clust.AR.b.45, name.clust.AR.b.50,
                               name.clust.AR.b.55, name.clust.AR.b.60, name.clust.AR.b.65, name.clust.AR.b.70,
                               name.clust.AR.b.75, name.clust.AR.b.80, name.clust.AR.b.85,
                               name.clust.AR.b.90, name.clust.AR.b.95)
  
  
  name.clust.AR.c.35 <- paste0("clust.AR.c.35.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.40 <- paste0("clust.AR.c.40.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.45 <- paste0("clust.AR.c.45.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.50 <- paste0("clust.AR.c.50.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.55 <- paste0("clust.AR.c.55.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.60 <- paste0("clust.AR.c.60.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.65 <- paste0("clust.AR.c.65.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.70 <- paste0("clust.AR.c.70.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.75 <- paste0("clust.AR.c.75.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.80 <- paste0("clust.AR.c.80.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.85 <- paste0("clust.AR.c.85.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.90 <- paste0("clust.AR.c.90.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  name.clust.AR.c.95 <- paste0("clust.AR.c.95.",c("av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size"))
  
  
  
  name.clust.AR.c.scenari <- c(name.clust.AR.c.35, name.clust.AR.c.40, name.clust.AR.c.45, name.clust.AR.c.50,
                               name.clust.AR.c.55, name.clust.AR.c.60, name.clust.AR.c.65, name.clust.AR.c.70,
                               name.clust.AR.c.75, name.clust.AR.c.80, name.clust.AR.c.85,
                               name.clust.AR.c.90, name.clust.AR.c.95)
  
  
  ## Names for stat in true transmission networks
  
  # MCAR
  
  
  
  name.true.transm.MCAR.35 <- paste0("true.transm.MCAR.35.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  name.true.transm.MCAR.40 <- paste0("true.transm.MCAR.40.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  name.true.transm.MCAR.45 <- paste0("true.transm.MCAR.45.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.MCAR.50 <- paste0("true.transm.MCAR.50.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.55 <- paste0("true.transm.MCAR.55.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.MCAR.60 <- paste0("true.transm.MCAR.60.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.65 <- paste0("true.transm.MCAR.65.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.70 <- paste0("true.transm.MCAR.70.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.75 <- paste0("true.transm.MCAR.75.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.80 <- paste0("true.transm.MCAR.80.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.85 <- paste0("true.transm.MCAR.85.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.MCAR.90 <- paste0("true.transm.MCAR.90.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.MCAR.95 <- paste0("true.transm.MCAR.95.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  
  name.true.transm.MCAR.100 <- paste0("true.transm.MCAR.100.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                                "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                                "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                                "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                                "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                                "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                                
                                                                "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.MCAR.scenari <- c(name.true.transm.MCAR.35, name.true.transm.MCAR.40, name.true.transm.MCAR.45, name.true.transm.MCAR.50,
                                     name.true.transm.MCAR.55, name.true.transm.MCAR.60, name.true.transm.MCAR.65, name.true.transm.MCAR.70,
                                     name.true.transm.MCAR.75, name.true.transm.MCAR.80, name.true.transm.MCAR.85,
                                     name.true.transm.MCAR.90, name.true.transm.MCAR.95)
  
  # MAR .a
  name.true.transm.AR.a.35 <- paste0("true.transm.AR.a.35.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.40 <- paste0("true.transm.AR.a.40.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.45 <- paste0("true.transm.AR.a.45.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.a.50 <- paste0("true.transm.AR.a.50.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.a.55 <- paste0("true.transm.AR.a.55.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.60 <- paste0("true.transm.AR.a.60.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.65 <- paste0("true.transm.AR.a.65.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.70 <- paste0("true.transm.AR.a.70.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.a.75 <- paste0("true.transm.AR.a.75.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.80 <- paste0("true.transm.AR.a.80.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.85 <- paste0("true.transm.AR.a.85.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.90 <- paste0("true.transm.AR.a.90.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.95 <- paste0("true.transm.AR.a.95.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.a.scenari <- c(name.true.transm.AR.a.35, name.true.transm.AR.a.40, name.true.transm.AR.a.45, name.true.transm.AR.a.50,
                                     name.true.transm.AR.a.55, name.true.transm.AR.a.60, name.true.transm.AR.a.65, name.true.transm.AR.a.70,
                                     name.true.transm.AR.a.75, name.true.transm.AR.a.80, name.true.transm.AR.a.85,
                                     name.true.transm.AR.a.90, name.true.transm.AR.a.95)
  
  
  # MAR .b
  
  name.true.transm.AR.b.35 <- paste0("true.transm.AR.b.35.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.b.40 <- paste0("true.transm.AR.b.40.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.45 <- paste0("true.transm.AR.b.45.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.50 <- paste0("true.transm.AR.b.50.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.55 <- paste0("true.transm.AR.b.55.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.b.60 <- paste0("true.transm.AR.b.60.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.65 <- paste0("true.transm.AR.b.65.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.70 <- paste0("true.transm.AR.b.70.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.b.75 <- paste0("true.transm.AR.b.75.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.80 <- paste0("true.transm.AR.b.80.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.85 <- paste0("true.transm.AR.b.85.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.b.90 <- paste0("true.transm.AR.b.90.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.95 <- paste0("true.transm.AR.b.95.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.b.scenari <- c(name.true.transm.AR.b.35, name.true.transm.AR.b.40, name.true.transm.AR.b.45, name.true.transm.AR.b.50,
                                     name.true.transm.AR.b.55, name.true.transm.AR.b.60, name.true.transm.AR.b.65, name.true.transm.AR.b.70,
                                     name.true.transm.AR.b.75, name.true.transm.AR.b.80, name.true.transm.AR.b.85,
                                     name.true.transm.AR.b.90, name.true.transm.AR.b.95)
  
  # MAR .c
  
  name.true.transm.AR.c.35 <- paste0("true.transm.AR.c.35.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.40 <- paste0("true.transm.AR.c.40.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.c.45 <- paste0("true.transm.AR.c.45.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.50 <- paste0("true.transm.AR.c.50.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.55 <- paste0("true.transm.AR.c.55.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.c.60 <- paste0("true.transm.AR.c.60.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.65 <- paste0("true.transm.AR.c.65.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.70 <- paste0("true.transm.AR.c.70.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.75 <- paste0("true.transm.AR.c.75.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.c.80 <- paste0("true.transm.AR.c.80.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.85 <- paste0("true.transm.AR.c.85.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  
  name.true.transm.AR.c.90 <- paste0("true.transm.AR.c.90.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male"))
  
  
  name.true.transm.AR.c.95 <- paste0("true.transm.AR.c.95.",c("num.men.15.25", "num.women.15.25", "num.men.25.40", "num.women.25.40", "num.men.40.50", 
                                                              "num.women.40.50", "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", 
                                                              "partners.men.15.25.w.40.50", "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", 
                                                              "partners.men.25.40.w.40.50", "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", 
                                                              "partners.men.40.50.w.40.50", "mean.AD", "median.AD", "sd.AD", "RT.AAD.male", "RT.SDAD.male", 
                                                              "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male", 
                                                              
                                                              "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male" ))
  
  name.true.transm.AR.c.scenari <- c(name.true.transm.AR.c.35, name.true.transm.AR.c.40, name.true.transm.AR.c.45, name.true.transm.AR.c.50,
                                     name.true.transm.AR.c.55, name.true.transm.AR.c.60, name.true.transm.AR.c.65, name.true.transm.AR.c.70,
                                     name.true.transm.AR.c.75, name.true.transm.AR.c.80, name.true.transm.AR.c.85,
                                     name.true.transm.AR.c.90, name.true.transm.AR.c.95)
  
  
  # ALL names together
  names.scenari <- c("flag.women", "flag.men", "Pop.mean.AD", "Pop.med.AD", "Pop.sd.AD", names(het.lme.val), "flag.lme", 
                     name.clust.MCAR.scenari,
                     name.clust.AR.a.scenari, name.clust.AR.b.scenari, 
                     name.clust.AR.c.scenari,
                     
                     name.clust.MCAR.100,
                     
                     names.epidemic.df,
                     
                     name.true.transm.MCAR.scenari,
                     name.true.transm.AR.a.scenari,
                     name.true.transm.AR.b.scenari,
                     name.true.transm.AR.c.scenari,
                     
                     name.true.transm.MCAR.100)
  
  
  outputvector <- c(flag.women, flag.men, mean.AD, med.AD, sd.AD, het.lme.val, flag.lme,
                    
                    transm.clust.MCAR.cov.35, transm.clust.MCAR.cov.40, transm.clust.MCAR.cov.45,
                    transm.clust.MCAR.cov.50, transm.clust.MCAR.cov.55, transm.clust.MCAR.cov.60,
                    transm.clust.MCAR.cov.65, transm.clust.MCAR.cov.70, transm.clust.MCAR.cov.75,
                    transm.clust.MCAR.cov.80, transm.clust.MCAR.cov.85, transm.clust.MCAR.cov.90,
                    transm.clust.MCAR.cov.95,
                    
                    transm.clust.AR.a.cov.35, transm.clust.AR.a.cov.40, transm.clust.AR.a.cov.45,
                    transm.clust.AR.a.cov.50, transm.clust.AR.a.cov.55, transm.clust.AR.a.cov.60,
                    transm.clust.AR.a.cov.65, transm.clust.AR.a.cov.70, transm.clust.AR.a.cov.75,
                    transm.clust.AR.a.cov.80, transm.clust.AR.a.cov.85, transm.clust.AR.a.cov.90,
                    transm.clust.AR.a.cov.95,
                    
                    transm.clust.AR.b.cov.35, transm.clust.AR.b.cov.40, transm.clust.AR.b.cov.45,
                    transm.clust.AR.b.cov.50, transm.clust.AR.b.cov.55, transm.clust.AR.b.cov.60,
                    transm.clust.AR.b.cov.65, transm.clust.AR.b.cov.70, transm.clust.AR.b.cov.75,
                    transm.clust.AR.b.cov.80, transm.clust.AR.b.cov.85, transm.clust.AR.b.cov.90,
                    transm.clust.AR.b.cov.95,
                    
                    transm.clust.AR.c.cov.35, transm.clust.AR.c.cov.40, transm.clust.AR.c.cov.45,
                    transm.clust.AR.c.cov.50, transm.clust.AR.c.cov.55, transm.clust.AR.c.cov.60,
                    transm.clust.AR.c.cov.65, transm.clust.AR.c.cov.70, transm.clust.AR.c.cov.75,
                    transm.clust.AR.c.cov.80, transm.clust.AR.c.cov.85, transm.clust.AR.c.cov.90,
                    transm.clust.AR.c.cov.95, 
                    
                    transm.clust.MCAR.cov.100,
                    
                    
                    summary.epidemic.df,
                    
                    true.transm.MCAR.cov.35, true.transm.MCAR.cov.40, true.transm.MCAR.cov.45,
                    true.transm.MCAR.cov.50, true.transm.MCAR.cov.55, true.transm.MCAR.cov.60,
                    true.transm.MCAR.cov.65, true.transm.MCAR.cov.70, true.transm.MCAR.cov.75,
                    true.transm.MCAR.cov.80, true.transm.MCAR.cov.85, true.transm.MCAR.cov.90,
                    true.transm.MCAR.cov.95,
                    
                    true.transm.AR.a.cov.35, true.transm.AR.a.cov.40, true.transm.AR.a.cov.45,
                    true.transm.AR.a.cov.50, true.transm.AR.a.cov.55, true.transm.AR.a.cov.60,
                    true.transm.AR.a.cov.65, true.transm.AR.a.cov.70, true.transm.AR.a.cov.75,
                    true.transm.AR.a.cov.80, true.transm.AR.a.cov.85, true.transm.AR.a.cov.90,
                    true.transm.AR.a.cov.95,
                    
                    true.transm.AR.b.cov.35, true.transm.AR.b.cov.40, true.transm.AR.b.cov.45,
                    true.transm.AR.b.cov.50, true.transm.AR.b.cov.55, true.transm.AR.b.cov.60,
                    true.transm.AR.b.cov.65, true.transm.AR.b.cov.70, true.transm.AR.b.cov.75,
                    true.transm.AR.b.cov.80, true.transm.AR.b.cov.85, true.transm.AR.b.cov.90,
                    true.transm.AR.b.cov.95,
                    
                    true.transm.AR.c.cov.35, true.transm.AR.c.cov.40, true.transm.AR.c.cov.45,
                    true.transm.AR.c.cov.50, true.transm.AR.c.cov.55, true.transm.AR.c.cov.60,
                    true.transm.AR.c.cov.65, true.transm.AR.c.cov.70, true.transm.AR.c.cov.75,
                    true.transm.AR.c.cov.80, true.transm.AR.c.cov.85, true.transm.AR.c.cov.90,
                    true.transm.AR.c.cov.95,
                    
                    true.transm.MCAR.cov.100)
  
  # Population level metrics: "flag.women", "flag.men", "Pop.mean.AD", "Pop.med.AD", "Pop.sd.AD", 
  # "av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "flag.lme", 
  
  # Sequence coverage metrics:
  # "av.age.male", "gendEffect.clust", "between.clust.var", "within.clust.var"
  
  # Flags:  
  # if we have at least one transmission from older men to younger women -->> flag.women
  # if we have at least one transmission from older women to younger men -->> flag.men
  # linear mixed effects models in transmission network given by gendEffect.clust > 5 -->> flag.lme
  
  outputvector <- as.numeric(outputvector)
  
  names(outputvector) <- names.scenari
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  
  return(outputvector)
  
  
  
}

# debugonce(LMEM.master.model.age.mixing.toy1)
# 
# LMEM.master.model.age.mixing.toy1(inputvector = inputvector)
# unlink(paste0("temp"), recursive = TRUE)

# 
# inputvector <- c(101,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
#                  -0.3, -0.3,
#                  -2.7, # conception
#                  -0.52, -0.05)
# 
# test.all <- wrapper.phylo.simpact.study.1(inputvector = inputvector) # L = 437


# Without age-mixing

inputvector <- c(-0.52, -0.05, 2, 0, 2, 0.25, -0.3, -0.1,
                 # 0.2,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18

# With age-mixing
# 
# inputvector <- c(-0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
#                  # 0.2,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18


# 
# 
# cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
# cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
# cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
# cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
# cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 

# 
# 
# 
# 
# 
# #####   PARALLELIZATION   ######
# ################################
# 
# 
# # replication number
# 

reps <- 2

# 
# 
# # Input parameters in matrix form reps times (rows).

inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)

# 
# 
# 
# sim.start.time <- proc.time()[3] # ! IDs.gender.men50.women50.age.group.features
# 
features.matrix <- simpact.parallel(model = LMEM.master.model.age.mixing.toy1,
                                    actual.input.matrix = inputmatrix,
                                    seed_count = 124,
                                    n_cluster = 8)



# View(features.matrix)

# 
# sim.end.time <- proc.time()[3] - sim.start.time
# 
# print(paste0("Simulation time: ", round(sim.end.time/60,2), " minutes"))
# 
# 
# 
# # save features in the working directory
# 

write.csv(features.matrix, file = paste0(work.dir,"/LMEM.master.model.age.mixing.toy1.csv"))

# 
# unlink(paste0("temp"), recursive = TRUE)
# 
# #######   CALIBRATION


