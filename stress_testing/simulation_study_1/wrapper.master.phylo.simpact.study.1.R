
# Define directory

work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)



# inputvector <- c(294052610,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
#                  -0.3, -0.3,
#                  -2.7, # conception
#                  -0.52, -0.05)


wrapper.master.phylo.simpact.study.1 <- function(inputvector = input.vector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  destDir <- "/home/david/Desktop/mastermodeltest/temp" # on laptop
  
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
                                   population.simtime = 50, #20, #40,  #25 for validation. 20 for calibration
                                   population.nummen = 400, #3000, #600, # 3800, #2500,
                                   population.numwomen = 400, # 3000, #600, #4200, #2500,
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
  
  # datalist.agemix <- get(load("datalist.agemix.RData"))
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  ############################ METRICS: TRANSMISSION NETWORK CHARACTERISTICS #####################
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/epi.metric.study.1.R")
  
  
  epid.metrics <- epi.metric.study.1(datalist.agemix = datalist.agemix)
  
  name.epid.metrics <- names(epid.metrics)
  
  names(epid.metrics) <- name.epid.metrics
    
  #################################### Features #############################
  
  # 2.1. Classic features
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  
  true.features.classic <- classic.features.study.1(datalist = datalist.agemix,
                                               work.dir = work.dir,
                                               sub.dir.rename = sub.dir.rename)
  
  name.true.features.classic <- names(true.features.classic)
  
  names(true.features.classic) <- name.true.features.classic

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
  
  
  
  #####################################################
  # Step 4: Construct time stamped phylogenetic trees #   + Sequence coverage Scenarios
  #####################################################
  
  
  dirfasttree <- work.dir
  
  
  if(file.exists(paste0(sub.dir.rename, "/C.Epidemic_seed.seq.bis.sim.nwk.fasta"))==TRUE){
    
    # check if the run has simulated sequence,
    # in other words: we have a transmission network with at least 6 individuals
    
    
    ### Featrures from all sequences ###
    ####################################
    
    ####### Features ALL sequences
    ### BEGIN
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
    
    
    
    # 2.2. Features from phylogenetic tree:
    
    # library(phytools)
    
    source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
    
    tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
    
    true.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal,
                                                               tree.calib.LTT = tree.calib.LTT,
                                                               work.dir = work.dir,
                                                               sub.dir.rename = sub.dir.rename,
                                                               simpact.trans.net = simpact.trans.net,
                                                               fasta.file = "C.Epidemic.fas",
                                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
    
    name.true.features.phylogenetic <- names(true.features.phylogenetic)
    
    names(true.features.phylogenetic) <- name.true.features.phylogenetic

    #### BEGIN Sequence Coverage Scenarios for MCAR ####
    
    cut.val <- 5 # consider at least 5 sequences
    
    ### 1st Scenario: 35 ###
    ########################
    
    cov.35.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 35,
                                 age.limit=65)
    
    
    if(length(cov.35.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.35.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.35.IDs.C.Epidemic.Fasta"))
      
 
      cov.35.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.35.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, "/cov.35.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.35.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.35.IDs,
                                                                 tree.calib.LTT = cov.35.IDs.tree.calib,
                                                                 work.dir = work.dir,
                                                                 sub.dir.rename = sub.dir.rename,
                                                                 simpact.trans.net = simpact.trans.net,
                                                                 fasta.file = "cov.35.IDs.C.Epidemic.Fasta",
                                                                 tree.file = "cov.35.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.35.IDs.features.phylogenetic <- rep(NA, length(cov.35.IDs.features.phylogenetic))
    }
    
    
    ### 2nd Scenario: 40 ###
    ########################
    
    cov.40.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 40,
                                 age.limit=65)
    
    
    if(length(cov.40.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.40.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.40.IDs.C.Epidemic.Fasta"))
      
      
      cov.40.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.40.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.40.IDs <- read.tree(paste0(sub.dir.rename, "/cov.40.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.40.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.40.IDs,
                                                                       tree.calib.LTT = cov.40.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.40.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.40.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.40.IDs.features.phylogenetic <- rep(NA, length(cov.40.IDs.features.phylogenetic))
    }
    
    
    ### 3rd Scenario: 45 ###
    ########################
    
    cov.45.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 45,
                                 age.limit=65)
    
    
    if(length(cov.45.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.45.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.45.IDs.C.Epidemic.Fasta"))
      
      
      cov.45.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.45.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.45.IDs <- read.tree(paste0(sub.dir.rename, "/cov.45.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.45.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.45.IDs,
                                                                       tree.calib.LTT = cov.45.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.45.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.45.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.45.IDs.features.phylogenetic <- rep(NA, length(cov.45.IDs.features.phylogenetic))
    }
    
    
    ### 4th Scenario: 50 ###
    ########################
    
    cov.50.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 50,
                                 age.limit=65)
    
    
    if(length(cov.50.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.50.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.50.IDs.C.Epidemic.Fasta"))
      
      
      cov.50.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.50.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.50.IDs <- read.tree(paste0(sub.dir.rename, "/cov.50.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.50.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.50.IDs,
                                                                       tree.calib.LTT = cov.50.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.50.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.50.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.50.IDs.features.phylogenetic <- rep(NA, length(cov.50.IDs.features.phylogenetic))
    }
    
    
    ### 5th Scenario: 55 ###
    ########################
    
    cov.55.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 55,
                                 age.limit=65)
    
    
    if(length(cov.55.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.55.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.55.IDs.C.Epidemic.Fasta"))
      
      
      cov.55.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.55.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.55.IDs <- read.tree(paste0(sub.dir.rename, "/cov.55.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.55.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.55.IDs,
                                                                       tree.calib.LTT = cov.55.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.55.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.55.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.55.IDs.features.phylogenetic <- rep(NA, length(cov.55.IDs.features.phylogenetic))
    }
    
    
    ### 6th Scenario: 60 ###
    ########################
    
    cov.60.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 60,
                                 age.limit=65)
    
    
    if(length(cov.60.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.60.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.60.IDs.C.Epidemic.Fasta"))
      
      
      cov.60.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.60.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.60.IDs <- read.tree(paste0(sub.dir.rename, "/cov.60.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.60.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.60.IDs,
                                                                       tree.calib.LTT = cov.60.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.60.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.60.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.60.IDs.features.phylogenetic <- rep(NA, length(cov.60.IDs.features.phylogenetic))
    }
    
    
    ### 7th Scenario: 65 ###
    ########################
    
    cov.65.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 65,
                                 age.limit=65)
    
    
    if(length(cov.65.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.65.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.65.IDs.C.Epidemic.Fasta"))
      
      
      cov.65.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.65.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.65.IDs <- read.tree(paste0(sub.dir.rename, "/cov.65.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.65.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.65.IDs,
                                                                       tree.calib.LTT = cov.65.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.65.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.65.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.65.IDs.features.phylogenetic <- rep(NA, length(cov.65.IDs.features.phylogenetic))
    }
    
    
    ### 8th Scenario: 70 ###
    ########################
    
    cov.70.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 70,
                                 age.limit=65)
    
    
    if(length(cov.70.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.70.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.70.IDs.C.Epidemic.Fasta"))
      
      
      cov.70.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.70.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.70.IDs <- read.tree(paste0(sub.dir.rename, "/cov.70.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.70.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.70.IDs,
                                                                       tree.calib.LTT = cov.70.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.70.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.70.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.70.IDs.features.phylogenetic <- rep(NA, length(cov.70.IDs.features.phylogenetic))
    }
    
    
    ### 9th Scenario: 75 ###
    ########################
    
    cov.75.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 75,
                                 age.limit=65)
    
    
    if(length(cov.75.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.75.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.75.IDs.C.Epidemic.Fasta"))
      
      
      cov.75.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.75.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.75.IDs <- read.tree(paste0(sub.dir.rename, "/cov.75.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.75.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.75.IDs,
                                                                       tree.calib.LTT = cov.75.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.75.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.75.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.75.IDs.features.phylogenetic <- rep(NA, length(cov.75.IDs.features.phylogenetic))
    }
    
    
    ### 10th Scenario: 80 ###
    ########################
    
    cov.80.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 80,
                                 age.limit=65)
    
    
    if(length(cov.80.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.80.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.80.IDs.C.Epidemic.Fasta"))
      
      
      cov.80.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.80.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.80.IDs <- read.tree(paste0(sub.dir.rename, "/cov.80.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.80.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.80.IDs,
                                                                       tree.calib.LTT = cov.80.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.80.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.80.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.80.IDs.features.phylogenetic <- rep(NA, length(cov.80.IDs.features.phylogenetic))
    }
    
    
    ### 11th Scenario: 85 ###
    ########################
    
    cov.85.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 85,
                                 age.limit=65)
    
    
    if(length(cov.85.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.85.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.85.IDs.C.Epidemic.Fasta"))
      
      
      cov.85.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.85.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.85.IDs <- read.tree(paste0(sub.dir.rename, "/cov.85.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.85.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.85.IDs,
                                                                       tree.calib.LTT = cov.85.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.85.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.85.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.85.IDs.features.phylogenetic <- rep(NA, length(cov.85.IDs.features.phylogenetic))
    }
    
    
    ### 12th Scenario: 90 ###
    ########################
    
    cov.90.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 90,
                                 age.limit=65)
    
    
    if(length(cov.90.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.90.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.90.IDs.C.Epidemic.Fasta"))
      
      
      cov.90.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.90.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.90.IDs <- read.tree(paste0(sub.dir.rename, "/cov.90.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.90.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.90.IDs,
                                                                       tree.calib.LTT = cov.90.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.90.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.90.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.90.IDs.features.phylogenetic <- rep(NA, length(cov.90.IDs.features.phylogenetic))
    }
    
    
    ### 13th Scenario: 95 ###
    ########################
    
    cov.95.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                                 limitTransmEvents = 7,
                                 timewindow = c(10,40),
                                 seq.cov = 95,
                                 age.limit=65)
    
    
    if(length(cov.95.IDs)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                          select.vec = cov.95.IDs, 
                          name.file = paste0(sub.dir.rename, "/cov.95.IDs.C.Epidemic.Fasta"))
      
      
      cov.95.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTree",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = "cov.95.IDs.C.Epidemic.Fasta",
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = FALSE)
      
      tree.cal.cov.95.IDs <- read.tree(paste0(sub.dir.rename, "/cov.95.IDs.C.Epidemic.Fasta.nwk"))
      
      cov.95.IDs.features.phylogenetic <- phylogenetic.features.study1(tree.topo=tree.cal.cov.95.IDs,
                                                                       tree.calib.LTT = cov.95.IDs.tree.calib,
                                                                       work.dir = work.dir,
                                                                       sub.dir.rename = sub.dir.rename,
                                                                       simpact.trans.net = simpact.trans.net,
                                                                       fasta.file = "cov.95.IDs.C.Epidemic.Fasta",
                                                                       tree.file = "cov.95.IDs.C.Epidemic.Fasta.nwk")
      
    }else{
      
      cov.95.IDs.features.phylogenetic <- rep(NA, length(cov.95.IDs.features.phylogenetic))
    }
    
    
    sim.features.phylogenetic <- c(cov.35.IDs.features.phylogenetic,
                                   cov.40.IDs.features.phylogenetic,
                                   cov.45.IDs.features.phylogenetic,
                                   cov.50.IDs.features.phylogenetic,
                                   cov.55.IDs.features.phylogenetic,
                                   cov.60.IDs.features.phylogenetic,
                                   cov.65.IDs.features.phylogenetic,
                                   cov.70.IDs.features.phylogenetic,
                                   cov.75.IDs.features.phylogenetic,
                                   cov.80.IDs.features.phylogenetic,
                                   cov.85.IDs.features.phylogenetic,
                                   cov.95.IDs.features.phylogenetic,
                                   cov.95.IDs.features.phylogenetic)
    
    name.cov.35.IDs.features.phylogenetic <- paste0("cov.35.",names(cov.35.IDs.features.phylogenetic))
    name.cov.40.IDs.features.phylogenetic <- paste0("cov.40.",names(cov.40.IDs.features.phylogenetic))
    name.cov.45.IDs.features.phylogenetic <- paste0("cov.45.",names(cov.45.IDs.features.phylogenetic))
    name.cov.50.IDs.features.phylogenetic <- paste0("cov.50.",names(cov.50.IDs.features.phylogenetic))
    name.cov.55.IDs.features.phylogenetic <- paste0("cov.55.",names(cov.55.IDs.features.phylogenetic))
    name.cov.60.IDs.features.phylogenetic <- paste0("cov.60.",names(cov.60.IDs.features.phylogenetic))
    name.cov.65.IDs.features.phylogenetic <- paste0("cov.65.",names(cov.65.IDs.features.phylogenetic))
    name.cov.70.IDs.features.phylogenetic <- paste0("cov.70.",names(cov.70.IDs.features.phylogenetic))
    name.cov.75.IDs.features.phylogenetic <- paste0("cov.75.",names(cov.75.IDs.features.phylogenetic))
    name.cov.80.IDs.features.phylogenetic <- paste0("cov.80.",names(cov.80.IDs.features.phylogenetic))
    name.cov.85.IDs.features.phylogenetic <- paste0("cov.85.",names(cov.85.IDs.features.phylogenetic))
    name.cov.90.IDs.features.phylogenetic <- paste0("cov.90.",names(cov.90.IDs.features.phylogenetic))
    name.cov.95.IDs.features.phylogenetic <- paste0("cov.95.",names(cov.95.IDs.features.phylogenetic))
    
    name.sim.features.phylogenetic <- c(name.cov.35.IDs.features.phylogenetic,
                                        name.cov.40.IDs.features.phylogenetic,
                                        name.cov.45.IDs.features.phylogenetic,
                                        name.cov.50.IDs.features.phylogenetic,
                                        name.cov.55.IDs.features.phylogenetic,
                                        name.cov.60.IDs.features.phylogenetic,
                                        name.cov.65.IDs.features.phylogenetic,
                                        name.cov.70.IDs.features.phylogenetic,
                                        name.cov.75.IDs.features.phylogenetic,
                                        name.cov.80.IDs.features.phylogenetic,
                                        name.cov.85.IDs.features.phylogenetic,
                                        name.cov.90.IDs.features.phylogenetic,
                                        name.cov.95.IDs.features.phylogenetic)
    
    names(sim.features.phylogenetic) <- name.sim.features.phylogenetic
    
  }else{
    
    
    # No phylogenetic features
    sim.features.phylogenetic <- rep(NA, 13*30)
  }
  #     c(name.epid.metrics, name.true.features.classic, name.true.features.phylogenetic)
  
  epi.metrics.features <- c(epid.metrics, true.features.classic, sim.features.phylogenetic)
  
  return(epi.metrics.features)
  
  unlink(paste0(sub.dir.rename), recursive = TRUE) # sub.dir.rename
  
  
}

# unlink(paste0("temp"), recursive = TRUE)

# 
# inputvector <- c(101,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
#                  -0.3, -0.3,
#                  -2.7, # conception
#                  -0.52, -0.05)
# 
# test.all <- wrapper.phylo.simpact.study.1(inputvector = inputvector) # L = 437

inputvector <- c(1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)
# 
# 
# 
# #####   PARALLELIZATION   ######
# ################################
# 
# 
# # replication number
# 

reps <- 4

# 
# 
# # Input parameters in matrix form reps times (rows).

inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)

# 
# 
# 
# sim.start.time <- proc.time()[3] # ! IDs.gender.men50.women50.age.group.features
# 
features.matrix <- phylo.simpact.parallel(model = wrapper.master.phylo.simpact.study.1,
                                          actual.input.matrix = inputmatrix,
                                          seed_count = 124,
                                          n_cluster = 4)
# 
# sim.end.time <- proc.time()[3] - sim.start.time
# 
# print(paste0("Simulation time: ", round(sim.end.time/60,2), " minutes"))
# 
# 
# 
# # save features in the working directory
# 

write.csv(features.matrix, file = paste0(work.dir,"/features.matrix.csv"))

# 
# unlink(paste0("temp"), recursive = TRUE)
# 
# #######   CALIBRATION
