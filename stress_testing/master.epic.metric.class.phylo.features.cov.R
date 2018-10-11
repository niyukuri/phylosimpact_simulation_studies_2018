

# Define directory

work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop


# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


# work.dir <- "~/Desktop/calibration/"

inputvector <- c(123, -0.52, -0.05, 2, 0, 2, 0.25, -0.3, -0.1,
                 # 0.2,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) # length(inputvector) = 18


master.epic.metric.class.phylo.features.cov <- function(inputvector){
  
  
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  # destDir <- "/home/david/Desktop/mastermodeltest/temp" # on laptop
  
  # destDir <- "/home/niyukuri/Desktop/mastermodeltest/temp" # on PC
  
  
  # library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  # library(ape)
  # library(expoTree)
  library(data.table)
  library(readr)
  # library(phangorn)
  library(lme4)
  library(nlme)
  library(minque) # with lmer
  library(dplyr)
  # library(adephylo)
  # library(treedater)
  # library(geiger)
  # library(picante)
  library(igraph)
  # library(phyloTop)
  # library(phytools)
  # library(Rsamtools) # select IDs sequences in a file
  library(robustbase) # colMedians
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 50, 
                                   population.nummen = 3000, 
                                   population.numwomen = 3000,
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
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  ###############################
  # Step 3: Sequence simulation #
  ###############################
  

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
  ### I. Compute transmission network characteristics #
  #####################################################
  
  
  # 1. Incidence trend #
  ######################
  
  # Overall incidence in 30 years
  incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25), timewindow = c(30, 40))
  
  METRICS.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
  
  METRICS.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
  METRICS.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]
  
  incidence.df.25.34 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 35), timewindow = c(30, 40))
  
  METRICS.incidence.df.25.34 <- incidence.df.25.34$incidence[3]
  
  METRICS.incidence.df.25.34.men <- incidence.df.25.34$incidence[1]
  METRICS.incidence.df.25.34.women <- incidence.df.25.34$incidence[2]
  
  incidence.df.35.44 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(35, 45), timewindow = c(30, 40))
  
  METRICS.incidence.df.35.44 <- incidence.df.35.44$incidence[3]
  
  METRICS.incidence.df.35.44.men <- incidence.df.35.44$incidence[1]
  METRICS.incidence.df.35.44.women <- incidence.df.35.44$incidence[2]
  
  # Incidence 35 - 40
  
  # # 35-36
  # incidence.df.15.24.int.35.36 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(15, 25), timewindow = c(35,36))
  # 
  # METRICS.incidence.df.15.24.int.35.36 <- incidence.df.15.24.int.35.36$incidence[3]
  # 
  # incidence.df.25.34.int.35.36 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(25, 35), timewindow = c(35,36))
  # 
  # METRICS.incidence.df.25.34.int.35.36 <- incidence.df.25.34.int.35.36$incidence[3]
  # 
  # 
  # incidence.df.35.44.int.35.36 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(35, 45), timewindow = c(35,36))
  # 
  # METRICS.incidence.df.35.44.int.35.36 <- incidence.df.35.44.int.35.36$incidence[3]
  # 
  # 
  # # 36-37
  # incidence.df.15.24.int.36.37 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(15, 25), timewindow = c(36,37))
  # 
  # METRICS.incidence.df.15.24.int.36.37 <- incidence.df.15.24.int.36.37$incidence[3]
  # 
  # incidence.df.25.34.int.36.37 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(25, 35), timewindow = c(36,37))
  # 
  # METRICS.incidence.df.25.34.int.36.37 <- incidence.df.25.34.int.36.37$incidence[3]
  # 
  # 
  # incidence.df.35.44.int.36.37 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(35, 45), timewindow = c(36,37))
  # 
  # METRICS.incidence.df.35.44.int.36.37 <- incidence.df.35.44.int.36.37$incidence[3]
  # 
  # # 37-38
  # incidence.df.15.24.int.37.38 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(15, 25), timewindow = c(37,38))
  # 
  # METRICS.incidence.df.15.24.int.37.38 <- incidence.df.15.24.int.37.38$incidence[3]
  # 
  # incidence.df.25.34.int.37.38 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(25, 35), timewindow = c(37,38))
  # 
  # METRICS.incidence.df.25.34.int.37.38 <- incidence.df.25.34.int.37.38$incidence[3]
  # 
  # 
  # incidence.df.35.44.int.37.38 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(35, 45), timewindow = c(37,38))
  # 
  # METRICS.incidence.df.35.44.int.37.38 <- incidence.df.35.44.int.37.38$incidence[3]
  # 
  # # 38-39
  # incidence.df.15.24.int.38.39 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(15, 25), timewindow = c(38,39))
  # 
  # METRICS.incidence.df.15.24.int.38.39 <- incidence.df.15.24.int.38.39$incidence[3]
  # 
  # incidence.df.25.34.int.38.39 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(25, 35), timewindow = c(38,39))
  # 
  # METRICS.incidence.df.25.34.int.38.39 <- incidence.df.25.34.int.38.39$incidence[3]
  # 
  # 
  # incidence.df.35.44.int.38.39 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(35, 45), timewindow = c(38,39))
  # 
  # METRICS.incidence.df.35.44.int.38.39 <- incidence.df.35.44.int.38.39$incidence[3]
  # 
  # # 39-40
  # incidence.df.15.24.int.39.40 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(15, 25), timewindow = c(39,40))
  # 
  # METRICS.incidence.df.15.24.int.39.40 <- incidence.df.15.24.int.39.40$incidence[3]
  # 
  # incidence.df.25.34.int.39.40 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(25, 35), timewindow = c(39,40))
  # 
  # METRICS.incidence.df.25.34.int.39.40 <- incidence.df.25.34.int.39.40$incidence[3]
  # 
  # 
  # incidence.df.35.44.int.39.40 <- incidence.calculator(datalist = datalist.agemix,
  #                                                      agegroup = c(35, 45), timewindow = c(39,40))
  # 
  # METRICS.incidence.df.35.44.int.39.40 <- incidence.df.35.44.int.39.40$incidence[3]
  # 
  
  # In Next 5 years
  
  # 40-41
  incidence.df.15.24.int.40.41 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 25), timewindow = c(40, 41))
  
  METRICS.incidence.df.15.24.int.40.41 <- incidence.df.15.24.int.40.41$incidence[3]
  
  incidence.df.25.34.int.40.41 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(25, 35), timewindow = c(40, 41))
  
  METRICS.incidence.df.25.34.int.40.41 <- incidence.df.25.34.int.40.41$incidence[3]
  
  
  incidence.df.35.44.int.40.41 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(35, 45), timewindow = c(40, 41))
  
  METRICS.incidence.df.35.44.int.40.41 <- incidence.df.35.44.int.40.41$incidence[3]
  
  
  # 41-42
  incidence.df.15.24.int.41.42 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 25), timewindow = c(41, 42))
  
  METRICS.incidence.df.15.24.int.41.42 <- incidence.df.15.24.int.41.42$incidence[3]
  
  
  incidence.df.25.34.int.41.42 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(25, 35), timewindow = c(41, 42))
  
  METRICS.incidence.df.25.34.int.41.42 <- incidence.df.25.34.int.41.42$incidence[3]
  
  
  incidence.df.35.44.int.41.42 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(35, 45), timewindow = c(41, 42))
  
  METRICS.incidence.df.35.44.int.41.42 <- incidence.df.35.44.int.41.42$incidence[3]
  
  
  # 42-43
  incidence.df.15.24.int.42.43 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 25), timewindow = c(42,43))
  
  METRICS.incidence.df.15.24.int.42.43 <- incidence.df.15.24.int.42.43$incidence[3]
  
  
  incidence.df.25.34.int.42.43 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(25, 35), timewindow = c(42,43))
  
  METRICS.incidence.df.25.34.int.42.43 <- incidence.df.25.34.int.42.43$incidence[3]
  
  
  incidence.df.35.44.int.42.43 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(35, 45), timewindow = c(42,43))
  
  METRICS.incidence.df.35.44.int.42.43 <- incidence.df.35.44.int.42.43$incidence[3]
  
  
  # 43-44
  incidence.df.15.24.int.43.44 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 25), timewindow = c(43,44))
  
  METRICS.incidence.df.15.24.int.43.44 <- incidence.df.15.24.int.43.44$incidence[3]
  
  
  incidence.df.25.34.int.43.44 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(25, 35), timewindow = c(43,44))
  
  METRICS.incidence.df.25.34.int.43.44 <- incidence.df.25.34.int.43.44$incidence[3]
  
  
  incidence.df.35.44.int.43.44 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(35, 45), timewindow = c(43,44))
  
  METRICS.incidence.df.35.44.int.43.44 <- incidence.df.35.44.int.43.44$incidence[3]
  
  # 44-45
  incidence.df.15.24.int.44.45 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 25), timewindow = c(44,45))
  
  METRICS.incidence.df.15.24.int.44.45 <- incidence.df.15.24.int.44.45$incidence[3]
  
  
  incidence.df.25.34.int.44.45 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(25, 35), timewindow = c(44,45))
  
  METRICS.incidence.df.25.34.int.44.45 <- incidence.df.25.34.int.44.45$incidence[3]
  
  
  incidence.df.35.44.int.44.45 <- incidence.calculator(datalist = datalist.agemix,
                                                       agegroup = c(35, 45), timewindow = c(44,45))
  
  METRICS.incidence.df.35.44.int.44.45 <- incidence.df.35.44.int.44.45$incidence[3]
  
  # 
  
  
  # 2. Age mixing #
  ################
  
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
    
    names(mix.rels.dat) <- c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
  }
  
  
  METRICS.LMEM.rels.age.mix <-  mix.rels.dat 
  
  
  
  # (ii) Age mixing in transmissions
  
  
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
  
  
  METRICS.LMEM.transm.age.mix <- het.lme.val 
  
  # 3. Onward transmissions #
  ###########################
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = simpact.trans.net)
  
  METRICS.transm.average <- mean(transm.count)
  
  METRICS.transm.median <- median(transm.count) # add
  
  METRICS.transm.sd <- sd(transm.count) # add
  
  
  # epi.Metrics <- c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
  #                  
  #                  METRICS.incidence.df.15.24.int.40.41, METRICS.incidence.df.25.34.int.40.41,
  #                  METRICS.incidence.df.35.44.int.40.41,
  #                  METRICS.incidence.df.15.24.int.41.42, METRICS.incidence.df.25.34.int.41.42,
  #                  METRICS.incidence.df.35.44.int.41.42,
  #                  METRICS.incidence.df.15.24.int.42.43, METRICS.incidence.df.25.34.int.42.43,
  #                  METRICS.incidence.df.35.44.int.42.43,
  #                  METRICS.incidence.df.15.24.int.43.44, METRICS.incidence.df.25.34.int.43.44,
  #                  METRICS.incidence.df.35.44.int.43.44,
  #                  METRICS.incidence.df.15.24.int.44.45, METRICS.incidence.df.25.34.int.44.45,
  #                  METRICS.incidence.df.35.44.int.44.45,
  #                  
  #                  METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average,
  #                  METRICS.transm.median, METRICS.transm.sd)
  
  epi.Metrics <- c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
                   
                   METRICS.incidence.df.15.24.int.40.41, METRICS.incidence.df.25.34.int.40.41,
                   METRICS.incidence.df.35.44.int.40.41,
                   METRICS.incidence.df.15.24.int.41.42, METRICS.incidence.df.25.34.int.41.42,
                   METRICS.incidence.df.35.44.int.41.42,
                   METRICS.incidence.df.15.24.int.42.43, METRICS.incidence.df.25.34.int.42.43,
                   METRICS.incidence.df.35.44.int.42.43,
                   METRICS.incidence.df.15.24.int.43.44, METRICS.incidence.df.25.34.int.43.44,
                   METRICS.incidence.df.35.44.int.43.44,
                   METRICS.incidence.df.15.24.int.44.45, METRICS.incidence.df.25.34.int.44.45,
                   METRICS.incidence.df.35.44.int.44.45,
                   
                   METRICS.LMEM.rels.age.mix,
                   
                   METRICS.LMEM.transm.age.mix,
                   
                   METRICS.transm.average, METRICS.transm.median, METRICS.transm.sd)
  
  
  
  metric.names <- c("METRICS.inc.df.15.24", "METRICS.incidence.df.25.34", "METRICS.incidence.df.35.44",
                    
                    "METRICS.inc.df.15.24.int.40.41", "METRICS.inc.df.25.34.int.40.41",
                    "METRICS.inc.df.35.44.int.40.41",
                    "METRICS.inc.df.15.24.int.41.42", "METRICS.inc.df.25.34.int.41.42",
                    "METRICS.inc.df.35.44.int.41.42",
                    "METRICS.inc.df.15.24.int.42.43", "METRICS.inc.df.25.34.int.42.43",
                    "METRICS.inc.df.35.44.int.42.43",
                    "METRICS.inc.df.15.24.int.43.44", "METRICS.inc.df.25.34.int.43.44",
                    "METRICS.inc.df.35.44.int.43.44",
                    "METRICS.inc.df.15.24.int.44.45", "METRICS.inc.df.25.34.int.44.45",
                    "METRICS.inc.df.35.44.int.44.45",
                    
                    paste0("METRICS.",names(METRICS.LMEM.rels.age.mix)), 
                    paste0("METRICS.",names(METRICS.LMEM.transm.age.mix)),
                    
                    "METRICS.transm.average",
                    "METRICS.transm.median", "METRICS.transm.sd")
  
  
  names(epi.Metrics) <- metric.names
  
  outputvector <- epi.Metrics
  
  
  
  ############################################################
  ### II. Compute features considering missingness scenarios #
  ############################################################
  
  
  
  
  
  
  
  
  
  
  return(outputvector)
}


# epi.metric.sim <- master.simulation.epi.metrics(inputvector = inputvector)



inputvector <- c(1,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)

reps <- 20


# Input parameters in matrix form reps times (rows).
inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)

sim.start.time <- proc.time()[3] 

epi.Metrics.Save <- phylo.simpact.parallel(model = master.simulation.epi.metrics,
                                           actual.input.matrix = inputmatrix,
                                           seed_count = 124,
                                           n_cluster = 4)

sim.end.time <- proc.time()[3] - sim.start.time

print(paste0("Simulation time: ", round(sim.end.time/60,2), " minutes"))

write.csv(epi.Metrics.Save, file = "epi.Metrics.Save.csv")


