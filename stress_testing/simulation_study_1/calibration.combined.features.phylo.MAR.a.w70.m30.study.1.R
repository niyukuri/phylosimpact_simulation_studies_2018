#######   CALIBRATION    WITH COMBINED CLASSIC  & PHYLO FEATURES  ##########


# I.a. MAR

# Missing at random


# Sequenced individulas are chosen based on age and gender
# Specifically: 30% men and 70 % women
# age.group.15.25 = c(15,25),
# age.group.25.40 = c(25,40),
# age.group.40.50 = c(40,50)



##############
### Models ###
##############



# 35%


simpact4ABC.classic.phylo.mAr.cov.35 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.35 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.35) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.35)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.35,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.35,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.35,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.35,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.35,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.35, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.35,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.35) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}



# 40%


simpact4ABC.classic.phylo.mAr.cov.40 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.40 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.40) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.40)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.40,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.40,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.40,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.40,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.40,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.40, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.40,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.40) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 45%


simpact4ABC.classic.phylo.mAr.cov.45 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.45 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.45) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.45)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.45,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.45,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.45,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.45,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.45,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.45, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.45,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.45) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 50%


simpact4ABC.classic.phylo.mAr.cov.50 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.50 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.50) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.50)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.50,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.50,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.50,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.50,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.50,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.50, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.50,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.50) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 55 %


simpact4ABC.classic.phylo.mAr.cov.55 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.55 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.55) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.55)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.55,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.55,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.55,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.55,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.55,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.55, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.55,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.55) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 60%


simpact4ABC.classic.phylo.mAr.cov.60 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.60 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.60) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.60)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.60,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.60,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.60,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.60,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.60,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.60, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.60,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.60) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 65%


simpact4ABC.classic.phylo.mAr.cov.65 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.65 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.65) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.65)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.65,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.65,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.65,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.65,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.65,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.65, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.65,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.65) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 70%


simpact4ABC.classic.phylo.mAr.cov.70 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.70 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.70) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.70)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.70,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.70,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.70,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.70,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.70,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.70, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.70,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.70) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 75%


simpact4ABC.classic.phylo.mAr.cov.75 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.75 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.75) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.75)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.75,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.75,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.75,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.75,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.75,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.75, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.75,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.75) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 80%


simpact4ABC.classic.phylo.mAr.cov.80 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.80 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.80) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.80)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.80,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.80,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.80,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.80,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.80,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.80, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.80,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.80) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 85%


simpact4ABC.classic.phylo.mAr.cov.85 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.85 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.85) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.85)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.85,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.85,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.85,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.85,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.85,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.85, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.85,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.85) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 90%


simpact4ABC.classic.phylo.mAr.cov.90 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.90 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.90) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.90)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.90,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.90,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.90,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.90,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.90,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.90, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.90,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.90) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# 95%


simpact4ABC.classic.phylo.mAr.cov.95 <- function(inputvector){
  
  
  library(EasyABC)
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  #library(readr)
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
  library(Rsamtools)
  library(robustbase)
  
  
  # work.dir <- "/home/david/Desktop/calibration" # on laptop
  
  work.dir <- "/home/niyukuri/Desktop/calibration" # on PC
  
  
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
  
  ABC_DestDir.classic.phylo.mAr.cov.95 <- paste0(work.dir,"/temp/",generate.filename(10))
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/classic.features.study.1.R")
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/phylogenetic.features.study.1.R")
  
  
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
  
  
  results <- simpact.run(cfg, ABC_DestDir.classic.phylo.mAr.cov.95) #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  
  datalist <- readthedata(results)
  
  
  classic.stat <- classic.features.study.1(datalist = datalist,
                                           work.dir = work.dir,
                                           sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.95)
  
  
  ## Phylodynamics components
  
  
  # Construct transmission networks #
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)
  
  # Transmission data frame
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    # Simulate sequences
    
    
    trans.net <- simpact.trans.net # all transmission networks
    
    
    dirseqgen <- work.dir
    
    seeds.num <- 123
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.95,
                                   simpact.trans.net = simpact.trans.net,
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = FALSE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(ABC_DestDir.classic.phylo.mAr.cov.95,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.95,"/C.Epidemic.fas") , format = "fasta")
    
    
    dirfasttree <- work.dir
    
    tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                 sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.95,
                                                 fasttree.tool = "FastTree",
                                                 calendar.dates = "samplingtimes.all.csv",
                                                 simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                                 count.start = 1977,
                                                 endsim = 40,
                                                 clust = FALSE)
    
    tree.calib.LTT <- tree.calib
    
    write.tree(tree.calib, file = paste0(ABC_DestDir.classic.phylo.mAr.cov.95,"/calibrated.tree.nwk"))
    
    
    N <- node.age(tree.calib)
    
    int.node.age <- N$Ti # internal nodes ages
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    tree.cal <- read.tree(paste0(ABC_DestDir.classic.phylo.mAr.cov.95, "/calibrated.tree.nwk"))
    
    
    
    # Compute targets
    
    phlyo.stat <- phylogenetic.features.study1(tree.topo=tree.cal,
                                               tree.calib.LTT = tree.calib.LTT,
                                               work.dir = work.dir,
                                               sub.dir.rename = ABC_DestDir.classic.phylo.mAr.cov.95,
                                               simpact.trans.net = simpact.trans.net,
                                               fasta.file = "C.Epidemic.fas",
                                               tree.file = "C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk")
    
  }else{
    
    phlyo.stat <- rep(NA,30) # length(phylo.target.mAr.cov.95) =30
    
  }
  
  
  summary.stat <- c(classic.stat, phlyo.stat)
  
  
  # relsperpersonperyear <- nrow(datalist$rtable) / (nrow(datalist$ptable)/2) / cfg$population.simtime
  # agegapsd <- sd(datalist$rtable$AgeGap)
  # outputvector <- c(relsperpersonperyear, agegapsd)
  
  outputvector <- as.numeric(summary.stat) # names.class , Length = 47
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}



##############
### Priors ###
##############


simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 1, 3), c("unif", -0.5, 0.5),
                      c("unif", 2, 4), c("unif", 0, 1), c("unif", -1, 0), c("unif", -0.9, 0), 
                      c("unif", 0, 0.5), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), 
                      c("unif", 0, 0.5), c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9),
                      c("unif", 10, 14),  c("unif", -3.5, -1.7))


##############################################
### Targets statistics: observed summaries ###
##############################################


source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")

targets.stat <- read.csv("~/Desktop/mastermodeltest/features.matrix.csv")

targets.stat <- targets.stat[,2:length(targets.stat)]

median.targets.stat <-  colMedians(targets.stat)  # Ok # library(robustbase)

# Classic stats

classic.target <- colMedians(targets.stat[,11:26])

classic.target <- as.numeric(classic.target)


# Phylo stats


phylo.target.mAr.cov.35 <-  colMedians(targets.stat[,417:446])

phylo.target.mAr.cov.40 <-  colMedians(targets.stat[,447:476])

phylo.target.mAr.cov.45 <-  colMedians(targets.stat[,477:506])

phylo.target.mAr.cov.50 <-  colMedians(targets.stat[,507:536])

phylo.target.mAr.cov.55 <-  colMedians(targets.stat[,537:566])

phylo.target.mAr.cov.60 <-  colMedians(targets.stat[,567:596])

phylo.target.mAr.cov.65 <-  colMedians(targets.stat[,597:626])

phylo.target.mAr.cov.70 <-  colMedians(targets.stat[,627:656])

phylo.target.mAr.cov.75 <-  colMedians(targets.stat[,657:686])

phylo.target.mAr.cov.80 <-  colMedians(targets.stat[,687:716])

phylo.target.mAr.cov.85 <-  colMedians(targets.stat[,717:746])

phylo.target.mAr.cov.90 <-  colMedians(targets.stat[,747:776])

phylo.target.mAr.cov.95 <-  colMedians(targets.stat[,777:806])




# Calibrations
###############

library(EasyABC)
library(data.table)

# 35%

sum_stat_obs.mAr.cov.35 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.35))


ABC_rej.classic.phylo.mAr.cov.35 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.35,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.35,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.35 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.35$param)

write.csv(output.params.classic.phylo.mAr.cov.35, file = "output.params.classic.phylo.mAr.cov.35.csv")

output.params.classic.phylo.mAr.cov.35 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.35))

median.targets.stat.phylo.mAr.cov.35 <-  colMedians((output.params.classic.phylo.mAr.cov.35))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.35 <- as.numeric(median.targets.stat.phylo.mAr.cov.35)

reps <- 4


inputmatrix.mAr.cov.35 <- matrix(rep(inputvector.mAr.cov.35, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.35 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.35,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.35, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.35.csv"))


# 40%


sum_stat_obs.mAr.cov.40 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.40))


ABC_rej.classic.phylo.mAr.cov.40 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.40,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.40,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.40 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.40$param)

write.csv(output.params.classic.phylo.mAr.cov.40, file = "output.params.classic.phylo.mAr.cov.40.csv")

output.params.classic.phylo.mAr.cov.40 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.40))

median.targets.stat.phylo.mAr.cov.40 <-  colMedians((output.params.classic.phylo.mAr.cov.40))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.40 <- as.numeric(median.targets.stat.phylo.mAr.cov.40)

reps <- 4


inputmatrix.mAr.cov.40 <- matrix(rep(inputvector.mAr.cov.40, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.40 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.40,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.40, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.40.csv"))


# 45%


sum_stat_obs.mAr.cov.45 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.45))


ABC_rej.classic.phylo.mAr.cov.45 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.45,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.45,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.45 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.45$param)

write.csv(output.params.classic.phylo.mAr.cov.45, file = "output.params.classic.phylo.mAr.cov.45.csv")

output.params.classic.phylo.mAr.cov.45 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.45))

median.targets.stat.phylo.mAr.cov.45 <-  colMedians((output.params.classic.phylo.mAr.cov.45))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.45 <- as.numeric(median.targets.stat.phylo.mAr.cov.45)

reps <- 4


inputmatrix.mAr.cov.45 <- matrix(rep(inputvector.mAr.cov.45, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.45 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.45,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.45, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.45.csv"))


# 50%


sum_stat_obs.mAr.cov.50 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.50))


ABC_rej.classic.phylo.mAr.cov.50 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.50,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.50,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.50 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.50$param)

write.csv(output.params.classic.phylo.mAr.cov.50, file = "output.params.classic.phylo.mAr.cov.50.csv")

output.params.classic.phylo.mAr.cov.50 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.50))

median.targets.stat.phylo.mAr.cov.50 <-  colMedians((output.params.classic.phylo.mAr.cov.50))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.50 <- as.numeric(median.targets.stat.phylo.mAr.cov.50)

reps <- 4


inputmatrix.mAr.cov.50 <- matrix(rep(inputvector.mAr.cov.50, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.50 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.50,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.50, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.50.csv"))


# 55%


sum_stat_obs.mAr.cov.55 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.55))


ABC_rej.classic.phylo.mAr.cov.55 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.55,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.55,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.55 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.55$param)

write.csv(output.params.classic.phylo.mAr.cov.55, file = "output.params.classic.phylo.mAr.cov.55.csv")

output.params.classic.phylo.mAr.cov.55 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.55))

median.targets.stat.phylo.mAr.cov.55 <-  colMedians((output.params.classic.phylo.mAr.cov.55))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.55 <- as.numeric(median.targets.stat.phylo.mAr.cov.55)

reps <- 4


inputmatrix.mAr.cov.55 <- matrix(rep(inputvector.mAr.cov.55, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.55 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.55,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.55, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.55.csv"))


# 60%


sum_stat_obs.mAr.cov.60 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.60))


ABC_rej.classic.phylo.mAr.cov.60 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.60,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.60,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.60 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.60$param)

write.csv(output.params.classic.phylo.mAr.cov.60, file = "output.params.classic.phylo.mAr.cov.60.csv")

output.params.classic.phylo.mAr.cov.60 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.60))

median.targets.stat.phylo.mAr.cov.60 <-  colMedians((output.params.classic.phylo.mAr.cov.60))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.60 <- as.numeric(median.targets.stat.phylo.mAr.cov.60)

reps <- 4


inputmatrix.mAr.cov.60 <- matrix(rep(inputvector.mAr.cov.60, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.60 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.60,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.60, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.60.csv"))


# 65%


sum_stat_obs.mAr.cov.65 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.65))


ABC_rej.classic.phylo.mAr.cov.65 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.65,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.65,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.65 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.65$param)

write.csv(output.params.classic.phylo.mAr.cov.65, file = "output.params.classic.phylo.mAr.cov.65.csv")

output.params.classic.phylo.mAr.cov.65 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.65))

median.targets.stat.phylo.mAr.cov.65 <-  colMedians((output.params.classic.phylo.mAr.cov.65))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.65 <- as.numeric(median.targets.stat.phylo.mAr.cov.65)

reps <- 4


inputmatrix.mAr.cov.65 <- matrix(rep(inputvector.mAr.cov.65, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.65 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.65,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.65, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.65.csv"))


# 70%


sum_stat_obs.mAr.cov.70 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.70))


ABC_rej.classic.phylo.mAr.cov.70 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.70,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.70,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.70 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.70$param)

write.csv(output.params.classic.phylo.mAr.cov.70, file = "output.params.classic.phylo.mAr.cov.70.csv")

output.params.classic.phylo.mAr.cov.70 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.70))

median.targets.stat.phylo.mAr.cov.70 <-  colMedians((output.params.classic.phylo.mAr.cov.70))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.70 <- as.numeric(median.targets.stat.phylo.mAr.cov.70)

reps <- 4


inputmatrix.mAr.cov.70 <- matrix(rep(inputvector.mAr.cov.70, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.70 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.70,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.70, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.70.csv"))


# 75%


sum_stat_obs.mAr.cov.75 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.75))


ABC_rej.classic.phylo.mAr.cov.75 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.75,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.75,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.75 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.75$param)

write.csv(output.params.classic.phylo.mAr.cov.75, file = "output.params.classic.phylo.mAr.cov.75.csv")

output.params.classic.phylo.mAr.cov.75 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.75))

median.targets.stat.phylo.mAr.cov.75 <-  colMedians((output.params.classic.phylo.mAr.cov.75))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.75 <- as.numeric(median.targets.stat.phylo.mAr.cov.75)

reps <- 4


inputmatrix.mAr.cov.75 <- matrix(rep(inputvector.mAr.cov.75, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.75 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.75,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.75, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.75.csv"))


# 80%


sum_stat_obs.mAr.cov.80 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.80))


ABC_rej.classic.phylo.mAr.cov.80 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.80,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.80,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.80 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.80$param)

write.csv(output.params.classic.phylo.mAr.cov.80, file = "output.params.classic.phylo.mAr.cov.80.csv")

output.params.classic.phylo.mAr.cov.80 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.80))

median.targets.stat.phylo.mAr.cov.80 <-  colMedians((output.params.classic.phylo.mAr.cov.80))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.80 <- as.numeric(median.targets.stat.phylo.mAr.cov.80)

reps <- 4


inputmatrix.mAr.cov.80 <- matrix(rep(inputvector.mAr.cov.80, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.80 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.80,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.80, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.80.csv"))


# 85%


sum_stat_obs.mAr.cov.85 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.85))


ABC_rej.classic.phylo.mAr.cov.85 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.85,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.85,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.85 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.85$param)

write.csv(output.params.classic.phylo.mAr.cov.85, file = "output.params.classic.phylo.mAr.cov.85.csv")

output.params.classic.phylo.mAr.cov.85 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.85))

median.targets.stat.phylo.mAr.cov.85 <-  colMedians((output.params.classic.phylo.mAr.cov.85))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.85 <- as.numeric(median.targets.stat.phylo.mAr.cov.85)

reps <- 4


inputmatrix.mAr.cov.85 <- matrix(rep(inputvector.mAr.cov.85, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.85 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.85,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.85, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.85.csv"))


# 90%


sum_stat_obs.mAr.cov.90 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.90))


ABC_rej.classic.phylo.mAr.cov.90 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.90,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.90,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.90 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.90$param)

write.csv(output.params.classic.phylo.mAr.cov.90, file = "output.params.classic.phylo.mAr.cov.90.csv")

output.params.classic.phylo.mAr.cov.90 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.90))

median.targets.stat.phylo.mAr.cov.90 <-  colMedians((output.params.classic.phylo.mAr.cov.90))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.90 <- as.numeric(median.targets.stat.phylo.mAr.cov.90)

reps <- 4


inputmatrix.mAr.cov.90 <- matrix(rep(inputvector.mAr.cov.90, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.90 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.90,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.90, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.90.csv"))


# 95%


sum_stat_obs.mAr.cov.95 <- c(as.numeric(classic.target), as.numeric(phylo.target.mAr.cov.95))


ABC_rej.classic.phylo.mAr.cov.95 <-  ABC_rejection(model = simpact4ABC.classic.phylo.mAr.cov.95,
                                                   prior = simpact_prior,
                                                   summary_stat_target = sum_stat_obs.mAr.cov.95,
                                                   nb_simul = 12,
                                                   use_seed = TRUE,
                                                   seed_count = 1,
                                                   n_cluster = 4,
                                                   tol = 2/12)



output.params.classic.phylo.mAr.cov.95 <- as.data.table(ABC_rej.classic.phylo.mAr.cov.95$param)

write.csv(output.params.classic.phylo.mAr.cov.95, file = "output.params.classic.phylo.mAr.cov.95.csv")

output.params.classic.phylo.mAr.cov.95 <- as.data.frame(as.matrix(output.params.classic.phylo.mAr.cov.95))

median.targets.stat.phylo.mAr.cov.95 <-  colMedians((output.params.classic.phylo.mAr.cov.95))  # Ok # library(robustbase)


# Run default model with parameters values from calibration

source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/wrapper.test.study.1.R")

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)

inputvector.mAr.cov.95 <- as.numeric(median.targets.stat.phylo.mAr.cov.95)

reps <- 4


inputmatrix.mAr.cov.95 <- matrix(rep(inputvector.mAr.cov.95, reps), byrow = TRUE, nrow = reps)


epic.metric.calibrates.mAr.cov.95 <- simpact.parallel(model = wrapper.test.study.1,
                                                      actual.input.matrix = inputmatrix.mAr.cov.95,
                                                      seed_count = 124,
                                                      n_cluster = 4)


write.csv(epic.metric.calibrates.mAr.cov.95, file = paste0(work.dir,"/epic.metric.calibrates.mAr.cov.95.csv"))



