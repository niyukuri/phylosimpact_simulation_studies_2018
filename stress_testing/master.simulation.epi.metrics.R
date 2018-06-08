

# Define directory

work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


# work.dir <- "~/Desktop/calibration/"


master.simulation.epi.metrics <- function(inputvector){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  destDir <- "/home/david/Desktop/mastermodeltest/temp" # on laptop
  
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
  
  
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                   # population.msm = "no",
                                   population.simtime = 50, #20, #40,  #25 for validation. 20 for calibration
                                   population.nummen = 300, #3000, #600, # 3800, #2500,
                                   population.numwomen = 300, # 3000, #600, #4200, #2500,
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
  
  seedid <- inputvector[1]
  
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
  
  simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  ############################ METRICS: TRANSMISSION NETWORK CHARACTERISTICS #####################
  
  # 1. Incidence trend #
  ######################
  
  # Overall incidence in 30 years
  incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25), timewindow = c(10, 40))
  
  METRICS.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
  
  METRICS.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
  METRICS.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]
  
  incidence.df.25.34 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 35), timewindow = c(10, 40))
  
  METRICS.incidence.df.25.34 <- incidence.df.25.34$incidence[3]
  
  METRICS.incidence.df.25.34.men <- incidence.df.25.34$incidence[1]
  METRICS.incidence.df.25.34.women <- incidence.df.25.34$incidence[2]
  
  incidence.df.35.44 <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(35, 45), timewindow = c(10, 40))
  
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
  
  
  # 2. Age mixing in transmissions #
  ##################################
  
  agemix.df <- agemixing.trans.df(trans.network = simpact.trans.net,
                                  limitTransmEvents = 7)
  
  agemix.fit <- lmer(AgeInfecDon ~ AgeInfecRec + (1 | DonId),
                     data = dplyr::filter(agemix.df, GenderDon =="0"),
                     REML = TRUE,
                     control=lmerControl(check.nobs.vs.nlev = "ignore",
                                         check.nobs.vs.rankZ = "ignore",
                                         check.nobs.vs.nRE="ignore"))
  
  # agemix.fit <- fit.agemix.trans(datatable = agemix.df)
  
  coef.inter <- fixef(agemix.fit)
  
  METRICS.age.mix.trans.interc <- coef.inter[[1]]
  METRICS.age.mix.slope <- coef.inter[2]
  
  # 3. Onward transmissions #
  ###########################
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = simpact.trans.net)
  
  METRICS.transm.average <- mean(transm.count)
  
  METRICS.transm.median <- median(transm.count) # add
  
  METRICS.transm.sd <- sd(transm.count) # add
  
  
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
                   
                   METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average,
                   METRICS.transm.median, METRICS.transm.sd)
  
  
  metric.names <- c("METRICS.incidence.df.15.24", "METRICS.incidence.df.25.34", "METRICS.incidence.df.35.44",
                    
                    "METRICS.incidence.df.15.24.int.40.41", "METRICS.incidence.df.25.34.int.40.41",
                    "METRICS.incidence.df.35.44.int.40.41",
                    "METRICS.incidence.df.15.24.int.41.42", "METRICS.incidence.df.25.34.int.41.42",
                    "METRICS.incidence.df.35.44.int.41.42",
                    "METRICS.incidence.df.15.24.int.42.43", "METRICS.incidence.df.25.34.int.42.43",
                    "METRICS.incidence.df.35.44.int.42.43",
                    "METRICS.incidence.df.15.24.int.43.44", "METRICS.incidence.df.25.34.int.43.44",
                    "METRICS.incidence.df.35.44.int.43.44",
                    "METRICS.incidence.df.15.24.int.44.45", "METRICS.incidence.df.25.34.int.44.45",
                    "METRICS.incidence.df.35.44.int.44.45",
                    
                    "METRICS.age.mix.trans.interc", "METRICS.age.mix.slope", "METRICS.transm.average",
                    "METRICS.transm.median", "METRICS.transm.sd")
  
  
  names(epi.Metrics) <- metric.names
  
  outputvector <- epi.Metrics
  
  # outputvector <- outputvector[1]
  
  return(outputvector)
}


# epi.metric.sim <- master.simulation.epi.metrics(inputvector = inputvector)



inputvector <- c(1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
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


