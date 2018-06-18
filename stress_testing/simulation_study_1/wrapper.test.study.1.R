


wrapper.test.study.1 <- function(inputvector){
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  destDir <- "/home/david/Desktop/mastermodeltest/temp" # on laptop
  
  # destDir <- "/home/niyukuri/Desktop/mastermodeltest/temp" # on PC
  
  # 
  # library(RSimpactCyan)
  # library(RSimpactHelper)
  # library(Rcpp)
  # library(ape)
  # library(expoTree)
  # library(data.table)
  # library(readr)
  # library(phangorn)
  # library(lme4)
  # library(nlme)
  # library(dplyr)
  # library(adephylo)
  # library(treedater)
  # library(geiger)
  # library(picante)
  # library(igraph)
  # library(phyloTop)
  # library(phytools)
  # library(Rsamtools)
  # library(robustbase)
  # 
  
  
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(expoTree)
  library(data.table)
  library(lme4)
  library(nlme)
  library(dplyr)
  library(geiger)
  library(igraph)
  library(robustbase)
  
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 50, 
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
  cfg.list["population.eyecap.fraction"] <- inputvector[10] # [9] # 0.2 c("unif", 0, 0.5)
  #
  # # HIV transmission
  # ###################
  #
  
  
  cfg.list["hivtransmission.param.a"] <- inputvector[11] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[12] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[13] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[14] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[15] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[16] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[17] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[18] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[19] # [18] # -2.7 c("unif", -3.5, -1.7)
  
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
  
  # datalist.agemix <- get(load("datalist.agemix.RData"))
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  # simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  ############################ METRICS: TRANSMISSION NETWORK CHARACTERISTICS #####################
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/simulation_study_1/epi.metric.study.1.R")
  
  
  epid.metrics <- epi.metric.study.1(datalist.agemix = datalist.agemix)
  
  name.epid.metrics <- names(epid.metrics)
  
  names(epid.metrics) <- name.epid.metrics
  
  
  return(epid.metrics)
  
  unlink(paste0(sub.dir.rename), recursive = TRUE) # sub.dir.rename
  
  
}

# 
# # d <- read.csv("inputvector.csv")
# # 
# # inputvector <- c(124,d$x)
# # f <- wrapper.test.study.1(inputvector = inputvector)
# 
# 
# df <- agemixing.trans.df(trans.network = trans.network,
#                          limitTransmEvents = 7)
# 
# 
# datatable <- df
# 
# datatable.men <- datatable
# datatable.men$AgeInfecDon <- datatable.men$AgeInfecDon-15
# 
# men.lmer <- lmer(AgeInfecDon ~ AgeInfecRec + (1 | DonId),
#                  data = dplyr::filter(datatable, GenderDon =="0"),
#                  REML = TRUE,
#                  control=lmerControl(check.nobs.vs.nlev = "ignore",
#                                      check.nobs.vs.rankZ = "ignore",
#                                      check.nobs.vs.nRE="ignore"))
# 
# coef.men <- coef(men.lmer)$DonId
# 
# coef.men.fix <- coef(summary(men.lmer))[, "Estimate"]
# 
# # data.men <- dplyr::filter(datatable, GenderDon =="0")
# 
# g.men <- ggplot(datatable.men, aes(AgeInfecDon, AgeInfecRec))
# 
# # Scatterplot
# g.men + geom_point() + 
#   geom_abline(slope = coef.men.fix[[2]], intercept = coef.men.fix[[1]], color="red") +
#   labs(subtitle="Heterosexual transmission: men to women", 
#        y="AgeInfecDon: Men", 
#        x="AgeInfecRec: Women", 
#        title="Age mixing in transmission")#, 
# #caption="Source: midwest")
# 
# 
# datatable.women <- datatable
# datatable.women$AgeInfecDon <- datatable.women$AgeInfecDon-15
# 
# women.lmer <- lmer(AgeInfecDon ~ AgeInfecRec + (1 | DonId),
#                    data = dplyr::filter(datatable, GenderDon =="1"),
#                    REML = TRUE,
#                    control=lmerControl(check.nobs.vs.nlev = "ignore",
#                                        check.nobs.vs.rankZ = "ignore",
#                                        check.nobs.vs.nRE="ignore"))
# 
# # data.women <- dplyr::filter(datatable, GenderDon =="1")
# 
# g.women <- ggplot(datatable.women, aes(AgeInfecDon, AgeInfecRec))
# 
# # Scatterplot
# g.women + geom_point() + 
#   geom_smooth(method="lm", se=F) +
#   labs(subtitle="Heterosexual transmission: women to men", 
#        y="AgeInfecDon: Women", 
#        x="AgeInfecRec: Men", 
#        title="Age mixing in transmission")#, 
# #caption="Source: midwest")
# 
# 
# # # Scatterplot
# # g.men + geom_point() + 
# #   geom_smooth(method="lm", se=F) +
# #   labs(subtitle="Heterosexual transmission: men to women", 
# #        y="AgeInfecDon: Men", 
# #        x="AgeInfecRec: Women", 
# #        title="Age mixing in transmission")#, 
# # #caption="Source: midwest")
# 
