age.distr <- agedistr.creator(shape = 5, scale = 65)


cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                 population.simtime = 40, #110, #calibration was based on 40, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 2500,#1000, # On the VSC it was 2500?,#2500,#2500, #2500,
                                 population.numwomen = 2500,#1000, # On the VSC it was 2500?,#2500,#2500, #2500,
                                 hivseed.time = 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 25, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 25,
                                 hivtransmission.param.a = -1,
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
cfg.list["monitoring.fraction.log_viralload"] <- 0.3
cfg.list["dropout.interval.dist.uniform.min"] <- 100
cfg.list["dropout.interval.dist.uniform.max"] <- 200

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 0
cfg.list["diagnosis.baseline"] <- -99
# cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # We should move this to the wrapper function so that it can be calibrated
cfg.list["population.msm"] = "no"


cfg.list["hivtransmission.param.f1"] = log(1.1) #log(inputvector[2])
cfg.list["hivtransmission.param.f2"] = log(log(sqrt(1.1)) / log(1.1)) / 5 #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[2] #inputvector[3]
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[2] #inputvector[3]
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[3] #inputvector[4]
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[3] #inputvector[4]
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[4] #inputvector[5]
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[4] #inputvector[5]
cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[5] #inputvector[6]
cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[6] #inputvector[7]
cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[7] #inputvector[8]
cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[8] #inputvector[9]




cfg <- cfg.list

cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
cfg["monitoring.fraction.log_viralload"] <- 0.3
cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
#cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
#cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[9] #inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[9] #inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.baseline"] <- inputvector[10] #inputvector[11]

cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[11] #inputvector[12]
cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[12] #inputvector[13]
cfg["conception.alpha_base"] <- inputvector[13] #inputvector[14] #is conception.alpha.base (higher up)
cfg["dissolution.alpha_0"] <- inputvector[14] #inputvector[15]
cfg["dissolution.alpha_4"] <- 0 #inputvector[16]



# Here we insert the ART acceptability paramter and the ART  interventions, so that we can calibrate them


cfg["person.art.accept.threshold.dist.fixed.value"] <- inputvector[15] # Let's search between 0.5 and 1


# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 25 #25
art.intro["diagnosis.baseline"] <- inputvector[16] #0#100 # We should move this to the wrapper function so that it can be calibrated
art.intro["monitoring.cd4.threshold"] <- 100 # 1200
#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

art.intro2 <- list()
art.intro2["time"] <- 25 + 5 #25 + 5 # inputvector[5] ######### 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 25 + 8 #25 + 8 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 25 + 11 #25 + 11 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 25 + 13 #25 + 13
art.intro5["monitoring.cd4.threshold"] <- 600 # This is equivalent to immediate access

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist # scenario(interventionlist, tasp.indicator)

cfg["hivtransmission.param.a"] = inputvector[17] #-1
cfg["hivtransmission.param.b"] = inputvector[18] #-90
cfg["hivtransmission.param.c"] = inputvector[19] #0.5


results <- simpact.run(configParams = cfg,
                                destDir = destDir,
                                agedist = age.distr,
                                seed = seedid, #, Introducing ART has helped to keep the prevalence high
                                intervention = intervention)