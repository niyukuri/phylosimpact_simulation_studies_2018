
#######   CALIBRATION

library(EasyABC)

library(RSimpactCyan)
library(RSimpactHelper)


ABC_DestDir <- "/home/david/Desktop/mastermodeltest/calibration"

inputvector <- c(123,1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)


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


simpact4ABC <- function(inputvector.cal){
  cfg <- cfg.list
  cfg["formation.hazard.agegapry.baseline"] <- inputvector.cal[1] # 2.8
  cfg["formation.hazard.agegapry.numrel_man"] <- inputvector.cal[2] # -0.3
  cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector.cal[2] # -0.3
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
n_init <- 100 #40
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

