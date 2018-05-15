
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


age.distr <- agedistr.creator(shape = 5, scale = 65)
#
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 40, 
                                 population.nummen = 100, 
                                 population.numwomen = 100,
                                 hivseed.time = 10, 
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, 
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 debut.debutage = 15
                                 )
#
# #

inputvector <- c(123, # seed.id
                 -052, # dissolution.alpha_0
                 -0.52, # dissolution.alpha_4
                 0, # person.agegap.man.dist.normal.mu
                 0, # person.agegap.woman.dist.normal.mu
                 3, # person.agegap.man.dist.normal.sigma
                 3, # person.agegap.woman.dist.normal.sigma
                 0.25, # formation.hazard.agegapry.gap_agescale_man
                 0.25, # formation.hazard.agegapry.gap_agescale_woman
                 -0.3, # formation.hazard.agegapry.numrel_man
                 -0.3, # formation.hazard.agegapry.numrel_woman
                 -0.1, # formation.hazard.agegapry.numrel_diff
                 0.2, # population.eyecap.fraction
                 -1, # hivtransmission.param.a
                 -90, # hivtransmission.param.b
                 0.5, # hivtransmission.param.c
                 0.04879016, # hivtransmission.param.f1
                 -0.1386294, # hivtransmission.param.f2
                 5, # person.vsp.toacute.x
                 7, # person.vsp.toaids.x
                 12, # person.vsp.tofinalaids.x
                 -2.7 # conception.alpha_base
                 )
#
#
# # Seed for reproducability
# ##########################
#
seedid <- inputvector[1]
#
#
# # Sexual behaviour
# ###################
#
cfg.list["dissolution.alpha_0"] <- inputvector[2] # -0.52
cfg.list["dissolution.alpha_4"] <- inputvector[3] # -0.05
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4] # 0
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # 0
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # 3
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[7] # 3
cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[8] # 0.25
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[9] # 0.25
cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[10] # -0.3
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[11] # -0.3
cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[12] # -0.1
cfg.list["population.eyecap.fraction"] <- inputvector[13] # 0.2
#
# # HIV transmission
# ###################
#


cfg.list["hivtransmission.param.a"] <- inputvector[14] # -1
cfg.list["hivtransmission.param.b"] <- inputvector[15] # -90
cfg.list["hivtransmission.param.c"] <- inputvector[16] # 0.5
cfg.list["hivtransmission.param.f1"] <- inputvector[17] # 0.04879016
cfg.list["hivtransmission.param.f2"] <- inputvector[18] # -0.1386294
cfg.list["person.vsp.toacute.x"] <- inputvector[19] # 5 
cfg.list["person.vsp.toaids.x"] <- inputvector[20] # 7
cfg.list["person.vsp.tofinalaids.x"] <- inputvector[21] # 12
#
# # Demographic
# ##############
#

cfg.list["conception.alpha_base"] <- inputvector[22] # -2.7

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
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.4
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

#
# Let's introduce ART,
art.intro <- list()
art.intro["time"] <- 25 #25
art.intro["diagnosis.baseline"] <- 100
art.intro["monitoring.cd4.threshold"] <- 100

# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

art.intro2 <- list()
art.intro2["time"] <- 25 + 5 #25 + 5 = 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 25 + 8 #25 + 8 = 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 25 + 11 #25 + 11 = 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 25 + 13 #25 + 13 = 38
art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access

interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist

# Events
cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3




