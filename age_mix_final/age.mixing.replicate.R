# Master model for simulation of age-mixing patterns


# Define directory

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/david/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)



inputvector <- c(123, -0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) 


source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/LMEMphylo.AR.groups.fun.agemix.R")
source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/CAR.groups.fun.agemixBIS.R")
source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/AR.groups.fun.agemixBIS.R")
source("~/phylosimpact_simulation_studies_2018/age_mixing_uncertainty/LMEMphylo.CAR.groups.fun.agemix.R") 


# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop

# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC

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

library(intergraph)
library(ggtree)
library(lubridate)
library(ggplot2)
library(ggnetwork)
library(geomnet)
library(tidyr)



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





#######################
# Step 1: Run Simpact #
#######################

## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
#
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 50, 
                                 population.nummen = 1000, 
                                 population.numwomen = 1000,
                                 hivseed.time = 10, 
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, 
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 debut.debutage = 15
)

# # Assumption of nature of sexual network
# #########################################
#
cfg.list["population.msm"] = "no"


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


# # Assumptions to avoid negative branch lengths
# ###############################################
# # + sampling == start ART
# # when someone start ART, he/she is sampled and becomes non-infectious

cfg.list["monitoring.fraction.log_viralload"] <- 0


#
# ## Add-ons
#
### BEGIN Add-on
cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
cfg.list["dropout.interval.dist.type"] <- "uniform"
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

#######################
# Step 1: Run Simpact #
#######################


# results <- tryCatch(simpact.run(configParams = cfg.list,
#                                 destDir = sub.dir.rename,
#                                 agedist = age.distr,
#                                 seed = seedid,
#                                 intervention = intervention),
#                     error = simpact.errFunction)

results <- simpact.run(configParams = cfg.list,
                       destDir = sub.dir.rename,
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)


DataListALL <- readthedata(results)


datalist.agemix <- DataListALL

# datalist.agemix <- get(load("datalist.agemix.RData"))




###########################################
# Step 2: Construct transmission networks #
###########################################


source("~/phylosimpact_simulation_studies_2018/age_mix_final/advanced.transmission.network.builder.R")



simpact.trans.net.adv <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)


simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)



# simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)

###############################
# Step 3: Sequence simulation #
###############################


trans.net <- simpact.trans.net.adv # simpact.trans.net # all transmission networks


dirseqgen <- work.dir

seeds.num <- inputvector[1]

# Sequence simulation is done for at least a transmission network with 6 individuals
# This means that limitTransmEvents equal at least 7

sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                               sub.dir.rename = sub.dir.rename,
                               simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net,
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

agemixing.df <- agemixing.trans.df(trans.network = simpact.trans.net.adv, # simpact.trans.net
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

CAR.100 <- CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net
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

age.group.40.50 = c(40, 50)
timewindow = c(30, 40)
seq.cov = 100

mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net 
                           limitTransmEvents = 7,
                           timewindow = timewindow, 
                           seq.cov = seq.cov, 
                           age.limit = age.group.40.50[2])


### Transmission network table

infectionTable <- vector("list", length(simpact.trans.net.adv))

for(j in 1:length(simpact.trans.net.adv)){
  
  p <- j
  
  trans.network.i <- as.data.frame(simpact.trans.net.adv[[p]])
  
  # trans.network.i <- trans.network.i[-1,]
  
  id.lab <- paste0(p,".",trans.network.i$id,".C")
  
  trans.network.i$id.lab <- id.lab
  trans.network.i$ageSampTimeRec <- trans.network.i$SampTime - trans.network.i$TOBRec
  
  infectionTable[[p]] <- trans.network.i
  
  
}


infecttable <- rbindlist(infectionTable)


table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)



# if(length(mCAr.IDs)>5){


choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                    select.vec = mCAr.IDs,
                    name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))


mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                      sub.dir.rename = sub.dir.rename,
                                                      fasttree.tool = "FastTree",
                                                      calendar.dates = "samplingtimes.all.csv",
                                                      simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                      count.start = 1977,
                                                      endsim = 40,
                                                      clust = FALSE)


tree <- mCAr.IDs.tree.calib
sim.start.year <- 1987
first.transmission <- min(mCAr.IDs.tree.calib$sts)
mrsd <- max(mCAr.IDs.tree.calib$sts)

class(tree) <- "phylo"

dates <- format(date_decimal(c(mrsd, first.transmission)), "%Y-%m-%d")
tree$root.edge <-  - sim.start.year
phylotree.plot <- ggtree(tree, mrsd = dates[1]) +
  theme_tree2() +
  theme_grey() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) + #)scales::pretty_breaks(n = 10))
  xlab("Time") +
  ylab("")
print(phylotree.plot)



sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv"))

## 
# Node age with picante package

N <- node.age(mCAr.IDs.tree.calib)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date

mrca.v <- mrca(mCAr.IDs.tree.calib, full = FALSE)


mrca.t <- mrca(mCAr.IDs.tree.calib, full = TRUE)


tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))


library(phylobase)

m <- MRCA(mCAr.IDs.tree.calib)

library(phytools)

findMRCA(tree.cal.cov.35.IDs, tips=NULL, type=c("node","height"))



M <-  findMRCA(tree.cal.cov.35.IDs, c("11.76.C" , "11.104.C" , "11.140.C" , "11.156.C"))

phy <- mCAr.IDs.tree.calib

node.depth(phy, method = 1)
node.depth.edgelength(phy)
node.height(phy, clado.style = FALSE)

library(phangorn)

tree <- mCAr.IDs.tree.calib # rtree(10)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()

Ancestors(tree)
Children(tree)
Descendants(tree)
Siblings(tree)



Ancestors(tree, 1:3, "all")
Children(tree, 11)
Descendants(tree, 11, "tips")
Siblings(tree, 3)
# Siblings of all nodes
Siblings(tree)
mrca.phylo(tree, 1:3)
mrca.phylo(tree, match(c("t1", "t2", "t3"), tree$tip))
mrca.phylo(tree)
# same as mrca(tree), but faster for large trees


ansestor <- Ancestors(mCAr.IDs.tree.calib)

ansestor.v <- vector()

for(i in 1:length(ansestor)){
  
  k <- ansestor[[i]]
  
  ansestor.v <- c(ansestor.v, unique(k))
  
}

sort.int.ansestor <- unique(sort(ansestor.v))
sort.int.node.age <- sort(int.node.age)

tip.names <- names(mrca.v[1,])

dates.tree.df <- dplyr::filter(sampling.dates, sampling.dates$V1%in%tip.names)

# rearrange dates in tips order
tip.names.f <- vector()
dates.tree.dat <- vector()
for(i in 1:nrow(dates.tree.df)){
  for(j in 1:length(tip.names)){
    if(tip.names[i] == dates.tree.df$V1[[j]]){
      tip.names.f <- c(tip.names.f, tip.names[i])
      dates.tree.dat <- c(dates.tree.dat, 1977+40-dates.tree.df$V2[[j]])
    }
  }
}


dates.tree.named <- dates.tree.dat
names(dates.tree.named) <- tip.names.f


# make mrca matrix diagonal 0 and other elements age of mrca

sort.int.ansestor <- unique(sort(ansestor.v))
sort.int.node.age <- sort(int.node.age)


mrca.v.age <- mrca.v

for(i in 1:nrow(mrca.v.age)){
  for(j in 1:nrow(mrca.v.age)){
    
    if(i==j){
      mrca.v.age[i,j] <- 0
    }else{
      
      if(mrca.v[i,j] %in% sort.int.ansestor){
        
        p.index <- which(sort.int.ansestor == mrca.v[i,j])
        
        mrca.v.age[i,j] <-  sort.int.node.age[p.index]
      }
      
    }
    
  }
}



# make mrca matrix elements: sampling date - age of mrca

# Fist contingency matrix

mrca.v.age.samp <- mrca.v.age

mrca.v.age.samp.cont1 <- mrca.v.age.samp

for(i in 1:nrow(mrca.v.age)){
  
  for(j in 1:nrow(mrca.v.age)){
    
    if(i!=j){
      
      i.dat <- tip.names.f[i]
      
      v.index <- which(tip.names.f == i.dat)
      
      samp.date.tip <- dates.tree.dat[v.index]
      
      mrca.v.age.samp.cont1[i,] <- samp.date.tip - mrca.v.age.samp[i,]
      
    }
    
  }
}



# Second contingency matrix

mrca.v.age.samp <- mrca.v.age

mrca.v.age.samp.cont2 <- mrca.v.age.samp

for(i in 1:nrow(mrca.v.age)){
  
  for(j in 1:nrow(mrca.v.age)){
    
    if(i!=j){
      
      i.dat <- tip.names.f[i]
      
      v.index <- which(tip.names.f == i.dat)
      
      samp.date.tip <- dates.tree.dat[v.index]
      
      mrca.v.age.samp.cont2[,i] <- samp.date.tip - mrca.v.age.samp[,i]
      
    }
    
  }
}


# Diagonal zero

for(i in 1:nrow(mrca.v.age.samp.cont1)){
  for(j in 1:nrow(mrca.v.age.samp.cont1)){
    
    if(i==j){
      mrca.v.age.samp.cont1[i,j] <- 0
    }
  }
}


for(i in 1:nrow(mrca.v.age.samp.cont2)){
  for(j in 1:nrow(mrca.v.age.samp.cont2)){
    
    if(i==j){
      mrca.v.age.samp.cont2[i,j] <- 0
    }
  }
}



attributes.table.simpact.trans.net.adv <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%tip.names)

V.gender <- vector()
V.cd4 <- vector()
V.vl <- vector()
V.x <- vector()
V.y <- vector()
iD <- vector()

for(i in 1:length(tip.names)){
  for(j in 1:nrow(attributes.table.simpact.trans.net.adv)){
    if(tip.names[i] == attributes.table.simpact.trans.net.adv$id.lab[j]){
      
      V.gender <- c(V.gender, attributes.table.simpact.trans.net.adv$GenderRec[j])
      V.cd4 <- c(V.cd4, attributes.table.simpact.trans.net.adv$cd4[j])
      V.vl <- c(V.vl, attributes.table.simpact.trans.net.adv$vl[j])
      V.x <- c(V.x, attributes.table.simpact.trans.net.adv$location.x[j])
      V.y <- c(V.y, attributes.table.simpact.trans.net.adv$location.y[j])
      iD <- c(iD, tip.names[i])
      
    }
    
  }
}


Node.gender.cd4.vl.x.y <- data.frame(V.gender,V.cd4, V.vl, V.x, V.y, iD)

mrca.times.final <- as.matrix(abs(mrca.v.age.samp.cont2))

library(igraph)

net=graph.adjacency(as.matrix(mrca.times.final),mode="undirected",weighted=T,diag=FALSE)

E(net)       # The edges of the "net" object

V(net)       # The vertices of the "net" object

V(net)$gender <- Node.gender.cd4.vl.x.y$V.gender
V(net)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
V(net)$vl <- Node.gender.cd4.vl.x.y$V.vl
V(net)$loc.x <- Node.gender.cd4.vl.x.y$V.x
V(net)$loc.y <- Node.gender.cd4.vl.x.y$V.y

# each tips is connected to another one
plot(net)


cut.off <- 20


E(net)$weight

net.sp <- delete_edges(net, E(net)[weight>=cut.off]) # remove link greater to the cuttoff

E(net.sp)$weight

plot(net.sp, layout=layout_with_kk) 




# Incompatible matrix

names.attributes.ngaha <- Node.gender.cd4.vl.x.y

names.matrix.contigency <- names(mrca.times.final[1,])

gender.l <- names.attributes.ngaha$V.gender

mrca.times.filter <- mrca.times.final

for (i in 1:length(names(mrca.times.final[1,]))) {
  
  name.col.i <- names.matrix.contigency[i]
  
  index.i <- which(names(mrca.times.final[1,]) == name.col.i)
  
  gender.i <- gender.l[index.i]
  
  for(j in 1:length(names(mrca.times.final[1,]))){
    
    if(i != j){
      
      name.col.j <- names.matrix.contigency[j]
      
      index.j <- which(names(mrca.times.final[1,]) == name.col.j)
      
      gender.j <- gender.l[index.j]
      
      if(gender.i == gender.j){
        
        mrca.times.filter[i,j] <- 0
        
      }
      
    }
    
  }
  
}




net=graph.adjacency(as.matrix(mrca.times.filter),mode="undirected",weighted=T,diag=FALSE)

E(net)       # The edges of the "net" object

V(net)       # The vertices of the "net" object


cut.off <- 20


E(net)$weight

net.sp <- delete_edges(net, E(net)[weight>=cut.off]) # remove link greater to the cuttoff

E(net.sp)$weight

plot(net.sp, layout=layout_with_kk) 


# Age groups filtering

transm.matrix <- as.data.table(get.edgelist(net.sp))

# table.simpact.trans.net.adv

# reduced transmission table

ids <-  unique(c(transm.matrix$V1, transm.matrix$V2))

table.simpact.trans.net.igraph <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%ids)



age.groups.filtering.network.fun <- function(table.simpact.trans.net.igraph, 
                                             transm.matrix,
                                             age.group.15.25 = c(15,25),
                                             age.group.25.40 = c(25,40),
                                             age.group.40.50 = c(40,50)){
  
  Age.groups.table <- NULL
  
  v1.dat <- vector()
  v2.dat <- vector()
  age1.dat <- vector()
  age2.dat <- vector()
  gender1.dat <- vector()
  gender2.dat <- vector()
  
  for(i in 1:nrow(transm.matrix)){
    
    v1 <- transm.matrix$V1[i]
    v2 <- transm.matrix$V2[i]
    
    index.v1 <- which(table.simpact.trans.net.igraph$id.lab == v1)
    index.v2 <- which(table.simpact.trans.net.igraph$id.lab == v2)
    
    age1 <- table.simpact.trans.net.igraph$ageSampTimeRec[index.v1]
    age2 <- table.simpact.trans.net.igraph$ageSampTimeRec[index.v2]
    
    gender1 <- table.simpact.trans.net.igraph$GenderRec[index.v1]
    gender2 <- table.simpact.trans.net.igraph$GenderRec[index.v2]
    
    v1.dat <- c(v1.dat, v1)
    v2.dat <- c(v2.dat, v2)
    age1.dat <- c(age1.dat, age1)
    age2.dat <- c(age2.dat, age2)
    gender1.dat <- c(gender1.dat, gender1)
    gender2.dat <- c(gender2.dat, gender2)
    
  }
  
  age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
  
  
  
  age.group.15.25 = c(15,25)
  age.group.25.40 = c(25,40)
  age.group.40.50 = c(40,50)
  
  # men
  men.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==0)
  
  # women
  women.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==1)
  
  
  # men 15.25 and women
  
  men.15.25.women.15.25.1 <- vector()
  men.15.25.women.25.40.1 <- vector()
  men.15.25.women.40.50.1 <- vector()
  
  for (j in 1:nrow(men.age.table.1)) {
    
    
    if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
      
      if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
        
      }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
        
      }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  # women 15.25 and men
  
  women.15.25.men.15.25.2 <- vector()
  women.15.25.men.25.40.2 <- vector()
  women.15.25.men.40.50.2 <- vector()
  
  for (j in 1:nrow(women.age.table.1)) {
    
    
    if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
      
      if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.2$age2.dat[j])
        
      }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
        
      }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  
  # men 25.40 and women
  
  men.25.40.women.15.25.1 <- vector()
  men.25.40.women.25.40.1 <- vector()
  men.25.40.women.40.50.1 <- vector()
  
  for (j in 1:nrow(men.age.table.1)) {
    
    
    if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
      
      if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
        
      }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
        
      }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  
  
  
  # women 25.40 and men
  
  women.25.40.men.15.25.2 <- vector()
  women.25.40.men.25.40.2 <- vector()
  women.25.40.men.40.50.2 <- vector()
  
  for (j in 1:nrow(women.age.table.1)) {
    
    
    if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
      
      if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
        
      }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
        
      }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  
  
  # men 40.50 and women
  
  men.40.50.women.15.25.1 <- vector()
  men.40.50.women.25.40.1 <- vector()
  men.40.50.women.40.50.1 <- vector()
  
  for (j in 1:nrow(men.age.table.1)) {
    
    
    if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
      
      if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
        
      }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
        
      }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  
  
  
  # women 40.50 and men
  
  women.40.50.men.15.25.2 <- vector()
  women.40.50.men.25.40.2 <- vector()
  women.40.50.men.40.50.2 <- vector()
  
  for (j in 1:nrow(women.age.table.1)) {
    
    
    if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
      
      if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
        
        women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
        
      }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
        
        women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
        
      }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
        
        women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
      }
      
    }
    
    
  }
  
  
  
  
  men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
  
  men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
  
  men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
  
  men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
  
  men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
  
  men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
  
  men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
  
  men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
  
  men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
  
  Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
                               length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
                               length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
                             ncol = 3,
                             byrow = TRUE)
  
  colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
  rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
  
  Age.groups.table <- as.table(Age.groups.table)
  
  
  
  # run ClusterPicker
  
  system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
  
  # Read clusters' files
  
  dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                   all.files = FALSE,
                   full.names = FALSE, recursive = FALSE)
  
  
  clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=dd, 
                                                             simpact.trans.net = simpact.trans.net,
                                                             limitTransmEvents = 7)
  
  
  # library(phytools)
  
  # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
  
  tree.cal <- tree.cal.cov.35.IDs # read.tree(paste0(tree.topo))
  
  
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
  
  # # Estimating confidence intervals for rates and dates using a parametric bootstrap
  # pb <- parboot.treedater(tree.calib.LTT) # Lineage Through Time
  # 
  # # Lineages through time
  # LTT <- plot.parboot.ltt.dat(pb)
  # 
  # lb.mean.feature <- mean(LTT$lb) # mean of low values of LTT
  # lb.median.feature <- median(LTT$lb) # median of low values of LTT
  # 
  # ub.mean.feature <- mean(LTT$ub) # mean of upper values of LTT
  # ub.median.feature <- median(LTT$ub) # median of upper values of LTT
  # 
  # median.mean.feature <- mean(LTT$median) # mean of medians of values of LTT
  # median.median.feature <- median(LTT$median) # median of medians of values of LTT
  
  phylo.features.summary <- c(mean.feature, colless.feature, sackin.feature, mean.tipsDepths.feature, mean.nodesDepths.feature,
                              maxHeight.feature)
  # lb.mean.feature, lb.median.feature, ub.mean.feature, ub.median.feature,
  # median.mean.feature, median.median.feature)
  
  features.names <- c("meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", "mean.nodesDepths.feature",
                      "maxHeight.feature")
  # , "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
  #                     "LTT.median.mean.feature", "LTT.median.median.feature")
  
  names(phylo.features.summary) <- features.names
  
  clust.phylo.fit.params <- c(clust.fit.params, phylo.features.summary)
  
}else{
  
  clust.phylo.fit.params <- rep(NA, 14)
}

