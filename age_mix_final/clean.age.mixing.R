# Master model for simulation of age-mixing patterns

# Make sure you have seq-gen, FastTree, and comandline ClusterPicker_1.2.3 in you working directory

# Define directory

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster

source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")

source("~/phylosimpact_simulation_studies_2018/age_mix_final/advanced.transmission.network.builder.R")

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




inputvector <- c(123, -0.52, -0.05, 5, 7, 3, 0.25, -0.3, -0.1, 
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -2.7) 



## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
#
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 50, 
                                 population.nummen = 7000, 
                                 population.numwomen = 7000,
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




# Running Simpact 
#################

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


simpact.trans.net.adv <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)


# simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)



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



###################################################################
# Step 4: Epidemic statistics and sexual behaviour: full data set #
###################################################################



# (i) Age mixing in relationships
##################################

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
  
  names(mix.rels.dat) <-  c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
  
}else{
  
  mix.rels.dat <- rep(NA, 6)
  
  names(mix.rels.dat) <-  c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
  
}

# age.scatter.df <- agemix.model[[1]]

#  (ii) Point 	prevalence of concurrency in the adult population
##################################################################


# Concurrency point prevalence 6 months before a survey, among men


pp.cp.6months.male.rels <- tryCatch(concurr.pointprev.calculator(datalist = datalist.agemix,
                                                                 timepoint = 40 - 0.5), error=function(e) return(NA))


# (iii) Prevalence
##################

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
#################

incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(15, 25), timewindow = c(30, 40))

epi.rels.incidence.df.15.24 <- incidence.df.15.24$incidence[3]

epi.rels.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
epi.rels.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]


incidence.df.25.39 <- incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(25, 40), timewindow = c(30, 40))

epi.rels.incidence.df.25.39 <- incidence.df.25.39$incidence[3]

epi.rels.incidence.df.25.39.men <- incidence.df.25.39$incidence[1]
epi.rels.incidence.df.25.39.women <- incidence.df.25.39$incidence[2]


incidence.df.40.49 <- incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(25, 40), timewindow = c(30, 40))

epi.rels.incidence.df.40.49 <- incidence.df.40.49$incidence[3]

epi.rels.incidence.df.40.49.men <- incidence.df.40.49$incidence[1] # res
epi.rels.incidence.df.40.49.women <- incidence.df.40.49$incidence[2] # res



summary.epidemic.rels.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                              hiv.prev.25.40.women, hiv.prev.25.40.men,
                              hiv.prev.40.50.women, hiv.prev.40.50.men, 
                              mix.rels.dat,
                              pp.cp.6months.male.rels,
                              
                              epi.rels.incidence.df.15.24.men, epi.rels.incidence.df.15.24.women, 
                              epi.rels.incidence.df.25.39.men, epi.rels.incidence.df.25.39.women,
                              epi.rels.incidence.df.40.49.men, epi.rels.incidence.df.40.49.women)



#################################################################################
# Step 5: Epidemic statistics and sexual behaviour: data set of infected people #
#################################################################################


# Data list of infected individuals


age.group.40.50 = c(40, 50)
timewindow = c(30, 40)
seq.cov = 100


# Select IDs in MCAR scenario

mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net 
                           limitTransmEvents = 7,
                           timewindow = timewindow, 
                           seq.cov = seq.cov, 
                           age.limit = age.group.40.50[2])


# Transmission network table as from transmission networks for further steps
############################################################################


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


Study.DataTable <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs) 


IDs.study <- Study.DataTable$RecId


transm.datalist.agemix <- datalist.agemix # assign full data set new age mix data set

# Transmission table of selected individuals
table.simpact.trans.net.cov <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs)

# Person table of selected individuals
transm.datalist.agemix$ptable <- dplyr::filter(transm.datalist.agemix$ptable, transm.datalist.agemix$ptable$ID%in%IDs.study)



# (i) Age mixing in relationships

# 

agemix.rels.transm.df <- agemix.df.maker(transm.datalist.agemix)

# 
agemix.model <- pattern.modeller(dataframe = agemix.rels.transm.df,
                                 agegroup = c(15, 50),
                                 timepoint = 40, # transm.datalist.agemix$itable$population.simtime[1],
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
  
  mix.rels.transm.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
  
  names(mix.rels.transm.dat) <-  c("T.AAD.male", "T.SDAD.male", "T.slope.male", "T.WSD.male", "T.BSD.male", "T.intercept.male")
  
}else{
  
  mix.rels.transm.dat <- rep(NA, 6)
  
  names(mix.rels.transm.dat) <-  c("T.AAD.male", "T.SDAD.male", "T.slope.male", "T.WSD.male", "T.BSD.male", "T.intercept.male")
  
}

# age.scatter.df <- agemix.model[[1]]

#  (ii) Point 	prevalence of concurrency in the adult population:

# Concurrency point prevalence 6 months before a survey, among men


pp.cp.6months.male.transm <- tryCatch(concurr.pointprev.calculator(datalist = transm.datalist.agemix,
                                                                   timepoint = 40 - 0.5), error=function(e) return(NA))


# (iii) Prevalence

hiv.prev.lt25.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = 40)$pointprevalence[2]
hiv.prev.lt25.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                           agegroup = c(15, 25),
                                           timepoint = 40)$pointprevalence[1]

hiv.prev.25.40.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                              agegroup = c(25, 40),
                                              timepoint = 40)$pointprevalence[2]
hiv.prev.25.40.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                            agegroup = c(25, 40),
                                            timepoint = 40)$pointprevalence[1]

hiv.prev.40.50.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                              agegroup = c(40, 50),
                                              timepoint = 40)$pointprevalence[2]
hiv.prev.40.50.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                            agegroup = c(40, 50),
                                            timepoint = 40)$pointprevalence[1]


# (iv) Incidence

incidence.df.15.24 <- incidence.calculator(datalist = transm.datalist.agemix,
                                           agegroup = c(15, 25), timewindow = c(30, 40))

epi.transm.incidence.df.15.24 <- incidence.df.15.24$incidence[3]

epi.transm.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
epi.transm.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]


incidence.df.25.39 <- incidence.calculator(datalist = transm.datalist.agemix,
                                           agegroup = c(25, 40), timewindow = c(30, 40))

epi.transm.incidence.df.25.39 <- incidence.df.25.39$incidence[3]

epi.transm.incidence.df.25.39.men <- incidence.df.25.39$incidence[1]
epi.transm.incidence.df.25.39.women <- incidence.df.25.39$incidence[2]


incidence.df.40.49 <- incidence.calculator(datalist = transm.datalist.agemix,
                                           agegroup = c(25, 40), timewindow = c(30, 40))

epi.transm.incidence.df.40.49 <- incidence.df.40.49$incidence[3]

epi.transm.incidence.df.40.49.men <- incidence.df.40.49$incidence[1] # res
epi.transm.incidence.df.40.49.women <- incidence.df.40.49$incidence[2] # res



summary.epidemic.transm.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                                hiv.prev.25.40.women, hiv.prev.25.40.men,
                                hiv.prev.40.50.women, hiv.prev.40.50.men, 
                                mix.rels.transm.dat,
                                pp.cp.6months.male.transm,
                                
                                epi.transm.incidence.df.15.24.men, epi.transm.incidence.df.15.24.women, 
                                epi.transm.incidence.df.25.39.men, epi.transm.incidence.df.25.39.women,
                                epi.transm.incidence.df.40.49.men, epi.transm.incidence.df.40.49.women)



######################################
# Step 5: Building phylogenetic tree #
######################################



dirfasttree <- work.dir



# if(length(mCAr.IDs)>5){


# Select sequences from the pool of alignment
##############################################


choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                    select.vec = mCAr.IDs,
                    name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))


# Build and calibrate the phylogenetic tree
############################################

mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                      sub.dir.rename = sub.dir.rename,
                                                      fasttree.tool = "FastTree",
                                                      calendar.dates = "samplingtimes.all.csv",
                                                      simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                      count.start = 1977,
                                                      endsim = 40,
                                                      clust = FALSE)



N <- node.age(mCAr.IDs.tree.calib)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


mrca.v <- mrca(mCAr.IDs.tree.calib, full = FALSE) # MRCA ids


sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv")) # sampling times

# 
# tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))
# 


# Compute transmission clusters
###############################

# run ClusterPicker

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))

# Read clusters' files

dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                 all.files = FALSE,
                 full.names = FALSE, recursive = FALSE)

# Transmission clusters.

d <- clust.names <- dd

data.list.simpact.trans.net.adv <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster


# Transmission table of individuals in the transmission clusters
#################################################################


# Binding all data tables of clusters as these information are captured in transmission networks

clust.size <- vector() # size of each cluster # table.simpact.trans.net.adv

transm.df <- table.simpact.trans.net.adv

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size <- c(clust.size, nrow(clus.read))
  
  data.table.simpact.trans.net.i <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
  
  data.table.simpact.trans.net.i$clust.ID <- rep(i, nrow(data.table.simpact.trans.net.i))
  
  data.list.simpact.trans.net.adv[[i]] <- as.data.frame(data.table.simpact.trans.net.i)
  
}


data.table.simpact.trans.clusts.net.adv <- as.data.frame(do.call(rbind, data.list.simpact.trans.net.adv)) # data.table & data.frame

# data.table.simpact.trans.clusts.net.adv <- data.table.simpact.trans.net.adv


## Aligning internal nodes IDs and their age: !they must get same length

ansestor <- Ancestors(mCAr.IDs.tree.calib) # ancestors of each tips and internal node
# All ancestors output are internal nodes

ansestor.v <- vector()

for(i in 1:length(ansestor)){
  
  k <- ansestor[[i]]
  
  ansestor.v <- c(ansestor.v, unique(k))
  
}

sort.int.ansestor <- unique(sort(ansestor.v))
sort.int.node.age <- sort(int.node.age)


tip.names <- names(mrca.v[1,])

dates.tree.df <- dplyr::filter(sampling.dates, sampling.dates$V1%in%tip.names) # dates of these tips

# rearrange dates in tips order as are displayed on the tree
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


# MRCA matrix
#############

# make mrca matrix diagonal 0 and other elements (internal nodes IDs) assign them the age of mrca

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


# Diagonal zero for mrca.v.age.samp.cont1 and mrca.v.age.samp.cont2

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


# filter table.simpact.trans.net.adv and remain with table of tips names (individulas in the tree)

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


Node.gender.cd4.vl.x.y <- data.table(V.gender,V.cd4, V.vl, V.x, V.y, iD)


# Adding clusters ID on the previous attributes table from attributes.table.simpact.trans.net.adv

clust.ID <- vector()

for(i in 1:nrow(Node.gender.cd4.vl.x.y)){ # attributes table ofr all tips on the tree: Node.gender.cd4.vl.x.y
  
  id <- Node.gender.cd4.vl.x.y$iD[i]
  
  if(id%in%data.table.simpact.trans.clusts.net.adv$id.lab){    # ID of tree which belongs to IDs of clusters
    # transmission table of individuls in the transmission clusters: data.table.simpact.trans.net.adv
    
    id.index <- which(data.table.simpact.trans.clusts.net.adv$id.lab == id)
    
    clust.ID <- c(clust.ID, data.table.simpact.trans.clusts.net.adv$clust.ID[id.index])
    
  }else{
    clust.ID <- c(clust.ID, 0) # tip ID which is not in any transmission cluster is assigned value 0
  }
  
}

Node.gender.cd4.vl.x.y$clust.ID <- clust.ID

Node.gender.cd4.vl.x.y.clusID <- Node.gender.cd4.vl.x.y # attributes table with clusters' IDs



## Building transmission network

# 1. consider contigency matrix 2

mrca.times.final <- as.matrix(abs(mrca.v.age.samp.cont2))


net <- graph.adjacency(as.matrix(mrca.times.final), mode="undirected",weighted=T,diag=FALSE)

# E(net)       # The edges of the "net" object
# 
# V(net)       # The vertices of the "net" object

V(net)$gender <- Node.gender.cd4.vl.x.y$V.gender
V(net)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
V(net)$vl <- Node.gender.cd4.vl.x.y$V.vl
V(net)$loc.x <- Node.gender.cd4.vl.x.y$V.x
V(net)$loc.y <- Node.gender.cd4.vl.x.y$V.y




## Filtering the network by breaking some edges due to conditions from individuals attributes:
##############################################################################################

# 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load

# Now considering 1 and 2

names.attributes.ngaha <- Node.gender.cd4.vl.x.y

names.matrix.contigency <- names(mrca.times.final[1,])

gender.l <- names.attributes.ngaha$V.gender


clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID


mrca.times.filter <- mrca.times.final


# 
# for (i in 1:length(names(mrca.times.final[1,]))) {
#   
#   name.col.i <- names.matrix.contigency[i]
#   
#   index.i <- which(names(mrca.times.final[1,]) == name.col.i)
#   
#   gender.i <- gender.l[index.i]
#   
#   cluster.i <- clusters.zose[index.i]
#   
#   
#   for(j in 1:length(names(mrca.times.final[1,]))){
#     
#     if(i != j){
#       
#       name.col.j <- names.matrix.contigency[j]
#       
#       index.j <- which(names(mrca.times.final[1,]) == name.col.j)
#       
#       gender.j <- gender.l[index.j]
#       
#       cluster.j <- clusters.zose[index.j]
#       
#       
#       if(gender.i == gender.j){ # if same gender break the link
#         
#         mrca.times.filter[i,j] <- 0
#         
#       }
#       
#       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
#         
#         mrca.times.filter[i,j] <- 0
#         
#       }
#       
#       
#     }
#     
#   }
#   
# }

# i. Gender 
############

for (i in 1:length(names(mrca.times.final[1,]))) {
  
  name.col.i <- names.matrix.contigency[i]
  
  index.i <- which(names(mrca.times.final[1,]) == name.col.i)
  
  gender.i <- gender.l[index.i]
  
  for(j in 1:length(names(mrca.times.final[1,]))){
    
    if(i != j){
      
      name.col.j <- names.matrix.contigency[j]
      
      index.j <- which(names(mrca.times.final[1,]) == name.col.j)
      
      gender.j <- gender.l[index.j]
      
      if(gender.i == gender.j){ # if same gender break the link
        
        mrca.times.filter[i,j] <- 0
        
      }
      
    }
    
  }
  
}

mrca.times.filter.gender <- mrca.times.filter


# ii. Cluster
#############

mrca.times.filter.gender.clust <- mrca.times.filter.gender

for (i in 1:length(names(mrca.times.final[1,]))) {
  
  name.col.i <- names.matrix.contigency[i]
  
  index.i <- which(names(mrca.times.final[1,]) == name.col.i)
  
  cluster.i <- clusters.zose[index.i]
  
  
  for(j in 1:length(names(mrca.times.final[1,]))){
    
    if(i != j){
      
      name.col.j <- names.matrix.contigency[j]
      
      index.j <- which(names(mrca.times.final[1,]) == name.col.j)
      
      cluster.j <- clusters.zose[index.j]
      
      
      if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
        
        mrca.times.filter.gender.clust[i,j] <- 0
        
      }
      
      
    }
    
  }
  
}


# iii. tMRCA
#############



net.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.gender.clust),mode="undirected",weighted=T,diag=FALSE)


# Consider plausible transmissions and difference between sampling time and tMRCA


cut.off <- 7 # years

# E(net.cont.1)$weight

net.cont.1 <- delete_edges(net.cont.1, E(net.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff

# E(net.cont.1)$weight

# plot(net.cont.1, layout=layout_with_kk) 



# Delete tips of the phylogenetic tree which are not part of transmission clusters:  they have clust.ID==0 >> deletes vertices 
###################################################################################


Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
Non.ids <- Non.ids.dat$iD

net.cleaned <- delete_vertices(net.cont.1, Non.ids)



# 
# # 2. consider contigency matrix 1
# 
# mrca.times.final.2 <- as.matrix(abs(mrca.v.age.samp.cont1))
# 
# 
# net.2 <- graph.adjacency(as.matrix(mrca.times.final.2), mode="undirected",weighted=T,diag=FALSE)
# 
# E(net.2)       # The edges of the "net.2" object
# 
# V(net.2)       # The vertices of the "net.2" object
# 
# V(net.2)$gender <- Node.gender.cd4.vl.x.y$V.gender
# V(net.2)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
# V(net.2)$vl <- Node.gender.cd4.vl.x.y$V.vl
# V(net.2)$loc.x <- Node.gender.cd4.vl.x.y$V.x
# V(net.2)$loc.y <- Node.gender.cd4.vl.x.y$V.y
# 
# 
# 
# 
# ## Filtering the net.2work by breaking some edges due to conditions from individuals attributes:
# 
# # 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load
# 
# 
# names.attributes.ngaha <- Node.gender.cd4.vl.x.y
# 
# names.matrix.contigency <- names(mrca.times.final.2[1,])
# 
# gender.l <- names.attributes.ngaha$V.gender
# 
# 
# clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID
# 
# 
# mrca.times.filter.2 <- mrca.times.final.2
# 
# 
# # 
# # for (i in 1:length(names(mrca.times.final.2[1,]))) {
# #   
# #   name.col.i <- names.matrix.contigency[i]
# #   
# #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
# #   
# #   gender.i <- gender.l[index.i]
# #   
# #   cluster.i <- clusters.zose[index.i]
# #   
# #   
# #   for(j in 1:length(names(mrca.times.final.2[1,]))){
# #     
# #     if(i != j){
# #       
# #       name.col.j <- names.matrix.contigency[j]
# #       
# #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
# #       
# #       gender.j <- gender.l[index.j]
# #       
# #       cluster.j <- clusters.zose[index.j]
# #       
# #       
# #       if(gender.i == gender.j){ # if same gender break the link
# #         
# #         mrca.times.filter.2[i,j] <- 0
# #         
# #       }
# #       
# #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
# #         
# #         mrca.times.filter.2[i,j] <- 0
# #         
# #       }
# #       
# #       
# #     }
# #     
# #   }
# #   
# # }
# 
# # i. Gender 
# 
# for (i in 1:length(names(mrca.times.final.2[1,]))) {
#   
#   name.col.i <- names.matrix.contigency[i]
#   
#   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
#   
#   gender.i <- gender.l[index.i]
#   
#   for(j in 1:length(names(mrca.times.final.2[1,]))){
#     
#     if(i != j){
#       
#       name.col.j <- names.matrix.contigency[j]
#       
#       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
#       
#       gender.j <- gender.l[index.j]
#       
#       if(gender.i == gender.j){ # if same gender break the link
#         
#         mrca.times.filter.2[i,j] <- 0
#         
#       }
#       
#     }
#     
#   }
#   
# }
# 
# mrca.times.filter.2.gender <- mrca.times.filter.2
# 
# 
# # ii. Cluster
# 
# mrca.times.filter.2.gender.clust <- mrca.times.filter.2.gender
# 
# for (i in 1:length(names(mrca.times.final.2[1,]))) {
#   
#   name.col.i <- names.matrix.contigency[i]
#   
#   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
#   
#   cluster.i <- clusters.zose[index.i]
#   
#   
#   for(j in 1:length(names(mrca.times.final.2[1,]))){
#     
#     if(i != j){
#       
#       name.col.j <- names.matrix.contigency[j]
#       
#       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
#       
#       cluster.j <- clusters.zose[index.j]
#       
#       
#       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
#         
#         mrca.times.filter.2.gender.clust[i,j] <- 0
#         
#       }
#       
#       
#     }
#     
#   }
#   
# }
# 
# 
# 
# net.2.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.2.gender.clust),mode="undirected",weighted=T,diag=FALSE)
# 
# 
# # Consider plausible transmissions and difference between sampling time and tMRCA
# 
# 
# cut.off <- 20
# 
# E(net.2.cont.1)$weight
# 
# net.2.cont.1 <- delete_edges(net.2.cont.1, E(net.2.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff
# 
# E(net.2.cont.1)$weight
# 
# plot(net.2.cont.1, layout=layout_with_kk) 
# 
# 
# # Delete tips which are not part of transmission clusters, they have clust.ID==0 >> deletes vertices 
# 
# Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
# Non.ids <- Non.ids.dat$iD
# 
# net.2.cleaned <- delete_vertices(net.2.cont.1, Non.ids)



# r=graph.union(net.cleaned, net.2.cleaned)




# Age structure in the transmission network built from phylogenetic tree
#########################################################################


# produce age table

net.sp <- net.cleaned


transm.matrix <- as.data.table(get.edgelist(net.sp)) # matrix of links of the ransmission network built from phylogenetic tree

# table.simpact.trans.net.adv

# reduced transmission table: table.simpact.trans.net.adv of ids in transmission clusters

ids <-  unique(c(transm.matrix$V1, transm.matrix$V2))


table.transm.clust.net.igraph <- dplyr::filter(data.table.simpact.trans.clusts.net.adv, data.table.simpact.trans.clusts.net.adv$id.lab%in%ids) 



# 1.

# True age structure in transmission clusters as observed from phylogenetic tree #
##################################################################################


age.groups.filtered.trans.clust.network.fun <- function(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                        transm.matrix = transm.matrix,
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
    
    index.v1 <- which(table.transm.clust.net.igraph$id.lab == v1)
    index.v2 <- which(table.transm.clust.net.igraph$id.lab == v2)
    
    age1 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v1]
    age2 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v2]
    
    gender1 <- table.transm.clust.net.igraph$GenderRec[index.v1]
    gender2 <- table.transm.clust.net.igraph$GenderRec[index.v2]
    
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
        
        women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
        
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
  
  return(Age.groups.table)
  
}


# 2.

# True age structure in transmission clusters as observed from transmission network #
#####################################################################################

age.groups.filtered.transmission.clust.fun <- function(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                       age.group.15.25 = c(15,25),
                                                       age.group.25.40 = c(25,40),
                                                       age.group.40.50 = c(40,50)){
  
  table.transm.clust.net.igraph$ageSampTimeDon <- table.transm.clust.net.igraph$SampTime - table.transm.clust.net.igraph$TOBDon
  
  Age.groups.table <- NULL
  
  v1.dat <- vector()
  v2.dat <- vector()
  age1.dat <- vector()
  age2.dat <- vector()
  gender1.dat <- vector()
  gender2.dat <- vector()
  
  for(i in 1:nrow(table.transm.clust.net.igraph)){
    
    v1 <- table.transm.clust.net.igraph$RecId[i]
    v2 <- table.transm.clust.net.igraph$DonId[i]
    
    index.v1 <- which(table.transm.clust.net.igraph$RecId == v1)
    
    age1 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v1]
    age2 <- table.transm.clust.net.igraph$ageSampTimeDon[index.v1]
    
    gender1 <- table.transm.clust.net.igraph$GenderRec[index.v1]
    gender2 <- table.transm.clust.net.igraph$GenderDon[index.v1]
    
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
        
        women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
        
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
  
  return(Age.groups.table)
  
}




# 3.

# True age structure in transmission transmission network for selected individuals #
#####################################################################################

age.groups.filtered.transmission.net.fun <- function(table.transmission.net.cov = table.simpact.trans.net.cov,
                                                     age.group.15.25 = c(15,25),
                                                     age.group.25.40 = c(25,40),
                                                     age.group.40.50 = c(40,50)){
  
  table.transmission.net.cov$ageSampTimeDon <- table.transmission.net.cov$SampTime - table.transmission.net.cov$TOBDon
  
  Age.groups.table <- NULL
  
  v1.dat <- vector()
  v2.dat <- vector()
  age1.dat <- vector()
  age2.dat <- vector()
  gender1.dat <- vector()
  gender2.dat <- vector()
  
  for(i in 1:nrow(table.transmission.net.cov)){
    
    v1 <- table.transmission.net.cov$RecId[i]
    v2 <- table.transmission.net.cov$DonId[i]
    
    index.v1 <- which(table.transmission.net.cov$RecId == v1)
    
    age1 <- table.transmission.net.cov$ageSampTimeRec[index.v1]
    age2 <- table.transmission.net.cov$ageSampTimeDon[index.v1]
    
    gender1 <- table.transmission.net.cov$GenderRec[index.v1]
    gender2 <- table.transmission.net.cov$GenderDon[index.v1]
    
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
        
        women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
        
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
  
  return(Age.groups.table)
  
}




age.structure.transm.clust <- age.groups.filtered.trans.clust.network.fun(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                          transm.matrix = transm.matrix,
                                                                          age.group.15.25 = c(15,25),
                                                                          age.group.25.40 = c(25,40),
                                                                          age.group.40.50 = c(40,50))

age.structure.transm.clust.true <- age.groups.filtered.transmission.clust.fun(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                              age.group.15.25 = c(15,25),
                                                                              age.group.25.40 = c(25,40),
                                                                              age.group.40.50 = c(40,50))



age.structure.transm.net.true <- age.groups.filtered.transmission.net.fun(table.transmission.net.cov = table.simpact.trans.net.cov,
                                                                          age.group.15.25 = c(15,25),
                                                                          age.group.25.40 = c(25,40),
                                                                          age.group.40.50 = c(40,50))

