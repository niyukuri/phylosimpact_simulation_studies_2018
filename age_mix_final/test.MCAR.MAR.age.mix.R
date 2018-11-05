# Master model for simulation of age-mixing patterns

# Make sure you have seq-gen, FastTree, and comandline ClusterPicker_1.2.3 in you working directory

# Define directory

# work.dir <- "/home/david/Desktop/mastermodeltest" # on laptop

# work.dir <- "/home/dniyukuri/lustre/agemix.25.10.2018.2" # on CHPC

# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


# setwd(paste0(work.dir))



test.MCAR.MAR.age.mix <- function(inputvector=inputvector){
  
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/advanced.transmission.network.builder.R")
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/age.mixing.MCAR.fun.R")
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/age.mixing.MAR.a.fun.R")
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/age.mixing.MAR.b.fun.R")
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/age.mixing.MAR.c.fun.R")
  
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/needed.functions.RSimpactHelp.R")
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/advanced.transmission.network.builder.R")
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/age.mixing.MCAR.fun.R")
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/age.mixing.MAR.a.fun.R")
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/age.mixing.MAR.b.fun.R")
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/age.mixing.MAR.c.fun.R")
  
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  work.dir <- "/home/dniyukuri/lustre/agemix.25.10.2018.2" # on PCHPC
  
  
  setwd(paste0(work.dir))
  
  
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
  
  library(lubridate)
  
  library(tidyr)
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  
  
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 5000, 
                                   population.numwomen = 5000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 40, 
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
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    results.mcar <- rep(NA, 487) # 37 + 82 + 33 + 1 + 1 + 53 = 207
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      results.mcar <- rep(NA, 487)
    }else{
      
      
      DataListALL <- readthedata(results)
      
      
      datalist.agemix <- DataListALL
      
      
      
      # datalist.agemix <- get(load("datalist.agemix.RData"))
      
      
      
      
      ###########################################
      # Step 2: Construct transmission networks #
      ###########################################
      
      
      simpact.trans.net.adv <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
      
      
      net.size.vector <- vector() # i_th seed in the list of seeds
      
      for(i in 1:length(simpact.trans.net.adv)){
        
        tree.n <- simpact.trans.net.adv[[i]] # transmission network for i^th seed
        
        net.size.vector <- c(net.size.vector, nrow(as.data.frame(tree.n)))
        
      }
      
      big.index <- which(net.size.vector>=100)
      
      if(length(big.index) >= 1){ 
        
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
        
        
        pp.cp.6months.male.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                                timepoint = 40 - 0.5) %>%
          dplyr::select(concurr.pointprev) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        pp.cp.6months.female.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                                  timepoint = 40 - 0.5) %>%
          dplyr::select(concurr.pointprev) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        
        # (iii) Prevalence
        ##################
        
        
        hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 40),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 40),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 50),
                                                      timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 50),
                                                    timepoint = 40) %>%
          dplyr::select(pointprevalence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        
        # (iv) Incidence
        #################
        
        
        
        epi.rels.incidence.df.15.24.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(15, 25),
                                                                timewindow = c(39, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.15.24.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(15, 25),
                                                                  timewindow = c(39, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        epi.rels.incidence.df.25.39.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(25, 40),
                                                                timewindow = c(39, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.25.39.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(25, 40),
                                                                  timewindow = c(39, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        epi.rels.incidence.df.40.49.men <- incidence.calculator(datalist = datalist.agemix,
                                                                agegroup = c(40, 50),
                                                                timewindow = c(39, 40),
                                                                only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(1) %>%
          as.numeric()
        
        epi.rels.incidence.df.40.49.women <- incidence.calculator(datalist = datalist.agemix,
                                                                  agegroup = c(40, 50),
                                                                  timewindow = c(39, 40),
                                                                  only.active = "No") %>%
          dplyr::select(incidence) %>%
          dplyr::slice(2) %>%
          as.numeric()
        
        
        
        summary.epidemic.rels.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                                      hiv.prev.25.40.women, hiv.prev.25.40.men,
                                      hiv.prev.40.50.women, hiv.prev.40.50.men, 
                                      mix.rels.dat,
                                      pp.cp.6months.male.rels, pp.cp.6months.female.rels, 
                                      
                                      epi.rels.incidence.df.15.24.men, epi.rels.incidence.df.15.24.women, 
                                      epi.rels.incidence.df.25.39.men, epi.rels.incidence.df.25.39.women,
                                      epi.rels.incidence.df.40.49.men, epi.rels.incidence.df.40.49.women)
        
        names(summary.epidemic.rels.df) <- c("R.prev.15.25.w", "R.prev.15.25.m", "R.prev.25.40.w", "R.prev.25.40.m", "R.prev.40.50.w", "R.prev.40.50.m",
                                             names(mix.rels.dat), 
                                             "R.p.prev.6months.m","R.p.prev.6months.f",
                                             "R.inc.15.25.w", "R.inc.15.25.m", "R.inc.25.40.w", "R.inc.25.40.m", "R.inc.40.50.w", "R.inc.40.50.m")
        
        
        # 
        # mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net.adv, # simpact.trans.net 
        #                            limitTransmEvents = 7,
        #                            timewindow = c(40,40), 
        #                            seq.cov = 50, 
        #                            age.limit = 50)
        # 
        
        
        # MCAR
        
        res.clust.MCAR.cov.35 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 35,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.40 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 40,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.45 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 45,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.50 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 50,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.55 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 55,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MCAR.cov.60 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 60,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MCAR.cov.65 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 65,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MCAR.cov.70 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 70,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        
        res.clust.MCAR.cov.75 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 75,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MCAR.cov.80 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 80,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.85 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 85,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MCAR.cov.90 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 90,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        res.clust.MCAR.cov.95 <- tryCatch(age.mixing.MCAR.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                              datalist.agemix = datalist.agemix, 
                                                              work.dir = work.dir,  
                                                              dirfasttree = work.dir, 
                                                              sub.dir.rename = sub.dir.rename,
                                                              limitTransmEvents = 7,
                                                              timewindow = c(30,40),
                                                              seq.cov = 95,
                                                              age.group.15.25 = c(15,25),
                                                              age.group.25.40 = c(25,40),
                                                              age.group.40.50 = c(40,50),
                                                              cut.off = 7),
                                          error=function(e) return(rep(NA, 133)))
        
        
        # MAR
        # (a) 0.7
        
        
        res.clust.MAR.a.cov.35 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 35,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.40 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 40,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.45 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 45,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.50 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 50,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.55 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 55,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.a.cov.60 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 60,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.a.cov.65 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 65,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.a.cov.70 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 70,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        res.clust.MAR.a.cov.75 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 75,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.a.cov.80 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 80,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.85 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 85,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.a.cov.90 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 90,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.a.cov.95 <- tryCatch(age.mixing.MAR.a.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 95,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        # MAR
        # (b) 0.3
        
        
        res.clust.MAR.b.cov.35 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 35,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.40 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 40,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.45 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 45,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.50 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 50,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.55 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 55,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.b.cov.60 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 60,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.b.cov.65 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 65,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.b.cov.70 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 70,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        res.clust.MAR.b.cov.75 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 75,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.b.cov.80 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 80,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.85 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 85,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.b.cov.90 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 90,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.b.cov.95 <- tryCatch(age.mixing.MAR.b.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 95,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        # MAR
        # (c) 0.5
        
        
        res.clust.MAR.c.cov.35 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 35,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.40 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 40,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.45 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 45,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.50 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 50,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.55 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 55,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.c.cov.60 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 60,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        res.clust.MAR.c.cov.65 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 65,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.c.cov.70 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 70,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        res.clust.MAR.c.cov.75 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 75,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.c.cov.80 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 80,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.85 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 85,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        
        
        res.clust.MAR.c.cov.90 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 90,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        res.clust.MAR.c.cov.95 <- tryCatch(age.mixing.MAR.c.fun(simpact.trans.net = simpact.trans.net.adv,  
                                                                datalist.agemix = datalist.agemix, 
                                                                work.dir = work.dir,  
                                                                dirfasttree = work.dir, 
                                                                sub.dir.rename = sub.dir.rename,
                                                                limitTransmEvents = 7,
                                                                timewindow = c(30,40),
                                                                seq.cov = 95,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50),
                                                                cut.off = 7),
                                           error=function(e) return(rep(NA, 133)))
        
        
        results.summary.epidemic.rels.df <- summary.epidemic.rels.df
        
        
        names(results.summary.epidemic.rels.df) <- paste0(names(summary.epidemic.rels.df))
        
        
        
        results.mcar <- c(res.clust.MCAR.cov.35, res.clust.MCAR.cov.40, 
                          res.clust.MCAR.cov.45, res.clust.MCAR.cov.50, 
                          res.clust.MCAR.cov.55, res.clust.MCAR.cov.60, 
                          res.clust.MCAR.cov.65, res.clust.MCAR.cov.70, 
                          res.clust.MCAR.cov.75, res.clust.MCAR.cov.80, 
                          res.clust.MCAR.cov.85,res.clust.MCAR.cov.90, 
                          res.clust.MCAR.cov.95)
        
        names(results.mcar) <- c(paste0("cov.MCAR.",35,".",paste0(names(res.clust.MCAR.cov.35))), paste0("cov.MCAR.",40,".",paste0(names(res.clust.MCAR.cov.40))),
                                 paste0("cov.MCAR.",45,".",paste0(names(res.clust.MCAR.cov.45))), paste0("cov.MCAR.",50,".",paste0(names(res.clust.MCAR.cov.50))),
                                 paste0("cov.MCAR.",55,".",paste0(names(res.clust.MCAR.cov.55))), paste0("cov.MCAR.",60,".",paste0(names(res.clust.MCAR.cov.60))),
                                 paste0("cov.MCAR.",65,".",paste0(names(res.clust.MCAR.cov.65))), paste0("cov.MCAR.",70,".",paste0(names(res.clust.MCAR.cov.70))),
                                 paste0("cov.MCAR.",75,".",paste0(names(res.clust.MCAR.cov.75))), paste0("cov.MCAR.",80,".",paste0(names(res.clust.MCAR.cov.80))),
                                 paste0("cov.MCAR.",85,".",paste0(names(res.clust.MCAR.cov.85))), paste0("cov.MCAR.",90,".",paste0(names(res.clust.MCAR.cov.90))),
                                 paste0("cov.MCAR.",95,".",paste0(names(res.clust.MCAR.cov.95))))
        
        
        
        results.mcar.a <- c(res.clust.MAR.a.cov.35, res.clust.MAR.a.cov.40, 
                            res.clust.MAR.a.cov.45, res.clust.MAR.a.cov.50, 
                            res.clust.MAR.a.cov.55, res.clust.MAR.a.cov.60, 
                            res.clust.MAR.a.cov.65, res.clust.MAR.a.cov.70, 
                            res.clust.MAR.a.cov.75, res.clust.MAR.a.cov.80, 
                            res.clust.MAR.a.cov.85,res.clust.MAR.a.cov.90, 
                            res.clust.MAR.a.cov.95)
        
        names(results.mcar.a) <- c(paste0("cov.MAR.a.",35,".",paste0(names(res.clust.MAR.a.cov.35))), paste0("cov.MAR.a.",40,".",paste0(names(res.clust.MAR.a.cov.40))),
                                   paste0("cov.MAR.a.",45,".",paste0(names(res.clust.MAR.a.cov.45))), paste0("cov.MAR.a.",50,".",paste0(names(res.clust.MAR.a.cov.50))),
                                   paste0("cov.MAR.a.",55,".",paste0(names(res.clust.MAR.a.cov.55))), paste0("cov.MAR.a.",60,".",paste0(names(res.clust.MAR.a.cov.60))),
                                   paste0("cov.MAR.a.",65,".",paste0(names(res.clust.MAR.a.cov.65))), paste0("cov.MAR.a.",70,".",paste0(names(res.clust.MAR.a.cov.70))),
                                   paste0("cov.MAR.a.",75,".",paste0(names(res.clust.MAR.a.cov.75))), paste0("cov.MAR.a.",80,".",paste0(names(res.clust.MAR.a.cov.80))),
                                   paste0("cov.MAR.a.",85,".",paste0(names(res.clust.MAR.a.cov.85))), paste0("cov.MAR.a.",90,".",paste0(names(res.clust.MAR.a.cov.90))),
                                   paste0("cov.MAR.a.",95,".",paste0(names(res.clust.MAR.a.cov.95))))
        
        results.mcar.b <- c(res.clust.MAR.b.cov.35, res.clust.MAR.b.cov.40, 
                            res.clust.MAR.b.cov.45, res.clust.MAR.b.cov.50, 
                            res.clust.MAR.b.cov.55, res.clust.MAR.b.cov.60, 
                            res.clust.MAR.b.cov.65, res.clust.MAR.b.cov.70, 
                            res.clust.MAR.b.cov.75, res.clust.MAR.b.cov.80, 
                            res.clust.MAR.b.cov.85,res.clust.MAR.b.cov.90, 
                            res.clust.MAR.b.cov.95)
        
        names(results.mcar.b) <- c(paste0("cov.MAR.b.",35,".",paste0(names(res.clust.MAR.a.cov.35))), paste0("cov.MAR.b.",40,".",paste0(names(res.clust.MAR.a.cov.40))),
                                   paste0("cov.MAR.b.",45,".",paste0(names(res.clust.MAR.a.cov.45))), paste0("cov.MAR.b.",50,".",paste0(names(res.clust.MAR.a.cov.50))),
                                   paste0("cov.MAR.b.",55,".",paste0(names(res.clust.MAR.a.cov.55))), paste0("cov.MAR.b.",60,".",paste0(names(res.clust.MAR.a.cov.60))),
                                   paste0("cov.MAR.b.",65,".",paste0(names(res.clust.MAR.a.cov.65))), paste0("cov.MAR.b.",70,".",paste0(names(res.clust.MAR.a.cov.70))),
                                   paste0("cov.MAR.b.",75,".",paste0(names(res.clust.MAR.a.cov.75))), paste0("cov.MAR.b.",80,".",paste0(names(res.clust.MAR.a.cov.80))),
                                   paste0("cov.MAR.b.",85,".",paste0(names(res.clust.MAR.a.cov.85))), paste0("cov.MAR.b.",90,".",paste0(names(res.clust.MAR.a.cov.90))),
                                   paste0("cov.MAR.b.",95,".",paste0(names(res.clust.MAR.a.cov.95))))
        
        
        results.mcar.c <- c(res.clust.MAR.c.cov.35, res.clust.MAR.c.cov.40, 
                            res.clust.MAR.c.cov.45, res.clust.MAR.c.cov.50, 
                            res.clust.MAR.c.cov.55, res.clust.MAR.c.cov.60, 
                            res.clust.MAR.c.cov.65, res.clust.MAR.c.cov.70, 
                            res.clust.MAR.c.cov.75, res.clust.MAR.c.cov.80, 
                            res.clust.MAR.c.cov.85,res.clust.MAR.c.cov.90, 
                            res.clust.MAR.c.cov.95)
        
        
        names(results.mcar.c) <- c(paste0("cov.MAR.c.",35,".",paste0(names(res.clust.MAR.a.cov.35))), paste0("cov.MAR.c.",40,".",paste0(names(res.clust.MAR.a.cov.40))),
                                   paste0("cov.MAR.c.",45,".",paste0(names(res.clust.MAR.a.cov.45))), paste0("cov.MAR.c.",50,".",paste0(names(res.clust.MAR.a.cov.50))),
                                   paste0("cov.MAR.c.",55,".",paste0(names(res.clust.MAR.a.cov.55))), paste0("cov.MAR.c.",60,".",paste0(names(res.clust.MAR.a.cov.60))),
                                   paste0("cov.MAR.c.",65,".",paste0(names(res.clust.MAR.a.cov.65))), paste0("cov.MAR.c.",70,".",paste0(names(res.clust.MAR.a.cov.70))),
                                   paste0("cov.MAR.c.",75,".",paste0(names(res.clust.MAR.a.cov.75))), paste0("cov.MAR.c.",80,".",paste0(names(res.clust.MAR.a.cov.80))),
                                   paste0("cov.MAR.c.",85,".",paste0(names(res.clust.MAR.a.cov.85))), paste0("cov.MAR.c.",90,".",paste0(names(res.clust.MAR.a.cov.90))),
                                   paste0("cov.MAR.c.",95,".",paste0(names(res.clust.MAR.a.cov.95))))
        
        
        results.outputvector <- c(results.summary.epidemic.rels.df, 
                                  results.mcar, results.mcar.a, results.mcar.b, results.mcar.c)
        
        
        
      }else{
        
        results.outputvector <- rep(NA, 6936) # 133*13*4 + 134+ 6
        
      }
      
    }
    
  }
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  
  return(results.outputvector)
  
  
}

