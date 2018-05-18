

age.group.gender.men30.women70.compute.summary.statistics.combined.95 <- function(datalist = datalist.agemix,
                                                                                  work.dir = work.dir,
                                                                                  sub.dir.rename = sub.dir.rename){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- paste0(work.dir)
  
  
  datalist.agemix <- datalist
  
  # datalist.agemix <- get(load("datalist.agemix.RData"))
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  #################################### Features #############################
  
  # 1.2. Features from sexual and transmission network
  
  
  # 
  # 1.2.1. Demographic feature:
  
  #   (i) Population growth rate (pop.growth.calculator function)
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = c(0,40)) # c(0, datalist.agemix$itable$population.simtime[1])
  
  
  # 1.2.2. Transmission features:	
  
  #   (i) Prevalence (prevalence.calculator function)
  # 
  #   hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
  #                                                agegroup = c(15, 25),
  #                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  #   
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = 40)$pointprevalence[2]
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = 40)$pointprevalence[1]
  hiv.prev.25.34.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 35),
                                                timepoint = 40)$pointprevalence[2]
  hiv.prev.25.34.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 35),
                                              timepoint = 40)$pointprevalence[1]
  hiv.prev.35.44.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(35, 45),
                                                timepoint = 40)$pointprevalence[2]
  hiv.prev.35.44.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(35, 45),
                                              timepoint = 40)$pointprevalence[1]
  
  
  # (ii) Transmission 	rate (transmission.rate.calculator function)
  
  # transm.rate <- transmission.rate.calculator(datalist = datalist.agemix,
  #                                             timewindow = c(10, 40), 
  #                                             int = FALSE, by=1)
  # 
  # 
  
  # (iii) ART coverage
  
  # cov.vector <- ART.coverage.vector.creator(datalist = datalist.agemix,
  #                                          agegroup = c(15, 50))
  # plot(cov.vector)
  
  # ART.coverage.vector.creator <- function(datalist = datalist,
  #          agegroup = c(15, 50)){
  #   ART.cov.eval.timepoints <- seq(from = datalist$itable$t[2],
  #                                  to = 40 #datalist$itable$population.simtime[1])
  #   ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  #   for (art.cov.index in 1:length(ART.cov.vector)){
  #     ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist,
  #                                                                  agegroup = agegroup,
  #                                                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
  #       sum(ART.coverage.calculator(datalist = datalist,
  #                                   agegroup = agegroup,
  #                                   timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  #   }
  #   return(ART.cov.vector)
  
  
  # 1.2.3. Sexual behaviour features:
  
  #  (i) Relationship 	rate (relationship.rate.calculator function) - rate of new relationship formation (partner turnover rate):  	
  
  # 
  # relas.rate <- relationship.rate.calculator(datalist = datalist.agemix,
  #                                            timewindow = c(0, 40), 
  #                                            int = FALSE, by=1)
  
  # (ii) Relationship per person per year
  
  relsperpersonperyear <- nrow(datalist.agemix$rtable) / (nrow(datalist.agemix$ptable)/2) / 40 #cfg.list$population.simtime
  
  # (iv) SD age gap between couples
  
  agegapsd <- sd(datalist.agemix$rtable$AgeGap)
  
  
  # (v) Age mixing in relationships
  
  # 
  # agemix.df <- agemix.df.maker(datalist.agemix)
  # 
  # agemix.model <- pattern.modeller(dataframe = agemix.df,
  #                                  agegroup = c(15, 60),
  #                                  timepoint = 40, # datalist.agemix$itable$population.simtime[1],
  #                                  timewindow = 5)#1)#3)
  # 
  # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  #
  # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
  
  # men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
  #                  data = dplyr::filter(agemix.model[[1]], Gender =="male"),
  #                  REML = TRUE,
  #                  control=lmerControl(check.nobs.vs.nlev = "ignore",
  #                                      check.nobs.vs.rankZ = "ignore",
  #                                      check.nobs.vs.nRE="ignore"))
  # #
  # bignumber <- NA # let's try if NA works (instead of 9999 for example)
  # #
  # AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  # SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
  # #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
  # slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
  # WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
  # 
  # BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
  # 
  # intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
  
  # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
  
  ## AAD: average age difference across all relationship
  ## VAD: variance of these age differences
  ## SDAD: standard deviation of age differences
  ## BSD: between-subject standard deviation of age differences
  
  
  # age.scatter.df <- agemix.model[[1]]
  
  #  (iii) Point 	prevalence of concurrency in the adult population:
  
  # Concurrency point prevalence 6 months before a survey, among men
  
  # pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
  #                                                    timepoint = datalist.agemix$itable$population.simtime[1] - 0.5)
  # 
  
  pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                     timepoint = 40 - 0.5)
  
  # 
  #   c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.34, METRICS.incidence.df.35.44,
  #     
  #     METRICS.incidence.df.15.24.int.40.41, METRICS.incidence.df.25.34.int.40.41,
  #     METRICS.incidence.df.35.44.int.40.41,
  #     METRICS.incidence.df.15.24.int.41.42, METRICS.incidence.df.25.34.int.41.42,
  #     METRICS.incidence.df.35.44.int.41.42,
  #     METRICS.incidence.df.15.24.int.42.43, METRICS.incidence.df.25.34.int.42.43,
  #     METRICS.incidence.df.35.44.int.42.43,
  #     METRICS.incidence.df.15.24.int.43.44, METRICS.incidence.df.25.34.int.43.44,
  #     METRICS.incidence.df.35.44.int.43.44,
  #     METRICS.incidence.df.15.24.int.44.45, METRICS.incidence.df.25.34.int.44.45,
  #     METRICS.incidence.df.35.44.int.44.45,
  #     
  #     METRICS.age.mix.trans.interc, METRICS.age.mix.slope, METRICS.transm.average,
  #     
  #     growthrate, hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
  #     hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, transm.rate, # cov.vector
  #     relas.rate,  relsperpersonperyear, agegapsd,
  #     
  #     # AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
  #     pp.cp.6months.male
  #   )
  
  
  
  ###############################
  # Step 3: Sequence simulation #
  ###############################
  
  
  trans.net <- simpact.trans.net # all transmission networks
  
  
  dirseqgen <- work.dir
  
  seeds.num <- 123
  
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
  
  
  
  #####################################################
  # Step 4: Construct time stamped phylogenetic trees #   + Sequence coverage Scenarios
  #####################################################
  
  
  dirfasttree <- work.dir
  
  
  if(file.exists(paste0(sub.dir.rename, "/C.Epidemic_seed.seq.bis.sim.nwk.fasta"))==TRUE){
    
    # check if the run has simulated sequence,
    # in other words: we have a transmission network with at least 6 individuals
    
    
    
    ### Sequence coverage sceanrios  ###: 35 - 45 - 55 - 65 - 75 - 85 - 95 - 100
    ####################################
    
    
    ## Transform the seuences in fasta format
    sequ.dna. <- read.dna(paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna., file = paste0(sub.dir.rename,"/C.Epidemic.Fasta"), format = "fasta")
    
    
    # Sequence coverage Scenarios
    #############################
    
    
    
    #### BEGIN Sequence Coverage Scenarios
    
    ### 7th Scenario: 95 ###
    ########################
    
    cov.95.gender.men30.women70.age.group <- IDs.indiv.seq.gender.age.group.fun(simpact.trans.net = simpact.trans.net,
                                                                                limitTransmEvents = 7, 
                                                                                perc.men = 30, 
                                                                                seq.cov = 95,
                                                                                age.men = c(15, 60),
                                                                                age.women = c(15, 40))
    cov.95.IDs.gender.men30.women70.age.group <- cov.95.gender.men30.women70.age.group$outputvector
    
    
    
    cut.val = 5
    
    source("~/phylosimpact_simulation_studies_2018/stress_testing/phylogenetic.features.fun.R")
    # men.50.women50
    if(length(cov.95.IDs.gender.men30.women70.age.group)>=cut.val){
      
      choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.Fasta"),
                          select.vec = cov.95.IDs.gender.men30.women70.age.group, 
                          name.file = paste0(sub.dir.rename,"/cov.95.IDs.gender.men30.women70.age.group.C.Epidemic.Fasta"))
      
      
      cov.95.IDs.gender.men30.women70.age.group.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                                                             sub.dir.rename = sub.dir.rename,
                                                                                             fasttree.tool = "FastTree",
                                                                                             calendar.dates = "samplingtimes.all.csv",
                                                                                             simseqfile = "cov.95.IDs.gender.men30.women70.age.group.C.Epidemic.Fasta",
                                                                                             count.start = 1977,
                                                                                             endsim = 40,
                                                                                             clust = FALSE)
      cov.95.IDs.gender.men30.women70.age.group.tree.calib.LTT <- cov.95.IDs.gender.men30.women70.age.group.tree.calib
      
      write.tree(cov.95.IDs.gender.men30.women70.age.group.tree.calib, file = paste0(sub.dir.rename,"/cov.95.IDs.gender.men30.women70.age.group.calibrated.tree.nwk"))
      
      cov.95.IDs.gender.men30.women70.age.group.features <- phylogenetic.features.fun(tree.topo=paste0(sub.dir.rename,"/cov.95.IDs.gender.men30.women70.age.group.calibrated.tree.nwk"))
      # ,
      #                                                                                 tree.calib.LTT = cov.95.IDs.gender.men30.women70.age.group.tree.calib.LTT)
      # 
      cov.95.IDs.gender.men30.women70.age.group.features <- as.numeric(cov.95.IDs.gender.men30.women70.age.group.features)
      cov.95.gender.men30.women70.age.group.ratio.seq <- cov.95.gender.men30.women70.age.group$ratio.seq
      cov.95.gender.men30.women70.age.group.ratio.emp <- cov.95.gender.men30.women70.age.group$ratio.emp
      
      cov.95.IDs.gender.men30.women70.age.group.features <- c(cov.95.IDs.gender.men30.women70.age.group.features,
                                                              cov.95.gender.men30.women70.age.group.ratio.seq,
                                                              cov.95.gender.men30.women70.age.group.ratio.emp)
      
    }else{
      cov.95.IDs.gender.men30.women70.age.group.features <- rep(NA, 8)
    }
    
    
    
    #### END Sequence Coverage Scenarios
    
    summary.df <- c(growthrate, 
                    
                    hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
                    hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, 
                    # transm.rate, # cov.vector
                    # relas.rate,  
                    relsperpersonperyear, agegapsd,
                    
                    # AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                    
                    pp.cp.6months.male,
                    
                    
                    as.numeric(cov.95.IDs.gender.men30.women70.age.group.features)
                    
                    # cov.95.gender.men30.women70.ratio.seq, cov.95.gender.men30.women70.ratio.emp
    )
    
    
    features.names <- c("Pop.growthrate", 
                        
                        "hiv.prev.lt25.women", "hiv.prev.lt25.men", "hiv.prev.25.34.women",
                        "hiv.prev.25.34.men", "hiv.prev.35.44.women", "hiv.prev.35.44.men",
                        # "transm.rate", # cov.vector
                        # "relas.rate",  
                        "relsperpersonperyear", "agegapsd",
                        
                        ## "AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        
                        "pp.cp.6months.male",
                        
                        "cov.95.IDs.gender.men30.women70.age.group.meanHeight.feature", "cov.95.IDs.gender.men30.women70.age.group.colless.feature", "cov.95.IDs.gender.men30.women70.age.group.sackin.feature", 
                        "cov.95.IDs.gender.men30.women70.age.group.mean.tipsDepths.feature", "cov.95.IDs.gender.men30.women70.age.group.mean.nodesDepths.feature",
                        "cov.95.IDs.gender.men30.women70.age.group.maxHeight.feature", 
                        
                        
                        "cov.95.gender.men30.women70.age.group.ratio.seq", "cov.95.gender.men30.women70.age.group.ratio.emp"
    )
    
    names(summary.df) <- features.names # > length(features.names) [1] 549
    
    
  }else{
    
    features.names <- c("Pop.growthrate", 
                        
                        "hiv.prev.lt25.women", "hiv.prev.lt25.men", "hiv.prev.25.34.women",
                        "hiv.prev.25.34.men", "hiv.prev.35.44.women", "hiv.prev.35.44.men",
                        # "transm.rate", # cov.vector
                        # "relas.rate",  
                        "relsperpersonperyear", "agegapsd",
                        
                        ## "AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        
                        "pp.cp.6months.male",
                        
                        "cov.95.IDs.gender.men30.women70.age.group.meanHeight.feature", "cov.95.IDs.gender.men30.women70.age.group.colless.feature", "cov.95.IDs.gender.men30.women70.age.group.sackin.feature", 
                        "cov.95.IDs.gender.men30.women70.age.group.mean.tipsDepths.feature", "cov.95.IDs.gender.men30.women70.age.group.mean.nodesDepths.feature",
                        "cov.95.IDs.gender.men30.women70.age.group.maxHeight.feature", 
                        
                        
                        "cov.95.gender.men30.women70.age.group.ratio.seq", "cov.95.gender.men30.women70.age.group.ratio.emp"
    )
    
    summary.NA <- rep(NA,8)
    
    summary.df.classic <- c(growthrate, 
                            
                            hiv.prev.lt25.women, hiv.prev.lt25.men, hiv.prev.25.34.women,
                            hiv.prev.25.34.men, hiv.prev.35.44.women, hiv.prev.35.44.men, 
                            # transm.rate, # cov.vector
                            # relas.rate,  
                            relsperpersonperyear, agegapsd,
                            
                            # AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male,
                            
                            pp.cp.6months.male)
    
    summary.df <- c(summary.df.classic,  summary.NA)
    
    names(summary.df) <- features.names
    
  }
  
  return(summary.df)
  
}


