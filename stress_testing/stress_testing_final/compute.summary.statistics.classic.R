



compute.summary.statistics.classic <- function(datalist = datalist.agemix,
                                               simpact.trans.net = simpact.trans.net,
                                               work.dir = work.dir,
                                               sub.dir.rename = sub.dir.rename,
                                               dirfasttree = work.dir,
                                               limitTransmEvents = 7,
                                               seq.cov = 100,
                                               age.group.15.25 = c(15,25),
                                               age.group.25.40 = c(25,40),
                                               age.group.40.50 = c(40,50),
                                               endpoint = 40,
                                               timewindow = c(30,40)){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- paste0(work.dir)
  
  
  datalist.agemix <- datalist
  
  
  # New datalist for relationships
  
  
  age.limit = age.group.40.50[2]
  
  
  datalist.new <- datalist
  
  
  person.datalist.new <- datalist.new$ptable
  
  perc.100 <- nrow(person.datalist.new)
  
  person.TOB.datalist.new <- person.datalist.new
  
  person.TOB.datalist.new$AgeDeath <- abs(person.TOB.datalist.new$TOB) + person.TOB.datalist.new$TOD
  
  person.TOB.datalist.new$AgEndpoint <- abs(person.TOB.datalist.new$TOB) + endpoint
  
  person.TOB.datalist.new$AgeLowWindow <- abs(person.TOB.datalist.new$TOB) + timewindow[1]
  
  person.TOB.datalist.new$AgeUppWindow <- abs(person.TOB.datalist.new$TOB) + timewindow[2]
  
  
  men.women.datalist.new.df.alive <- dplyr::filter(person.TOB.datalist.new, 
                                                   person.TOB.datalist.new$AgeDeath=="Inf" & person.TOB.datalist.new$AgeLowWindow >= age.group.15.25[1] & person.TOB.datalist.new$AgeUppWindow < age.group.40.50[2])
  
  men.women.datalist.new.df.died <- dplyr::filter(person.TOB.datalist.new, 
                                                  person.TOB.datalist.new$AgeDeath!="Inf" & person.TOB.datalist.new$AgeDeath >= age.group.15.25[1] & person.TOB.datalist.new$AgeDeath < age.group.40.50[2])
  
  
  men.women.datalist.new.df <- rbind(men.women.datalist.new.df.alive, men.women.datalist.new.df.died)
  
  
  perc.100.limit.window <- nrow(men.women.datalist.new.df)
  
  
  perc.seq.coverage <- round(perc.100.limit.window*seq.cov/100) # total number of wanted individuals at seq.cov sequence coverage
  
  
  samp.IDs <- sample(men.women.datalist.new.df$ID, perc.seq.coverage)
  
  
  
  
  
  # Persons' table of selected individuals within the time winedow
  
  pers.datalist.selected <- dplyr::filter(person.datalist.new, person.datalist.new$ID%in%samp.IDs)
  
  
  # Data list of infected individuals within the time window
  
  datalist.selec <- datalist
  
  datalist.selec$ptable <- pers.datalist.selected
  
  
  # Using new datalist for next steps
  
  datalist.agemix <- datalist.selec
  
  
  ########################################
  # I. Behavioural and epidemic features #
  ########################################
  
  
  # 1.2. Features from sexual and transmission network
  
  # 
  # 1.2.1. Demographic feature:
  
  #   (i) Population growth rate (pop.growth.calculator function)
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = timewindow) # c(0, datalist.agemix$itable$population.simtime[1])
  
  
  # 1.2.2. Transmission features:	
  
  #   (i) Prevalence (prevalence.calculator function)
  # 
  #   hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
  #                                                agegroup = c(15, 25),
  #                                                timepoint = datalist.agemix$itable$population.simtime[1])$pointprevalence[2]
  #   
  
  hiv.prev.15.25.women <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                                agegroup = age.group.15.25,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.15.25.men <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                              agegroup = age.group.15.25,
                                              timepoint = endpoint)$pointprevalence[1]
  hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                                agegroup = age.group.25.40,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                              agegroup = age.group.25.40,
                                              timepoint = endpoint)$pointprevalence[1]
  hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                                agegroup = age.group.40.50,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist, # datalist.agemix,
                                              agegroup = age.group.40.50,
                                              timepoint = endpoint)$pointprevalence[1]
  
  
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
  
  
  agemix.rels.df <- agemix.df.maker(datalist.agemix)
  
  # 
  agemix.model <- pattern.modeller(dataframe = agemix.rels.df,
                                   agegroup = c(15, 50),
                                   timepoint = endpoint, # datalist.agemix$itable$population.simtime[1],
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
    
    names(mix.rels.dat) <- c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
    names(mix.rels.dat) <- c("R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male")
    
  }
  
  
  pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                     timepoint = endpoint - 0.5)
  
  
  classic.features <- c(growthrate, 
                        
                        hiv.prev.15.25.women, hiv.prev.15.25.men, hiv.prev.25.40.women,
                        hiv.prev.25.40.men, hiv.prev.40.50.women, hiv.prev.40.50.men, 
                        
                        relsperpersonperyear, agegapsd,
                        
                        mix.rels.dat,
                        
                        pp.cp.6months.male)
  
  
  classic.features.names <- c("Pop.growthrate", 
                              
                              "hiv.prev.15.25.women", "hiv.prev.15.25.men", "hiv.prev.25.40.women",
                              "hiv.prev.25.40.men", "hiv.prev.40.50.women", "hiv.prev.40.50.men",
                              
                              "relsperpersonperyear", "agegapsd",
                              
                              "R.AAD.male", "R.SDAD.male", "R.slope.male", "R.WSD.male", "R.BSD.male", "R.intercept.male",
                              
                              "pp.cp.6months.male")
  
  names(classic.features) <- classic.features.names 
  
  
  outputvector.epi.rels.behav  <- classic.features
  
  
  #############################
  # II. Transmission features #
  #############################
  
  
  # New datalist for transmissions
  
  # Data table of infected individuals with at least limitTransmEvents transmissions events
  
  # extract sampling time for further computations
  
  seeds.id <- length(simpact.trans.net)
  
  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)
  
  for(i in 1:seeds.id){
    
    transm.age.df.ic <- as.data.frame(simpact.trans.net[[i]])
    
    age.samp.Rec <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBRec
    age.samp.Don <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBDon
    
    transm.age.i <- cbind(transm.age.df.ic, age.samp.Rec, age.samp.Don)
    
    new.transm.tab[[i]] <- transm.age.i
    
  }
  
  # id of people who got infection by seed event: seeds.id
  trans.network <- new.transm.tab
  
  seeds.id <- length(trans.network)
  
  
  ID.select <- vector() # ID of selected transmission network
  ID.select.count <- vector() # number of individuals in these networks
  
  for (i in 1: seeds.id) {
    
    
    trans.network.i <- as.data.frame(trans.network[[i]])
    
    if(nrow(trans.network.i)>=limitTransmEvents){
      
      
      ID.select <- c(ID.select, i)
      ID.select.count <- c(ID.select.count, nrow(trans.network.i))
      
    } # X if
    
  } # Y for
  
  
  infectionTable <- vector("list", length(ID.select))
  
  for(j in 1:length(ID.select)){
    
    p <- ID.select[j]
    
    trans.network.i <- as.data.frame(trans.network[[p]])
    
    trans.network.i <- trans.network.i[-1,]
    
    id.lab <- paste0(p,".",trans.network.i$id,".C")
    
    trans.network.i$id.lab <- id.lab
    
    infectionTable[[p]] <- trans.network.i
  }
  
  
  infecttable <- rbindlist(infectionTable)
  
  
  # IDs fro completely random subset of sequences within the given time window
  
  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net, 
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow, 
                             seq.cov = seq.cov, 
                             age.limit = age.group.40.50[2])
  
  
  
  # Data table of infected individuals within the time window
  
  data.transm.agemix <- dplyr::filter(infecttable, infecttable$id.lab%in%mCAr.IDs) 
  
  
  # Data list of infected individuals within the time window
  
  datalist.agemix.transm <- datalist
  
  datalist.agemix.transm$ptable <- dplyr::filter(datalist.agemix.transm$ptable, datalist.agemix.transm$ptable$ID%in%data.transm.agemix$RecId)
  
  
  # (i) Age mixing in transmissions
  
  agemix.transm.df <- agemix.df.maker(datalist.agemix.transm)
  
  # 
  agemix.model <- pattern.modeller(dataframe = agemix.transm.df,
                                   agegroup = c(age.group.15.25[1], age.group.40.50[2]),
                                   timepoint = 40, # datalist.agemix$itable$population.simtime[1],
                                   timewindow = 10)#1)#3)
  
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
    
    mix.rels.dat.names <- c("TR.AAD.male", "TR.SDAD.male", "TR.slope.male", "TR.WSD.male", "TR.BSD.male", "TR.intercept.male")

    
    names(mix.rels.dat) <- mix.rels.dat.names 
    
    
    outputvector.epi.rels.transm  <- mix.rels.dat
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
    mix.rels.dat.names <- c("TR.AAD.male", "TR.SDAD.male", "TR.slope.male", "TR.WSD.male", "TR.BSD.male", "TR.intercept.male")
    
    names(mix.rels.dat) <- mix.rels.dat.names 
    
    outputvector.epi.rels.transm  <- mix.rels.dat
    
  }
  
  
  # (ii) Fiting age mixing transmission table with mixed-effect linear model
  
  
  het.fit.lme.agemixing <- lme(age.samp.Rec ~ GenderRec, data = data.transm.agemix, random = ~ 1|DonId,
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
  
  names(het.lme.val) <-  c("het.av.age.male", "het.gendEffect", "het.between.transm.var", "het.within.transm.var", "het.SD.female", "het.SD.male")
  
  
  outputvector.lme.transm <- het.lme.val
  
  
  features.here <- c(outputvector.epi.rels.behav, outputvector.epi.rels.transm, outputvector.lme.transm)
  
  
  return(features.here)
  
}


