



MAR.compute.classic.phylo.features.cov <- function(datalist = datalist.agemix,
                                                   simpact.trans.net = simpact.trans.net, 
                                                   work.dir = work.dir,
                                                   sub.dir.rename = sub.dir.rename,
                                                   dirfasttree = work.dir,
                                                   limitTransmEvents = 7,
                                                   seq.cov = 35,
                                                   seq.gender.ratio = 0.7,
                                                   age.group.15.25 = c(15,25),
                                                   age.group.25.40 = c(25,40),
                                                   age.group.40.50 = c(40,50),
                                                   endpoint = 40,
                                                   timewindow = c(30,40)){
  
  
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- paste0(work.dir)
  
  
  datalist.agemix <- datalist
  
  
  
  # New datalist for relationships
  
  
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
  
  
  men.datalist.new.df <- dplyr::filter(men.women.datalist.new.df, men.women.datalist.new.df$Gender=="0")
  
  women.datalist.new.df <- dplyr::filter(men.women.datalist.new.df, men.women.datalist.new.df$Gender=="1")
  
  
  # Group 15 - 25
  ###############
  
  sum.men.alive.15.25 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="0", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.15.25[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.15.25[2])
  sum.women.alive.15.25 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="1", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.15.25[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.15.25[2])
  
  sum.men.died.15.25 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="0", men.women.datalist.new.df.died$AgeDeath>=age.group.15.25[1] & men.women.datalist.new.df.died$AgeDeath<age.group.15.25[2])
  sum.women.died.15.25 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="1", men.women.datalist.new.df.died$AgeDeath>=age.group.15.25[1] & men.women.datalist.new.df.died$AgeDeath<age.group.15.25[2])
  
  sum.men.15.25 <- rbind(sum.men.alive.15.25, sum.men.died.15.25)
  sum.women.15.25 <- rbind(sum.women.alive.15.25, sum.women.died.15.25)
  
  
  perc.100.15.25 <- nrow(sum.men.alive.15.25) + nrow(sum.women.alive.15.25) + nrow(sum.men.died.15.25) + nrow(sum.women.died.15.25)
  
  
  perc.women.15.25 <- round(perc.100.15.25 * seq.cov * seq.gender.ratio/100)
  
  
  perc.men.15.25 <- round(perc.100.15.25 * seq.cov * (1-seq.gender.ratio)/100)
  
  
  if(perc.men.15.25 <= length(sum.men.15.25$ID)){
    
    x.id.15.25 <- sample(sum.men.15.25$ID, perc.men.15.25)
    
  }else{
    
    x.id.15.25 <- sample(sum.men.15.25$ID, length(sum.men.15.25$ID))
    
  }
  
  
  if(perc.women.15.25 <= length(sum.women.15.25$ID)){
    
    y.id.15.25 <- sample(sum.women.15.25$ID, perc.women.15.25)
    
  }else{
    
    y.id.15.25 <- sample(sum.women.15.25$ID, length(sum.women.15.25$ID))
    
  }
  
  
  samp.all.15.25 <- c(x.id.15.25, y.id.15.25)
  
  
  
  # Group 25 - 40
  ###############
  
  
  sum.men.alive.25.40 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="0", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.25.40[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.25.40[2])
  sum.women.alive.25.40 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="1", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.25.40[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.25.40[2])
  
  sum.men.died.25.40 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="0", men.women.datalist.new.df.died$AgeDeath>=age.group.25.40[1] & men.women.datalist.new.df.died$AgeDeath<age.group.25.40[2])
  sum.women.died.25.40 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="1", men.women.datalist.new.df.died$AgeDeath>=age.group.25.40[1] & men.women.datalist.new.df.died$AgeDeath<age.group.25.40[2])
  
  sum.men.25.40 <- rbind(sum.men.alive.25.40, sum.men.died.25.40)
  sum.women.25.40 <- rbind(sum.women.alive.25.40, sum.women.died.25.40)
  
  
  perc.100.25.40 <- nrow(sum.men.alive.25.40) + nrow(sum.women.alive.25.40) + nrow(sum.men.died.25.40) + nrow(sum.women.died.25.40)
  
  
  perc.women.25.40 <- round(perc.100.25.40 * seq.cov * seq.gender.ratio/100)
  
  
  perc.men.25.40 <- round(perc.100.25.40 * seq.cov * (1-seq.gender.ratio)/100)
  
  
  if(perc.men.25.40 <= length(sum.men.25.40$ID)){
    
    x.id.25.40 <- sample(sum.men.25.40$ID, perc.men.25.40)
    
  }else{
    
    x.id.25.40 <- sample(sum.men.25.40$ID, length(sum.men.25.40$ID))
    
  }
  
  
  if(perc.women.25.40 <= length(sum.women.25.40$ID)){
    
    y.id.25.40 <- sample(sum.women.25.40$ID, perc.women.25.40)
    
  }else{
    
    y.id.25.40 <- sample(sum.women.25.40$ID, length(sum.women.25.40$ID))
    
  }
  
  
  samp.all.25.40 <- c(x.id.25.40, y.id.25.40)
  
  
  # Group 40 - 50
  ###############
  
  
  sum.men.alive.40.50 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="0", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.40.50[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.40.50[2])
  sum.women.alive.40.50 <- dplyr::filter(men.women.datalist.new.df.alive, men.women.datalist.new.df.alive$Gender=="1", men.women.datalist.new.df.alive$AgeLowWindow>=age.group.40.50[1] & men.women.datalist.new.df.alive$AgeUppWindow<age.group.40.50[2])
  
  sum.men.died.40.50 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="0", men.women.datalist.new.df.died$AgeDeath>=age.group.40.50[1] & men.women.datalist.new.df.died$AgeDeath<age.group.40.50[2])
  sum.women.died.40.50 <- dplyr::filter(men.women.datalist.new.df.died, men.women.datalist.new.df.died$Gender=="1", men.women.datalist.new.df.died$AgeDeath>=age.group.40.50[1] & men.women.datalist.new.df.died$AgeDeath<age.group.40.50[2])
  
  sum.men.40.50 <- rbind(sum.men.alive.40.50, sum.men.died.40.50)
  sum.women.40.50 <- rbind(sum.women.alive.40.50, sum.women.died.40.50)
  
  
  perc.100.40.50 <- nrow(sum.men.alive.40.50) + nrow(sum.women.alive.40.50) + nrow(sum.men.died.40.50) + nrow(sum.women.died.40.50)
  
  
  perc.women.40.50 <- round(perc.100.40.50 * seq.cov * seq.gender.ratio/100)
  
  
  perc.men.40.50 <- round(perc.100.40.50 * seq.cov * (1-seq.gender.ratio)/100)
  
  
  if(perc.men.40.50 <= length(sum.men.40.50$ID)){
    
    x.id.40.50 <- sample(sum.men.40.50$ID, perc.men.40.50)
    
  }else{
    
    x.id.40.50 <- sample(sum.men.40.50$ID, length(sum.men.40.50$ID))
    
  }
  
  
  if(perc.women.40.50 <= length(sum.women.40.50$ID)){
    
    y.id.40.50 <- sample(sum.women.40.50$ID, perc.women.40.50)
    
  }else{
    
    y.id.40.50 <- sample(sum.women.40.50$ID, length(sum.women.40.50$ID))
    
  }
  
  
  samp.all.40.50 <- c(x.id.40.50, y.id.40.50)
  
  
  
  samp.all <- c(samp.all.15.25, samp.all.25.40, samp.all.40.50)
  
  
  samp.IDs <- samp.all
  
  
  
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
  hiv.prev.15.25.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = age.group.15.25,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.15.25.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = age.group.15.25,
                                              timepoint = endpoint)$pointprevalence[1]
  hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = age.group.25.40,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = age.group.25.40,
                                              timepoint = endpoint)$pointprevalence[1]
  hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = age.group.40.50,
                                                timepoint = endpoint)$pointprevalence[2]
  hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist.agemix,
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
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
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
  
  
  outputvector.epi.behav  <- classic.features
  
  
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
  
  
  
  mAr.IDs <- IDs.Seq.Age.Groups(simpact.trans.net = simpact.trans.net,
                                limitTransmEvents = limitTransmEvents,
                                timewindow = timewindow,
                                seq.cov = seq.cov,
                                seq.gender.ratio = seq.gender.ratio,
                                age.group.15.25 = age.group.15.25,
                                age.group.25.40 = age.group.25.40,
                                age.group.40.50 = age.group.40.50)
  
  
  
  
  # Data table of infected individuals within the time window
  
  data.transm.agemix <- dplyr::filter(infecttable, infecttable$id.lab%in%mAr.IDs) 
  
  
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
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
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
  
  
  outputvector.transmission <- het.lme.val
  
  ##############################
  # III. Phylogenetic features #
  ##############################
  
  
  
  
  # Function for linear mixed effect models in transmission clusters
  ###################################################################
  
  
  mixed.effect.fit.transmission.clusters <- function(clust.names=clust.names,
                                                     simpact.trans.net = simpact.trans.net,
                                                     limitTransmEvents = 7){
    
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
    
    
    
    transm.df <- infecttable
    
    
    
    
    # Transmission clusters.
    
    d <- clust.names
    
    transmission.clust.list <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster
    
    
    # Binding all data tables of clusters (age, gender, number of cluster)
    
    clust.size <- vector() # size of each cluster
    
    for (i in 1:length(d)) {
      
      transm.df.cl.dat <- NULL
      
      clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
      
      clust.size <- c(clust.size, nrow(clus.read))
      
      transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
      
      transm.df.cl.dat$age <- transm.df.cl$age.samp.Rec
      transm.df.cl.dat$gender <- transm.df.cl$GenderRec
      transm.df.cl.dat$clust.id <- as.factor(rep(i, nrow(transm.df.cl)))
      
      transmission.clust.list[[i]] <- as.data.frame(transm.df.cl.dat)
      
    }
    
    
    clust.table.df <- as.data.frame(do.call(rbind, transmission.clust.list)) # data.table & data.frame
    
    
    het.fit.lme.agemixing <- lme(age ~ gender, data = clust.table.df, random = ~ 1|clust.id,
                                 weights = varIdent( c("1" = 0.5), ~ 1 |gender))
    
    
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
    
    
    clust.lme.val <- c(het.a, het.beta, het.b1, het.b2, SD.female, SD.male)
    
    
    Num.Clus <- length(d)
    
    av.Clust.size <- sum(clust.size)/length(d)
    
    ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
    
    names(ouptuvector.clust) <- c("clust.av.age.male", "clust.gendEffect", "clust.between.transm.var", "clust.within.transm.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size")
    
    return(ouptuvector.clust)
    
    
  }
  
  
  
  mAr.IDs <- IDs.Seq.Age.Groups(simpact.trans.net = simpact.trans.net,
                                limitTransmEvents = limitTransmEvents,
                                timewindow = timewindow,
                                seq.cov = seq.cov,
                                seq.gender.ratio = seq.gender.ratio,
                                age.group.15.25 = age.group.15.25,
                                age.group.25.40 = age.group.25.40,
                                age.group.40.50 = age.group.40.50)
  
  # nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 
  
  dirfasttree <- dirfasttree
  
  
  
  if(length(mAr.IDs)>5){
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")))
    
    
    mAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                         sub.dir.rename = sub.dir.rename,
                                                         fasttree.tool = "FastTree",
                                                         calendar.dates = "samplingtimes.all.csv",
                                                         simseqfile = paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),
                                                         count.start = 1977,
                                                         endsim = endpoint,
                                                         clust = FALSE)
    
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.tree")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
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
  
  
  outputvector.phylo <-  clust.phylo.fit.params
  
  
  outputvector <- c(outputvector.epi.behav, outputvector.transmission, outputvector.phylo)
  
  
  return(outputvector)
  
}

# 
# w <- MAR.compute.summary.statistics.phylo (datalist = datalist.agemix,
#                                            work.dir = work.dir,
#                                            sub.dir.rename = sub.dir.rename,
#                                            dirfasttree = work.dir,
#                                            limitTransmEvents = 7,
#                                            seq.cov = 35,
#                                            seq.gender.ratio = 0.7,
#                                            age.group.15.25 = c(15,25),
#                                            age.group.25.40 = c(25,40),
#                                            age.group.40.50 = c(40,50),
#                                            endpoint = 40,
#                                            timewindow = c(30,40))
# 
