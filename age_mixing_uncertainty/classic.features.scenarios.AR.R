



classic.features.scenarios.AR <- function(datalist = datalist.agemix,
                                          limitTransmEvents = limitTransmEvents,
                                          work.dir = work.dir,
                                          sub.dir.rename = sub.dir.rename,
                                          time.window = c(10,30),
                                          seq.cov = 70,
                                          seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                          age.group.15.25 = c(15,25),
                                          age.group.25.40 = c(25,40),
                                          age.group.40.50 = c(40,50)){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  work.dir <- paste0(work.dir)
  
  
  datalist.agemix <- datalist
  
  # datalist.agemix <- get(load("datalist.agemix.RData"))
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  #################################### Features on condition of transmission network #############################
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = limitTransmEvents)
  
  #  Selected individuals under missingness scenario
  
  IDs.study <- IDs.Seq.Age.Groups(simpact.trans.net = simpact.trans.net,
                                  limitTransmEvents = limitTransmEvents,
                                  timewindow = time.window,
                                  seq.cov = seq.cov,
                                  seq.gender.ratio = seq.gender.ratio,
                                  age.group.15.25 = age.group.15.25,
                                  age.group.25.40 = age.group.25.40,
                                  age.group.40.50 = age.group.40.50)
  
  
  
  agemixing.df.IDs <- dplyr::filter(agemixing.df, agemixing.df$id.lab%in%IDs.study)
  
  Rec.Id <- agemixing.df.IDs$RecId
  
  person.table <- datalist.agemix
  
  # data list of selected individuals under missingness scenario
  
  person.table$ptable <- dplyr::filter(datalist.agemix$ptable, datalist.agemix$ptable$ID%in%Rec.Id)
  
  
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
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
    
    
    
    # (v) Age mixing in relationships
    
    # 

    agemix.rels.df <- agemix.df.maker(person.table)
    
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
    
    if(nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
      
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
      
      mix.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
      
    }else{
      
      mix.dat <- rep(NA, 6)
      
    }
    
    # age.scatter.df <- agemix.model[[1]]
    
    #  (iii) Point 	prevalence of concurrency in the adult population:
    
    # Concurrency point prevalence 6 months before a survey, among men
    
    # pp.cp.6months.male <- concurr.pointprev.calculator(datalist = datalist.agemix,
    #                                                    timepoint = datalist.agemix$itable$population.simtime[1] - 0.5)
    # 
    
    
    
    pp.cp.6months.male <- tryCatch(concurr.pointprev.calculator(datalist = datalist.agemix,
                                                                timepoint = 40 - 0.5), error=function(e) return(NA))
    
    
    
    summary.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                    hiv.prev.25.40.women, hiv.prev.25.40.men,
                    hiv.prev.40.50.women, hiv.prev.40.50.men, 
                    
                    mix.dat,
                    
                    pp.cp.6months.male)
    
    # name.cov.vector <- paste0(seq(from=1, to=length(cov.vector)),paste0(".cov.vector")) 
    
    features.names <- c("hiv.prev.lt25.women", "hiv.prev.lt25.men", 
                        "hiv.prev.25.40.women", "hiv.prev.25.40.men",
                        "hiv.prev.40.50.women", "hiv.prev.40.50.men", 
                        
                        "AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        
                        "pp.cp.6months.male")
    
    names(summary.df) <- features.names # > length(features.names) [1] 16
    
    
  }else{
    
    features.names <- c("hiv.prev.lt25.women", "hiv.prev.lt25.men", 
                        "hiv.prev.25.40.women", "hiv.prev.25.40.men",
                        "hiv.prev.40.50.women", "hiv.prev.40.50.men", 
                        
                        "AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male",
                        
                        "pp.cp.6months.male")
    
    summary.df <- rep(NA, length(features.names)) # NA -> 0
    
    names(summary.df) <- features.names # > length(features.names) [1] 16
    
  }
  
  
  return(summary.df)
  
}

# 
# x <- classic.features.scenarios.AR(datalist = datalist.agemix,
#                                    limitTransmEvents = 7,
#                                    work.dir = work.dir,
#                                    sub.dir.rename = sub.dir.rename,
#                                    time.window = c(30,40),
#                                    seq.cov = 70,
#                                    seq.gender.ratio = 0.3, 
#                                    age.group.15.25 = c(15,25),
#                                    age.group.25.40 = c(25,40),
#                                    age.group.40.50 = c(40,50))
                                   