

epi.metric.study.1 <- function(datalist.agemix = datalist){
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  agemix.transm.df <- agemixing.trans.df(trans.network = simpact.trans.net, 
                                         limitTransmEvents = 7)
  
  if(!is.null(agemix.transm.df) == TRUE){
    
    
    
    # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
    
    
    ############################ METRICS: TRANSMISSION NETWORK CHARACTERISTICS #####################
    
    # 1. Incidence trend: past five years
    ########################################
    
    # age.group.25 <- 25
    # age.group.25.40 <- c(25,40)
    # age.group.40.50 <- c(40,50)
    
    incidence.df.15.24 <- incidence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25), timewindow = c(35, 40))
    
    METRICS.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
    
    METRICS.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
    METRICS.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]
    
    
    incidence.df.25.39 <- incidence.calculator(datalist = datalist.agemix,
                                               agegroup = c(25, 40), timewindow = c(35, 40))
    
    METRICS.incidence.df.25.39 <- incidence.df.25.39$incidence[3]
    
    METRICS.incidence.df.25.39.men <- incidence.df.25.39$incidence[1]
    METRICS.incidence.df.25.39.women <- incidence.df.25.39$incidence[2]
    
    
    incidence.df.40.49 <- incidence.calculator(datalist = datalist.agemix,
                                               agegroup = c(25, 40), timewindow = c(35, 40))
    
    METRICS.incidence.df.40.49 <- incidence.df.40.49$incidence[3]
    
    METRICS.incidence.df.40.49.men <- incidence.df.40.49$incidence[1] # res
    METRICS.incidence.df.40.49.women <- incidence.df.40.49$incidence[2] # res
    
    
    # 2. Age mixing in transmissions #
    ##################################
    
    
    agemix.fit.men <- fit.agemix.trans.men(datatable = agemix.transm.df)
    
    coef.inter.men <- fixef(agemix.fit.men)
    
    METRICS.age.mix.trans.interc.men <- coef.inter.men[[1]]
    METRICS.age.mix.slope.men <- coef.inter.men[[2]]
    
    
    agemix.fit.women <- fit.agemix.trans.women(datatable = agemix.transm.df)
    
    coef.inter.women <- fixef(agemix.fit.women)
    
    METRICS.age.mix.trans.interc.women <- coef.inter.women[[1]]
    METRICS.age.mix.slope.women <- coef.inter.women[[2]]
    
    # XXXXXXXXX
    age.group.25 <- 25 # < 25
    age.group.25.40 <- c(25,40) # >=25, <40
    age.group.40.50 <- c(40,55)  # >=40, <55
    
    IDs.new.infec <- new.transmissions.dat(datalist = datalist.agemix, 
                                           time.window=c(30,40))
    
    transm.dfdat <- subset(agemix.transm.df, agemix.transm.df$RecId%in%IDs.new.infec) # transmission table of IDs of new infections
    
    transm.dfdat.men <- dplyr::filter(transm.dfdat, transm.dfdat$GenderRec==0) # transmission table of men IDs
    
    transm.dfdat.women <- dplyr::filter(transm.dfdat, transm.dfdat$GenderRec==1) # transmission table of women IDs
    
    
    age.group.25.df.men <- dplyr::filter(transm.dfdat.men, transm.dfdat.men$AgeInfecRec < age.group.25)
    age.group.25.40.df.men <- dplyr::filter(transm.dfdat.men, transm.dfdat.men$AgeInfecRec >= age.group.25.40[1] & transm.dfdat.men$AgeInfecRec < age.group.25.40[2])
    age.group.40.50.df.men <- dplyr::filter(transm.dfdat.men, transm.dfdat.men$AgeInfecRec >= age.group.40.50[1] & transm.dfdat.men$AgeInfecRec < age.group.40.50[2])
    
    # age.group.25.df.women <- dplyr::filter(transm.dfdat.women, transm.dfdat.women$AgeInfecRec < age.group.25)
    # age.group.25.40.df.women <- dplyr::filter(transm.dfdat.women, transm.dfdat.women$AgeInfecRec >= age.group.25.40[1] & transm.dfdat.women$AgeInfecRec < age.group.25.40[2])
    # age.group.40.50.df.women <- dplyr::filter(transm.dfdat.women, transm.dfdat.women$AgeInfecRec >= age.group.40.50[1] & transm.dfdat.women$AgeInfecRec < age.group.40.50[2])
    
    
    # XXXXXXXXXXXXXXX
    
    # 3. Onward transmissions #
    ###########################
    
    transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                            trans.network = simpact.trans.net,
                                            limitTransmEvents = 1,
                                            time.window=c(10,40))
    
    METRICS.transm.mean <- mean(transm.count)
    
    METRICS.transm.median<- median(transm.count)
    
    METRICS.transm.sd <- sd(transm.count)
    
    
    metrics.df <- c(METRICS.incidence.df.15.24, METRICS.incidence.df.25.39, METRICS.incidence.df.40.49,
                    METRICS.age.mix.trans.interc.men, METRICS.age.mix.slope.men,
                    METRICS.age.mix.trans.interc.women, METRICS.age.mix.slope.women,
                    METRICS.transm.mean, METRICS.transm.median, METRICS.transm.sd)
    
    names.metrics.df <- c("METRICS.incidence.df.15", "24METRICS.incidence.df.25.39", "METRICS.incidence.df.40.49",
                          "METRICS.age.mix.trans.interc.men", "METRICS.age.mix.slope.men",
                          "METRICS.age.mix.trans.interc.women", "METRICS.age.mix.slope.women",
                          "METRICS.transm.mean", "METRICS.transm.median", "METRICS.transm.sd")
    
    
    names(metrics.df) <- names.metrics.df
    
  }else{
    
    names.metrics.df <- c("METRICS.incidence.df.15", "24METRICS.incidence.df.25.39", "METRICS.incidence.df.40.49",
                          "METRICS.age.mix.trans.interc.men", "METRICS.age.mix.slope.men",
                          "METRICS.age.mix.trans.interc.women", "METRICS.age.mix.slope.women",
                          "METRICS.transm.mean", "METRICS.transm.median", "METRICS.transm.sd")
    
    metrics.df <- rep(NA, length(names.metrics.df))
    
    names(metrics.df) <- names.metrics.df
  }
  
  return(metrics.df)
  
}


# test.metric <- epi.metric.study.1(datalist.agemix = datalist)

