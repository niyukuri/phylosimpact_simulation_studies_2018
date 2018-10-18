#' A function that returns age mixing patterns quantities in transmission networks
#' in scenarios where individuals are missing at completly at random
#' Quantities from fitting mixed effect models to true relationships and transmission records
#' @param simpact.trans.net a list of transmission networks produced by \code{\link{transm.network.builder}}
#' @param datalist a list of data frame produced by simpact
#' @param limitTransmEvents Number of minimum transmission events to be considered in each transmission networks
#' @param timewindow Time interval
#' @param seq.cov Percentage of individulas considered for this transmission pattern scenario
#' @param age.group.15.25 age group between 15 and 25 years old
#' @param age.group.25.40 age group between 25 and 40 years old
#' @param age.group.40.50 age group between 40 and 50 years old
#' @return a vector of number of men and women in different age group, number of transmissions within all age groups, and mean and SD of age different between infectors and infectees
#' @examples
#' w <- CAR.groups.fun.agemixBIS(simpact.trans.net = simpact.trans.net,
#'                            limitTransmEvents = 7,
#'                            timewindow = c(30,40),
#'                            seq.cov = 70,
#'                            age.group.15.25 = c(15,25),
#'                            age.group.25.40 = c(25,40),
#'                            age.group.40.50 = c(40,50))

#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
#'


# true we record of infection - infection time


# true we record of infection - infection time

CAR.groups.fun.agemixBIS <- function(simpact.trans.net = simpact.trans.net,
                                     datalist = datalist.agemix,
                                     limitTransmEvents = 7,
                                     timewindow = c(30,40),
                                     seq.cov = 70,
                                     age.group.15.25 = c(15,25),
                                     age.group.25.40 = c(25,40),
                                     age.group.40.50 = c(40,50)){
  
  
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
  
  
  # (ii) Fiting age mixing transmission table with mixed-effect linear model

  
  # SD for the two strata
  
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
  
  
  
  
  # (iii) Age mixing in relationships
  
  # Function to sort transmissions pairings between different age groups
  
  sort.partners.fun <- function(partner.table = partner.table){ # for receivers
    
    # age and gender structured receiver individuals 
    
    num.15.25.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                              error=function(e) return(NULL))
    
    num.15.25.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                error=function(e) return(NULL))
    
    
    num.25.40.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                              error=function(e) return(NULL))
    
    num.25.40.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                error=function(e) return(NULL))
    
    
    num.40.50.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                              error=function(e) return(NULL))
    
    num.40.50.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                error=function(e) return(NULL))
    
    # consider filter == men
    
    part.men.15.25.women.15.25 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.15.25[1] & num.15.25.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL)) # table of women partners of men between 15 and 25 years old
    
    part.men.15.25.women.25.40 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.25.40[1] & num.15.25.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL)) 
    
    part.men.15.25.women.40.50 <- tryCatch(dplyr::filter(num.15.25.men, num.15.25.men$age.samp.Don >= age.group.40.50[1] & num.15.25.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL)) 
    
    
    part.men.25.40.women.15.25 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.15.25[1] & num.25.40.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL)) 
    
    part.men.25.40.women.25.40 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.25.40[1] & num.25.40.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL)) 
    
    part.men.25.40.women.40.50 <- tryCatch(dplyr::filter(num.25.40.men, num.25.40.men$age.samp.Don >= age.group.40.50[1] & num.25.40.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL)) 
    
    
    part.men.40.50.women.15.25 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.15.25[1] & num.40.50.men$age.samp.Don < age.group.15.25[2]),
                                           error=function(e) return(NULL)) 
    
    part.men.40.50.women.25.40 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.25.40[1] & num.40.50.men$age.samp.Don < age.group.25.40[2]),
                                           error=function(e) return(NULL)) 
    
    part.men.40.50.women.40.50 <- tryCatch(dplyr::filter(num.40.50.men, num.40.50.men$age.samp.Don >= age.group.40.50[1] & num.40.50.men$age.samp.Don < age.group.40.50[2]),
                                           error=function(e) return(NULL)) 
    
    
    N.partners <- c(nrow(part.men.15.25.women.15.25), nrow(part.men.15.25.women.25.40), nrow(part.men.15.25.women.40.50),
                    nrow(part.men.25.40.women.15.25), nrow(part.men.25.40.women.25.40), nrow(part.men.25.40.women.40.50),
                    nrow(part.men.40.50.women.15.25), nrow(part.men.40.50.women.25.40), nrow(part.men.40.50.women.40.50))
    
    return(N.partners)
    
  }
  
  
  trans.sum.age.limit <- data.transm.agemix
  
  # Group 15 - 25
  ###############
  
  trans.sum.men.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.samp.Rec < age.group.15.25[2])
  
  trans.sum.women.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.samp.Rec < age.group.15.25[2])
  
  
  # Group 25 - 40
  ###############
  
  trans.sum.men.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.samp.Rec < age.group.25.40[2])
  
  trans.sum.women.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.samp.Rec < age.group.25.40[2])
  
  
  
  # Group 40 - 50
  ###############
  
  trans.sum.men.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.samp.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.samp.Rec < age.group.40.50[2])
  
  trans.sum.women.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.samp.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.samp.Rec < age.group.40.50[2])
  
  
  
  partners.age.str <- sort.partners.fun(data.transm.agemix) # 154
  
  ouput.transm.dat <- c(nrow(trans.sum.men.15.25), nrow(trans.sum.women.15.25),
                        nrow(trans.sum.men.25.40), nrow(trans.sum.women.25.40),
                        nrow(trans.sum.men.40.50), nrow(trans.sum.women.40.50),
                        
                        partners.age.str)
  
  
  # Age difference statistics #
  #############################
  AD <- abs(abs(data.transm.agemix$TOBDon) - abs(data.transm.agemix$TOBRec))
  mean.AD <- mean(AD)
  med.AD <- median(AD)
  sd.AD <- sd(AD)
  
  
  
  
  ouput.transm.dat.AD <- c(ouput.transm.dat, mean.AD, med.AD, sd.AD,
                           
                           mix.rels.dat,
                           
                           as.numeric(het.lme.val)) 
  
  
  val.names <- c("num.men.15.25", "num.women.15.25",
                 "num.men.25.40", "num.women.25.40",
                 "num.men.40.50", "num.women.40.50",
                 
                 "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                 "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                 "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                 
                 "mean.AD", "median.AD", "sd.AD",
                 
                 "RT.AAD.male", "RT.SDAD.male", "RT.slope.male", "RT.WSD.male", "RT.BSD.male", "RT.intercept.male",
                 
                 "T.het.av.age.male", "T.het.gendEffect", "T.het.between.transm.var", "T.het.within.transm.var", "T.het.SD.female", "T.het.SD.male")
  
  
  names(ouput.transm.dat.AD) <- val.names
  
  
  return(ouput.transm.dat.AD)
}


