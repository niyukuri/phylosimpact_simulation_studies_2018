



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
  
  
  
  
  # Sexual behaviour
  
  
  #  Point 	prevalence of concurrency in the adult population
  
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
  
  
  
  # (ii) Relationship per person per year ??
  
  relsperpersonperyear <- nrow(datalist.agemix$rtable) / (nrow(datalist.agemix$ptable)/2) / (timewindow[2] - timewindow[1])
  
  # (iv) SD age gap between couples
  
  agegap.mean <- mean(datalist.agemix$rtable$AgeGap)
  
  agegap.med <- median(datalist.agemix$rtable$AgeGap)
  
  agegap.sd <- sd(datalist.agemix$rtable$AgeGap)
  
  
  classic.features <-   c(growthrate, 
                          
                          hiv.prev.lt25.women, hiv.prev.lt25.men, 
                          hiv.prev.25.40.women, hiv.prev.25.40.men,
                          hiv.prev.40.50.women, hiv.prev.40.50.men,
                          
                          pp.cp.6months.male.rels, pp.cp.6months.female.rels,
                          
                          relsperpersonperyear, 
                          agegap.mean, agegap.med, agegap.sd)
  
  
  classic.features.names <- c("Pop.growthrate", 
                              
                              "hiv.prev.15.25.women", "hiv.prev.15.25.men", 
                              "hiv.prev.25.40.women", "hiv.prev.25.40.men",
                              "hiv.prev.40.50.women", "hiv.prev.40.50.men",
                              
                              "relsperpersonperyear", 
                              "agegap.mean", "agegap.med", "agegap.sd")
  
  names(classic.features) <- classic.features.names 
  
  
  return(classic.features)
  
}


