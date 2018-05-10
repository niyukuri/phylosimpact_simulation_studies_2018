


# Calculate the transmission rate
###################################


transmission.rate.calculator <- function(datalist = datalist,
                                         timewindow = c(20, 40), int = FALSE, by=1){
  
  if(int==FALSE){
    
    Infec.pop.table <- datalist$ptable[InfectType==1]
    
    numb.infec.pop <- nrow(Infec.pop.table %>%
                             subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))
    
    
    transm.rate <- numb.infec.pop / diff(timewindow)
    
    return(transm.rate)
  }
  
  if(int==TRUE){
    
    Infec.pop.table <- datalist$ptable[InfectType==1]
    
    upper.limit <- ceiling(diff(timewindow)/by)
    
    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
    
    infec.pop.int <- vector()
    trans.rate.int <- vector()
    
    
    for(i in 0:(upper.limit-2)){
      
      timewindow.int <- c(interv.time[1+i], interv.time[2+i])
      
      infec.pop.numb <- nrow(Infec.pop.table %>%
                               subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))
      
      
      infec.pop.int <- c(infec.pop.int,  infec.pop.numb)
      
      
      trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))
      
    }
    
    return(trans.rate.int)
  }
  
}


# Relationship rate
#####################


relationship.rate.calculator <- function(datalist = datalist,
                                         timewindow = c(20, 40), int = FALSE, by=1){
  
  if(int==FALSE){
    Rels.table <- datalist$rtable
    
    Rels.table.window <- Rels.table %>%
      subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])
    
    numb.rels <- nrow(Rels.table.window)
    
    numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women
    
    numb.rels.women <- length(unique(Rels.table.window$ID2))
    
    rels.rate <- (numb.rels)/ diff(timewindow)
    
    return(rels.rate)
  }
  
  if(int==TRUE){
    
    Rels.table <- datalist$rtable
    
    upper.limit <- ceiling(diff(timewindow)/by)
    
    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
    
    rels.int <- vector()
    rels.rate.int <- vector()
    
    
    for(i in 0:(upper.limit-2)){
      
      timewindow.int <- c(interv.time[1+i], interv.time[2+i])
      
      rels.numb <- nrow(Rels.table %>%
                          subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))
      
      
      rels.int <- c(rels.int,  rels.numb)
      
      
      rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))
      
    }
    
    return(rels.rate.int)
  }
  
  
}


# age mixing pattern model
###########################

ampmodel <- function(data = dplyr::filter(agemix.model[[1]], Gender =="male")) {
  lmer(pagerelform ~ agerelform0 + (1 | ID),
       data = data,
       REML = TRUE)
}


# Concurrency point prevalence 6 months before a survey, among men
####################################################################


concurr.pointprev.calculator <- function(datalist,
                                         timepoint = datalist$itable$population.simtime[1] - 0.5){
  
  output <- data.table()
  DTalive.infected <- alive.infected(datalist = datalist,
                                     timepoint = timepoint, site = "All") # First we only take the data of people who were alive at time_i
  agemix.df <- agemix.df.maker(datalist)
  degrees.df <- degree.df.maker(dataframe.df = agemix.df,
                                agegroup = c(15, 50),
                                hivstatus = 2,
                                survey.time = timepoint,
                                window.width = 0,
                                gender.degree = "male",
                                only.new = FALSE)
  
  number.people.with.cps <- sum(degrees.df$Degree > 1)
  popsize <- nrow(DTalive.infected)
  concurr.pointprevalence <- number.people.with.cps / popsize
  
  return(concurr.pointprevalence)
}



# ART coverage
###############

ART.coverage.vector.creator <- function(datalist = datalist,
                                        agegroup = c(15, 50)){
  ART.cov.eval.timepoints <- seq(from = datalist$itable$t[2],
                                 to = datalist$itable$population.simtime[1])
  ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  for (art.cov.index in 1:length(ART.cov.vector)){
    ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist,
                                                                 agegroup = agegroup,
                                                                 timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
      sum(ART.coverage.calculator(datalist = datalist,
                                  agegroup = agegroup,
                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  }
  return(ART.cov.vector)
}


# Function to compute lineages through time with confidence intervals
# ######################################################################

plot.parboot.ltt.dat <- function (pbtd, t0 = NA, res = 100, ...)
{
  t1 <- max(pbtd$td$sts, na.rm = T)
  if (is.na(t0))
    t0 <- min(sapply(pbtd$trees, function(tr) tr$timeOf))
  times <- seq(t0, t1, l = res)
  ltt <- cbind(times = times, t(sapply(times, function(t) {
    c(pml = sum(pbtd$td$sts > t) - sum(pbtd$td$Ti > t), setNames(quantile(sapply(pbtd$trees,
                                                                                 function(tre) sum(tre$sts > t) - sum(tre$Ti > t)),
                                                                          probs = c(0.025, 0.5, 0.975)), c("lb", "median",
                                                                                                           "ub")))
  })))
  pl.df <- as.data.frame(ltt)
  return(pl.df)
  # p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
  #                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
  # p <- p + geom_path(aes(x = times, y = pml))
  # (p <- p + ylab("Lineages through time") + xlab("Time"))
}

# LTT <- plot.parboot.ltt.dat(pb)

# names(LTT)
# # [1] "times"  "pml"    "lb"     "median" "ub" 
# p <- ggplot(LTT) + geom_ribbon(aes(x = times, ymin = lb,
#                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
# p <- p + geom_path(aes(x = times, y = pml))
# (p <- p + ylab("Lineages through time") + xlab("Time"))

