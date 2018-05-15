


bvar <- function(model) {
  
  # Outputs a df with between-subject variance, upr & lwr limits
  
  # Must take an merMod  object
  
  bsd <- as.numeric(as.data.frame(VarCorr(model))[1,5])
  
  
}


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


# Age mixing in transmission
#############################

agemixing.trans.df <- function(datalist = datalist, 
                               trans.network = trans.network){
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= datalist$itable$population.simtime[1]),]
  
  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
  
  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do
  
  infectionTable <- vector("list", length(seeds.id))
  
  for (i in 1: length(seeds.id)) {
    
    trans.network.i <- as.data.frame(trans.network[[i]])
    
    if(nrow(trans.network.i) >= 3){
      
      trans.network.i <- trans.network.i[-1,]
      
      rrtable <- as.data.frame(cbind(trans.network.i$RecId, trans.network.i$DonId, 
                                     trans.network.i$GenderRec, trans.network.i$GenderDon,
                                     trans.network.i$TOBRec, trans.network.i$TOBDon,
                                     trans.network.i$InfecTime, trans.network.i$SampTime))
      
      names(rrtable) <- c("RecId", "DonId", "GenderRec",
                          "GenderDon", "TOBRec", "TOBDon", "InfecTime", "SampTime")
      
      rrtable.men <- subset(rrtable, rrtable$GenderDon=="0")
      rrtable.women <- subset(rrtable, rrtable$GenderDon=="1")
      
      ids.men <- rrtable.men$DonId
      ids.men.part.w <- rrtable.men$RecId
      age.gap.ID1 <- abs(rrtable.men$TOBDon) - abs(rrtable.men$TOBRec) # men are donors
      tob.men.ID1 <- rrtable.men$TOBDon
      tob.women.ID1 <- rrtable.men$TOBRec
      
      age.men.ID1 <- abs(rrtable.men$TOBDon) + rrtable.men$InfecTime
      age.women.ID1 <- abs(rrtable.men$TOBRec) + rrtable.men$InfecTime
      
      infectime.m <- rrtable.men$InfecTime
      samptime.m <- rrtable.men$SampTime
      
      ID1.m <- ids.men
      ID2.m <- ids.men.part.w
      age.gap.m <- age.gap.ID1
      infectime.m <- infectime.m 
      
      infectable.m <- cbind(ID1.m, ID2.m, tob.men.ID1, tob.women.ID1, age.men.ID1, age.women.ID1, age.gap.m, infectime.m, samptime.m)
      
      ids.women <- rrtable.women$DonId
      ids.women.part.m <- rrtable.women$RecId
      age.gap.ID2 <- abs(rrtable.women$TOBRec) - abs(rrtable.women$TOBDon) # men are receptors
      tob.men.ID2 <- rrtable.women$TOBRec
      tob.women.ID2 <- rrtable.women$TOBDon
      
      age.men.ID2 <- abs(rrtable.women$TOBDon) + rrtable.women$InfecTime
      age.women.ID2 <- abs(rrtable.women$TOBRec) + rrtable.women$InfecTime
      
      infectime.w <- rrtable.women$InfecTime
      samptime.w <- rrtable.women$SampTime
      
      ID1.w <- ids.women.part.m
      ID2.w <- ids.women
      age.gap.w <- age.gap.ID2
      infectime.w <- infectime.w
      
      infectable.w <- cbind(ID1.w, ID2.w, tob.men.ID2, tob.women.ID2, age.men.ID2, age.women.ID2, age.gap.w, infectime.w, samptime.w)
      
      infectable.i <- as.data.frame(rbind(infectable.m, infectable.w))
      
      names(infectable.i) <- c("ID1", "ID2", "TOBID1", "TOBID2", "AgeID1", "AgeID2", "AgeGap", "infecttime", "samptime")
      infectionTable[[i]] <- infectable.i
    }
    
    
  }
  
  
  infecttable <- rbindlist(infectionTable) 
  
  return(infecttable)
  
}


# Fit age mixing in transmission
################################

fit.agemix.trans <- function(datatable = agemix.df){
  
  datatable <- datatable
  
  agemix.inter <- lmer(AgeID2 ~ AgeID1 + (1|ID1), data = datatable)
  
  return(agemix.inter)
  
}


# Onward transmissions
######################


onwardtransmissions.dat <- function(datalist = datalist, 
                                    trans.network = trans.network){
  
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= datalist$itable$population.simtime[1]),]
  
  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)
  
  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do
  
  
  # Onward transmissions in each transmission network
  
  onwardtransm <- vector("list", length(seeds.id))
  
  for (j in 1: length(seeds.id)) {
    
    trans.network.j <- as.data.frame(trans.network[[j]])
    
    trans.network.j <- trans.network.j[-1,] # remove the universal infector
    
    if(nrow(trans.network.j) > 1){ # consider transmission networks with at least one onward transmission
      
      d.j <- table(trans.network.j$DonId) # in the transmission table, the number of times DonId appears is the number of Onward transmissions after acuiring the infection 
      num.j <- as.data.frame(as.numeric(d.j))
      names(num.j) <- c("TransCount")
      onwardtransm[[j]] <- num.j
      
    }
    
  }
  
  onwardtransmissions <- rbindlist(onwardtransm) 
  
  count.dat <- onwardtransmissions$TransCount
  
  return(count.dat) # count.dat = all infections - seeds which didn;t produce at least one transmission
  
}



# compute phylogenetic features





phylogenetic.features.fun <- function(tree.topo=tree,
                                      tree.calib.LTT = tree.calib.LTT){
  
  ########################################
  ### FEATURES FROM PHYLOGENETIC TREE ####
  ########################################
  
  
  # 1.3. Features from phylogenetic tree:
  
  # library(phytools)
  
  # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
  
  tree.cal <- read.tree(paste0(tree.topo))
  
  # tree <- read.tree("calibrated.tree.save.nwk")
  
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
  
  # Estimating confidence intervals for rates and dates using a parametric bootstrap
  pb <- parboot.treedater(tree.calib.LTT) # Lineage Through Time
  
  # Lineages through time
  LTT <- plot.parboot.ltt.dat(pb)
  
  lb.mean.feature <- mean(LTT$lb) # mean of low values of LTT
  lb.median.feature <- median(LTT$lb) # median of low values of LTT
  
  ub.mean.feature <- mean(LTT$ub) # mean of upper values of LTT
  ub.median.feature <- median(LTT$ub) # median of upper values of LTT
  
  median.mean.feature <- mean(LTT$median) # mean of medians of values of LTT
  median.median.feature <- median(LTT$median) # median of medians of values of LTT
  
  phylo.faetures.summary <- c(mean.feature, colless.feature, sackin.feature, mean.tipsDepths.feature, mean.nodesDepths.feature,
                              maxHeight.feature, lb.mean.feature, lb.median.feature, ub.mean.feature, ub.median.feature,
                              median.mean.feature, median.median.feature)
  
  features.names <- c("meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", "mean.nodesDepths.feature",
                      "maxHeight.feature", "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
                      "LTT.median.mean.feature", "LTT.median.median.feature")
  
  names(phylo.faetures.summary) <- features.names
  
  return(phylo.faetures.summary)
  
}

