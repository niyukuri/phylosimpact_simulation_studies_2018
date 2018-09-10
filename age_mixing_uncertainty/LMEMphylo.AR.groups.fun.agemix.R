#' A function that returns age mixing patterns quantities based on fitting mixed effect models to  transmission clusters
#' in scenarios where individuals are missing at random
#' @param simpact.trans.net a list of transmission networks produced by \code{\link{transm.network.builder}}
#' @param work.dir Working directory
#' @param dirfasttree Directory where FastTree is located
#' @param sub.dir.rename Subdirectory each run when we perform parallelization
#' @param limitTransmEvents Number of minimum transmission events to be considered in each transmission networks
#' @param timewindow Time interval
#' @param seq.cov Percentage of individulas considered for this transmission pattern scenario
#' @param age.group.15.25 age group between 15 and 25 years old
#' @param age.group.25.40 age group between 25 and 40 years old
#' @param age.group.40.50 age group between 40 and 50 years old
#' @return a vector of number of men and women in different age group, number of transmissions within all age groups, and mean and SD of age different between infectors and infectees
#' @examples
#' w <- LMEMphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#'                                      work.dir = work.dir,
#'                                      dirfasttree = dirfasttree,
#'                                      sub.dir.rename = sub.dir.rename,
#'                                      limitTransmEvents = 7,
#'                                      timewindow = c(30,40),
#'                                      seq.cov = 70,
#'                                      age.group.15.25 = c(15,25),
#'                                      age.group.25.40 = c(25,40),
#'                                      age.group.40.50 = c(40,50))

#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
#'


# true we record of infection - infection time


# true we record of infection - infection time

LMEMphylo.AR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
                                           work.dir = work.dir,
                                           dirfasttree = dirfasttree,
                                           sub.dir.rename = sub.dir.rename,
                                           limitTransmEvents = 7,
                                           timewindow = c(30,40),
                                           seq.cov = 70,
                                           seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                           age.group.15.25 = c(15,25),
                                           age.group.25.40 = c(25,40),
                                           age.group.40.50 = c(40,50)){
  
  
  
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
  
  
  ##
  
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
  # transm.df <- agemixing.trans.df(trans.network = simpact.trans.net,
  #                                 limitTransmEvents = 7)
  
  
  # Function for linear mixed effect models
  ##########################################
  
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
    
    
    fit.lme.transm.clust <- lme(age ~ gender, data = clust.table.df, random = ~ 1|clust.id)
    
    
    
    a <- coef(summary(fit.lme.transm.clust))[1] # average age in transmission clusters
    
    beta.va <- coef(summary(fit.lme.transm.clust))[2] # average age difference in transmission clusters: 
    # seen as bridge width which shows potential cross-generation transmission
    
    
    b1 <- as.numeric(VarCorr(fit.lme.transm.clust)[3]) # between cluster variation
    
    b2 <- as.numeric(VarCorr(fit.lme.transm.clust)[4]) # within cluster variation
    
    clust.lme.val <- c(a, beta.va, b1, b2)
    
    
    Num.Clus <- length(d)
    
    av.Clust.size <- sum(clust.size)/length(d)
    
    ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
    
    names(ouptuvector.clust) <- c("av.age.clust.male", "gendEffect.clust", "between.clust.var", "within.clust.var", "Num.Clusters", "av.Clust.Size")
    
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
                                                         endsim = 40,
                                                         clust = FALSE)
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.nwk")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    d <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                    all.files = FALSE,
                    full.names = FALSE, recursive = FALSE)
    

    
    
    clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=d, 
                                                               simpact.trans.net = simpact.trans.net,
                                                               limitTransmEvents = 7)
    
    
  }
  
  return(clust.fit.params)
  
}


