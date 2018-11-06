
mixed.effect.fit.transmission.clusters <- function(clust.names=clust.names,
                                                   simpact.trans.net = simpact.trans.net,
                                                   limitTransmEvents = limitTransmEvents,
                                                   select.IDs = select.IDs){
  
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
  
  
  table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)
  
  
  Study.DataTable <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%select.IDs) 
  
  
  IDs.study <- Study.DataTable$RecId
  
  
  transm.datalist.agemix <- datalist.agemix # assign full data set new age mix data set
  
  # Transmission table of selected individuals
  table.simpact.trans.net.cov <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%select.IDs)
  
  # Person table of selected individuals
  transm.datalist.agemix$ptable <- dplyr::filter(transm.datalist.agemix$ptable, transm.datalist.agemix$ptable$ID%in%IDs.study)
  
  
  
  dirfasttree <- work.dir
  
  
  
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
  
  if(nrow(clust.table.df) >= length(d)){
    
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
    
    mean.Num.Clus <- mean(clust.size)
    med.Num.Clus <- median(clust.size)
    SD.Num.Clus <- sd(clust.size)
    
    ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size, mean.Num.Clus, med.Num.Clus, SD.Num.Clus)
    
    names(ouptuvector.clust) <- c("ME.cl.av.age.male", "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var", "ME.cl.SD.female", "ME.cl.SD.male", 
                                  "Num.Clusters", "av.Clust.Size", "mean.Cl.size", "med.Cl.size", "SD.Cl.size")
    
  }else{
    
    ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
    
    clust.stat.names <-  c("ME.cl.av.age.male", "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var", "ME.cl.SD.female", "ME.cl.SD.male", 
                           "Num.Clusters", "av.Clust.Size", "mean.Cl.size", "med.Cl.size", "SD.Cl.size")
    
    ouptuvector.clust <- rep(NA, length(clust.stat.names))
    
    names(ouptuvector.clust) <- clust.stat.names
    
    
  }
  
  return(ouptuvector.clust)
  
  
}