
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
  
  
  # Check how many indiv in a cluster
  
  # for (i in 1:length(d)) {
  #   
  #   clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  #   
  #   t.v <- c(t.v, nrow(clus.read))
  #   
  # }
  
  for (i in 1:length(d)) {
    
    transm.df.cl.dat <- NULL
    
    clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
    #size <- c(size, nrow(clus.read))
    
    transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
    
    transm.df.cl.dat$age <- transm.df.cl$age.samp.Rec
    transm.df.cl.dat$gender <- transm.df.cl$GenderRec
    transm.df.cl.dat$clust.id <- as.factor(rep(i, nrow(transm.df.cl)))
    
    transmission.clust.list[[i]] <- as.data.frame(transm.df.cl.dat)
    
  }
  
  
  clust.table.df <- as.data.frame(do.call(rbind, transmission.clust.list)) # data.table & data.frame
  
  
  fit.lme.transm.clust <- lme(age ~ gender, data = clust.table.df, random = ~ 1|clust.id)
  
  
  
  a <- coef(summary(fit.lme.transm.clust))[1] # average age in transmission clusters
  
  beta <- coef(summary(fit.lme.transm.clust))[2] # average age difference in transmission clusters: 
  # seen as bridge width which shows potential cross-generation transmission
  
  
  b1 <- as.numeric(VarCorr(fit.lme.transm.clust)[3]) # between cluster variation
  
  b2 <- as.numeric(VarCorr(fit.lme.transm.clust)[4]) # within cluster variation
  
  clust.lme.val <- c(a, beta, b1, b2)
  
  names(clust.lme.val) <- c("av.age.male", "av.age.diff", "between.clust.var", "within.clust.var")
  
  return(clust.lme.val)
  
  
}
