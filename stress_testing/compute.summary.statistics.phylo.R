



compute.summary.statistics.phylo <- function(simpact.trans.net = simpact.trans.net,
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
  
  dirseqgen <- work.dir
  
  seeds.num <- floor(runif(1, min=0, max=1000))
  
  # Sequence simulation is done for at least a transmission network with 6 individuals
  # This means that limitTransmEvents equal at least 7
  
  sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                 sub.dir.rename = sub.dir.rename,
                                 simpact.trans.net = simpact.trans.net,
                                 seq.gen.tool = "seq-gen",
                                 seeds.num = seeds.num,
                                 endpoint = 40,
                                 limitTransmEvents = 7, # no less than 7
                                 hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                 clust = FALSE) # hiv.seq.file lodged in work.dir
  
  # Transform the sequence format to be handled by ClusterPicker
  sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
  write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")
  
  
  
  
  # Function for linear mixed effect models in transmission clusters
  ###################################################################
  
  
  mixed.effect.fit.transmission.clusters <- function(clust.names=clust.names,
                                                     simpact.trans.net = simpact.trans.net,
                                                     limitTransmEvents = limitTransmEvents){
    
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
    
    ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
    
    names(ouptuvector.clust) <- c("clust.av.age.male", "clust.gendEffect", "clust.between.transm.var", "clust.within.transm.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size")
    
    
    }else{
      
      ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
      
      clust.stat.names <- c("clust.av.age.male", "clust.gendEffect", "clust.between.transm.var", "clust.within.transm.var", "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size")
      
      ouptuvector.clust <- rep(NA, length(clust.stat.names))
      
      names(ouptuvector.clust) <- clust.stat.names
      
      
    }
    
    return(ouptuvector.clust)
    
    
  }
  
  
  
  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow,
                             seq.cov = seq.cov,
                             age.limit = age.group.40.50[2])
  
  # nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 
  
  dirfasttree <- dirfasttree
  
  
  
  if(length(mCAr.IDs)>5){
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mCAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))
    
    
    mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                          sub.dir.rename = sub.dir.rename,
                                                          fasttree.tool = "FastTree",
                                                          calendar.dates = "samplingtimes.all.csv",
                                                          simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                          count.start = 1977,
                                                          endsim = endpoint,
                                                          clust = FALSE)
    
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                     all.files = FALSE,
                     full.names = FALSE, recursive = FALSE)
    
    
    clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=dd, 
                                                               simpact.trans.net = simpact.trans.net,
                                                               limitTransmEvents = limitTransmEvents)
    
    
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
    
    names(clust.phylo.fit.params) <- c("clust.av.age.male", "clust.gendEffect", "clust.between.transm.var", "clust.within.transm.var", 
                                       "clust.SD.female", "clust.SD.male", "Num.Clusters", "av.Clust.Size",
                                       "meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", 
                                       "mean.nodesDepths.feature", "maxHeight.feature")
  }
  
  return(clust.phylo.fit.params)
  
}
