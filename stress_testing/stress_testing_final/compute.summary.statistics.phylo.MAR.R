

compute.summary.statistics.phylo.MAR <- function(simpact.trans.net = simpact.trans.net,
                                                 datalist.agemix = datalist.agemix,
                                                 work.dir = work.dir,
                                                 sub.dir.rename = sub.dir.rename,
                                                 dirfasttree = work.dir,
                                                 limitTransmEvents = 7,
                                                 seq.cov = 100,
                                                 seq.gender.ratio = 0.7,
                                                 age.group.15.25 = c(15,25),
                                                 age.group.25.40 = c(25,40),
                                                 age.group.40.50 = c(40,50),
                                                 endpoint = 40,
                                                 timewindow = c(30,40),
                                                 cut.off = 7){
  
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/mixed.effect.fit.transmission.clusters.R")
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/stats.age.groups.trans.clust.network.fun.R")
  
  
  work.dir <- paste0(work.dir)
  
  dirseqgen <- work.dir
  
  
  # Function for linear mixed effect models in transmission clusters
  ###################################################################
  
  
  mAr.IDs <- IDs.Seq.Age.Groups(simpact.trans.net = simpact.trans.net,
                                limitTransmEvents = limitTransmEvents, 
                                timewindow = timewindow, 
                                seq.cov = seq.cov,
                                seq.gender.ratio = seq.gender.ratio, 
                                age.group.15.25 = age.group.15.25,
                                age.group.25.40 = age.group.25.40, 
                                age.group.40.50 = age.group.40.50)
  
  # nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 
  
  dirfasttree <- dirfasttree
  
  
  
  if(length(mAr.IDs)>10){
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")))
    
    
    mAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                         sub.dir.rename = sub.dir.rename,
                                                         fasttree.tool = "FastTree",
                                                         calendar.dates = "samplingtimes.all.csv",
                                                         simseqfile = paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"), # paste0("cov.",seq.cov,".", seq.gender.ratio, ".mAr.IDs.C.Epidemic.Fasta")
                                                         count.start = 1977,
                                                         endsim = endpoint,
                                                         clust = FALSE)
    
    N <- node.age(mAr.IDs.tree.calib)
    
    # Time to MRCA: internal nodes ages
    
    int.node.age <- N$Ti
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    # Fitting internal nodes: Poisson regression
    
    # Number of internal nodes aggregated
    
    
    min.node <- round(min(int.node.age))
    max.node <- round(max(int.node.age))
    
    intervals <- max.node - min.node
    
    dt.node.age.dt <- int.node.age
    
    
    i.vec <- vector() # initialise time intervals
    int.node.vec <- vector() # initialize internal nodes
    
    for (i in 1:intervals) {
      inf <- min.node - 1 + i
      sup <- min.node+i
      i.vec <- c(i.vec, sup)
      int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
      int.node.vec <- c(int.node.vec,length(int.node.age.i))
    }
    
    
    node.df <- data.frame(i.vec, int.node.vec)
    
    output <-glm(formula = int.node.vec ~ i.vec, data = node.df, family = poisson)
    
    intercept.fit <- output$coefficients[[1]]
    
    slope.fit <- output$coefficients[[2]]
    
    
    names(intercept.fit) <-     "node.interc.fit"
    names(slope.fit) <- "node.slope.fit"
    
    
    mrca.v <- mrca(mAr.IDs.tree.calib, full = FALSE) # MRCA ids
    
    
    sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv")) # sampling times
    
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.tree")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                     all.files = FALSE,
                     full.names = FALSE, recursive = FALSE)
    
    
    
    
    clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=dd, 
                                                               simpact.trans.net = simpact.trans.net,
                                                               limitTransmEvents = limitTransmEvents,
                                                               select.IDs = mAr.IDs)
    
    
    
    stats.AG <- stats.age.groups.trans.clust.network.fun(simpact.trans.net = simpact.trans.net, 
                                                         datalist.agemix = datalist.agemix,
                                                         work.dir = work.dir,  
                                                         dirfasttree = dirfasttree, 
                                                         sub.dir.rename = sub.dir.rename,
                                                         limitTransmEvents = limitTransmEvents,
                                                         timewindow = timewindow,
                                                         seq.cov = seq.cov,
                                                         select.IDs = mAr.IDs,
                                                         age.group.15.25 = age.group.15.25,
                                                         age.group.25.40 = age.group.25.40,
                                                         age.group.40.50 = age.group.40.50,
                                                         cut.off = cut.off)
    
    
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
    
    features.names <- c("meanHeight", "colless", "sackin", "mean.tipsDepths", "mean.nodesDepths",
                        "maxHeight")
    # , "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
    #                     "LTT.median.mean.feature", "LTT.median.median.feature")
    
    names(phylo.features.summary) <- features.names
    
    
    clust.phylo.fit.params <- c(intercept.fit,
                                slope.fit,
                                
                                clust.fit.params, 
                                phylo.features.summary,
                                stats.AG)
    
    
  }else{
    
    clust.phylo.fit.params <- rep(NA, 37)
    
    names(clust.phylo.fit.params) <- c("node.interc.fit",
                                       "node.slope.fit",
                                       
                                       "ME.cl.av.age.male",  "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var",  "ME.cl.SD.female",          
                                       "ME.cl.SD.male", "Num.Clusters",  "av.Clust.Size",   "mean.Num.Clus",  "med.Num.Clus",   "SD.Num.Clus" , "meanHeight",
                                       "colless",  "sackin", "mean.tipsDepths", "mean.nodesDepths", "maxHeight",  "cl.prop.men15.25.F.15.25",   "cl.prop.men25.40.F.15.25",
                                       "cl.prop.men40.50.F.15.25", "cl.prop.men15.25.F.25.40",   "cl.prop.men25.40.F.25.40",   "cl.prop.men40.50.F.25.40",  
                                       "cl.prop.men15.25.F.40.50" ,  "cl.prop.men25.40.F.40.50" , "cl.prop.men40.50.F.40.50",   "cl.prop.women15.25.M.15.25", 
                                       "cl.prop.women25.40.M.15.25", "cl.prop.women40.50.M.15.25" ,"cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40", 
                                       "cl.prop.women40.50.M.25.40", "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50")
    
    
    
    
  }
  
  return(clust.phylo.fit.params)
  
}
