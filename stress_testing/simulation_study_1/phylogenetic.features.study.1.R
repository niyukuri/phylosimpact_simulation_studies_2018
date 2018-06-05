



phylogenetic.features.study1 <- function(tree.topo=tree,
                                         tree.calib.LTT = tree.calib.LTT,
                                         work.dir = work.dir,
                                         sub.dir.rename = sub.dir.rename)
{
  
  ########################################
  ### FEATURES FROM PHYLOGENETIC TREE ####
  ########################################
  
  
  # 1. Vector of internal nodes
  #############################
  
  
  N <- node.age(tree.calib.LTT)
  
  int.node.age <- N$Ti # internal nodes ages
  
  latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
  
  
  dt.node.age.dt <- int.node.age
  
  min.val = 1985
  max.val = 2017
  
  step.int=1
  d <- (max.val-min.val)/step.int
  
  int.node.vec <- vector()
  years <- vector()
  
  for (i in 1:d) {
    inf <- 1985+i*step.int
    sup <- 1987+i*step.int
    int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
    int.node.vec <- c(int.node.vec,length(int.node.age.i))
    years <- c(years, inf)
  }
  
  
  poisson.reg.df <- as.data.frame(cbind(years, int.node.vec))
  
  # plot(years, int.node.vec, xlab = "Calendar time", ylab = "Branching events", pch = 16)
  
  count.poisson.fit <- glm(int.node.vec ~ years, data = poisson.reg.df,  family=poisson())
  
  count.coeff <-  count.poisson.fit[[1]][[1]]  # Ok
  time.coeff <-  count.poisson.fit[[1]][[2]]  # Ok
  
  
  # 2. Mean, median, and Sd of weithed age-structured of number of individuals in transmission clusters
  #####################################################################################################
  
  # run ClusterPicker
  
  system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/C.Epidemic.fas"), paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta.nwk"),  paste0("0.9 0.9 0.045 2 gap"))))
  
  # Read clusters' files
  
  d <- list.files(path = paste0(sub.dir.rename), pattern = "C.Epidemic_C.Epidemic_seed.seq.bis.simta_clusterPicks_cluster", 
                  all.files = FALSE,
                  full.names = FALSE, recursive = FALSE)
  
  
  # age groups: < 25 years, 25 - 40 years, 40 - 50 years
  
  age.group.25 <- 25
  age.group.25.40 <- c(25,40)
  age.group.40.50 <- c(40,70)
  
  
  # transmission table
  
  transm.df <- agemixing.trans.df(trans.network = simpact.trans.net,
                                  limitTransmEvents = 7)
  
  stat.clust <- list() # initialise age-structured weithed number of female/male in each transission cluster
  
  count.clust <- vector() # initialise vector of size of transmission clusters
  
  for (i in 1:length(d)) {
    
    clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
    #size <- c(size, nrow(clus.read))
    
    transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of a cluster
    
    transm.df.cl.men <- dplyr::filter(transm.df.cl, transm.df.cl$GenderRec==0) # transmission data table of men IDs of a cluster
    
    transm.df.cl.women <- dplyr::filter(transm.df.cl, transm.df.cl$GenderRec==1) # transmission data table of women IDs of a cluster
    
    # vector of sizes of clusters
    count.clust <- c(count.clust, nrow(transm.df.cl)) 
    
    # age-structured data table
    age.group.25.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$AgeInfecRec < age.group.25)
    age.group.25.40.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$AgeInfecRec >= age.group.25.40[1] & transm.df.cl.men$AgeInfecRec < age.group.25.40[2])
    age.group.40.50.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$AgeInfecRec >= age.group.40.50[1] & transm.df.cl.men$AgeInfecRec < age.group.40.50[2])
    
    age.group.25.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$AgeInfecRec < age.group.25)
    age.group.25.40.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$AgeInfecRec >= age.group.25.40[1] & transm.df.cl.women$AgeInfecRec < age.group.25.40[2])
    age.group.40.50.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$AgeInfecRec >= age.group.40.50[1] & transm.df.cl.women$AgeInfecRec < age.group.40.50[2])
    
    stat.clust.i <- list()
    
    men.25 <- nrow(age.group.25.df.men)
    men.25.40 <- nrow(age.group.25.40.df.men)
    men.40.50 <- nrow(age.group.40.50.df.men)
    
    women.25 <- nrow(age.group.25.df.women)
    women.25.40 <- nrow(age.group.25.40.df.women)
    women.40.50 <- nrow(age.group.40.50.df.women)
    
    pop.count <- nrow(transm.df.cl)
    
    
    stat.clust.i$men.25 <- men.25*pop.count
    stat.clust.i$men.25.40 <- men.25.40*pop.count
    stat.clust.i$men.40.50 <- men.40.50*pop.count
    
    stat.clust.i$women.25 <- women.25*pop.count
    stat.clust.i$women.25.40 <- women.25.40*pop.count
    stat.clust.i$women.40.50 <- women.40.50*pop.count
    
    stat.clust.i$pop.clust <- pop.count
    
    stat.clust[[i]] <- stat.clust.i
    
  }
  
  
  stat.clust.list <- list() # initialise a list age-structured average weithed number of female/male for all transission cluster
  
  stat.clust.list.j <- list()# initialise a list age-structured average weithed number of female/male in each transission cluster
  
  
  for(j in 1:length(stat.clust)){
    
    stat.clust.j <- stat.clust[[j]]
    
    av.men.25 <- stat.clust.j$men.25/sum(count.clust) # sum(count.clust) is \sum_i n_i
    av.men.25.40 <- stat.clust.j$men.25.40/sum(count.clust)
    av.men.40.50 <- stat.clust.j$men.40.50/sum(count.clust)
    
    av.women.25 <- stat.clust.j$women.25/sum(count.clust)
    av.women.25.40 <- stat.clust.j$women.25.40/sum(count.clust)
    av.women.40.50 <- stat.clust.j$women.40.50/sum(count.clust)
    
    stat.clust.list.j$av.men.25 <- av.men.25
    stat.clust.list.j$av.men.25.40 <- av.men.25.40
    stat.clust.list.j$av.men.40.50 <- av.men.40.50
    stat.clust.list.j$av.women.25 <- av.women.25
    stat.clust.list.j$av.women.25.40 <- av.women.25.40
    stat.clust.list.j$av.women.40.50 <- av.women.40.50
    
    stat.clust.list[[j]] <- stat.clust.list.j
    
  }
  
  
  clust.stat.table <- rbindlist(stat.clust.list)
  
  # sum.clust.stat <- sapply(clust.stat.table, sum)
  mean.clust.stat <- sapply(clust.stat.table, mean)  # Ok
  median.clust.stat <- colMedians(as.matrix(clust.stat.table))  # Ok # library(robustbase)
  sd.clust.stat <- sapply(clust.stat.table, sd)  # Ok
  
  
  
  # 3. Mean, median, and SD of size of transmission custers
  ###########################################################
  
  # count.clust
  mean.count.clust <- mean(count.clust) # Ok
  median.count.clust <- median(count.clust)  # Ok
  sd.count.clust <- sd(count.clust)  # Ok
  
  
  # 4. Colless index, Sackin's index, mean & median for mean.nodesDepths, mean branch length
  ###########################################################################################
  
  # library(phytools)
  
  # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
  
  tree.cal <- tree.topo # read.tree(paste0(tree.topo))
  
  
  # library(phyloTop)
  
  colless.feature <- colless.phylo(tree.cal, normalise = TRUE) # Ok
  
  sackin.feature <- sackin.phylo(tree.cal, normalise = TRUE) # Ok
  
  
  # Mean height of internal nodes
  
  H <- nodeHeights(tree.cal) # similar to node.depth.edgelength(tree)
  
  # It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
  # is represented twice and only twice. Thus, if we exclude the root node (zero height), 
  # we can just take the mean of H[,1].
  
  mean.height.internal.nodes <- mean(sort(H[,1])[3:nrow(H)]) # Ok
  
  median.height.internal.nodes <- median(sort(H[,1])[3:nrow(H)]) # Ok
  
  Depths <- getDepths(tree.cal) # depth of tips and nodes
  
  # mean.tipsDepths.feature <- mean(Depths$tipDepths)
  
  mean.nodesDepths.feature <- mean(Depths$nodeDepths) # Ok
  
  median.nodesDepths.feature <- median(Depths$nodeDepths) # Ok
  
  maxHeight.feature <- maxHeight(tree.cal, normalise = TRUE) # Ok
  
  
  phylo.features.summary <- c(count.coeff, time.coeff,
                              mean.clust.stat, median.clust.stat, sd.clust.stat,
                              mean.count.clust, median.count.clust, sd.count.clust,
                              colless.feature, sackin.feature, mean.height.internal.nodes,
                              median.height.internal.nodes, mean.nodesDepths.feature, 
                              median.nodesDepths.feature, maxHeight.feature)
  
  
  name.mean.clust.stat <- names(mean.clust.stat)
  name.median.clust.stat <- names(median.clust.stat)
  name.sd.clust.stat <- names(sd.clust.stat)
  
  features.names <- c("count.coeff", "time.coeff",
                      name.mean.clust.stat, name.median.clust.stat,name.sd.clust.stat,
                      "mean.count.clust", "median.count.clust", "sd.count.clust",
                      "colless.feature", "sackin.feature", "mean.height.internal.nodes",
                      "median.height.internal.nodes", "mean.nodesDepths.feature", 
                      "median.nodesDepths.feature", "maxHeight.feature")
  
  names(phylo.features.summary) <- features.names
  
  return(phylo.features.summary)
  
}


