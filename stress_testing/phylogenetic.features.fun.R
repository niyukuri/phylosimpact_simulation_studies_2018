



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
