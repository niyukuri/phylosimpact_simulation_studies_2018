#' A function that returns age mixing patterns quantities in transmission clusters
#' in scenarios where individuals are missing at completly at random
#' @param simpact.trans.net a list of transmission networks produced by \code{\link{transm.network.builder}}
#' @param work.dir working directory
#' @param dirfasttree directory where is the fastTree tool
#' @param sub.dir.rename subdurectory required when we have to run more than one simulations
#' @param limitTransmEvents Number of minimum transmission events to be considered in each transmission networks
#' @param timewindow Time interval
#' @param seq.cov Percentage of individulas considered for this transmission pattern scenario
#' @param age.group.15.25 age group between 15 and 25 years old
#' @param age.group.25.40 age group between 25 and 40 years old
#' @param age.group.40.50 age group between 40 and 50 years old
#' @return a vector of number of men and women in different age group, number of transmissions within all age groups, and mean and SD of age different between infectors and infectees
#' @examples
#' w <- MLphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#'                                work.dir = work.dir,
#'                                dirfasttree = dirfasttree,
#'                                sub.dir.rename = sub.dir.rename,
#'                                limitTransmEvents = 7,
#'                                timewindow = c(30,40),
#'                                seq.cov = 70,
#'                                age.group.15.25 = c(15,25),
#'                                age.group.25.40 = c(25,40),
#'                                age.group.40.50 = c(40,50))

#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export


# infection with phylo - sampling time 

LMEMphylo.CAR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
                                            work.dir = work.dir,
                                            dirfasttree = dirfasttree,
                                            sub.dir.rename = sub.dir.rename,
                                            limitTransmEvents = 7,
                                            timewindow = c(30,40),
                                            seq.cov = 70,
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
  
  
  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow,
                             seq.cov = seq.cov,
                             age.limit = age.group.40.50[2])
  
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
                                                          endsim = 40,
                                                          clust = FALSE)
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    d <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                    all.files = FALSE,
                    full.names = FALSE, recursive = FALSE)
    
    
    clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=d, 
                                                               simpact.trans.net = simpact.trans.net,
                                                               limitTransmEvents = 7)
    
    
  }
  
  return(clust.fit.params)
  
}


