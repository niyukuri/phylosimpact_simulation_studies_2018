#' A function that returns age mixing patterns quantities in transmission clusters and phylogenetic features
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

MLphylo.CAR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
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
    
    
    # Function to count number of men and women in clusters and possible pairings within each cluster
    
    sort.partners.fun.phylo <- function(partner.table = partner.table){ # for receivers
      
      # age and gender structured receiver individuals
      
      num.15.25.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                error=function(e) return(NULL))
      
      num.15.25.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.15.25[1] & partner.table$age.samp.Rec < age.group.15.25[2]),
                                  error=function(e) return(NULL))
      
      
      num.25.40.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                error=function(e) return(NULL))
      
      num.25.40.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.25.40[1] & partner.table$age.samp.Rec < age.group.25.40[2]),
                                  error=function(e) return(NULL))
      
      
      num.40.50.men <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="0" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                error=function(e) return(NULL))
      
      num.40.50.women <- tryCatch(dplyr::filter(partner.table, partner.table$GenderRec=="1" & partner.table$age.samp.Rec >= age.group.40.50[1] & partner.table$age.samp.Rec < age.group.40.50[2]),
                                  error=function(e) return(NULL))
      
      comb = function(n, x) {
        
        if(n>=x){
          va <- factorial(n) / factorial(n-x) / factorial(x)
        }else{
          va <- factorial(x) / factorial(x-n) / factorial(n)
        }
        return(va)
      }
      
      
      # Possibles pairings
      
      
      pairs.15.25.men.women.15.25 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.15.25.women)), # C(1,n)
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.15.25 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.15.25.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.40.50.men.women.15.25 <- comb(1, nrow(num.40.50.men)) * comb(1, nrow(num.15.25.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.15.25.women)),
      #                                       error=function(e) return(NA))
      #
      
      pairs.15.25.men.women.25.40 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.25.40 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.40.50.men.women.25.40 <- comb(1, nrow(num.40.50.men)) * comb(1, nrow(num.25.40.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.25.40.women)),
      #                                       error=function(e) return(NA))
      #
      
      
      pairs.15.25.men.women.40.50 <- comb(1, nrow(num.15.25.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.15.25.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #
      pairs.25.40.men.women.40.50 <- comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.25.40.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #z
      pairs.40.50.men.women.40.50 <-  comb(1, nrow(num.25.40.men)) * comb(1, nrow(num.40.50.women))
      # tryCatch(comb(nrow(num.40.50.men), nrow(num.40.50.women)),
      #                                       error=function(e) return(NA))
      #
      
      
      pairings.al <- c(pairs.15.25.men.women.15.25, pairs.15.25.men.women.25.40, pairs.15.25.men.women.40.50,
                       pairs.25.40.men.women.15.25, pairs.25.40.men.women.25.40, pairs.25.40.men.women.40.50,
                       pairs.40.50.men.women.15.25, pairs.40.50.men.women.25.40, pairs.40.50.men.women.40.50)
      
      men.women  <- c(nrow(num.15.25.men), nrow(num.15.25.women),
                      nrow(num.25.40.men), nrow(num.25.40.women),
                      nrow(num.40.50.men), nrow(num.40.50.women))
      
      
      val.names <- c("num.men.15.25", "num.women.15.25",
                     "num.men.25.40", "num.women.25.40",
                     "num.men.40.50", "num.women.40.50",
                     
                     "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                     "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                     "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50")
      
      num.ind.pairings.al <- c(men.women, pairings.al)
      
      
      names(num.ind.pairings.al) <- val.names
      
      
      return(num.ind.pairings.al)
      
      
    }
    
    
    if(length(d)>=1){
      
      pairings.clust.tab.list <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster
      
      
      for (i in 1:length(d)) {
        
        clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
        #size <- c(size, nrow(clus.read))
        
        transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
        
        pairings.clust.tab <- sort.partners.fun.phylo(partner.table = transm.df.cl)
        
        
        # Age difference statistics #
        #############################
        
        AD <- abs(abs(transm.df.cl$TOBDon) - abs(transm.df.cl$TOBRec))
        mean.AD <- mean(AD)
        med.AD <- median(AD)
        sd.AD <- sd(AD)

        
        AD.stat <- c(mean.AD, med.AD, sd.AD)
        
        pairings.clust.tab.AD <- c(pairings.clust.tab, AD.stat)
        
        pairings.clust.tab.AD <- as.numeric(pairings.clust.tab.AD)
        
        val.names <- c("num.men.15.25", "num.women.15.25",
                       "num.men.25.40", "num.women.25.40",
                       "num.men.40.50", "num.women.40.50",
                       
                       "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                       "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                       "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                       
                       "mean.AD", "median.AD", "sd.AD")
        
        names(pairings.clust.tab.AD) <- val.names
        
        
        pairings.clust.tab.list[[i]] <- pairings.clust.tab.AD
        
      }
      
      
      clust.stat.table <- as.data.frame(do.call(rbind, pairings.clust.tab.list)) # data.table & data.frame
      
      clust.stat.table.up <- clust.stat.table
      
      
    }else{
      
      val.names <- c("num.men.15.25", "num.women.15.25",
                     "num.men.25.40", "num.women.25.40",
                     "num.men.40.50", "num.women.40.50",
                     
                     "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                     "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                     "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                     
                     "mean.AD", "median.AD", "sd.AD")
      
      clust.stat.table <- rep(NA, length(val.names))
      
      names(clust.stat.table) <- val.names
      
      clust.stat.table.up <- clust.stat.table
      
      
    }
    
    
    # clust.stat.table <- as.data.frame(do.call(rbind, pairings.clust.tab.list)) # data.table & data.frame
    
  }else{
    
    val.names <- c("num.men.15.25", "num.women.15.25",
                   "num.men.25.40", "num.women.25.40",
                   "num.men.40.50", "num.women.40.50",
                   
                   "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                   "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                   "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                   
                   "mean.AD", "median.AD", "sd.AD")
    
    clust.stat.table <- rep(NA, length(val.names))
    
    names(clust.stat.table) <- val.names
    
    clust.stat.table.up <- clust.stat.table
    
  }
  
  
  
  ################## ADD PHYLOGENETIC FEATURES ############
  
  # This is for machine learning purpose
  
  # 2. Mean, median, and Sd of weithed age-structured of number of individuals in transmission clusters
  #####################################################################################################
  
  # # run ClusterPicker
  # 
  # system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0(fasta.file)), paste0(sub.dir.rename,"/", paste0(tree.file)),  paste0("0.9 0.9 0.045 2 gap"))))
  # 
  # # Read clusters' files
  # 
  # d <- list.files(path = paste0(sub.dir.rename), pattern = "C.Epidemic_C.Epidemic_seed.seq.bis.simta_clusterPicks_cluster", 
  #                 all.files = FALSE,
  #                 full.names = FALSE, recursive = FALSE)
  
  
  # age groups: < 25 years, 25 - 40 years, 40 - 50 years
  
  age.group.25 <- age.group.15.25[1]
  age.group.25.40 <- age.group.25.40
  age.group.40.50 <- age.group.40.50
  
  # 
  #   # transmission table
  # 
  #   
  #   # id of people who got infection by seed event: seeds.id
  #   trans.network <- simpact.trans.net
  #   seeds.id <- length(trans.network)
  #   limitTransmEvents <- limitTransmEvents
  #   
  #   ID.select <- vector() # ID of selected transmission network
  #   ID.select.count <- vector() # number of individuals in these networks
  #   
  #   for (i in 1: seeds.id) {
  #     
  #     
  #     trans.network.i <- as.data.frame(trans.network[[i]])
  #     
  #     if(nrow(trans.network.i)>=limitTransmEvents){
  #       
  #       
  #       ID.select <- c(ID.select, i)
  #       ID.select.count <- c(ID.select.count, nrow(trans.network.i))
  #       
  #     } # X if
  #     
  #   } # Y for
  #   
  #   
  #   infectionTable <- vector("list", length(ID.select))
  #   
  #   for(j in 1:length(ID.select)){
  #     
  #     p <- ID.select[j]
  #     
  #     trans.network.i <- as.data.frame(trans.network[[p]])
  #     
  #     trans.network.i <- trans.network.i[-1,]
  #     
  #     trans.network.i$AgeInfecDon <- abs(trans.network.i$TOBDon) + trans.network.i$InfecTime
  #     trans.network.i$AgeInfecRec <- abs(trans.network.i$TOBRec) + trans.network.i$InfecTime
  #     
  #     id.lab <- paste0(p,".",trans.network.i$id,".C")
  #     
  #     trans.network.i$id.lab <- id.lab
  #     
  #     infectionTable[[p]] <- trans.network.i
  #   }
  #   
  #   
  #   infecttable <- rbindlist(infectionTable)
  #   
  #   
  #   
  #   transm.df <- infecttable
  #   # transm.df <- agemixing.trans.df(trans.network = simpact.trans.net,
  #   #                                 limitTransmEvents = 7)
  
  stat.clust <- list() # initialise age-structured weithed number of female/male in each transission cluster
  
  count.clust <- vector() # initialise vector of size of transmission clusters
  
  
  #    # Change AgeInfecRec with age.samp.Rec
  
  for (i in 1:length(d)) {
    
    clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
    #size <- c(size, nrow(clus.read))
    
    transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of a cluster
    
    transm.df.cl.men <- dplyr::filter(transm.df.cl, transm.df.cl$GenderRec==0) # transmission data table of men IDs of a cluster
    
    transm.df.cl.women <- dplyr::filter(transm.df.cl, transm.df.cl$GenderRec==1) # transmission data table of women IDs of a cluster
    
    # vector of sizes of clusters
    count.clust <- c(count.clust, nrow(transm.df.cl)) 
    
    
    # age-structured data table
    age.group.25.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$age.samp.Rec < age.group.25)
    age.group.25.40.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$age.samp.Rec >= age.group.25.40[1] & transm.df.cl.men$age.samp.Rec < age.group.25.40[2])
    age.group.40.50.df.men <- dplyr::filter(transm.df.cl.men, transm.df.cl.men$age.samp.Rec >= age.group.40.50[1] & transm.df.cl.men$age.samp.Rec < age.group.40.50[2])
    
    age.group.25.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$age.samp.Rec < age.group.25)
    age.group.25.40.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$age.samp.Rec >= age.group.25.40[1] & transm.df.cl.women$age.samp.Rec < age.group.25.40[2])
    age.group.40.50.df.women <- dplyr::filter(transm.df.cl.women, transm.df.cl.women$age.samp.Rec >= age.group.40.50[1] & transm.df.cl.women$age.samp.Rec < age.group.40.50[2])
    
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
  
  
  clust.stat.table <- rbindlist(stat.clust.list) # data.table & data.frame
  
  clust.stat.table <- as.data.frame(as.matrix(clust.stat.table)) # make it data.frame ONLY to be handled by colMedians
  
  # sum.clust.stat <- sapply(clust.stat.table, sum)
  mean.clust.stat <- sapply(clust.stat.table, mean)  # Ok
  median.clust.stat <- colMedians(clust.stat.table)  # Ok # library(robustbase)
  sd.clust.stat <- sapply(clust.stat.table, sd)  # Ok
  
  
  
  # 3. Mean, median, and SD of size of transmission clusters
  ###########################################################
  
  # count.clust
  mean.count.clust <- mean(count.clust) # Ok
  median.count.clust <- median(count.clust)  # Ok
  sd.count.clust <- sd(count.clust)  # Ok
  
  
  # 4. Colless index, Sackin's index, mean & median for mean.nodesDepths, mean branch length
  ###########################################################################################
  
  # library(phytools)
  
  # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
  
  tree.cal <- multi2di(tree.cal.cov.35.IDs) # read.tree(paste0(tree.topo))
  
  
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
  
  
  phylo.features.summary <- c(mean.clust.stat, median.clust.stat, sd.clust.stat,
                              mean.count.clust, median.count.clust, sd.count.clust,
                              colless.feature, sackin.feature, mean.height.internal.nodes,
                              median.height.internal.nodes, mean.nodesDepths.feature, 
                              median.nodesDepths.feature, maxHeight.feature)
  
  
  name.mean.clust.stat <- paste0("clust.mean.", names(mean.clust.stat)) # average number of men and women in different age groups
  name.median.clust.stat <- paste0("clust.median.", names(median.clust.stat))
  name.sd.clust.stat <- paste0("clust.SD.", names(sd.clust.stat))
  
  #   name.clust.AR.c.35 <- paste0("clust.AR.c.35.",names(transm.clust.AR.c.cov.35.val))
  
  features.names <- c(name.mean.clust.stat, name.median.clust.stat,name.sd.clust.stat,
                      "mean.Sizes.clust", "median.Sizes.clust", "sd.Sizes.clust",
                      "colless.feature", "sackin.feature", "mean.height.internal.nodes",
                      "median.height.internal.nodes", "mean.nodesDepths.feature", 
                      "median.nodesDepths.feature", "maxHeight.feature")
  
  names(phylo.features.summary) <- features.names
  
  # return(phylo.features.summary)
  # 
  # 
  # return(clust.stat.table) # sapply(transm.clust.AR.c.cov.90, mean)
  
  clust.stat.table.mean <- sapply(clust.stat.table.up, mean)
  
  ALL.summary <- c(clust.stat.table.mean, phylo.features.summary)
  
  
  return(ALL.summary)
  
}

# 
# d <- MLphylo.CAR.groups.fun.agemix(simpact.trans.net = simpact.trans.net,
#                                    work.dir = work.dir,
#                                    dirfasttree = dirfasttree,
#                                    sub.dir.rename = sub.dir.rename,
#                                    limitTransmEvents = 7,
#                                    timewindow = c(30,40),
#                                    seq.cov = 70,
#                                    age.group.15.25 = c(15,25),
#                                    age.group.25.40 = c(25,40),
#                                    age.group.40.50 = c(40,50))





