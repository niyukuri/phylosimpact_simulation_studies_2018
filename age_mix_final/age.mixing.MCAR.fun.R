




age.mixing.MCAR.fun <- function(simpact.trans.net = simpact.trans.net.adv, 
                                datalist.agemix = datalist.agemix,
                                work.dir = work.dir,  
                                dirfasttree = dirfasttree, 
                                sub.dir.rename = sub.dir.rename,
                                limitTransmEvents = 7,
                                timewindow = c(30,40),
                                seq.cov = 35,
                                age.group.15.25 = c(15,25),
                                age.group.25.40 = c(25,40),
                                age.group.40.50 = c(40,50),
                                cut.off = 7){
  
  
  
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/needed.functions.RSimpactHelp.R")
  
  #################################################################################
  # Step 5: Epidemic statistics and sexual behaviour: data set of infected people #
  #################################################################################
  
  
  # Data list of infected individuals
  
  # Select IDs in MCAR scenario
  
  simpact.trans.net <- simpact.trans.net
  limitTransmEvents <- limitTransmEvents
  timewindow <- timewindow
  seq.cov <- seq.cov
  age.group.40.50 <- age.group.40.50
  
  mCAr.IDs.100 <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net, # simpact.trans.net 
                                 limitTransmEvents = limitTransmEvents,
                                 timewindow = timewindow, 
                                 seq.cov = 100, 
                                 age.limit = age.group.40.50[2])
  
  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net, # simpact.trans.net 
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow, 
                             seq.cov = seq.cov, 
                             age.limit = age.group.40.50[2])
  
  
  if(length(mCAr.IDs) >= 20){
    
    
    simpact.trans.net.adv <- simpact.trans.net
    
    
    # Transmission network table as from transmission networks for further steps
    ############################################################################
    
    
    infectionTable <- vector("list", length(simpact.trans.net.adv))
    
    for(j in 1:length(simpact.trans.net.adv)){
      
      p <- j
      
      trans.network.i <- as.data.frame(simpact.trans.net.adv[[p]])
      
      # trans.network.i <- trans.network.i[-1,]
      
      id.lab <- paste0(p,".",trans.network.i$id,".C")
      
      trans.network.i$id.lab <- id.lab
      trans.network.i$ageSampTimeRec <- trans.network.i$SampTime - trans.network.i$TOBRec
      
      infectionTable[[p]] <- trans.network.i
      
      
    }
    
    
    infecttable <- rbindlist(infectionTable)
    
    
    table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)
    
    
    Study.DataTable <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs) 
    
    
    IDs.study <- Study.DataTable$RecId
    
    
    transm.datalist.agemix <- datalist.agemix # assign full data set new age mix data set
    
    # Transmission table of selected individuals
    table.simpact.trans.net.cov <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs)
    
    table.simpact.trans.net.cov.100 <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs.100)
    
    # Person table of selected individuals
    transm.datalist.agemix$ptable <- dplyr::filter(transm.datalist.agemix$ptable, transm.datalist.agemix$ptable$ID%in%IDs.study)
    
    
    
    # (i) Age mixing in relationships
    
    # 
    
    agemix.rels.transm.df <- agemix.df.maker(transm.datalist.agemix)
    
    # 
    agemix.model <- pattern.modeller(dataframe = agemix.rels.transm.df,
                                     agegroup = c(15, 50),
                                     timepoint = 40, # transm.datalist.agemix$itable$population.simtime[1],
                                     timewindow = 10)#1)#3)
    # 
    # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
    # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
    #
    # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
    
    data = dplyr::filter(agemix.model[[1]], Gender =="male")
    
    if( nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
      
      men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                       data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                       REML = TRUE,
                       control=lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))
      
      bignumber <- NA # let's try if NA works (instead of 9999 for example)
      AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
      SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
      #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
      slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
      WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
      
      BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
      
      intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
      
      # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
      
      ## AAD: average age difference across all relationship
      ## VAD: variance of these age differences
      ## SDAD: standard deviation of age differences
      ## BSD: between-subject standard deviation of age differences
      
      mix.rels.transm.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
      
      names(mix.rels.transm.dat) <-  c("T.AAD.male", "T.SDAD.male", "T.slope.male", "T.WSD.male", "T.BSD.male", "T.intercept.male")
      
    }else{
      
      mix.rels.transm.dat <- rep(NA, 6)
      
      names(mix.rels.transm.dat) <-  c("T.AAD.male", "T.SDAD.male", "T.slope.male", "T.WSD.male", "T.BSD.male", "T.intercept.male")
      
    }
    
    # age.scatter.df <- agemix.model[[1]]
    
    #  (ii) Point 	prevalence of concurrency in the adult population:
    
    # Concurrency point prevalence 6 months before a survey, among men
    
    
    pp.cp.6months.male.transm <- tryCatch(concurr.pointprev.calculator(datalist = transm.datalist.agemix,
                                                                       timepoint = 40 - 0.5), error=function(e) return(NA))
    
    
    # (iii) Prevalence
    
    hiv.prev.lt25.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                                 agegroup = c(15, 25),
                                                 timepoint = 40)$pointprevalence[2]
    hiv.prev.lt25.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = 40)$pointprevalence[1]
    
    hiv.prev.25.40.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                                  agegroup = c(25, 40),
                                                  timepoint = 40)$pointprevalence[2]
    hiv.prev.25.40.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                                agegroup = c(25, 40),
                                                timepoint = 40)$pointprevalence[1]
    
    hiv.prev.40.50.women <- prevalence.calculator(datalist = transm.datalist.agemix,
                                                  agegroup = c(40, 50),
                                                  timepoint = 40)$pointprevalence[2]
    hiv.prev.40.50.men <- prevalence.calculator(datalist = transm.datalist.agemix,
                                                agegroup = c(40, 50),
                                                timepoint = 40)$pointprevalence[1]
    
    
    # (iv) Incidence
    
    incidence.df.15.24 <- incidence.calculator(datalist = transm.datalist.agemix,
                                               agegroup = c(15, 25), timewindow = c(30, 40))
    
    epi.transm.incidence.df.15.24 <- incidence.df.15.24$incidence[3]
    
    epi.transm.incidence.df.15.24.men <- incidence.df.15.24$incidence[1]
    epi.transm.incidence.df.15.24.women <- incidence.df.15.24$incidence[2]
    
    
    incidence.df.25.39 <- incidence.calculator(datalist = transm.datalist.agemix,
                                               agegroup = c(25, 40), timewindow = c(30, 40))
    
    epi.transm.incidence.df.25.39 <- incidence.df.25.39$incidence[3]
    
    epi.transm.incidence.df.25.39.men <- incidence.df.25.39$incidence[1]
    epi.transm.incidence.df.25.39.women <- incidence.df.25.39$incidence[2]
    
    
    incidence.df.40.49 <- incidence.calculator(datalist = transm.datalist.agemix,
                                               agegroup = c(25, 40), timewindow = c(30, 40))
    
    epi.transm.incidence.df.40.49 <- incidence.df.40.49$incidence[3]
    
    epi.transm.incidence.df.40.49.men <- incidence.df.40.49$incidence[1] # res
    epi.transm.incidence.df.40.49.women <- incidence.df.40.49$incidence[2] # res
    
    
    
    summary.epidemic.transm.df <- c(hiv.prev.lt25.women, hiv.prev.lt25.men, 
                                    hiv.prev.25.40.women, hiv.prev.25.40.men,
                                    hiv.prev.40.50.women, hiv.prev.40.50.men, 
                                    mix.rels.transm.dat,
                                    pp.cp.6months.male.transm,
                                    
                                    epi.transm.incidence.df.15.24.men, epi.transm.incidence.df.15.24.women, 
                                    epi.transm.incidence.df.25.39.men, epi.transm.incidence.df.25.39.women,
                                    epi.transm.incidence.df.40.49.men, epi.transm.incidence.df.40.49.women)
    
    
    names(summary.epidemic.transm.df) <- c("T.prev.15.25.w", "T.prev.15.25.m", "T.prev.25.40.w", "T.prev.25.40.m", "T.prev.40.50.w", "T.prev.40.50.m",
                                           names(mix.rels.transm.dat), 
                                           "T.p.prev.6months.m",
                                           "T.inc.15.25.w", "T.inc.15.25.m", "T.inc.25.40.w", "T.inc.25.40.m", "T.inc.40.50.w", "T.inc.40.50.m")
    
    
    ######################################
    # Step 5: Building phylogenetic tree #
    ######################################
    
    
    
    dirfasttree <- work.dir
    
    
    
    # Select sequences from the pool of alignment
    ##############################################
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mCAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))
    
    
    # Build and calibrate the phylogenetic tree
    ############################################
    
    mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                          sub.dir.rename = sub.dir.rename,
                                                          fasttree.tool = "FastTree",
                                                          calendar.dates = "samplingtimes.all.csv",
                                                          simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                          count.start = 1977,
                                                          endsim = 40,
                                                          clust = FALSE)
    
    
    
    N <- node.age(mCAr.IDs.tree.calib)
    
    # Time to MRCA: internal nodes ages
    
    int.node.age <- N$Ti
    
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    mrca.v <- mrca(mCAr.IDs.tree.calib, full = FALSE) # MRCA ids
    
    
    sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv")) # sampling times
    
    # 
    # tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))
    # 
    
    
    # Compute transmission clusters
    ###############################
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                     all.files = FALSE,
                     full.names = FALSE, recursive = FALSE)
    
    # Transmission clusters.
    
    d <- clust.names <- dd
    
    data.list.simpact.trans.net.adv <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster
    
    
    # Transmission table of individuals in the transmission clusters
    #################################################################
    
    
    # Binding all data tables of clusters as these information are captured in transmission networks
    
    clust.size <- vector() # size of each cluster # table.simpact.trans.net.adv
    
    transm.df <- table.simpact.trans.net.adv
    
    for (i in 1:length(d)) {
      
      transm.df.cl.dat <- NULL
      
      clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
      
      clust.size <- c(clust.size, nrow(clus.read))
      
      data.table.simpact.trans.net.i <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
      
      data.table.simpact.trans.net.i$clust.ID <- rep(i, nrow(data.table.simpact.trans.net.i))
      
      data.list.simpact.trans.net.adv[[i]] <- as.data.frame(data.table.simpact.trans.net.i)
      
    }
    
    
    data.table.simpact.trans.clusts.net.adv <- as.data.frame(do.call(rbind, data.list.simpact.trans.net.adv)) # data.table & data.frame
    
    data.table.simpact.trans.clusts.net.adv <- data.table.simpact.trans.clusts.net.adv[!duplicated(data.table.simpact.trans.clusts.net.adv[c("id","id.lab")]),] # remove duplicate id.lab
    # t may happen that one seq.ID appear in more than one cluster
    
    # data.table.simpact.trans.clusts.net.adv <- data.table.simpact.trans.net.adv
    
    
    ## Aligning internal nodes IDs and their age: !they must get same length
    
    ancestor <- Ancestors(mCAr.IDs.tree.calib) # ancestors of each tips and internal node
    # All ancestors output are internal nodes
    
    ancestor.v <- vector()
    
    for(i in 1:length(ancestor)){
      
      k <- ancestor[[i]]
      
      ancestor.v <- c(ancestor.v, unique(k))
      
    }
    
    sort.int.ancestor <- unique(sort(ancestor.v))
    sort.int.node.age <- sort(int.node.age)
    
    
    tip.names <- names(mrca.v[1,])
    
    dates.tree.df <- dplyr::filter(sampling.dates, sampling.dates$V1%in%tip.names) # dates of these tips
    
    # rearrange dates in tips order as are displayed on the tree
    tip.names.f <- vector()
    dates.tree.dat <- vector()
    for(i in 1:nrow(dates.tree.df)){
      for(j in 1:length(tip.names)){
        if(tip.names[i] == dates.tree.df$V1[[j]]){
          tip.names.f <- c(tip.names.f, tip.names[i])
          dates.tree.dat <- c(dates.tree.dat, 1977+40-dates.tree.df$V2[[j]])
        }
      }
    }
    
    
    dates.tree.named <- dates.tree.dat
    names(dates.tree.named) <- tip.names.f
    
    
    # MRCA matrix
    #############
    
    # make mrca matrix diagonal 0 and other elements (internal nodes IDs) assign them the age of mrca
    
    mrca.v.age <- mrca.v
    
    for(i in 1:nrow(mrca.v.age)){
      for(j in 1:nrow(mrca.v.age)){
        
        if(i==j){
          mrca.v.age[i,j] <- 0
        }else{
          
          if(mrca.v[i,j] %in% sort.int.ancestor){
            
            p.index <- which(sort.int.ancestor == mrca.v[i,j])
            
            mrca.v.age[i,j] <-  sort.int.node.age[p.index]
          }
          
        }
        
      }
    }
    
    
    
    # make mrca matrix elements: sampling date - age of mrca
    
    # Fist contingency matrix
    
    mrca.v.age.samp <- mrca.v.age
    
    mrca.v.age.samp.cont1 <- mrca.v.age.samp
    
    for(i in 1:nrow(mrca.v.age)){
      
      for(j in 1:nrow(mrca.v.age)){
        
        if(i!=j){
          
          i.dat <- tip.names.f[i]
          
          v.index <- which(tip.names.f == i.dat)
          
          samp.date.tip <- dates.tree.dat[v.index]
          
          mrca.v.age.samp.cont1[i,] <- samp.date.tip - mrca.v.age.samp[i,]
          
        }
        
      }
    }
    
    
    
    # Second contingency matrix
    
    mrca.v.age.samp <- mrca.v.age
    
    mrca.v.age.samp.cont2 <- mrca.v.age.samp
    
    for(i in 1:nrow(mrca.v.age)){
      
      for(j in 1:nrow(mrca.v.age)){
        
        if(i!=j){
          
          i.dat <- tip.names.f[i]
          
          v.index <- which(tip.names.f == i.dat)
          
          samp.date.tip <- dates.tree.dat[v.index]
          
          mrca.v.age.samp.cont2[,i] <- samp.date.tip - mrca.v.age.samp[,i]
          
        }
        
      }
    }
    
    
    # Diagonal zero for mrca.v.age.samp.cont1 and mrca.v.age.samp.cont2
    
    for(i in 1:nrow(mrca.v.age.samp.cont1)){
      for(j in 1:nrow(mrca.v.age.samp.cont1)){
        
        if(i==j){
          mrca.v.age.samp.cont1[i,j] <- 0
        }
      }
    }
    
    
    for(i in 1:nrow(mrca.v.age.samp.cont2)){
      for(j in 1:nrow(mrca.v.age.samp.cont2)){
        
        if(i==j){
          mrca.v.age.samp.cont2[i,j] <- 0
        }
      }
    }
    
    
    # filter table.simpact.trans.net.adv and remain with table of tips names (individulas in the tree)
    
    attributes.table.simpact.trans.net.adv <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%tip.names)
    
    V.gender <- vector()
    V.cd4 <- vector()
    V.vl <- vector()
    V.x <- vector()
    V.y <- vector()
    iD <- vector()
    
    for(i in 1:length(tip.names)){
      for(j in 1:nrow(attributes.table.simpact.trans.net.adv)){
        if(tip.names[i] == attributes.table.simpact.trans.net.adv$id.lab[j]){
          
          V.gender <- c(V.gender, attributes.table.simpact.trans.net.adv$GenderRec[j])
          V.cd4 <- c(V.cd4, attributes.table.simpact.trans.net.adv$cd4[j])
          V.vl <- c(V.vl, attributes.table.simpact.trans.net.adv$vl[j])
          V.x <- c(V.x, attributes.table.simpact.trans.net.adv$location.x[j])
          V.y <- c(V.y, attributes.table.simpact.trans.net.adv$location.y[j])
          iD <- c(iD, tip.names[i])
        }
        
      }
    }
    
    
    Node.gender.cd4.vl.x.y <- data.table(V.gender,V.cd4, V.vl, V.x, V.y, iD)
    
    
    # Adding clusters ID on the previous attributes table from attributes.table.simpact.trans.net.adv
    
    clust.ID.vec <- vector()
    id.vec <- vector()
    
    for(k in 1:nrow(Node.gender.cd4.vl.x.y)){ # attributes table ofr all tips on the tree: Node.gender.cd4.vl.x.y
      
      id <- Node.gender.cd4.vl.x.y$iD[k]
      
      
      if(id%in%data.table.simpact.trans.clusts.net.adv$id.lab){    # ID of tree which belongs to IDs of clusters
        # transmission table of individuls in the transmission clusters: data.table.simpact.trans.net.adv
        
        id.index <- which(data.table.simpact.trans.clusts.net.adv$id.lab == id)
        
        clust.ID.vec.i <- data.table.simpact.trans.clusts.net.adv$clust.ID[id.index]
        
      }else{
        
        clust.ID.vec.i <-  0 # tip ID which is not in any transmission cluster is assigned value 0
        
      }
      
      clust.ID.vec <-  c(clust.ID.vec, clust.ID.vec.i)
      id.vec <- c(id.vec, id)
      
    }
    
    Node.gender.cd4.vl.x.y$clust.ID <- clust.ID.vec
    
    Node.gender.cd4.vl.x.y.clusID <- Node.gender.cd4.vl.x.y # attributes table with clusters' IDs
    
    
    
    ## Building transmission network
    
    # 1. consider contigency matrix 2
    
    mrca.times.final <- as.matrix(abs(mrca.v.age.samp.cont2))
    
    
    net <- graph.adjacency(as.matrix(mrca.times.final), mode="undirected",weighted=T,diag=FALSE)
    
    # E(net)       # The edges of the "net" object
    # 
    # V(net)       # The vertices of the "net" object
    
    V(net)$gender <- Node.gender.cd4.vl.x.y$V.gender
    V(net)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
    V(net)$vl <- Node.gender.cd4.vl.x.y$V.vl
    V(net)$loc.x <- Node.gender.cd4.vl.x.y$V.x
    V(net)$loc.y <- Node.gender.cd4.vl.x.y$V.y
    
    
    
    
    ## Filtering the network by breaking some edges due to conditions from individuals attributes:
    ##############################################################################################
    
    # 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load
    
    # Now considering 1 and 2
    
    names.attributes.ngaha <- Node.gender.cd4.vl.x.y
    
    names.matrix.contigency <- names(mrca.times.final[1,])
    
    gender.l <- names.attributes.ngaha$V.gender
    
    
    clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID
    
    
    mrca.times.filter <- mrca.times.final
    
    
    # 
    # for (i in 1:length(names(mrca.times.final[1,]))) {
    #   
    #   name.col.i <- names.matrix.contigency[i]
    #   
    #   index.i <- which(names(mrca.times.final[1,]) == name.col.i)
    #   
    #   gender.i <- gender.l[index.i]
    #   
    #   cluster.i <- clusters.zose[index.i]
    #   
    #   
    #   for(j in 1:length(names(mrca.times.final[1,]))){
    #     
    #     if(i != j){
    #       
    #       name.col.j <- names.matrix.contigency[j]
    #       
    #       index.j <- which(names(mrca.times.final[1,]) == name.col.j)
    #       
    #       gender.j <- gender.l[index.j]
    #       
    #       cluster.j <- clusters.zose[index.j]
    #       
    #       
    #       if(gender.i == gender.j){ # if same gender break the link
    #         
    #         mrca.times.filter[i,j] <- 0
    #         
    #       }
    #       
    #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
    #         
    #         mrca.times.filter[i,j] <- 0
    #         
    #       }
    #       
    #       
    #     }
    #     
    #   }
    #   
    # }
    
    # i. Gender 
    ############
    
    for (i in 1:length(names(mrca.times.final[1,]))) {
      
      name.col.i <- names.matrix.contigency[i]
      
      index.i <- which(names(mrca.times.final[1,]) == name.col.i)
      
      gender.i <- gender.l[index.i]
      
      for(j in 1:length(names(mrca.times.final[1,]))){
        
        if(i != j){
          
          name.col.j <- names.matrix.contigency[j]
          
          index.j <- which(names(mrca.times.final[1,]) == name.col.j)
          
          gender.j <- gender.l[index.j]
          
          if(gender.i == gender.j){ # if same gender break the link
            
            mrca.times.filter[i,j] <- 0
            
          }
          
        }
        
      }
      
    }
    
    mrca.times.filter.gender <- mrca.times.filter
    
    
    # ii. Cluster
    #############
    
    mrca.times.filter.gender.clust <- mrca.times.filter.gender
    
    for (i in 1:length(names(mrca.times.final[1,]))) {
      
      name.col.i <- names.matrix.contigency[i]
      
      index.i <- which(names(mrca.times.final[1,]) == name.col.i)
      
      cluster.i <- clusters.zose[index.i]
      
      
      for(j in 1:length(names(mrca.times.final[1,]))){
        
        if(i != j){
          
          name.col.j <- names.matrix.contigency[j]
          
          index.j <- which(names(mrca.times.final[1,]) == name.col.j)
          
          cluster.j <- clusters.zose[index.j]
          
          
          if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
            
            mrca.times.filter.gender.clust[i,j] <- 0
            
          }
          
          
        }
        
      }
      
    }
    
    
    # iii. tMRCA
    #############
    
    
    
    net.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.gender.clust),mode="undirected",weighted=T,diag=FALSE)
    
    
    # Consider plausible transmissions and difference between sampling time and tMRCA
    
    
    cut.off <- cut.off # years
    
    # E(net.cont.1)$weight
    
    net.cont.1 <- delete_edges(net.cont.1, E(net.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff
    
    # E(net.cont.1)$weight
    
    # plot(net.cont.1, layout=layout_with_kk) 
    
    
    
    # Delete tips of the phylogenetic tree which are not part of transmission clusters:  they have clust.ID==0 >> deletes vertices 
    ###################################################################################
    
    
    Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
    Non.ids <- Non.ids.dat$iD
    
    net.cleaned <- delete_vertices(net.cont.1, Non.ids)
    
    
    
    # 
    # # 2. consider contigency matrix 1
    # 
    # mrca.times.final.2 <- as.matrix(abs(mrca.v.age.samp.cont1))
    # 
    # 
    # net.2 <- graph.adjacency(as.matrix(mrca.times.final.2), mode="undirected",weighted=T,diag=FALSE)
    # 
    # E(net.2)       # The edges of the "net.2" object
    # 
    # V(net.2)       # The vertices of the "net.2" object
    # 
    # V(net.2)$gender <- Node.gender.cd4.vl.x.y$V.gender
    # V(net.2)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
    # V(net.2)$vl <- Node.gender.cd4.vl.x.y$V.vl
    # V(net.2)$loc.x <- Node.gender.cd4.vl.x.y$V.x
    # V(net.2)$loc.y <- Node.gender.cd4.vl.x.y$V.y
    # 
    # 
    # 
    # 
    # ## Filtering the net.2work by breaking some edges due to conditions from individuals attributes:
    # 
    # # 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load
    # 
    # 
    # names.attributes.ngaha <- Node.gender.cd4.vl.x.y
    # 
    # names.matrix.contigency <- names(mrca.times.final.2[1,])
    # 
    # gender.l <- names.attributes.ngaha$V.gender
    # 
    # 
    # clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID
    # 
    # 
    # mrca.times.filter.2 <- mrca.times.final.2
    # 
    # 
    # # 
    # # for (i in 1:length(names(mrca.times.final.2[1,]))) {
    # #   
    # #   name.col.i <- names.matrix.contigency[i]
    # #   
    # #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
    # #   
    # #   gender.i <- gender.l[index.i]
    # #   
    # #   cluster.i <- clusters.zose[index.i]
    # #   
    # #   
    # #   for(j in 1:length(names(mrca.times.final.2[1,]))){
    # #     
    # #     if(i != j){
    # #       
    # #       name.col.j <- names.matrix.contigency[j]
    # #       
    # #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
    # #       
    # #       gender.j <- gender.l[index.j]
    # #       
    # #       cluster.j <- clusters.zose[index.j]
    # #       
    # #       
    # #       if(gender.i == gender.j){ # if same gender break the link
    # #         
    # #         mrca.times.filter.2[i,j] <- 0
    # #         
    # #       }
    # #       
    # #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
    # #         
    # #         mrca.times.filter.2[i,j] <- 0
    # #         
    # #       }
    # #       
    # #       
    # #     }
    # #     
    # #   }
    # #   
    # # }
    # 
    # # i. Gender 
    # 
    # for (i in 1:length(names(mrca.times.final.2[1,]))) {
    #   
    #   name.col.i <- names.matrix.contigency[i]
    #   
    #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
    #   
    #   gender.i <- gender.l[index.i]
    #   
    #   for(j in 1:length(names(mrca.times.final.2[1,]))){
    #     
    #     if(i != j){
    #       
    #       name.col.j <- names.matrix.contigency[j]
    #       
    #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
    #       
    #       gender.j <- gender.l[index.j]
    #       
    #       if(gender.i == gender.j){ # if same gender break the link
    #         
    #         mrca.times.filter.2[i,j] <- 0
    #         
    #       }
    #       
    #     }
    #     
    #   }
    #   
    # }
    # 
    # mrca.times.filter.2.gender <- mrca.times.filter.2
    # 
    # 
    # # ii. Cluster
    # 
    # mrca.times.filter.2.gender.clust <- mrca.times.filter.2.gender
    # 
    # for (i in 1:length(names(mrca.times.final.2[1,]))) {
    #   
    #   name.col.i <- names.matrix.contigency[i]
    #   
    #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
    #   
    #   cluster.i <- clusters.zose[index.i]
    #   
    #   
    #   for(j in 1:length(names(mrca.times.final.2[1,]))){
    #     
    #     if(i != j){
    #       
    #       name.col.j <- names.matrix.contigency[j]
    #       
    #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
    #       
    #       cluster.j <- clusters.zose[index.j]
    #       
    #       
    #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
    #         
    #         mrca.times.filter.2.gender.clust[i,j] <- 0
    #         
    #       }
    #       
    #       
    #     }
    #     
    #   }
    #   
    # }
    # 
    # 
    # 
    # net.2.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.2.gender.clust),mode="undirected",weighted=T,diag=FALSE)
    # 
    # 
    # # Consider plausible transmissions and difference between sampling time and tMRCA
    # 
    # 
    # cut.off <- 20
    # 
    # E(net.2.cont.1)$weight
    # 
    # net.2.cont.1 <- delete_edges(net.2.cont.1, E(net.2.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff
    # 
    # E(net.2.cont.1)$weight
    # 
    # plot(net.2.cont.1, layout=layout_with_kk) 
    # 
    # 
    # # Delete tips which are not part of transmission clusters, they have clust.ID==0 >> deletes vertices 
    # 
    # Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
    # Non.ids <- Non.ids.dat$iD
    # 
    # net.2.cleaned <- delete_vertices(net.2.cont.1, Non.ids)
    
    
    
    # r=graph.union(net.cleaned, net.2.cleaned)
    
    
    
    
    # Age structure in the transmission network built from phylogenetic tree
    #########################################################################
    
    
    # produce age table
    
    net.sp <- net.cleaned
    
    
    transm.matrix <- as.data.table(get.edgelist(net.sp)) # matrix of links of the transmission network built from phylogenetic tree
    
    # table.simpact.trans.net.adv
    
    # reduced transmission table: table.simpact.trans.net.adv of ids in transmission clusters
    
    ids <-  unique(c(transm.matrix$V1, transm.matrix$V2))
    
    
    table.transm.clust.net.igraph <- dplyr::filter(data.table.simpact.trans.clusts.net.adv, data.table.simpact.trans.clusts.net.adv$id.lab%in%ids) 
    
    
    
    # 1.
    
    # True age structure in transmission clusters as observed from phylogenetic tree #
    ##################################################################################
    
    
    age.groups.filtered.trans.clust.network.fun <- function(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                            transm.matrix = transm.matrix,
                                                            age.group.15.25 = c(15,25),
                                                            age.group.25.40 = c(25,40),
                                                            age.group.40.50 = c(40,50)){
      
      Age.groups.table <- NULL
      
      v1.dat <- vector()
      v2.dat <- vector()
      age1.dat <- vector()
      age2.dat <- vector()
      gender1.dat <- vector()
      gender2.dat <- vector()
      
      for(i in 1:nrow(transm.matrix)){
        
        v1 <- transm.matrix$V1[i]
        v2 <- transm.matrix$V2[i]
        
        index.v1 <- which(table.transm.clust.net.igraph$id.lab == v1)
        index.v2 <- which(table.transm.clust.net.igraph$id.lab == v2)
        
        age1 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v1]
        age2 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v2]
        
        gender1 <- table.transm.clust.net.igraph$GenderRec[index.v1]
        gender2 <- table.transm.clust.net.igraph$GenderRec[index.v2]
        
        v1.dat <- c(v1.dat, v1)
        v2.dat <- c(v2.dat, v2)
        age1.dat <- c(age1.dat, age1)
        age2.dat <- c(age2.dat, age2)
        gender1.dat <- c(gender1.dat, gender1)
        gender2.dat <- c(gender2.dat, gender2)
        
      }
      
      age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
      
      
      
      # men
      men.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==0)
      
      # women
      women.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==1)
      
      
      # men 15.25 and women
      
      men.15.25.women.15.25.1 <- vector()
      men.15.25.women.25.40.1 <- vector()
      men.15.25.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1)>1){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      # women 15.25 and men
      
      women.15.25.men.15.25.2 <- vector()
      women.15.25.men.25.40.2 <- vector()
      women.15.25.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1)>1){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
      }
      
      
      
      
      
      
      # men 25.40 and women
      
      men.25.40.women.15.25.1 <- vector()
      men.25.40.women.25.40.1 <- vector()
      men.25.40.women.40.50.1 <- vector()
      
      
      if(nrow(men.age.table.1) >1 ){
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
      }
      
      
      
      
      
      
      # women 25.40 and men
      
      women.25.40.men.15.25.2 <- vector()
      women.25.40.men.25.40.2 <- vector()
      women.25.40.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
      }
      
      
      
      
      
      # men 40.50 and women
      
      men.40.50.women.15.25.1 <- vector()
      men.40.50.women.25.40.1 <- vector()
      men.40.50.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) >1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      
      # women 40.50 and men
      
      women.40.50.men.15.25.2 <- vector()
      women.40.50.men.25.40.2 <- vector()
      women.40.50.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
      }
      
      
      men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
      
      men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
      
      men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
      
      men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
      
      men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
      
      men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
      
      men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
      
      men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
      
      men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
      
      Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
                                   length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
                                   length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
                                 ncol = 3,
                                 byrow = TRUE)
      
      colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      
      Age.groups.table <- as.table(Age.groups.table)
      
      
      men.15.25.T <- sum(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50))
      men.25.40.T <- sum(length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50))
      men.40.50.T <- sum(length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50))
      
      prop.men.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/men.15.25.T, length(men.15.25.women.25.40)/men.15.25.T, length(men.15.25.women.40.50)/men.15.25.T,
                                            length(men.25.40.women.15.25)/men.25.40.T, length(men.25.40.women.25.40)/men.25.40.T, length(men.25.40.women.40.50)/men.25.40.T,
                                            length(men.40.50.women.15.25)/men.40.50.T, length(men.40.50.women.25.40)/men.40.50.T, length(men.40.50.women.40.50)/men.40.50.T),
                                          ncol = 3,
                                          byrow = TRUE)
      
      colnames(prop.men.age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(prop.men.age.groups.table) <- c("prop.Male.15.25", "prop.Male.25.40", "prop.Male.40.50")
      
      
      
      
      women.15.25.T <- sum(length(men.15.25.women.15.25), length(men.25.40.women.15.25), length(men.40.50.women.15.25))
      women.25.40.T <- sum(length(men.15.25.women.25.40), length(men.25.40.women.25.40), length(men.40.50.women.25.40))
      women.40.50.T <- sum(length(men.15.25.women.40.50), length(men.25.40.women.40.50), length(men.40.50.women.40.50))
      
      prop.women.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/women.15.25.T, length(men.25.40.women.15.25)/women.15.25.T, length(men.40.50.women.15.25)/women.15.25.T,
                                              length(men.15.25.women.25.40)/women.25.40.T, length(men.25.40.women.25.40)/women.25.40.T, length(men.40.50.women.25.40)/women.25.40.T,
                                              length(men.15.25.women.40.50)/women.40.50.T, length(men.25.40.women.40.50)/women.40.50.T, length(men.40.50.women.40.50)/women.40.50.T),
                                            ncol = 3,
                                            byrow = TRUE)
      
      colnames(prop.women.age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      rownames(prop.women.age.groups.table) <- c("prop.Female.15.25", "prop.Female.25.40", "prop.Female.40.50")
      
      outputlist <- NULL
      outputlist$Age.groups.table <- Age.groups.table
      outputlist$prop.men.age.groups.table <- prop.men.age.groups.table
      outputlist$prop.women.age.groups.table <- prop.women.age.groups.table
      
      
      return(outputlist)
      
    }
    
    
    
    
    
    # 2.
    
    # True age structure in transmission clusters as observed from transmission network #
    #####################################################################################
    
    age.groups.filtered.transmission.clust.fun <- function(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                           age.group.15.25 = c(15,25),
                                                           age.group.25.40 = c(25,40),
                                                           age.group.40.50 = c(40,50)){
      
      table.transm.clust.net.igraph$ageSampTimeDon <- table.transm.clust.net.igraph$SampTime - table.transm.clust.net.igraph$TOBDon
      
      Age.groups.table <- NULL
      
      v1.dat <- vector()
      v2.dat <- vector()
      age1.dat <- vector()
      age2.dat <- vector()
      gender1.dat <- vector()
      gender2.dat <- vector()
      
      for(i in 1:nrow(table.transm.clust.net.igraph)){
        
        v1 <- table.transm.clust.net.igraph$RecId[i]
        v2 <- table.transm.clust.net.igraph$DonId[i]
        
        index.v1 <- which(table.transm.clust.net.igraph$RecId == v1)
        
        age1 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v1]
        age2 <- table.transm.clust.net.igraph$ageSampTimeDon[index.v1]
        
        gender1 <- table.transm.clust.net.igraph$GenderRec[index.v1]
        gender2 <- table.transm.clust.net.igraph$GenderDon[index.v1]
        
        v1.dat <- c(v1.dat, v1)
        v2.dat <- c(v2.dat, v2)
        age1.dat <- c(age1.dat, age1)
        age2.dat <- c(age2.dat, age2)
        gender1.dat <- c(gender1.dat, gender1)
        gender2.dat <- c(gender2.dat, gender2)
        
      }
      
      age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
      
      
      
      age.group.15.25 = c(15,25)
      age.group.25.40 = c(25,40)
      age.group.40.50 = c(40,50)
      
      # men
      men.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==0)
      
      # women
      women.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==1)
      
      
      # men 15.25 and women
      
      men.15.25.women.15.25.1 <- vector()
      men.15.25.women.25.40.1 <- vector()
      men.15.25.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) >1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      # women 15.25 and men
      
      women.15.25.men.15.25.2 <- vector()
      women.15.25.men.25.40.2 <- vector()
      women.15.25.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1 ){
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        } 
      }
      
      
      
      # men 25.40 and women
      
      men.25.40.women.15.25.1 <- vector()
      men.25.40.women.25.40.1 <- vector()
      men.25.40.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) >1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }  
      }
      
      
      
      
      
      # women 25.40 and men
      
      women.25.40.men.15.25.2 <- vector()
      women.25.40.men.25.40.2 <- vector()
      women.25.40.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1 ){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
        }
      }
      
      
      
      
      # men 40.50 and women
      
      men.40.50.women.15.25.1 <- vector()
      men.40.50.women.25.40.1 <- vector()
      men.40.50.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) > 1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      
      # women 40.50 and men
      
      women.40.50.men.15.25.2 <- vector()
      women.40.50.men.25.40.2 <- vector()
      women.40.50.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1 ){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      
      men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
      
      men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
      
      men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
      
      men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
      
      men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
      
      men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
      
      men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
      
      men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
      
      men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
      
      Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
                                   length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
                                   length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
                                 ncol = 3,
                                 byrow = TRUE)
      
      colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      
      Age.groups.table <- as.table(Age.groups.table)
      
      
      men.15.25.T <- sum(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50))
      men.25.40.T <- sum(length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50))
      men.40.50.T <- sum(length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50))
      
      prop.men.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/men.15.25.T, length(men.15.25.women.25.40)/men.15.25.T, length(men.15.25.women.40.50)/men.15.25.T,
                                            length(men.25.40.women.15.25)/men.25.40.T, length(men.25.40.women.25.40)/men.25.40.T, length(men.25.40.women.40.50)/men.25.40.T,
                                            length(men.40.50.women.15.25)/men.40.50.T, length(men.40.50.women.25.40)/men.40.50.T, length(men.40.50.women.40.50)/men.40.50.T),
                                          ncol = 3,
                                          byrow = TRUE)
      
      colnames(prop.men.age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(prop.men.age.groups.table) <- c("prop.Male.15.25", "prop.Male.25.40", "prop.Male.40.50")
      
      
      
      
      women.15.25.T <- sum(length(men.15.25.women.15.25), length(men.25.40.women.15.25), length(men.40.50.women.15.25))
      women.25.40.T <- sum(length(men.15.25.women.25.40), length(men.25.40.women.25.40), length(men.40.50.women.25.40))
      women.40.50.T <- sum(length(men.15.25.women.40.50), length(men.25.40.women.40.50), length(men.40.50.women.40.50))
      
      prop.women.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/women.15.25.T, length(men.25.40.women.15.25)/women.15.25.T, length(men.40.50.women.15.25)/women.15.25.T,
                                              length(men.15.25.women.25.40)/women.25.40.T, length(men.25.40.women.25.40)/women.25.40.T, length(men.40.50.women.25.40)/women.25.40.T,
                                              length(men.15.25.women.40.50)/women.40.50.T, length(men.25.40.women.40.50)/women.40.50.T, length(men.40.50.women.40.50)/women.40.50.T),
                                            ncol = 3,
                                            byrow = TRUE)
      
      colnames(prop.women.age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      rownames(prop.women.age.groups.table) <- c("prop.Female.15.25", "prop.Female.25.40", "prop.Female.40.50")
      
      outputlist <- NULL
      outputlist$Age.groups.table <- Age.groups.table
      outputlist$prop.men.age.groups.table <- prop.men.age.groups.table
      outputlist$prop.women.age.groups.table <- prop.women.age.groups.table
      
      
      return(outputlist)
      
      
    }
    
    
    
    
    # 3.
    
    # True age structure in transmission transmission network for selected individuals #
    #####################################################################################
    
    age.groups.filtered.transmission.net.fun <- function(table.transmission.net.cov = table.simpact.trans.net.cov,
                                                         age.group.15.25 = c(15,25),
                                                         age.group.25.40 = c(25,40),
                                                         age.group.40.50 = c(40,50)){
      
      table.transmission.net.cov$ageSampTimeDon <- table.transmission.net.cov$SampTime - table.transmission.net.cov$TOBDon
      
      Age.groups.table <- NULL
      
      v1.dat <- vector()
      v2.dat <- vector()
      age1.dat <- vector()
      age2.dat <- vector()
      gender1.dat <- vector()
      gender2.dat <- vector()
      
      for(i in 1:nrow(table.transmission.net.cov)){
        
        v1 <- table.transmission.net.cov$RecId[i]
        v2 <- table.transmission.net.cov$DonId[i]
        
        index.v1 <- which(table.transmission.net.cov$RecId == v1)
        
        age1 <- table.transmission.net.cov$ageSampTimeRec[index.v1]
        age2 <- table.transmission.net.cov$ageSampTimeDon[index.v1]
        
        gender1 <- table.transmission.net.cov$GenderRec[index.v1]
        gender2 <- table.transmission.net.cov$GenderDon[index.v1]
        
        v1.dat <- c(v1.dat, v1)
        v2.dat <- c(v2.dat, v2)
        age1.dat <- c(age1.dat, age1)
        age2.dat <- c(age2.dat, age2)
        gender1.dat <- c(gender1.dat, gender1)
        gender2.dat <- c(gender2.dat, gender2)
        
      }
      
      age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
      
      
      
      age.group.15.25 = c(15,25)
      age.group.25.40 = c(25,40)
      age.group.40.50 = c(40,50)
      
      # men
      men.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==0)
      
      # women
      women.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==1)
      
      
      # men 15.25 and women
      
      men.15.25.women.15.25.1 <- vector()
      men.15.25.women.25.40.1 <- vector()
      men.15.25.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) >1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      # women 15.25 and men
      
      women.15.25.men.15.25.2 <- vector()
      women.15.25.men.25.40.2 <- vector()
      women.15.25.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1 ){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      # men 25.40 and women
      
      men.25.40.women.15.25.1 <- vector()
      men.25.40.women.25.40.1 <- vector()
      men.25.40.women.40.50.1 <- vector()
      
      if(nrow(men.age.table.1) > 1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
        }
        
      }
      
      
      
      
      
      # women 25.40 and men
      
      women.25.40.men.15.25.2 <- vector()
      women.25.40.men.25.40.2 <- vector()
      women.25.40.men.40.50.2 <- vector()
      
      if(nrow(women.age.table.1) >1 ){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      # men 40.50 and women
      
      men.40.50.women.15.25.1 <- vector()
      men.40.50.women.25.40.1 <- vector()
      men.40.50.women.40.50.1 <- vector()
      
      
      if(nrow(men.age.table.1) >1 ){
        
        for (j in 1:nrow(men.age.table.1)) {
          
          
          if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
              
            }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
              
            }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
      }
      
      
      
      
      
      # women 40.50 and men
      
      women.40.50.men.15.25.2 <- vector()
      women.40.50.men.25.40.2 <- vector()
      women.40.50.men.40.50.2 <- vector()
      
      
      if(nrow(women.age.table.1) > 1 ){
        
        for (j in 1:nrow(women.age.table.1)) {
          
          
          if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
            
            if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
              
              women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
              
            }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
              
              women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
              
            }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
              
              women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
            }
            
          }
          
          
        }
        
      }
      
      
      
      
      
      men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
      
      men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
      
      men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
      
      men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
      
      men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
      
      men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
      
      men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
      
      men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
      
      men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
      
      Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
                                   length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
                                   length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
                                 ncol = 3,
                                 byrow = TRUE)
      
      colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      
      Age.groups.table <- as.table(Age.groups.table)
      
      
      men.15.25.T <- sum(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50))
      men.25.40.T <- sum(length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50))
      men.40.50.T <- sum(length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50))
      
      prop.men.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/men.15.25.T, length(men.15.25.women.25.40)/men.15.25.T, length(men.15.25.women.40.50)/men.15.25.T,
                                            length(men.25.40.women.15.25)/men.25.40.T, length(men.25.40.women.25.40)/men.25.40.T, length(men.25.40.women.40.50)/men.25.40.T,
                                            length(men.40.50.women.15.25)/men.40.50.T, length(men.40.50.women.25.40)/men.40.50.T, length(men.40.50.women.40.50)/men.40.50.T),
                                          ncol = 3,
                                          byrow = TRUE)
      
      colnames(prop.men.age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
      rownames(prop.men.age.groups.table) <- c("prop.Male.15.25", "prop.Male.25.40", "prop.Male.40.50")
      
      
      
      
      women.15.25.T <- sum(length(men.15.25.women.15.25), length(men.25.40.women.15.25), length(men.40.50.women.15.25))
      women.25.40.T <- sum(length(men.15.25.women.25.40), length(men.25.40.women.25.40), length(men.40.50.women.25.40))
      women.40.50.T <- sum(length(men.15.25.women.40.50), length(men.25.40.women.40.50), length(men.40.50.women.40.50))
      
      prop.women.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/women.15.25.T, length(men.25.40.women.15.25)/women.15.25.T, length(men.40.50.women.15.25)/women.15.25.T,
                                              length(men.15.25.women.25.40)/women.25.40.T, length(men.25.40.women.25.40)/women.25.40.T, length(men.40.50.women.25.40)/women.25.40.T,
                                              length(men.15.25.women.40.50)/women.40.50.T, length(men.25.40.women.40.50)/women.40.50.T, length(men.40.50.women.40.50)/women.40.50.T),
                                            ncol = 3,
                                            byrow = TRUE)
      
      colnames(prop.women.age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
      rownames(prop.women.age.groups.table) <- c("prop.Female.15.25", "prop.Female.25.40", "prop.Female.40.50")
      
      outputlist <- NULL
      outputlist$Age.groups.table <- Age.groups.table
      outputlist$prop.men.age.groups.table <- prop.men.age.groups.table
      outputlist$prop.women.age.groups.table <- prop.women.age.groups.table
      
      
      return(outputlist)
      
    }
    
    
    
    
    # Results
    ###########
    
    # Function to handle NAs
    
    
    NA.handle.fun <- function(input=input){
      
      v.names <- names(input)
      
      v <- as.numeric(input)
      
      v.vec <- vector()
      
      for(i in 1:length(v)){
        
        v.i <- v[i]
        
        if(is.na(v.i)==TRUE){
          v.j <- 0
        }else{
          v.j <- v.i
        }
        v.vec <- c(v.vec, v.j)
      }
      
      names(v.vec) <- v.names
      return(v.vec)
    }
    
    # 1. Age structure table obtained from transmission network built from transmission clusters
    
    age.structure.transm.clust.List <- age.groups.filtered.trans.clust.network.fun(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                                   transm.matrix = transm.matrix,
                                                                                   age.group.15.25 = c(15,25),
                                                                                   age.group.25.40 = c(25,40),
                                                                                   age.group.40.50 = c(40,50))
    age.structure.transm.clust <- age.structure.transm.clust.List$Age.groups.table
    
    cl.age.str.M.15.25.F.15.25 <- age.structure.transm.clust[1,][1]
    cl.age.str.M.25.40.F.15.25 <- age.structure.transm.clust[2,][1]
    cl.age.str.M.40.50.F.15.25 <- age.structure.transm.clust[3,][1]
    
    cl.age.str.M.15.25.F.25.40 <- age.structure.transm.clust[1,][2]
    cl.age.str.M.25.40.F.25.40 <- age.structure.transm.clust[2,][2]
    cl.age.str.M.40.50.F.25.40 <- age.structure.transm.clust[3,][2]
    
    cl.age.str.M.15.25.F.40.50 <- age.structure.transm.clust[1,][3]
    cl.age.str.M.25.40.F.40.50 <- age.structure.transm.clust[2,][3]
    cl.age.str.M.40.50.F.40.50 <- age.structure.transm.clust[3,][3]
    
    table.cl.age.str <- c(cl.age.str.M.15.25.F.15.25, cl.age.str.M.25.40.F.15.25, cl.age.str.M.40.50.F.15.25,
                          cl.age.str.M.15.25.F.25.40, cl.age.str.M.25.40.F.25.40, cl.age.str.M.40.50.F.25.40,
                          cl.age.str.M.15.25.F.40.50, cl.age.str.M.25.40.F.40.50, cl.age.str.M.40.50.F.40.50)
    
    names(table.cl.age.str) <- c("cl.M.15.25.F.15.25", "cl.M.25.40.F.15.25", "cl.M.40.50.F.15.25",
                                 "cl.M.15.25.F.25.40", "cl.M.25.40.F.25.40", "cl.M.40.50.F.25.40",
                                 "cl.M.15.25.F.40.50", "cl.M.25.40.F.40.50", "cl.M.40.50.F.40.50")
    
    
    # Men prop
    
    age.structure.transm.clust.prop.men <- age.structure.transm.clust.List$prop.men.age.groups.table
    
    cl.age.str.prop.men.15.25.F.15.25 <- age.structure.transm.clust.prop.men[1,][1]
    cl.age.str.prop.men.25.40.F.15.25 <- age.structure.transm.clust.prop.men[2,][1]
    cl.age.str.prop.men.40.50.F.15.25 <- age.structure.transm.clust.prop.men[3,][1]
    
    cl.age.str.prop.men.15.25.F.25.40 <- age.structure.transm.clust.prop.men[1,][2]
    cl.age.str.prop.men.25.40.F.25.40 <- age.structure.transm.clust.prop.men[2,][2]
    cl.age.str.prop.men.40.50.F.25.40 <- age.structure.transm.clust.prop.men[3,][2]
    
    cl.age.str.prop.men.15.25.F.40.50 <- age.structure.transm.clust.prop.men[1,][3]
    cl.age.str.prop.men.25.40.F.40.50 <- age.structure.transm.clust.prop.men[2,][3]
    cl.age.str.prop.men.40.50.F.40.50 <- age.structure.transm.clust.prop.men[3,][3]
    
    table.cl.age.str.prop.men <- c(cl.age.str.prop.men.15.25.F.15.25, cl.age.str.prop.men.25.40.F.15.25, cl.age.str.prop.men.40.50.F.15.25,
                                   cl.age.str.prop.men.15.25.F.25.40, cl.age.str.prop.men.25.40.F.25.40, cl.age.str.prop.men.40.50.F.25.40,
                                   cl.age.str.prop.men.15.25.F.40.50, cl.age.str.prop.men.25.40.F.40.50, cl.age.str.prop.men.40.50.F.40.50)
    
    names(table.cl.age.str.prop.men) <- c("cl.prop.men15.25.F.15.25", "cl.prop.men25.40.F.15.25", "cl.prop.men40.50.F.15.25",
                                          "cl.prop.men15.25.F.25.40", "cl.prop.men25.40.F.25.40", "cl.prop.men40.50.F.25.40",
                                          "cl.prop.men15.25.F.40.50", "cl.prop.men25.40.F.40.50", "cl.prop.men40.50.F.40.50")
    
    table.cl.age.str.prop.men <- NA.handle.fun(input = table.cl.age.str.prop.men)
    
    
    # Women prop
    
    age.structure.transm.clust.prop.women <- age.structure.transm.clust.List$prop.women.age.groups.table
    
    cl.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.prop.women[1,][1]
    cl.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.prop.women[2,][1]
    cl.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.prop.women[3,][1]
    
    cl.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.prop.women[1,][2]
    cl.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.prop.women[2,][2]
    cl.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.prop.women[3,][2]
    
    cl.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.prop.women[1,][3]
    cl.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.prop.women[2,][3]
    cl.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.prop.women[3,][3]
    
    table.cl.age.str.prop.women <- c(cl.age.str.prop.women.15.25.M.15.25, cl.age.str.prop.women.25.40.M.15.25, cl.age.str.prop.women.40.50.M.15.25,
                                     cl.age.str.prop.women.15.25.M.25.40, cl.age.str.prop.women.25.40.M.25.40, cl.age.str.prop.women.40.50.M.25.40,
                                     cl.age.str.prop.women.15.25.M.40.50, cl.age.str.prop.women.25.40.M.40.50, cl.age.str.prop.women.40.50.M.40.50)
    
    names(table.cl.age.str.prop.women) <- c("cl.prop.women15.25.M.15.25", "cl.prop.women25.40.M.15.25", "cl.prop.women40.50.M.15.25",
                                            "cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40", "cl.prop.women40.50.M.25.40",
                                            "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50")
    
    table.cl.age.str.prop.women <- NA.handle.fun(input = table.cl.age.str.prop.women)
    
    
    
    # 2. True age structure table from transmission network of these individuals in transmission clusters
    
    age.structure.transm.clus.true.List <- age.groups.filtered.transmission.clust.fun(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                                      age.group.15.25 = c(15,25),
                                                                                      age.group.25.40 = c(25,40),
                                                                                      age.group.40.50 = c(40,50))
    
    
    age.structure.transm.clust.true <- age.structure.transm.clus.true.List$Age.groups.table
    
    cl.true.age.str.M.15.25.F.15.25 <- age.structure.transm.clust.true[1,][1]
    cl.true.age.str.M.25.40.F.15.25 <- age.structure.transm.clust.true[2,][1]
    cl.true.age.str.M.40.50.F.15.25 <- age.structure.transm.clust.true[3,][1]
    
    cl.true.age.str.M.15.25.F.25.40 <- age.structure.transm.clust.true[1,][2]
    cl.true.age.str.M.25.40.F.25.40 <- age.structure.transm.clust.true[2,][2]
    cl.true.age.str.M.40.50.F.25.40 <- age.structure.transm.clust.true[3,][2]
    
    cl.true.age.str.M.15.25.F.40.50 <- age.structure.transm.clust.true[1,][3]
    cl.true.age.str.M.25.40.F.40.50 <- age.structure.transm.clust.true[2,][3]
    cl.true.age.str.M.40.50.F.40.50 <- age.structure.transm.clust.true[3,][3]
    
    table.cl.true.age.str <- c(cl.true.age.str.M.15.25.F.15.25, cl.true.age.str.M.25.40.F.15.25, cl.true.age.str.M.40.50.F.15.25,
                               cl.true.age.str.M.15.25.F.25.40, cl.true.age.str.M.25.40.F.25.40, cl.true.age.str.M.40.50.F.25.40,
                               cl.true.age.str.M.15.25.F.40.50, cl.true.age.str.M.25.40.F.40.50, cl.true.age.str.M.40.50.F.40.50)
    
    
    names(table.cl.true.age.str) <- c("cl.true.M.15.25.F.15.25", "cl.true.M.25.40.F.15.25", "cl.true.M.40.50.F.15.25",
                                      "cl.true.M.15.25.F.25.40", "cl.true.M.25.40.F.25.40", "cl.str.M.40.50.F.25.40",
                                      "cl.true.M.15.25.F.40.50", "cl.true..M.25.40.F.40.50", "cl.M.40.50.F.40.50")
    
    
    # Men prop
    
    age.structure.transm.clus.true.prop.men <- age.structure.transm.clus.true.List$prop.men.age.groups.table
    
    cl.true.age.str.prop.men.15.25.F.15.25 <- age.structure.transm.clus.true.prop.men[1,][1]
    cl.true.age.str.prop.men.25.40.F.15.25 <- age.structure.transm.clus.true.prop.men[2,][1]
    cl.true.age.str.prop.men.40.50.F.15.25 <- age.structure.transm.clus.true.prop.men[3,][1]
    
    cl.true.age.str.prop.men.15.25.F.25.40 <- age.structure.transm.clus.true.prop.men[1,][2]
    cl.true.age.str.prop.men.25.40.F.25.40 <- age.structure.transm.clus.true.prop.men[2,][2]
    cl.true.age.str.prop.men.40.50.F.25.40 <- age.structure.transm.clus.true.prop.men[3,][2]
    
    cl.true.age.str.prop.men.15.25.F.40.50 <- age.structure.transm.clus.true.prop.men[1,][3]
    cl.true.age.str.prop.men.25.40.F.40.50 <- age.structure.transm.clus.true.prop.men[2,][3]
    cl.true.age.str.prop.men.40.50.F.40.50 <- age.structure.transm.clus.true.prop.men[3,][3]
    
    table.cl.true.age.str.prop.men <- c(cl.true.age.str.prop.men.15.25.F.15.25, cl.true.age.str.prop.men.25.40.F.15.25, cl.true.age.str.prop.men.40.50.F.15.25,
                                        cl.true.age.str.prop.men.15.25.F.25.40, cl.true.age.str.prop.men.25.40.F.25.40, cl.true.age.str.prop.men.40.50.F.25.40,
                                        cl.true.age.str.prop.men.15.25.F.40.50, cl.true.age.str.prop.men.25.40.F.40.50, cl.true.age.str.prop.men.40.50.F.40.50)
    
    names(table.cl.true.age.str.prop.men) <- c("cl.true.prop.men15.25.F.15.25", "cl.true.prop.men25.40.F.15.25", "cl.true.prop.men40.50.F.15.25",
                                               "cl.true.prop.men15.25.F.25.40", "cl.true.prop.men25.40.F.25.40", "cl.true.prop.men40.50.F.25.40",
                                               "cl.true.prop.men15.25.F.40.50", "cl.true.prop.men25.40.F.40.50", "cl.true.prop.men40.50.F.40.50")
    
    table.cl.true.age.str.prop.men <- NA.handle.fun(input = table.cl.true.age.str.prop.men)
    
    
    
    # Women prop
    
    age.structure.transm.clust.true.prop.women <- age.structure.transm.clust.List$prop.women.age.groups.table
    
    cl.true.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.true.prop.women[1,][1]
    cl.true.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.true.prop.women[2,][1]
    cl.true.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.true.prop.women[3,][1]
    
    cl.true.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.true.prop.women[1,][2]
    cl.true.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.true.prop.women[2,][2]
    cl.true.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.true.prop.women[3,][2]
    
    cl.true.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.true.prop.women[1,][3]
    cl.true.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.true.prop.women[2,][3]
    cl.true.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.true.prop.women[3,][3]
    
    table.cl.true.age.str.prop.women <- c(cl.true.age.str.prop.women.15.25.M.15.25, cl.true.age.str.prop.women.25.40.M.15.25, cl.true.age.str.prop.women.40.50.M.15.25,
                                          cl.true.age.str.prop.women.15.25.M.25.40, cl.true.age.str.prop.women.25.40.M.25.40, cl.true.age.str.prop.women.40.50.M.25.40,
                                          cl.true.age.str.prop.women.15.25.M.40.50, cl.true.age.str.prop.women.25.40.M.40.50, cl.true.age.str.prop.women.40.50.M.40.50)
    
    names(table.cl.true.age.str.prop.women) <- c("cl.true.prop.women15.25.M.15.25", "cl.true.prop.women25.40.M.15.25", "cl.true.prop.women40.50.M.15.25",
                                                 "cl.true.prop.women15.25.M.25.40", "cl.true.prop.women25.40.M.25.40", "cl.true.prop.women40.50.M.25.40",
                                                 "cl.true.prop.women15.25.M.40.50", "cl.true.prop.women25.40.M.40.50", "cl.true.prop.women40.50.M.40.50")
    
    table.cl.true.age.str.prop.women <- NA.handle.fun(input = table.cl.true.age.str.prop.women)
    
    
    
    
    # 3. True age structure table from transmission network of all selected individuals (people in the phylogenetic tree)
    
    age.structure.transm.net.true.List <- age.groups.filtered.transmission.net.fun(table.transmission.net.cov = table.simpact.trans.net.cov,
                                                                                   age.group.15.25 = c(15,25),
                                                                                   age.group.25.40 = c(25,40),
                                                                                   age.group.40.50 = c(40,50))
    
    age.structure.transm.net.true <- age.structure.transm.net.true.List$Age.groups.table
    
    tree.tra.age.str.M.15.25.F.15.25 <- age.structure.transm.net.true[1,][1]
    tree.tra.age.str.M.25.40.F.15.25 <- age.structure.transm.net.true[2,][1]
    tree.tra.age.str.M.40.50.F.15.25 <- age.structure.transm.net.true[3,][1]
    
    tree.tra.age.str.M.15.25.F.25.40 <- age.structure.transm.net.true[1,][2]
    tree.tra.age.str.M.25.40.F.25.40 <- age.structure.transm.net.true[2,][2]
    tree.tra.age.str.M.40.50.F.25.40 <- age.structure.transm.net.true[3,][2]
    
    tree.tra.age.str.M.15.25.F.40.50 <- age.structure.transm.net.true[1,][3]
    tree.tra.age.str.M.25.40.F.40.50 <- age.structure.transm.net.true[2,][3]
    tree.tra.age.str.M.40.50.F.40.50 <- age.structure.transm.net.true[3,][3]
    
    table.tree.tra.age.str <- c(tree.tra.age.str.M.15.25.F.15.25, tree.tra.age.str.M.25.40.F.15.25, tree.tra.age.str.M.40.50.F.15.25,
                                tree.tra.age.str.M.15.25.F.25.40, tree.tra.age.str.M.25.40.F.25.40, tree.tra.age.str.M.40.50.F.25.40,
                                tree.tra.age.str.M.15.25.F.40.50, tree.tra.age.str.M.25.40.F.40.50, tree.tra.age.str.M.40.50.F.40.50)
    
    
    names(table.tree.tra.age.str) <- c("tree.tra.M.15.25.F.15.25", "tree.tra.M.25.40.F.15.25", "tree.tra.M.40.50.F.15.25",
                                       "tree.tra.M.15.25.F.25.40", "tree.tra.M.25.40.F.25.40", "tree.tra.M.40.50.F.25.40",
                                       "tree.tra.M.15.25.F.40.50", "tree.tra.M.25.40.F.40.50", "tree.tra.M.40.50.F.40.50")
    
    
    # Men prop
    
    age.structure.transm.net.true.prop.men <- age.structure.transm.net.true.List$prop.men.age.groups.table
    
    tree.trans.true.age.str.prop.men.15.25.F.15.25 <- age.structure.transm.net.true.prop.men[1,][1]
    tree.trans.true.age.str.prop.men.25.40.F.15.25 <- age.structure.transm.net.true.prop.men[2,][1]
    tree.trans.true.age.str.prop.men.40.50.F.15.25 <- age.structure.transm.net.true.prop.men[3,][1]
    
    tree.trans.true.age.str.prop.men.15.25.F.25.40 <- age.structure.transm.net.true.prop.men[1,][2]
    tree.trans.true.age.str.prop.men.25.40.F.25.40 <- age.structure.transm.net.true.prop.men[2,][2]
    tree.trans.true.age.str.prop.men.40.50.F.25.40 <- age.structure.transm.net.true.prop.men[3,][2]
    
    tree.trans.true.age.str.prop.men.15.25.F.40.50 <- age.structure.transm.net.true.prop.men[1,][3]
    tree.trans.true.age.str.prop.men.25.40.F.40.50 <- age.structure.transm.net.true.prop.men[2,][3]
    tree.trans.true.age.str.prop.men.40.50.F.40.50 <- age.structure.transm.net.true.prop.men[3,][3]
    
    table.tree.trans.true.age.str.prop.men <- c(tree.trans.true.age.str.prop.men.15.25.F.15.25, tree.trans.true.age.str.prop.men.25.40.F.15.25, tree.trans.true.age.str.prop.men.40.50.F.15.25,
                                                tree.trans.true.age.str.prop.men.15.25.F.25.40, tree.trans.true.age.str.prop.men.25.40.F.25.40, tree.trans.true.age.str.prop.men.40.50.F.25.40,
                                                tree.trans.true.age.str.prop.men.15.25.F.40.50, tree.trans.true.age.str.prop.men.25.40.F.40.50, tree.trans.true.age.str.prop.men.40.50.F.40.50)
    
    names(table.tree.trans.true.age.str.prop.men) <- c("tree.trans.true.prop.men15.25.F.15.25", "tree.trans.true.prop.men25.40.F.15.25", "tree.trans.true.prop.men40.50.F.15.25",
                                                       "tree.trans.true.prop.men15.25.F.25.40", "tree.trans.true.prop.men25.40.F.25.40", "tree.trans.true.prop.men40.50.F.25.40",
                                                       "tree.trans.true.prop.men15.25.F.40.50", "tree.trans.true.prop.men25.40.F.40.50", "tree.trans.true.prop.men40.50.F.40.50")
    
    
    table.tree.trans.true.age.str.prop.men <- NA.handle.fun(input = table.tree.trans.true.age.str.prop.men)
    
    
    
    # Women prop
    
    age.structure.transm.clust.true.prop.women <- age.structure.transm.net.true.List$prop.women.age.groups.table
    
    tree.trans.true.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.true.prop.women[1,][1]
    tree.trans.true.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.true.prop.women[2,][1]
    tree.trans.true.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.true.prop.women[3,][1]
    
    tree.trans.true.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.true.prop.women[1,][2]
    tree.trans.true.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.true.prop.women[2,][2]
    tree.trans.true.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.true.prop.women[3,][2]
    
    tree.trans.true.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.true.prop.women[1,][3]
    tree.trans.true.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.true.prop.women[2,][3]
    tree.trans.true.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.true.prop.women[3,][3]
    
    table.tree.trans.true.age.str.prop.women <- c(tree.trans.true.age.str.prop.women.15.25.M.15.25, tree.trans.true.age.str.prop.women.25.40.M.15.25, tree.trans.true.age.str.prop.women.40.50.M.15.25,
                                                  tree.trans.true.age.str.prop.women.15.25.M.25.40, tree.trans.true.age.str.prop.women.25.40.M.25.40, tree.trans.true.age.str.prop.women.40.50.M.25.40,
                                                  tree.trans.true.age.str.prop.women.15.25.M.40.50, tree.trans.true.age.str.prop.women.25.40.M.40.50, tree.trans.true.age.str.prop.women.40.50.M.40.50)
    
    names(table.tree.trans.true.age.str.prop.women) <- c("tree.trans.true.prop.women15.25.M.15.25", "tree.trans.true.prop.women25.40.M.15.25", "tree.trans.true.prop.women40.50.M.15.25",
                                                         "tree.trans.true.prop.women15.25.M.25.40", "tree.trans.true.prop.women25.40.M.25.40", "tree.trans.true.prop.women40.50.M.25.40",
                                                         "tree.trans.true.prop.women15.25.M.40.50", "tree.trans.true.prop.women25.40.M.40.50", "tree.trans.true.prop.women40.50.M.40.50")
    
    
    
    table.tree.trans.true.age.str.prop.women <- NA.handle.fun(input = table.tree.trans.true.age.str.prop.women)
    
    
    
    
    # 4. True age structure table from transmission network of all (coverage 100%) selected individuals
    
    # table.simpact.trans.net.cov.100
    
    age.structure.transm.net.true.100.List <- age.groups.filtered.transmission.net.fun(table.transmission.net.cov = table.simpact.trans.net.cov.100,
                                                                                       age.group.15.25 = c(15,25),
                                                                                       age.group.25.40 = c(25,40),
                                                                                       age.group.40.50 = c(40,50))
    
    
    
    age.struc.trans.net.true.cov.100 <- age.structure.transm.net.true.100.List$Age.groups.table
    
    cov.100.age.str.M.15.25.F.15.25 <- age.struc.trans.net.true.cov.100[1,][1]
    cov.100.age.str.M.25.40.F.15.25 <- age.struc.trans.net.true.cov.100[2,][1]
    cov.100.age.str.M.40.50.F.15.25 <- age.struc.trans.net.true.cov.100[3,][1]
    
    cov.100.age.str.M.15.25.F.25.40 <- age.struc.trans.net.true.cov.100[1,][2]
    cov.100.age.str.M.25.40.F.25.40 <- age.struc.trans.net.true.cov.100[2,][2]
    cov.100.age.str.M.40.50.F.25.40 <- age.struc.trans.net.true.cov.100[3,][2]
    
    cov.100.age.str.M.15.25.F.40.50 <- age.struc.trans.net.true.cov.100[1,][3]
    cov.100.age.str.M.25.40.F.40.50 <- age.struc.trans.net.true.cov.100[2,][3]
    cov.100.age.str.M.40.50.F.40.50 <- age.struc.trans.net.true.cov.100[3,][3]
    
    table.cov.100.age.str <- c(cov.100.age.str.M.15.25.F.15.25, cov.100.age.str.M.25.40.F.15.25, cov.100.age.str.M.40.50.F.15.25,
                               cov.100.age.str.M.15.25.F.25.40, cov.100.age.str.M.25.40.F.25.40, cov.100.age.str.M.40.50.F.25.40,
                               cov.100.age.str.M.15.25.F.40.50, cov.100.age.str.M.25.40.F.40.50, cov.100.age.str.M.40.50.F.40.50)
    
    
    names(table.cov.100.age.str) <- c("cov.100.M.15.25.F.15.25", "cov.100.M.25.40.F.15.25", "cov.100.M.40.50.F.15.25",
                                      "cov.100.M.15.25.F.25.40", "cov.100.M.25.40.F.25.40", "cov.100.M.40.50.F.25.40",
                                      "cov.100.M.15.25.F.40.50", "cov.100.M.25.40.F.40.50", "cov.100.M.40.50.F.40.50")
    
    
    # Men prop
    
    age.struc.trans.net.true.cov.100.prop.men <- age.structure.transm.net.true.100.List$prop.men.age.groups.table
    
    cov.100.true.age.str.prop.men.15.25.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[1,][1]
    cov.100.true.age.str.prop.men.25.40.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[2,][1]
    cov.100.true.age.str.prop.men.40.50.F.15.25 <- age.struc.trans.net.true.cov.100.prop.men[3,][1]
    
    cov.100.true.age.str.prop.men.15.25.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[1,][2]
    cov.100.true.age.str.prop.men.25.40.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[2,][2]
    cov.100.true.age.str.prop.men.40.50.F.25.40 <- age.struc.trans.net.true.cov.100.prop.men[3,][2]
    
    cov.100.true.age.str.prop.men.15.25.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[1,][3]
    cov.100.true.age.str.prop.men.25.40.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[2,][3]
    cov.100.true.age.str.prop.men.40.50.F.40.50 <- age.struc.trans.net.true.cov.100.prop.men[3,][3]
    
    table.cov.100.true.age.str.prop.men <- c(cov.100.true.age.str.prop.men.15.25.F.15.25, cov.100.true.age.str.prop.men.25.40.F.15.25, cov.100.true.age.str.prop.men.40.50.F.15.25,
                                             cov.100.true.age.str.prop.men.15.25.F.25.40, cov.100.true.age.str.prop.men.25.40.F.25.40, cov.100.true.age.str.prop.men.40.50.F.25.40,
                                             cov.100.true.age.str.prop.men.15.25.F.40.50, cov.100.true.age.str.prop.men.25.40.F.40.50, cov.100.true.age.str.prop.men.40.50.F.40.50)
    
    names(table.cov.100.true.age.str.prop.men) <- c("cov.100.true.prop.men15.25.F.15.25", "cov.100.true.prop.men25.40.F.15.25", "cov.100.true.prop.men40.50.F.15.25",
                                                    "cov.100.true.prop.men15.25.F.25.40", "cov.100.true.prop.men25.40.F.25.40", "cov.100.true.prop.men40.50.F.25.40",
                                                    "cov.100.true.prop.men15.25.F.40.50", "cov.100.true.prop.men25.40.F.40.50", "cov.100.true.prop.men40.50.F.40.50")
    
    
    table.cov.100.true.age.str.prop.men <- NA.handle.fun(input = table.cov.100.true.age.str.prop.men)
    
    
    
    # Women prop
    
    age.structure.transm.clust.true.prop.women <- age.structure.transm.net.true.100.List$prop.women.age.groups.table
    
    cov.100.true.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.true.prop.women[1,][1]
    cov.100.true.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.true.prop.women[2,][1]
    cov.100.true.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.true.prop.women[3,][1]
    
    cov.100.true.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.true.prop.women[1,][2]
    cov.100.true.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.true.prop.women[2,][2]
    cov.100.true.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.true.prop.women[3,][2]
    
    cov.100.true.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.true.prop.women[1,][3]
    cov.100.true.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.true.prop.women[2,][3]
    cov.100.true.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.true.prop.women[3,][3]
    
    table.cov.100.true.age.str.prop.women <- c(cov.100.true.age.str.prop.women.15.25.M.15.25, cov.100.true.age.str.prop.women.25.40.M.15.25, cov.100.true.age.str.prop.women.40.50.M.15.25,
                                               cov.100.true.age.str.prop.women.15.25.M.25.40, cov.100.true.age.str.prop.women.25.40.M.25.40, cov.100.true.age.str.prop.women.40.50.M.25.40,
                                               cov.100.true.age.str.prop.women.15.25.M.40.50, cov.100.true.age.str.prop.women.25.40.M.40.50, cov.100.true.age.str.prop.women.40.50.M.40.50)
    
    names(table.cov.100.true.age.str.prop.women) <- c("cov.100.true.prop.women15.25.M.15.25", "cov.100.true.prop.women25.40.M.15.25", "cov.100.true.prop.women40.50.M.15.25",
                                                      "cov.100.true.prop.women15.25.M.25.40", "cov.100.true.prop.women25.40.M.25.40", "cov.100.true.prop.women40.50.M.25.40",
                                                      "cov.100.true.prop.women15.25.M.40.50", "cov.100.true.prop.women25.40.M.40.50", "cov.100.true.prop.women40.50.M.40.50")
    
    
    
    table.cov.100.true.age.str.prop.women <- NA.handle.fun(input = table.cov.100.true.age.str.prop.women)
    
    
    
    table.statistics <- c(table.cl.age.str, table.cl.age.str.prop.men, table.cl.age.str.prop.women, 
                          table.cl.true.age.str, table.cl.true.age.str.prop.men, table.cl.true.age.str.prop.women, 
                          table.tree.tra.age.str, table.tree.trans.true.age.str.prop.men, table.tree.trans.true.age.str.prop.women)
    
    cov.100.stats <- c(table.cov.100.age.str, table.cov.100.true.age.str.prop.men, table.cov.100.true.age.str.prop.women)
    
    
    # Age mixing metrics in the data sets of selected individuals
    
    age.mixing.transm <- mix.rels.transm.dat
    
    
    # Clusters statistics
    
    mean.clust.size <- mean(clust.size)
    median.clust.size <- median(clust.size)
    sd.clust.size <- sd(clust.size)
    
    clust.size.stat <- c(mean.clust.size, median.clust.size, sd.clust.size)
    
    names(clust.size.stat) <- c("mean.cl.size", "med.cl.size", "sd.cl.size")
    
    
    
    # table.statistics <- c(table.cl.age.str, table.cl.age.str.prop.men, table.cl.age.str.prop.women, 
    #                       table.cl.true.age.str, table.cl.true.age.str.prop.men, table.cl.true.age.str.prop.women, 
    #                       table.tree.tra.age.str, table.tree.trans.true.age.str.prop.men, table.tree.trans.true.age.str.prop.women)
    
    men.prop.clust <- table.statistics[10:18]
    women.prop.clust <- table.statistics[19:27]
    men.prop.true.clust <- table.statistics[37:45]
    women.prop.true.clust <- table.statistics[46:54]
    men.prop.tree.true <- table.statistics[64:72]
    women.prop.tree.true <- table.statistics[73:81]
    
    men.prop.cov.100 <- cov.100.stats[10:18]
    women.prop.cov.100 <- cov.100.stats[19:27]
    
    
    
    # Function for Root Mean Squared Error
    RMSE <- function(error) {
      rmse <- sqrt(mean(error^2)) 
      return(rmse)}
    
    # Function for Mean Absolute Error
    MAE <- function(error) { 
      mae <- mean(abs(error)) 
      return(mae)
    }
    
    
    # Difference between clusters' inference and true values in 100% coverage
    error.infer.clust.cov.100.men <- as.numeric(men.prop.clust) - as.numeric(men.prop.cov.100)
    RMSE.error.infer.clust.cov.100.men <- RMSE(error.infer.clust.cov.100.men)
    MAE.error.infer.clust.cov.100.men <- MAE(error.infer.clust.cov.100.men)
    
    error.infer.clust.cov.100.women <- as.numeric(women.prop.clust) - as.numeric(women.prop.cov.100)
    RMSE.infer.clust.cov.100.women <- RMSE(error.infer.clust.cov.100.women)
    MAE.infer.clust.cov.100.women <- MAE(error.infer.clust.cov.100.women)
    
    
    
    # Difference between clusters' inference and true values in the transmission tree
    error.infer.clust.tree.men <- as.numeric(men.prop.clust) - as.numeric(men.prop.tree.true)
    RMSE.error.infer.clust.tree.men <- RMSE(error.infer.clust.tree.men)
    MAE.error.infer.clust.tree.men <- MAE(error.infer.clust.tree.men)
    
    error.infer.clust.tree.women <- as.numeric(women.prop.clust) - as.numeric(women.prop.tree.true)
    RMSE.error.infer.clust.tree.women <- RMSE(error.infer.clust.tree.women)
    MAE.error.infer.clust.tree.women <- MAE(error.infer.clust.tree.women)
    
    
    # Difference between clusters' inference and true values in the transmission clusters
    error.infer.clust.true.clust.men <- as.numeric(men.prop.clust) - as.numeric(men.prop.true.clust)
    RMSE.error.infer.clust.true.clust.men <- RMSE(error.infer.clust.true.clust.men)
    MAE.error.infer.clust.true.clust.men <- MAE(error.infer.clust.true.clust.men)
    
    error.infer.clust.true.clust.women <- as.numeric(women.prop.clust) - as.numeric(women.prop.true.clust)
    RMSE.error.infer.clust.true.clust.women <- RMSE(error.infer.clust.true.clust.women)
    MAE.error.infer.clust.true.clust.women <- MAE(error.infer.clust.true.clust.women)
    
    
    # Difference between truth in clusters and truth in the transmission tree
    error.clust.true.tree.men <- as.numeric(men.prop.true.clust) - as.numeric(men.prop.tree.true)
    RMSE.error.clust.true.tree.men <- RMSE(error.clust.true.tree.men)
    MAE.error.clust.true.tree.men <- MAE(error.clust.true.tree.men)
    
    error.clust.true.tree.women <- as.numeric(women.prop.true.clust) - as.numeric(women.prop.tree.true)
    RMSE.error.clust.true.tree.women <- RMSE(error.clust.true.tree.women)
    MAE.error.clust.true.tree.women <- MAE(error.clust.true.tree.women)
    
    GoF.errors <- c(RMSE.error.infer.clust.cov.100.men, MAE.error.infer.clust.cov.100.men,
                    RMSE.infer.clust.cov.100.women, MAE.infer.clust.cov.100.women,
                    RMSE.error.infer.clust.tree.men, MAE.error.infer.clust.tree.men,
                    RMSE.error.infer.clust.tree.women, MAE.error.infer.clust.tree.women,
                    RMSE.error.infer.clust.true.clust.men, MAE.error.infer.clust.true.clust.men,
                    RMSE.error.infer.clust.true.clust.women, MAE.error.infer.clust.true.clust.women,
                    RMSE.error.clust.true.tree.men, MAE.error.clust.true.tree.men,
                    RMSE.error.clust.true.tree.women, MAE.error.clust.true.tree.women)
    
    Names.GoF.errors <- c("RMSE.error.infer.clust.cov.100.men", "MAE.error.infer.clust.cov.100.men",
                          "RMSE.infer.clust.cov.100.women", "MAE.infer.clust.cov.100.women",
                          "RMSE.error.infer.clust.tree.men", "MAE.error.infer.clust.tree.men",
                          "RMSE.error.infer.clust.tree.women", "MAE.error.infer.clust.tree.women",
                          "RMSE.error.infer.clust.true.clust.men", "MAE.error.infer.clust.true.clust.men",
                          "RMSE.error.infer.clust.true.clust.women", "MAE.error.infer.clust.true.clust.women",
                          "RMSE.error.clust.true.tree.men", "MAE.error.clust.true.tree.men",
                          "RMSE.error.clust.true.tree.women", "MAE.error.clust.true.tree.women")
    
    names(GoF.errors) <- Names.GoF.errors
    
    
    output.results.vector.error <- c(table.statistics, age.mixing.transm, clust.size.stat,
                                     cov.100.stats,
                                     GoF.errors)
    
  }else{
    
    output.results.vectot <- rep(NA, 133)
    
    names(output.results.vector.error) <-  c("cl.M.15.25.F.15.25", "cl.M.25.40.F.15.25",  "cl.M.40.50.F.15.25", "cl.M.15.25.F.25.40",                     
                                             "cl.M.25.40.F.25.40", "cl.M.40.50.F.25.40", "cl.M.15.25.F.40.50", "cl.M.25.40.F.40.50",                     
                                             "cl.M.40.50.F.40.50", "cl.prop.men15.25.F.15.25", "cl.prop.men25.40.F.15.25", "cl.prop.men40.50.F.15.25",               
                                             "cl.prop.men15.25.F.25.40", "cl.prop.men25.40.F.25.40", "cl.prop.men40.50.F.25.40", "cl.prop.men15.25.F.40.50",               
                                             "cl.prop.men25.40.F.40.50", "cl.prop.men40.50.F.40.50", "cl.prop.women15.25.M.15.25", "cl.prop.women25.40.M.15.25",             
                                             "cl.prop.women40.50.M.15.25", "cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40" ,  "cl.prop.women40.50.M.25.40",             
                                             "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50", "cl.true.M.15.25.F.15.25",                
                                             "cl.true.M.25.40.F.15.25", "cl.true.M.40.50.F.15.25", "cl.true.M.15.25.F.25.40", "cl.true.M.25.40.F.25.40",                
                                             "cl.str.M.40.50.F.25.40", "cl.true.M.15.25.F.40.50", "cl.true..M.25.40.F.40.50", "cl.M.40.50.F.40.50",                     
                                             "cl.true.prop.men15.25.F.15.25", "cl.true.prop.men25.40.F.15.25", "cl.true.prop.men40.50.F.15.25", "cl.true.prop.men15.25.F.25.40",          
                                             "cl.true.prop.men25.40.F.25.40", "cl.true.prop.men40.50.F.25.40", "cl.true.prop.men15.25.F.40.50", "cl.true.prop.men25.40.F.40.50",          
                                             "cl.true.prop.men40.50.F.40.50", "cl.true.prop.women15.25.M.15.25", "cl.true.prop.women25.40.M.15.25", "cl.true.prop.women40.50.M.15.25",        
                                             "cl.true.prop.women15.25.M.25.40", "cl.true.prop.women25.40.M.25.40", "cl.true.prop.women40.50.M.25.40", "cl.true.prop.women15.25.M.40.50",        
                                             "cl.true.prop.women25.40.M.40.50", "cl.true.prop.women40.50.M.40.50", "tree.tra.M.15.25.F.15.25", "tree.tra.M.25.40.F.15.25",               
                                             "tree.tra.M.40.50.F.15.25", "tree.tra.M.15.25.F.25.40", "tree.tra.M.25.40.F.25.40", "tree.tra.M.40.50.F.25.40",               
                                             "tree.tra.M.15.25.F.40.50", "tree.tra.M.25.40.F.40.50", "tree.tra.M.40.50.F.40.50", "tree.trans.true.prop.men15.25.F.15.25",  
                                             "tree.trans.true.prop.men25.40.F.15.25", "tree.trans.true.prop.men40.50.F.15.25", "tree.trans.true.prop.men15.25.F.25.40",   "tree.trans.true.prop.men25.40.F.25.40",  
                                             "tree.trans.true.prop.men40.50.F.25.40", "tree.trans.true.prop.men15.25.F.40.50", "tree.trans.true.prop.men25.40.F.40.50",   "tree.trans.true.prop.men40.50.F.40.50",  
                                             "tree.trans.true.prop.women15.25.M.15.25", "tree.trans.true.prop.women25.40.M.15.25", "tree.trans.true.prop.women40.50.M.15.25", "tree.trans.true.prop.women15.25.M.25.40",
                                             "tree.trans.true.prop.women25.40.M.25.40", "tree.trans.true.prop.women40.50.M.25.40", "tree.trans.true.prop.women15.25.M.40.50", "tree.trans.true.prop.women25.40.M.40.50",
                                             "tree.trans.true.prop.women40.50.M.40.50", "T.AAD.male", "T.SDAD.male",  "T.slope.male",                           
                                             "T.WSD.male", "T.BSD.male",  "T.intercept.male", "mean.cl.size", "med.cl.size", "sd.cl.size",
                                             
                                             "cov.100.M.15.25.F.15.25", "cov.100.M.25.40.F.15.25", "cov.100.M.40.50.F.15.25", "cov.100.M.15.25.F.25.40", "cov.100.M.25.40.F.25.40",             
                                             "cov.100.M.40.50.F.25.40",  "cov.100.M.15.25.F.40.50", "cov.100.M.25.40.F.40.50", "cov.100.M.40.50.F.40.50", "cov.100.true.prop.men15.25.F.15.25",  
                                             "cov.100.true.prop.men25.40.F.15.25", "cov.100.true.prop.men40.50.F.15.25", "cov.100.true.prop.men15.25.F.25.40", "cov.100.true.prop.men25.40.F.25.40",   
                                             "cov.100.true.prop.men40.50.F.25.40" , "cov.100.true.prop.men15.25.F.40.50", "cov.100.true.prop.men25.40.F.40.50", "cov.100.true.prop.men40.50.F.40.50", 
                                             "cov.100.true.prop.women15.25.M.15.25", "cov.100.true.prop.women25.40.M.15.25", "cov.100.true.prop.women40.50.M.15.25", "cov.100.true.prop.women15.25.M.25.40", 
                                             "cov.100.true.prop.women25.40.M.25.40", "cov.100.true.prop.women40.50.M.25.40", "cov.100.true.prop.women15.25.M.40.50", "cov.100.true.prop.women25.40.M.40.50", 
                                             "cov.100.true.prop.women40.50.M.40.50",
                                             
                                             "RMSE.error.infer.clust.cov.100.men", "MAE.error.infer.clust.cov.100.men",
                                             "RMSE.infer.clust.cov.100.women", "MAE.infer.clust.cov.100.women",
                                             "RMSE.error.infer.clust.tree.men", "MAE.error.infer.clust.tree.men",
                                             "RMSE.error.infer.clust.tree.women", "MAE.error.infer.clust.tree.women",
                                             "RMSE.error.infer.clust.true.clust.men", "MAE.error.infer.clust.true.clust.men",
                                             "RMSE.error.infer.clust.true.clust.women", "MAE.error.infer.clust.true.clust.women",
                                             "RMSE.error.clust.true.tree.men", "MAE.error.clust.true.tree.men",
                                             "RMSE.error.clust.true.tree.women", "MAE.error.clust.true.tree.women")
    
  }
  
  
  return(output.results.vector.error)
  
  
}
