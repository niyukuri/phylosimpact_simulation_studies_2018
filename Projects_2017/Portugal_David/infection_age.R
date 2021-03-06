#' Compute data list of infection age
#'
#' infage <- infection_age(datalist = datalist, endpoint = 40)

infection_age <- function(datalist = datalist, endpoint = 40){

  # HIV seed time
  hivseed.time <- datalist$etable[eventname=="HIV seeding"]$eventtime


  # 1. Table of donors and recipients and time of infection

  ## person table of all infected people
  infec.all.raw <- as.data.frame(datalist$ptable[InfectTime !="Inf"])

  infec.all <- infec.all.raw[which(infec.all.raw$InfectTime <= endpoint), ]

  ## recipients
  id.infec.all <- infec.all$ID

  ## donors
  id.infec.origin.all <- infec.all$InfectOrigID

  ## time of infection
  time.infec.all <- infec.all$InfectTime

  ## table of donors, recipients,  and time of infection including the eternal donor -1

  recDonTime <- cbind(id.infec.origin.all, id.infec.all, time.infec.all)
  recDonTime.data <- as.data.frame(recDonTime)

  ## order the table by time of infection (increasing)
  recDonTimeOrd <- recDonTime.data[order(recDonTime.data$time.infec.all),]


  # 2. Among donors and recipients who got infection by seed/transmission event

  #### Consider the infection types

  # -1: non infected people
  # 0: infected due to seed event >> this help to identifie the seeds IDs
  # 1: infected due to a transmission event

  # person table with all infected people (by seed or transmission event),
  # we then remove those who are not infected
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])

  pers.infec <- pers.infec.raw[which(pers.infec.raw$InfectTime <= endpoint),]

  # person table of infected individuals by seed event
  pers.table.seed <- subset(pers.infec, pers.infec$InfectType==0)

  # id of people who got infection by seed event: seeds.id
  seeds.id <- pers.table.seed$ID # do

  # person table of infected individuals by transmission event
  pers.table.trans <- subset(pers.infec,pers.infec$InfectType==1)



  # 3. One Network per Seed ID: now automated for all the seeds


  ## Now as we have the table of donors, recipients and time of infection,
  ## as we known the seed IDs, let them trace pathways of each seed
  ## bear in mind that each individual can be infected one time only
  # this means that a recipient is unique in the column of recipients

  ## To get the pathways of each seed, let use igraph function to filter these pathways

  # Take the table of donors, recipients and time of infection
  # which is NewrecDonTimeOrd and make a graph object

  library(igraph)

  ga <- recDonTimeOrd

  ga[,1] <- as.character(recDonTimeOrd[,1]) # donors
  ga[,2] <- as.character(recDonTimeOrd[,2]) # recipients
  gag = as.matrix(ga)
  ga.graph = graph.edgelist(gag[,1:2])


  # Let know extract from the above big network
  # the subcomponent for each seed (transmission pathways due to each seed)

  subcomps <- vector("list", length(seeds.id))
  for (i in 1: length(seeds.id)){
    subcomps[[i]] <- subcomponent(ga.graph, paste(seeds.id[i]), "out")
  }

  # table of donors, recipients,  and time of infection of infected people due to each seed
  net.recdontime.raw <- vector("list", length(seeds.id)) # with eternal infector -1
  init.recdontime <- vector("list", length(seeds.id)) # take the eternal infector -1
  net.recdontime <- vector("list", length(seeds.id)) # remove init.recdontime in net.recdontime.raw
  for(i in 1:length(seeds.id)){
    for(j in 1:length(subcomps)){
      if(seeds.id[i] == as.numeric(names(subcomps[[j]]))[1]){
        net.recdontime.raw[[i]] <- subset(recDonTimeOrd,recDonTimeOrd$id.infec.origin.all%in%as.numeric(names(subcomps[[j]])))
        init.recdontime[[i]] <- c(-1,seeds.id[i], unique(as.numeric(datalist$itable$hivseed.time)))
        net.recdontime[[i]] <- rbind(init.recdontime[[i]],net.recdontime.raw[[i]]) # 1st col donors, 2nd recipients and 3rd time of infection
      }
    }
  }



  ################################################################################
  ### itimes, dtimes, id and parent  for epi2tree function of expoTRee package ###
  ################################################################################

  # Note that for dtimes are times of sampling or removing,
  # they are diagnosis time or death time.
  # Thus, we have
  # diagnosed people nd died >> save death time
  # diagnosed people and alive >> save diagnosed time
  # non diagnosed people and alive >> save endpoint time
  # non diagnosed people and died >> save death time


  # Recode the ids of donors and recipients
  # Recipients help to recorde since we know that they are unique in their column

  # ids (infections ids = recipients ids)

  # table of donors, recipients,  and time of infection of infected people due to see 920
  # augmented with new recorded ids

  id1 <- vector("list", length(seeds.id))
  new.recdontime <- vector("list", length(seeds.id))
  # new.parent <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id1[[i]] <- as.vector(seq(from = 0, to = length(net.recdontime[[i]][,2])-1, by = 1))
    new.recdontime[[i]] <- cbind(id1[[i]],net.recdontime[[i]][,2],net.recdontime[[i]][,1],net.recdontime[[i]][,3])
  }

  # parents ids (donors)
  # the first donor is the universal donor -1
  NewParent <- vector("list", length(seeds.id))
  dat.recdontime <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    NewParent1 <- c(-1)
    for(j in 1:length(new.recdontime[[i]][,3])){
      for(k in 1:length(new.recdontime[[i]][,2])){
        if(new.recdontime[[i]][j,3] == new.recdontime[[i]][k,2]){
          NewParent1 <- c(NewParent1, new.recdontime[[i]][k,1])

        }
        else{
          NewParent1 <- c(NewParent1)
        }
      }
    }
    NewParent[[i]] <- c(NewParent1)

    # per seed we have a contact matrix with
    # c("NewID", "OldID", "OldParent", "NewParent", "InfecTime")
    dat.recdontime[[i]] <- cbind(new.recdontime[[i]][,1],new.recdontime[[i]][,2],new.recdontime[[i]][,3],NewParent[[i]], new.recdontime[[i]][,4])
  }

  # Searching sampling and removal times

  diagnosis.event <- as.data.frame(datalist$etable[eventname == "diagnosis"])
  death.hiv.event <- as.data.frame(datalist$etable[eventname == "aidsmortality"])
  death.norm.event <- as.data.frame(datalist$etable[eventname == "normalmortality"])
  keeps <- c("p1ID","eventtime")

  diagnosis <- diagnosis.event[keeps]
  aids.death <- death.hiv.event[keeps]
  normal.death <- death.norm.event[keeps]
  deaths.all <- as.data.frame(rbind(aids.death,normal.death))
  deaths.all.ord <- deaths.all[order(deaths.all$eventtime),]

  # since people can be diagnosed more than one time remove duplicates!
  diagnosis.unique <-  diagnosis[!duplicated(diagnosis$p1ID), ]
  # or subset(diagnosis, !duplicated(p1ID))

  # people alive at the end of simulation >> use of alive.infected() FUNCTION
  alive.at.endpoint <- alive.infected(datalist = datalist, timepoint = endpoint,
                                      site = "All")

  # person table of all people alive at the end of simulation but are HIV positive
  alive.hiv <- subset(alive.at.endpoint, alive.at.endpoint$InfectType !=-1)

  # person table of people alive at the end of simulation but are HIV positive
  # due to only each seed
  transm.dat <- vector("list", length(seeds.id))
  transm.dat.id <- vector("list", length(seeds.id))
  alive.hiv.pers <- vector("list", length(seeds.id))

  for(i in 1:length(seeds.id)){
    transm.dat[[i]] <- dat.recdontime[[i]]
    transm.dat.id[[i]] <- transm.dat[[i]][,2]
    alive.hiv.pers[[i]] <- subset(alive.hiv, alive.hiv$ID%in%transm.dat.id[[i]]) # person table of people alive at the end of simulation but are HIV positive
  }


  # person table of dead people
  pers <- vector("list", length(seeds.id))
  pers.alive.dat <- vector("list", length(seeds.id))
  pers.alive <- vector("list", length(seeds.id))
  pers.death.dat <- vector("list", length(seeds.id))
  pers.death <- vector("list", length(seeds.id))

  for(i in 1:length(seeds.id)){
    pers.alive.dat[[i]] <- subset(alive.hiv.pers[[i]], alive.hiv.pers[[i]]$ID %in% transm.dat.id[[i]])
    pers.alive[[i]] <- pers.alive.dat[[i]]$ID # alive people in the transmission network

    pers.death.dat[[i]] <- subset(deaths.all.ord, deaths.all.ord$p1ID%in%transm.dat.id[[i]]) # pers replaced by transm.dat.id[[i]]
    pers.death[[i]] <- pers.death.dat[[i]]$p1ID # dead people in the transmission network

  }


  # id of dead people in the transmission network = pers - alive.hiv
  id.dead <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id.dead[[i]] <- subset(transm.dat.id[[i]], !transm.dat.id[[i]]%in%alive.hiv.pers[[i]]$ID)
  }

  # id of alive people

  id.alive <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    id.alive[[i]] <- setdiff(transm.dat.id[[i]],id.dead[[i]])
  }


  # split times
  ##############

  # 1. alive and diagnosed people > id & time of diagnosis
  alive.diag <- vector("list", length(seeds.id))
  d1 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    alive.diag[[i]] <- subset(diagnosis.unique, diagnosis.unique$p1ID%in%pers.alive[[i]]) # OK
    d1[[i]] <- as.data.frame(cbind(alive.diag[[i]]$p1ID, alive.diag[[i]]$eventtime))
  }


  # 2. died and diagnosed people > id & time of diagnosis
  died.diag <- vector("list", length(seeds.id))
  d2 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    died.diag[[i]] <- subset(diagnosis.unique, diagnosis.unique$p1ID%in%id.dead[[i]]) # OK
    d2[[i]] <- as.data.frame(cbind(died.diag[[i]]$p1ID, died.diag[[i]]$eventtime))
  }


  # 3. alive non diagnosed people, below person table  > Id & time = endpoint
  id.non.diag.all <- vector("list", length(seeds.id))
  alive.non.diag <- vector("list", length(seeds.id))
  d3 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    # Ids of all no diagnosed
    id.non.diag.all[[i]] <- setdiff(transm.dat.id[[i]],c(alive.diag[[i]]$p1ID, died.diag[[i]]$p1ID))
    alive.non.diag[[i]] <- subset(alive.hiv.pers[[i]], alive.hiv.pers[[i]]$ID%in%id.non.diag.all[[i]])

    d3[[i]] <- as.data.frame(cbind(alive.non.diag[[i]]$ID, rep(endpoint, length(alive.non.diag[[i]]$ID))))
  }


  # 4. died and non diagnosed > Id & death time
  died.non.diag <- vector("list", length(seeds.id))
  d4 <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    died.non.diag[[i]] <- subset(deaths.all.ord, deaths.all.ord$p1ID%in%id.non.diag.all[[i]])
    d4[[i]] <- as.data.frame(cbind(died.non.diag[[i]]$p1ID, died.non.diag[[i]]$eventtime))
  }


  # All together: dtimes

  # id and sampling times together

  times.raw <- vector("list", length(seeds.id))
  id.times.raw <- vector("list", length(seeds.id))
  for(i in 1:length(seeds.id)){
    times.raw[[i]] <- do.call("rbind", list(d1[[i]], d2[[i]], d3[[i]], d4[[i]]))

    id.times.raw[[i]] <- cbind(times.raw[[i]]) # , transm.dat.id[[i]] are infectees for i^th seed
  }


  # rearange the above table according to the order of infection time which is the same
  # order in which we have old ids for each seed
  time.id.dt <- vector("list", length(seeds.id))
  time.dt.all <- vector("list", length(seeds.id))
  time.dt <- vector("list", length(seeds.id))
  time.id <-vector("list", length(seeds.id))


  # define the function to rearange the order anmed times.pers.sort()
  times.pers.sort <- function(id.times.raw,transm.dat.id){
    for(k in 1:length(seeds.id)){
      id.time.ord <- as.data.frame(id.times.raw[[k]]) # old ids and their sampling times
      id.old <- transm.dat.id[[k]] # old ids of infectees for each seed in the infection time order
      time.id.dt <- NULL
      time.dt.init <- vector()
      time.id.init <- vector()
      for(i in 1:nrow(id.time.ord)){
        for(j in 1:length(id.old)){
          if(id.old[i] == id.time.ord$V1[j]){
            # time.id.dt <- cbind(c(time.id, id.times.raw$pers[i]), c(time.dt, id.times.raw$V2[j]))

            time.dt.init <- c(time.dt.init, id.time.ord$V2[j])

            time.id.init <- c(time.id.init, id.old[i])

          }
        }
        time.id.dt <- cbind(time.dt.init, time.id.init)
      }
      time.id[[k]] <- time.id.dt
    }
    return(time.id)
  }

  b = times.pers.sort(id.times.raw, transm.dat.id)

  # to build the epi object handled by epi2tree function
  # to build a transmission tree, we reverse the time for infections time
  # itimes show how old is a given transmitted infection


  # initialize the list of of epi object (one per seed)

  infec.age <- vector("list", length(seeds.id))

  infec <- list()
  #transNet <- list()
  for(i in 1:length(seeds.id)){
    infecage <- (endpoint - (dat.recdontime[[i]][,5]))-hivseed.time
    transNet.id <- dat.recdontime[[i]][,1]
    transNet.parent <- dat.recdontime[[i]][,4]
    all.infage <- cbind(transNet.id, transNet.parent, infecage)
    infec.age[[i]] <- as.data.frame(all.infage)
  }

  return(infec.age)

}
