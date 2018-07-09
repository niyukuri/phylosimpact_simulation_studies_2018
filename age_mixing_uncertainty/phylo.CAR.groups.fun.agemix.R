
phylo.CAR.groups.fun.agemix <- function(simpact.trans.net = simpact.trans.net,
                                  limitTransmEvents = 7,
                                  timewindow = c(30,40),
                                  seq.cov = 70,
                                  #    seq.gender.ratio = 0.7, # within same age group women have 70% of being sampled & men have only 30%
                                  age.group.15.25 = c(15,25),
                                  age.group.25.40 = c(25,40),
                                  age.group.40.50 = c(40,50)){
  
  seeds.id <- length(simpact.trans.net)
  
  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)
  
  for(i in 1:seeds.id){
    
    transm.age.df.ic <- as.data.frame(simpact.trans.net[[i]])
    
    age.inf.Rec <- transm.age.df.ic$InfecTime - transm.age.df.ic$TOBRec
    age.inf.Don <- transm.age.df.ic$InfecTime - transm.age.df.ic$TOBDon
    
    transm.age.i <- cbind(transm.age.df.ic, age.inf.Rec, age.inf.Don)
    
    new.transm.tab[[i]] <- transm.age.i
    
  }
  
  # ID numbers of Selected networks with at least limitTransmEvents + 1 indiviuals
  
  IDs.transm <- vector()
  
  TransmEventsCountVector <- vector()
  
  for(k in 1:seeds.id){
    trans.net.i.check <- as.data.frame(new.transm.tab[[k]])
    
    if(nrow(trans.net.i.check)>=limitTransmEvents){
      
      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
      
      IDs.transm <- c(IDs.transm, k)
    }
  }
  
  if(length(IDs.transm)>=1){
    
    ## Binding together all selected transmission transmission networks ##
    
    for (q in 1:length(IDs.transm)){
      
      if(q==1){
        p <- IDs.transm[q]
        trans.sum <- new.transm.tab[[p]]
        rename.id <- paste0(p,".",trans.sum$id,".C")
        trans.sum$id <- rename.id
        trans.sum.rename.id <- trans.sum
      }
      else{
        
        p <- IDs.transm[q]
        
        read.trans.sum <- new.transm.tab[[p]]
        rename.id.read <- paste0(p,".",read.trans.sum$id,".C")
        read.trans.sum$id <- rename.id.read
        trans.sum.rename.id <- rbind(trans.sum.rename.id, read.trans.sum)
      }
      
    }
    
    # trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$age.inf.Rec<=age.group.40.50[2])
    
    trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$InfecTime >= timewindow[1] & trans.sum.rename.id$InfecTime <= timewindow[2])
    
    
    sort.partners.fun <- function(partner.table = partner.table){
      
      part.15.25 <- dplyr::filter(partner.table, partner.table$age.inf.Don >= age.group.15.25[1] & partner.table$age.inf.Don < age.group.15.25[2])
      part.25.40 <- dplyr::filter(partner.table, partner.table$age.inf.Don >= age.group.25.40[1] & partner.table$age.inf.Don < age.group.25.40[2])
      part.40.50 <- dplyr::filter(partner.table, partner.table$age.inf.Don >= age.group.40.50[1] & partner.table$age.inf.Don < age.group.40.50[2])
      
      N.part.15.25 <- nrow(part.15.25)
      N.part.25.40 <- nrow(part.25.40)
      N.part.40.50 <- nrow(part.40.50)
      
      return(c(N.part.15.25, N.part.25.40, N.part.40.50))
      
    }
    
    # Group 15 - 25
    ###############
    
    trans.sum.men.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.inf.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.inf.Rec < age.group.15.25[2])
    
    trans.sum.women.15.25 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.inf.Rec >= age.group.15.25[1] & trans.sum.age.limit$age.inf.Rec < age.group.15.25[2])
    
    perc.100.15.25 <- nrow(trans.sum.men.15.25) + nrow(trans.sum.women.15.25) # total number of individuals with age limit
    
    trans.sum.men.women.15.25 <- rbind(trans.sum.men.15.25, trans.sum.women.15.25)
    
    
    perc.men.women.15.25 <- round(perc.100.15.25 * seq.cov /100)
    
    
    if(perc.men.women.15.25 <= length(trans.sum.men.women.15.25$id)){
      
      x.y.id.15.25 <- sample(trans.sum.men.women.15.25$id, perc.men.women.15.25)
      
    }else{
      
      x.y.id.15.25 <- sample(trans.sum.men.women.15.25$id, length(trans.sum.men.women.15.25$id))
      
    }
    
    
    samp.all.15.25 <- c(x.y.id.15.25)
    
    
    partners.men.rec.15.25 <- sort.partners.fun(trans.sum.men.15.25) # partners of men when they are recipients
    partners.women.rec.15.25 <- sort.partners.fun(trans.sum.women.15.25) # partners of women when they are recipients
    
    # Group 25 - 40
    ###############
    
    trans.sum.men.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.inf.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.inf.Rec < age.group.25.40[2])
    
    trans.sum.women.25.40 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.inf.Rec >= age.group.25.40[1] & trans.sum.age.limit$age.inf.Rec < age.group.25.40[2])
    
    perc.100.25.40 <- nrow(trans.sum.men.25.40) + nrow(trans.sum.women.25.40) # total number of individuals with age limit
    
    trans.sum.men.women.25.40 <- rbind(trans.sum.men.25.40, trans.sum.women.25.40)
    
    
    perc.men.women.25.40 <- round(perc.100.25.40 * seq.cov /100)
    
    
    if(perc.men.women.25.40 <= length(trans.sum.men.women.25.40$id)){
      
      x.y.id.25.40 <- sample(trans.sum.men.women.25.40$id, perc.men.women.25.40)
      
    }else{
      
      x.y.id.25.40 <- sample(trans.sum.men.women.25.40$id, length(trans.sum.men.women.25.40$id))
      
    }
    
    
    samp.all.25.40 <- c(x.y.id.25.40)
    
    
    
    partners.men.rec.25.40 <- sort.partners.fun(trans.sum.men.25.40) # partners of men when they are recipients
    partners.women.rec.25.40 <- sort.partners.fun(trans.sum.women.25.40) # partners of women when they are recipients
    
    
    # Group 40 - 50
    ###############
    
    trans.sum.men.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$age.inf.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.inf.Rec < age.group.40.50[2])
    
    trans.sum.women.40.50 <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$age.inf.Rec >= age.group.40.50[1] & trans.sum.age.limit$age.inf.Rec < age.group.40.50[2])
    
    perc.100.40.50 <- nrow(trans.sum.men.40.50) + nrow(trans.sum.women.40.50) # total number of individuals with age limit
    
    trans.sum.men.women.40.50 <- rbind(trans.sum.men.40.50, trans.sum.women.40.50)
    
    
    perc.men.women.40.50 <- round(perc.100.40.50 * seq.cov /100)
    
    
    if(perc.men.women.40.50 <= length(trans.sum.men.women.40.50$id)){
      
      x.y.id.40.50 <- sample(trans.sum.men.women.40.50$id, perc.men.women.40.50)
      
    }else{
      
      x.y.id.40.50 <- sample(trans.sum.men.women.40.50$id, length(trans.sum.men.women.40.50$id))
      
    }
    
    
    samp.all.40.50 <- c(x.y.id.40.50)
    
    
    partners.men.rec.40.50 <- sort.partners.fun(trans.sum.men.40.50) # partners of men when they are recipients
    partners.women.rec.40.50 <- sort.partners.fun(trans.sum.women.40.50) # partners of women when they are recipients
    
    
    # samp.all <- c(samp.all.15.25, samp.all.25.40, samp.all.40.50)
    
    ouput.transm.dat <- c(nrow(trans.sum.men.15.25), nrow(trans.sum.women.15.25),
                          nrow(trans.sum.men.25.40), nrow(trans.sum.women.25.40),
                          nrow(trans.sum.men.40.50), nrow(trans.sum.women.40.50),
                          
                          partners.men.rec.15.25, partners.men.rec.25.40, partners.men.rec.40.50,
                          partners.women.rec.15.25, partners.women.rec.25.40, partners.women.rec.40.50)
    
    val.names <- c("num.men.15.25", "num.women.15.25",
                   "num.men.25.40", "num.women.25.40",
                   "num.men.40.50", "num.women.40.50",
                   
                   "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                   "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                   "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                   
                   "partners.women.15.25.w.15.25", "partners.women.15.25.w.25.40", "partners.women.15.25.w.40.50",
                   "partners.women.25.40.w.15.25", "partners.women.25.40.w.25.40", "partners.women.25.40.w.40.50",
                   "partners.women.40.50.w.15.25", "partners.women.40.50.w.25.40", "partners.women.40.50.w.40.50")
    
    names(ouput.transm.dat) <- val.names
    
  }else{
    
    # samp.all <- NA
    
    val.names <- c("num.men.15.25", "num.women.15.25",
                   "num.men.25.40", "num.women.25.40",
                   "num.men.40.50", "num.women.40.50",
                   
                   "partners.men.15.25.w.15.25", "partners.men.15.25.w.25.40", "partners.men.15.25.w.40.50",
                   "partners.men.25.40.w.15.25", "partners.men.25.40.w.25.40", "partners.men.25.40.w.40.50",
                   "partners.men.40.50.w.15.25", "partners.men.40.50.w.25.40", "partners.men.40.50.w.40.50",
                   
                   "partners.women.15.25.w.15.25", "partners.women.15.25.w.25.40", "partners.women.15.25.w.40.50",
                   "partners.women.25.40.w.15.25", "partners.women.25.40.w.25.40", "partners.women.25.40.w.40.50",
                   "partners.women.40.50.w.15.25", "partners.women.40.50.w.25.40", "partners.women.40.50.w.40.50")
    
    ouput.transm.dat <- rep(NA, length(val.names))
    
    names(ouput.transm.dat) <- val.names
    
    
  }
  
  # return(samp.all)
  
  return(ouput.transm.dat)
}





