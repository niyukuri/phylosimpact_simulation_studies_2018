
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


