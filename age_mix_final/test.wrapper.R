# Wrapper function

# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC


work.dir <- "/home/dniyukuri/lustre/agemix.25.10.2018.2" # on PCHPC


setwd(paste0(work.dir))


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


test.wrapper <- function(inputvector=inputvector){
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on PC
  
  work.dir <- "/home/dniyukuri/lustre/agemix.25.10.2018.2" # on CHPC
  
  
  setwd(paste0(work.dir))
  
  # source("~/phylosimpact_simulation_studies_2018/age_mix_final/test.MCAR.MAR.age.mix.R")
  
  source("/home/dniyukuri/lustre/agemix.25.10.2018.2/test.MCAR.MAR.age.mix.R")
  
  
  results.f <- tryCatch(test.MCAR.age.mix(inputvector = inputvector),
                        error=function(e) return(rep(NA, 6936)))
  
  return(results.f)
  
}




inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 


reps <- 1000


inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)

mcar.sim1 <- simpact.parallel(model = wrapper.test,
                              actual.input.matrix = inputmatrix,
                              seed_count = 100,
                              n_cluster = 24)

write.csv(mcar.sim1, file = "Results.mcar.sim.complete.csv")




# 
# inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)
# 
# mcar.sim2 <- simpact.parallel(model = wrapper.test,
#                               actual.input.matrix = inputmatrix,
#                               seed_count = 200,
#                               n_cluster = 20)
# 
# write.csv(mcar.sim2, file = "Results.mcar.sim2.csv")
# 
# 
# 

