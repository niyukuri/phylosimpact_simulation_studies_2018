

# Define directory


work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop


# work.dir <- "/home/dniyukuri/lustre/stress_testing_master_model" # on CHPC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


# work.dir <- "~/Desktop/calibration/"



wrapper.stress.testing <- function(inputvector){
  
  
  
  work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop
  
  
  # work.dir <- "/home/dniyukuri/lustre/stress_testing_master_model" # on CHPC
  
  
  setwd(paste0(work.dir))
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(readr)
  library(phangorn)
  library(lme4)
  library(nlme)
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  library(Rsamtools)
  library(robustbase)
  library(intergraph)
  library(lubridate)
  library(tidyr)
  
  source("/home/niyukuri/phylosimpact_simulation_studies_2018/age_mix_final/advanced.transmission.network.builder.R")
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MAR.R")
  
  source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metric.class.phylo.features.cov.R")
  
  
  # source("/home/dniyukuri/lustre/stress_testing_master_model/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/compute.summary.statistics.phylo.MAR.R")
  # 
  # source("/home/dniyukuri/lustre/stress_testing_master_model/complete.master.epic.metric.class.phylo.features.cov.R")
  
  
  results.f <- tryCatch(complete.master.epic.metric.class.phylo.features.cov(inputvector = inputvector),
                        error=function(e) return(rep(NA, 1990)))
  
  
  return(results.f)
  
  
}


reps <- 4

inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 




inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)


epi.mm.stats <- simpact.parallel(model = wrapper.stress.testing,
                                 actual.input.matrix = inputmatrix,
                                 seed_count = 124,
                                 n_cluster = 8)


write.csv(epi.mm.stats, file = "Results.epi.mm.stats.csv")



inputvector <- c(1000, -0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 

source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.classic.R")

source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MCAR.R")

source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MAR.R")

source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/calibration.ABC.R")

epi.metric.behav.stats.stat.missingness


median.epi.metric.behav.stats.stat.missingness <- median(epi.metric.behav.stats.stat.missingness)


epidemic.metrics # len = 39
epi.behav.stats  # len = 27
missingness # len = 37 * 13 * 4 [MCAR, MAR.a.0.7, MAR.b.0.3, MAR.c.0.5]


epidemic.metrics <- median.epi.metric.behav.stats.stat.missingness[1:39]

epi.behav.stats <- median.epi.metric.behav.stats.stat.missingness[40:66]

missingness.cov


