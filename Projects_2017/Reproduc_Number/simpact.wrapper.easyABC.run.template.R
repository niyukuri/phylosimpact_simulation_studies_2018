#!/usr/bin/env/ Rscript
#get the necessary libraries
pacman::p_load(dplyr, EasyABC, RSimpactCyan, RSimpactHelper, lhs)
#data file to read

comp <- "lin" #lin #mac #chpc #gent

if(comp == "win"){dirname <- "~/MaxART/RSimpactHelp"}else if(comp=="lin"){
    dirname <- "~/Documents/GIT_Projects/RSimpactHelp"}else if(comp=="chpc"){
    dirname <- "/mnt/lustre/users/tchibawara/MaxART/data"}else if(comp=="gent"){
    dirname <- "/user/data/gent/vsc400/vsc40070/simpact-test/data"}else{
    dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

all.sim.start <- proc.time()

set.new.seed <- 1
init.design.points <- 1 #set the initial design points

###### Generate the input parameters for the simulation ###############################################
inPUT.df.complete <- simpact.config.inputs(design.points = init.design.points, resample.count = 1,
                                           conception.alpha_base = c(-17.79167, -17.79167), #c(-4, -1.5)
                                           person.art.accept.threshold.dist.fixed.value = c(-6.301956, -6.301956), #good
                                           person.eagerness.man.dist.gamma.a = c(1.958693, 1.958693), #good
                                           person.eagerness.man.dist.gamma.b = c(-237.3288, -237.3288), #c(10,50)
                                           person.eagerness.woman.dist.gamma.a = c(9.434333, 9.434333), #c(0.3, 1.5)
                                           formation.hazard.agegapry.eagerness_diff =  c(-0.8333111, -0.8333111), #c(-0.1, 0),
                                           person.eagerness.woman.dist.gamma.b = c(-126.5245, -126.5245), #good
                                           formation.hazard.agegapry.numrel_man = c(26.39768, 26.39768), #c(-2, -0.1)
                                           formation.hazard.agegapry.numrel_woman = c(75.05623, 75.05623), #c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_man_exp = c(111.7771, 111.7771),#c(-2, -0.1),
                                           formation.hazard.agegapry.gap_factor_woman_exp = c(43.89739, 43.89739),#c(-2, -0.1),
                                           person.agegap.man.dist.normal.mu = c(460.1552, 460.1552), #c(1, 5),
                                           person.agegap.woman.dist.normal.mu = c(396.9776, 396.9776), #c(1, 5),
                                           person.agegap.man.dist.normal.sigma = c(507.0523, 507.0523), #c(0.5, 3.5),
                                           person.agegap.woman.dist.normal.sigma = c( 45.31334,  45.31334), #c(0.5, 3.5),
                                           hivtransmission.param.f1 = c(34.84124, 34.84124)
                                           )
############### Argument sample above ##############################

#Select a chunk to send to process
min.chunk <- 1
max.chunk <- 3

if(max.chunk > nrow(inPUT.df.complete)){max.chunk <- nrow(inPUT.df.complete)}
if(min.chunk > nrow(inPUT.df.complete) || min.chunk < 1){min.chunk <- max.chunk}

inANDout.df.chunk <- inPUT.df.complete[min.chunk:max.chunk,]

#make sure there are no empty rows
inANDout.df.chunk <- inANDout.df.chunk[!is.na(inANDout.df.chunk$sim.id),]

#set how many time the single row will be repeated
sim_repeat <- 1

# number of cores per node
ncluster.use <- 1

#indicate the target statitics that you want to hit
target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "inc.men.30.35", "inc.wom.30.35",
                      "prev.men.18.20", "prev.wom.18.20", "prev.men.25.30", "prev.wom.25.30","prev.men.35.40",
                      "prev.wom.35.40", "ART.cov.men.18.50", "ART.cov.wom.18.50", "median.wom.18.50.AD")

##Each of these should be calculated after each run, else we give an NA

#set the prior names - varied parameters
preprior.chunk <- names(dplyr::select(inANDout.df.chunk, contains(".")))
preprior.names.chunk <- preprior.chunk[2:length(preprior.chunk)]

#rbind all the results for this chunk to be merged after
#Create a dataframe with NA for the summary statistics Will collect all the chunks with the sim.id to link back
chunk.summary.stats.df <- data.frame(matrix(NA, nrow = 0, ncol = length(target.variables)+2))
names(chunk.summary.stats.df) <- c(target.variables, "sim.id")

############   MAIN Simulation is here #######################

simpact4ABC.chunk.wrapper <- function(simpact.chunk.prior){

  #This needs to be read by each processor
  pacman::p_load(RSimpactHelper)
  target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "inc.men.30.35", "inc.wom.30.35",
                        "prev.men.18.20", "prev.wom.18.20", "prev.men.25.30", "prev.wom.25.30","prev.men.35.40",
                        "prev.wom.35.40", "ART.cov.men.18.50", "ART.cov.wom.18.50", "median.wom.18.50.AD")
  err.functionGEN <- function(e){
    if (length(grep("MAXEVENTS",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    if (length(grep("internal event time",e$message)) != 0)
      return(chunk.summary.stats = rep(NA,length(target.variables)))
    stop(e)
  }

  simpact.chunk.run <- function(input.chunk.params){

    pacman::p_load(RSimpactCyan, RSimpactHelper, dplyr, data.table, magrittr, exactci, tidyr)

    ## Run preprior.names.chunk and copy the results here.
    input.varied.params.plus <- c("conception.alpha_base", "person.art.accept.threshold.dist.fixed.value",
                                  "person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b",
                                  "person.eagerness.woman.dist.gamma.a", "formation.hazard.agegapry.eagerness_diff",
                                  "person.eagerness.woman.dist.gamma.b", "formation.hazard.agegapry.numrel_man",
                                  "formation.hazard.agegapry.numrel_woman", "formation.hazard.agegapry.gap_factor_man_exp",
                                  "formation.hazard.agegapry.gap_factor_woman_exp", "person.agegap.man.dist.normal.mu",
                                  "person.agegap.woman.dist.normal.mu", "person.agegap.man.dist.normal.sigma",
                                  "person.agegap.woman.dist.normal.sigma", "hivtransmission.param.f1")

    target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "inc.men.30.35", "inc.wom.30.35",
                          "prev.men.18.20", "prev.wom.18.20", "prev.men.25.30", "prev.wom.25.30","prev.men.35.40",
                          "prev.wom.35.40", "ART.cov.men.18.50", "ART.cov.wom.18.50", "median.wom.18.50.AD")

    simulation.type <- ("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
    simpact.set.simulation(simulation.type)
    agedist.chunk.data.frame <- agedistr.creator(shape = 5, scale = 65)

    #### Set input params
    ##Specifying the initially chosen values for the simulation.
    cfg.chunk <- input.params.creator(population.simtime = 40, population.numwomen = 1000, population.nummen = 1000,
                                      simulation.type = simulation.type)

    #intervention introduced See the intervention.introduced
    # Simulation starts in 1977. After 27 years (in 2004), ART is introduced.
    iv.chunk <- intervention.introduced(simulation.type = simulation.type)

    #The first parameter is set to be the seed value
    seed.chunk.id <- input.chunk.params[1]

    #set up the parameters to be varied in the model starting from 2:length of the varied params.
    j <- 1
    for (cfg.chunk.par in input.varied.params.plus){
      j <- j + 1
      assign.chunk.cfg.value <- input.chunk.params[j]
      cfg.chunk[cfg.chunk.par][[1]] <- assign.chunk.cfg.value
      #setting up a value that is depended on the other input (we can do this for many other as needed)
      if(cfg.chunk.par == "hivtransmission.param.f1"){
        f2.num <- log((1+assign.chunk.cfg.value)/2)
        f2.den <- log(assign.chunk.cfg.value)
        cfg.chunk["hivtransmission.param.f2"][[1]] <- log(f2.num / f2.den)/5
      }
    }

    ## Keep the files produced in subfolders
    generate.filename <- function(how.long){
      chars <- c(letters, LETTERS)
      paste0(sample(chars,how.long), collapse = "")
    }

    sub.dir.sim.id <- generate.filename(8)
    sub.dir.rename <- paste0("temp/",sub.dir.sim.id,"/")

    testoutput <- simpact.run(configParams = cfg.chunk,
                              destDir = sub.dir.rename,
                              agedist = agedist.chunk.data.frame,
                              intervention = iv.chunk,
                              identifierFormat = paste0("%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r_",
                                                        sub.dir.sim.id,"-"),
                              seed = seed.chunk.id)

    if(testoutput$simulationtime < cfg.chunk$population.simtime)
    {
      if (testoutput$eventsexecuted >= cfg.chunk$population.maxevents-1)
      {
        stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
      }
      else
      {
        stop("Simulation stopped prematurely, probably ran out of events")
      }
    }
    chunk.datalist.test <- readthedata(testoutput)
    #save each of the run output.
    #save(chunk.datalist.test, file = paste0("temp/","chunk.datalist.",sub.dir.sim.id,".rda"))

    #delete all the file created during the current simulation
    unlink(paste0("temp/",sub.dir.sim.id), recursive = TRUE, force = TRUE)

    if(length(chunk.datalist.test)>1){
      #get the summary statistics for each run
      end.time.wind <- unique(chunk.datalist.test$itable$population.simtime)
      growth.rate <- pop.growth.calculator(datalist = chunk.datalist.test,
                                          timewindow = c(0, timewindow.max = end.time.wind))

      inc.20.25 <- incidence.calculator(datalist = chunk.datalist.test, agegroup = c(20, 25),
                                        timewindow = c(32, 34), only.active = "No")
      inc.men.20.25 <- inc.20.25$incidence[1]
      inc.wom.20.25 <- inc.20.25$incidence[2]

      inc.30.35 <- incidence.calculator(datalist = chunk.datalist.test, agegroup = c(30, 35),
                                        timewindow = c(32, 34), only.active = "No")
      inc.men.30.35 <- inc.30.35$incidence[1]
      inc.wom.30.35 <- inc.30.35$incidence[2]

      prev.18.20 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(18, 20),
                                         timepoint = 34)
      prev.men.18.20 = prev.18.20$pointprevalence[1]
      prev.wom.18.20 = prev.18.20$pointprevalence[2]

      prev.25.30 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(25, 30),
                                         timepoint = 34)
      prev.men.25.30 = prev.25.30$pointprevalence[1]
      prev.wom.25.30 = prev.25.30$pointprevalence[2]

      prev.35.40 = prevalence.calculator(datalist = chunk.datalist.test, agegroup = c(35, 40),
                                         timepoint = 34)
      prev.men.35.40 = prev.35.40$pointprevalence[1]
      prev.wom.35.40 = prev.35.40$pointprevalence[2]

      ARTcov <- ART.coverage.calculator(datalist = chunk.datalist.test, agegroup = c(18, 50),
                                        timepoint = 34, site="All")
      ART.cov.men.18.50 <- ARTcov$ART.coverage[1]
      ART.cov.wom.18.50 <- ARTcov$ART.coverage[2]

      agemix.df <- agemix.df.maker(chunk.datalist.test)
      pattern <- pattern.modeller(dataframe = agemix.df, agegroup = c(18, 50),
                                  timepoint = 34, timewindow = 1, start = FALSE)
      median.wom.18.50.AD <- as.numeric(median(pattern[[1]]$AgeGap[pattern[[1]]$Gender == "female"]))

      ##get the summary statistics as indicated by target.variables
      out.statistic <- c(growth.rate,
                         inc.men.20.25, inc.wom.20.25, inc.men.30.35, inc.wom.30.35,
                         prev.men.18.20, prev.wom.18.20, prev.men.25.30, prev.wom.25.30,
                         prev.men.35.40, prev.wom.35.40,
                         ART.cov.men.18.50, ART.cov.wom.18.50,
                         median.wom.18.50.AD)
      ##out.test.degree <- out.statistic[[2]]
    }else{
      out.statistic <- rep(NA,length(target.variables))
      ##out.statistic.degree <- NA
    }

    chunk.summary.stats <- out.statistic

    return(chunk.summary.stats)
  }

  chunk.summary.stats <- tryCatch(simpact.chunk.run(simpact.chunk.prior),
                                  error = err.functionGEN)
}


start.chunk.time <- proc.time()
for (chunk.sim.id in inANDout.df.chunk$sim.id){

  simpact.chunk.prior = list()

  for (i in preprior.names.chunk){

    #col.index <- which(colnames(preprior.names.chunk)==i)

    prior.chunk.val <- list(c("runif",1,as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i]),
                              as.numeric(inANDout.df.chunk[inANDout.df.chunk$sim.id==chunk.sim.id,i])), c("dunif",0,1))
    simpact.chunk.prior[[length(simpact.chunk.prior)+1]] <- prior.chunk.val
  }

  print(paste("Working on simulation number: ", chunk.sim.id, sep=" "))
  #invoke the ABC_rejection method repeating the number of simulation X* for each chunk row.
  ABC.chunk.result <- ABC_rejection(model = simpact4ABC.chunk.wrapper,
                                        prior = simpact.chunk.prior,
                                        nb_simul= sim_repeat,
                                        use_seed = TRUE,
                                        seed_count = 0,
                                        n_cluster = ncluster.use)


  #Save the statistics results with the chunk row sim.id repeated X* from the ABC_rejection method
  ABC.results.chunk.statistics <- data.frame(ABC.chunk.result$stats)
  names(ABC.results.chunk.statistics) <- target.variables
  ABC.results.chunk.statistics$sim.id <- chunk.sim.id

  chunk.summary.stats.df <- rbind(chunk.summary.stats.df, ABC.results.chunk.statistics)

}

# target.stats <- c(0.015, #gr
#                   0.016, 0.043, 0.031, 0.027,  #inc
#                   0.008, 0.143, 0.21, 0.47, 0.47, 0.538, #prev
#                   0.33, 0.34, #art
#                   5) #ad

inputANDoutput.chunk.df  <- left_join(chunk.summary.stats.df, inANDout.df.chunk, by = "sim.id")

rand.string <- paste0(sample(c(LETTERS,letters), 10), collapse = "")

filename.run <- paste0(dirname,"/","SummaryOutPut-df-",min.chunk,"-",max.chunk,".csv")

write.csv(inputANDoutput.chunk.df, file = filename.run, row.names = FALSE)

end.chunk.time <- proc.time() - start.chunk.time
print(paste0("Total time to compute simulation: ", end.chunk.time))

all.sim.end <- proc.time() - all.sim.start
print(paste0("Total time to finish simulation: ", all.sim.end))








