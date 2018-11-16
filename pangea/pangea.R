
library(PANGEA.HIV.sim)

outdir          <- getwd()                                          #set to a new empty directory; dir name must not contain whitespace, brackets, etc

pipeline.args	<- sim.regional.args( 	seed=42,                    #random number seed for reproducibility
                                     yr.end=2020,				#end of simulation
                                     s.PREV.max.n=1600,          #number of sequences
                                     s.INTERVENTION.prop=0.5,    #proportion of sampled sequences after intervention start in 2015
                                     epi.acute='high',           #frequency of early infections (high or low)
                                     epi.intervention='fast',    #intervention scale-up (none, slow or high)
                                     epi.import=0.05 )			#proportion of transmissions from outside the regional population


cat(sim.regional(outdir, pipeline.args=pipeline.args))              #run this script from the command line
