#######################################################
#######################################################
#######################################################
#
# start: run sim.regional
#
#######################################################
#######################################################
#######################################################
mkdir -p /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim
cd /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim
#######################################################
# start: run HPTN071.simulator
####################################################### 
echo 'run popart-highacute'
cp -R /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/ext/PangeaParamsHighAcute PangeaParamsHighAcute 
 /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/ext/popart-highacute PangeaParamsHighAcute 42 1 1
mv phylogenetic_individualdata* /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim/150129_HPTN071_scHF_IND.csv
mv phylogenetic_transmission* /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim/150129_HPTN071_scHF_TRM.csv
mv Annual_outputs* /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim/150129_HPTN071_scHF_EPI.csv
rm -rf PangeaParamsHighAcute 
 echo 'end popart-highacute'
#######################################################
# end: run HPTN071.simulator
#######################################################
cd ..
mkdir -p /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains
#######################################################
# start: run HPTN071.input.parser.v4
####################################################### 
echo 'run Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/HPTN071.input.parser.v4.Rscript'
 Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/HPTN071.input.parser.v4.Rscript -indir=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim -infile.args=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/150129_HPTN071_scHF_PipeArgs.R -outdir=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains -infile.ind=150129_HPTN071_scHF_IND.csv -infile.trm=150129_HPTN071_scHF_TRM.csv -outfile.ind=150129_HPTN071_scHF_IND.csv -outfile.trm=150129_HPTN071_scHF_TRM.csv 
echo 'end Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/HPTN071.input.parser.v4.Rscript'
#######################################################
# end: run HPTN071.input.parser.v4
#######################################################
mkdir -p /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/VirusTreeSimulator
#######################################################
# start: run VirusTreeSimulator
####################################################### 
echo 'run /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/ext/VirusTreeSimulator.jar'
 java -Xms64m -Xmx400m -jar /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/ext/VirusTreeSimulator.jar -seed 42 -demoModel Logistic -N0 1 -growthRate 2.851904 -t50 -2  /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains/150129_HPTN071_scHF_TRM.csv /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains/150129_HPTN071_scHF_IND.csv /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/VirusTreeSimulator/150129_HPTN071_scHF_
find /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/VirusTreeSimulator -name "*simple*" -delete
 echo 'end /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/ext/VirusTreeSimulator.jar'
#######################################################
# end: run VirusTreeSimulator
#######################################################
mkdir -p /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/SeqGen
#######################################################
# start: run SeqGen.createInputFile 
####################################################### 
echo 'run Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.createInputFile.Rscript'
 Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.createInputFile.Rscript -indir.epi=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains -infile.epi=150129_HPTN071_scHF_SAVE.R -indir.vts=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/VirusTreeSimulator -infile.vts=150129_HPTN071_scHF_ -infile.args=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/150129_HPTN071_scHF_PipeArgs.R -outdir=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/SeqGen 
 echo 'end Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.createInputFile.Rscript'
#######################################################
# end: run SeqGen.createInputFile
#######################################################
#######################################################
# start: run SeqGen.run 
####################################################### 
echo 'run Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.run.v4.Rscript'
 Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.run.v4.Rscript -indir.epi=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains -infile.epi=150129_HPTN071_scHF_SAVE.R -indir.sg=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/SeqGen -infile.sg=150129_HPTN071_scHF_ -infile.args=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea/150129_HPTN071_scHF_PipeArgs.R -outdir=/home/niyukuri/phylosimpact_simulation_studies_2018/pangea 
 echo 'end Rscript /home/niyukuri/R/x86_64-pc-linux-gnu-library/3.4/PANGEA.HIV.sim/PANGEA.SeqGen.run.v4.Rscript'
#######################################################
# end: run SeqGen.run
#######################################################
rm -rf /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/EpiSim /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/TrChains /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/VirusTreeSimulator /home/niyukuri/phylosimpact_simulation_studies_2018/pangea/SeqGen
#######################################################
#######################################################
#######################################################
#
# end: run sim.regional
#
#######################################################
#######################################################
#######################################################
