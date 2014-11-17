# Set up script

setwd("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder")
load("FinalGeneLists.Rda")
outdir = ("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/enrichedRegions/")
setwd("/home/vsochat/SCRIPT/R/DisorderEnrichment")

# Submission down here
for (f in 1:length(genes)){
  disorder = names(genes[f])
  dis = gsub(" ","_",strsplit(disorder,"[(]")[[1]][1])
  # Output file with fisher test results
  outfile = paste(outdir,dis,"_regionEnrichment.Rda",sep="")
  if (!file.exists(outfile)) {
    cat(disorder,"\n")
    jobby = paste(dis,".job",sep="")
    sink(paste(".job/",jobby,sep=""))
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=",jobby,"\n",sep="")  
    cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
    cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
    cat("#SBATCH --time=0-03:00\n",sep="")
    cat("#SBATCH --mem=8000\n",sep="")
    cat("Rscript /home/vsochat/SCRIPT/R/DisorderEnrichment/disorderEnrichment.R",dis,outfile,"\n")
    sink()
  
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    system(paste("sbatch -p dpwall ",paste(".job/",jobby,sep="")))
  }
}


