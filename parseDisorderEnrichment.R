# This script will parse the disorderEnrichment results, and
# calculate corrected p values (in two ways - correcting for each disorder,
# and corrected for ALL tests.)

# disorderEnrichment Parsing
setwd("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/enrichedRegions/")
files = list.files(pattern="_regionEnrichment.Rda")
raw = c()

# Save all Pvalues
for (f in 1:length(files)){
  cat("Loading file",f,"of",length(files),"\n")
  load(files[f])
  raw = rbind(raw,fisher)
}

FDR = p.adjust(raw$PVALUE)
raw = cbind(raw,FDR)

save(raw,file="/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/results/raw65Disorders.Rda")

sig = raw[which(raw$FDR<=0.05),]
write.table(sig,file="/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/results/sigResultFDR001.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)
save(sig,file="/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/results/sigResultFDR001.Rda")

# How many unique disorders?
unique(sig$DISORDER)

# ASSIGN TO SAMPLES - what regions of the brain?
# TABLE:
# DISORDER	PVALUE	CORRECTED	REGIONID	REGIONINFO ...
# Now load the region sample info, and create result file with significant results
samples = read.csv("/scratch/users/vsochat/DATA/ALLEN/Samples.csv",head=TRUE,sep=",")

# Let's save the output files here'
outdir = "/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/results/"

# For each disorder, we should save a table of regions
uniquedisorder = unique(sig$DISORDER)
enrichedRegions = c()
for (s in 1:length(uniquedisorder)){
  cat("Processing",s,"of",length(uniquedisorder),"\n")
  disorder = uniquedisorder[s]
  tmp = sig[which(sig$DISORDER == disorder),]
  # Now grab the region IDs
  sampleSubset = samples[as.numeric(tmp$SAMPLEID),]
  # Now append the p value to each sampleID
  enrichedRegions = rbind(enrichedRegions,cbind(sampleSubset,tmp))
}
save(enrichedRegions,file=paste(outdir,"disorder16_sigresult_fdr05.Rda",sep=""))

# For each disorder, we should now make a brain map that shows where the enrichment is
setwd("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/results")
load("disorder16_sigresult_fdr05.Rda")
library("scatterplot3d")

# For each disorder
uniquedisorders = as.character(unique(enrichedRegions$DISORDER))
uniqueregions = as.character(unique(enrichedRegions$structure_name))

# First let's make a matrix of disorder by region
regionmatrix = array(0,dim=c(length(uniquedisorders),length(uniqueregions)))
colnames(regionmatrix) = sort(uniqueregions)
rownames(regionmatrix) = sort(uniquedisorders)

# Parse the result
for (r in 1:nrow(enrichedRegions)){
  cat("Processing",r,"of",nrow(enrichedRegions),"\n")
  region = as.character(enrichedRegions[r,"structure_name"])
  disorder = as.character(enrichedRegions[r,"DISORDER"])
  regionmatrix[disorder,region] = 1
}

save(regionmatrix,file="Disorder16_sigresult_fdr05_regionmatrix.Rda")

# Cluster disorders based on the regionmatrix
disty = dist(regionmatrix)
hc = hclust(disty)

pdf(file=paste('Disorder16_sigresult_fdr05_Brainmaps.pdf',sep=""),onefile=TRUE)
plot(hc,main="Clustering Disorders Based on Similar Enriched Regions",xaxt="n",sub="")
for (s in 1:length(uniquedisorders)){
  cat("Processing",s,"of",nrow(enrichedRegions),"\n")
  disorder = uniquedisorders[s]
  tmp = enrichedRegions[which(enrichedRegions$DISORDER==disorder),]
  # Now let's visualize the coordinates
  coords = cbind(tmp$mni_x,tmp$mni_y,tmp$mni_z)
  scd = scatterplot3d(coords,pch=19,color=sample(colors(),1),main=paste("Enriched Regions for",disorder),xlab="MNIx",ylab="MNIy",zlab="MNIz",ylim = c(-100,70),xlim = c(-70,70),zlim = c(-70,80))
  # We need to convert to 2d point labels
  twodcoords = scd$xyz.convert(coords)
  text(twodcoords$x,twodcoords$y,labels=tmp$structure_name,cex=.5, pos=4)  
}
dev.off()  # Close file

# Now let's look at the disorders for each brain region
uniqueregions = unique(enrichedRegions$REGION)
filey="result/regionDisorderList.txt"
for (u in uniqueregions){
  cat(file=filey,u,"\n",append=TRUE)
  tmp = enrichedRegions$DISORDER[which(enrichedRegions$REGION == u)]
  cat(file=filey,tmp,sep="\n",append=TRUE)
  cat(file=filey,"\n",append=TRUE)
}
