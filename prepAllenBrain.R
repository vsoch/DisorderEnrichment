# This script will prepare the data objects with up and down regulated probes for each region

# Load the Allen Brain Atlas data
# This is probes (columns) by samples (rows)
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/AllenBrain58692.Rda")
colnames(allen$probes) = gsub("pid_","",colnames(allen$probes))

# Here is probe meta info
meta = read.csv("/scratch/users/vsochat/DATA/ALLEN/Probes.csv",head=FALSE,sep=",")
# 1058685,"A_23_P20713",729,"C8G","complement component 8, gamma polypeptide",733,9
colnames(meta) = c("ID","PROBE","PLATE","GENE","DESCRIPTION","NOTSURE1","NOTSURE2")

# We also need the sample information
samples = read.csv("/scratch/users/vsochat/DATA/ALLEN/Samples.csv",head=TRUE,sep=",")

# Take mean across locations
uniqueregions = unique(samples$structure_name)
probes = array(dim=c(length(uniqueregions),ncol(allen$probes)))
colnames(probes) = colnames(allen$probes)
rownames(probes) = uniqueregions

for (r in 1:length(uniqueregions)){
  cat("processing",r,"of",length(uniqueregions),"\n")
  region = uniqueregions[r]
  idx = which(samples$structure_name == region)
  if (length(idx)!=1){
    probes[region,] = colMeans(allen$probes[idx,])
  } else {
    probes[region,] = allen$probes[idx,]
  }
}
allen = list(probes=probes,samples=samples,meta=meta)
save(allen,file="AllenProbesRegionMeans.Rda")

# For each region, extract top up and down genes (later we will want ranking)
regionsup = list()
regionsdown = list()
data = probes
means = colMeans(data)
stds = apply(data,2,sd)
upperlimit = means + 3*stds
lowerlimit = means - 3*stds
for (g in 1:nrow(probes)){
  region = rownames(probes)[g]
  cat("Processing",g,"of",nrow(probes),"\n")
  # We want to find probes that are 3sd greater than their means
  upset = which(data[g,]>=upperlimit)
  downset = which(data[g,]<=lowerlimit)
  regionsup[[region]] = names(upset)
  regionsdown[[region]] = names(downset)
}

readme = "probes in list are <>= 3sd the mean expression for the probe across the entire brain"
regions = list(upprobes=regionsup,downprobes=regionsdown,regions=samples,readme=readme)
save(regions,file="3SDProbesAllenBrainRegions.Rda")

# In the Allen Brain Atlas, for each gene save a list of upregulated regions
upreg = list()
# And down regulated regions
downreg = list()
for (g in 1:ncol(probes)){
  cat("Processing",g,"of",ncol(probes),"\n")
  gen = colnames(probes)[g]
  # Calculate the mean and standard deviation
  meany = mean(probes[,g])
  stdy = sd(probes[,g])
  # Which probes are greater than 3std?
  upregulated = probes[which(probes[,g]>= meany+3*stdy),g]  
  # Which regions are less than 3std?
  downregulated = probes[which(probes[,g]<= meany-3*stdy),g]  
  upreg[[gen]] = upregulated
  downreg[[gen]] = downregulated
}

# Now let's make a matrix of complete region info for each
regionup = list()
regiondown = list()
for (g in 1:length(upreg)){
  cat("Processing",g,"of",length(upreg),"\n")
  gen = names(upreg)[g]
  up = upreg[[g]]
  down = downreg[[g]]
  # Look up probe info
  upidx = names(up)
  downidx = names(down)
  upmeta = samples[which(samples$structure_name==upidx),]
  downmeta = samples[which(samples$structure_name==downidx),]
  regionup[[gen]] = upmeta
  regiondown[[gen]] = downmeta
}

regulation = list(upreg=upreg,downreg=downreg,regionup=regionup,regiondown=regiondown,samples=samples,probes=meta,rawdata=probes)
save(regulation,file="AllenRegionNormalUpDownRegulatedGenes3std.Rda")

