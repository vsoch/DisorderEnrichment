# Get command args
args <- commandArgs(TRUE)
disorder = args[1]
resultfile = args[2]

load("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/AllenProbesRegionMeans.Rda")
load("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/3SDProbesAllenBrainRegions.Rda")
load("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/AllenRegionNormalUpDownRegulatedGenes3std.Rda")
load("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/FinalGeneLists.Rda")

# Here is probe meta info
meta = read.csv("/scratch/users/vsochat/DATA/ALLEN/Probes.csv",head=FALSE,sep=",")
# 1058685,"A_23_P20713",729,"C8G","complement component 8, gamma polypeptide",733,9
colnames(meta) = c("ID","PROBE","PLATE","GENE","DESCRIPTION","NOTSURE1","NOTSURE2")
# We also need the sample information
samples = read.csv("/scratch/users/vsochat/DATA/ALLEN/Samples.csv",head=TRUE,sep=",")

# REGION ENRICHMENT - I need to do this for each region - for each disorder
# Save output file, and after we will correct p values for multiple testing

# Find the disorder in the list
tmp = names(genes)
dis = c()
for (d in tmp){
  dis = c(dis,gsub(" ","_",strsplit(d,"[(]")[[1]][1]))
}

savename = disorder
disorder = names(genes)[which(dis==disorder)]

# Remove numerical id
tmp = genes[[disorder]]
gen = c()
for (g in tmp){
  gen = c(gen,strsplit(g,"[(]")[[1]][1])
}

# Filter down to genes that overlap
gen = gen[which(gen %in% meta$GENE)]
# Save the list of genes in case we need
#save(gen,file=paste("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/genelists/",savename,".Rda",sep=""))

fisher = array(dim=c(length(regions$upprobes),4))
# For each region,  
for (u in 1:length(regions$upprobes)){
  region = names(regions$upprobes[u])
  cat("Region",u,"of",length(regions$upprobes),"\n")
  up = regions$upprobes[[region]]
  down = regions$downprobes[[region]]
  # Get gene names for up
  geneup = as.character(unique(meta$GENE[which(meta$ID %in% up)]))
  genedown = as.character(unique(meta$GENE[which(meta$ID %in% down)]))
  weirdexpression = as.character(unique(c(geneup,genedown)))
  normalexpression = as.character(unique(meta$GENE[-which(meta$GENE %in% weirdexpression)]))

  # Make a 2x2 table for the regions/disorder
  # Weird expression and disorder gene
  a = length(which(gen %in% weirdexpression))
  # Not weird expression and disorder gene
  b = length(which(gen %in% normalexpression))
  # Weird expression and non-disorder gene
  c = length(weirdexpression) - a
  # Not weird expression and non-disorder gene
  d = length(normalexpression) - b

  tabley = matrix(c(a,b,c,d),nrow=2,dimnames=list(c("DisorderGenes","~DisorderGenes"),c("WeirdExpression","~WeirdExpression")))

  ft = fisher.test(tabley, conf.level = 0.95)
  fisher[u,] = c(disorder,u,ft$p.value,region)
}
fisher = as.data.frame(fisher,stringsAsFactors=FALSE)
colnames(fisher) = c("DISORDER","SAMPLEID","PVALUE","REGION")
fisher$PVALUE = as.numeric(fisher$PVALUE)
# Save fisher result to outfile
save(fisher,file=resultfile)
