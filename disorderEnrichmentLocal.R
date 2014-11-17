# This script will look for enrichment of disorder genes in the Allen Brain Atlas - meant
# to be run locally to assess significant findings

# I have my disorder genes - we want to know how important are they in normal brain function?
# importance == highly up or down regulated
# determined by a fisher test with the following table:
#                 regionshighlyupdown  ~regionshighlyupdown
# asdGenes
# ~asdGenes
# a significant result means that there is a high count in the top left cell

# Here is the data we need
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/AllenRegionNormalUpDownRegulatedGenes3std.Rda")
load("/home/vanessa/Documents/Work/ALLEN/AllenProbesRegionMeans.Rda")
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/3SDProbesAllenBrainRegions.Rda")
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/FinalGeneLists.Rda")

# Here is probe meta info
meta = read.csv("/home/vanessa/Documents/Work/ALLEN/Probes.csv",head=FALSE,sep=",")
# 1058685,"A_23_P20713",729,"C8G","complement component 8, gamma polypeptide",733,9
colnames(meta) = c("ID","PROBE","PLATE","GENE","DESCRIPTION","NOTSURE1","NOTSURE2")
samples = read.csv("/home/vanessa/Documents/Work/ALLEN/Samples.csv",head=TRUE,sep=",")

# Get the unique disorders with significant result
disorders = as.character(unique(enrichedRegions$DISORDER))

# First let's plot the "weird expression" genes in each region.  We want to be
# convinced, visually, that these genes are relatively up and down.
all_regions = as.character(unique(samples$structure_name))
probeMeans = colMeans(allen$probes)
probestds = apply(allen$probes,2,sd)

for (r in 1:length(all_regions)){
  region = all_regions[r]
  cat("Region",r,"of",length(all_regions),"\n")
  up = regions$upprobes[[region]]
  down = regions$downprobes[[region]]
  if (length(up)+length(down) > 0) {
    # Let's create a vector of the mean expression vs. the actual expression values
    actual_expression_up = allen$probes[region,up]
    actual_expression_down = allen$probes[region,down]
    actual_expression_all = allen$probes[region,c(up,down)]
    # Now get the mean values
    if (length(actual_expression_all) > 1) {
      mean_expression_values = probeMeans[names(actual_expression_all)]
      sd_expression_values = probestds[names(actual_expression_all)]
      # Change name of probes to gene names
      names(actual_expression_all) = as.character(meta$GENE[which(meta$ID %in% names(actual_expression_all))])
    } else {
      mean_expression_values = probeMeans[c(up,down)]
      sd_expression_values = probestds[c(up,down)]
      # Change name of probes to gene names
      names(actual_expression_all) = as.character(meta$GENE[which(meta$ID %in% c(up,down))])
    }
    # calculate differences between actual and mean
    differences = array(dim=length(actual_expression_all))
    higher = which(actual_expression_all >= mean_expression_values)
    lower = which(actual_expression_all < mean_expression_values)
    differences[higher] = actual_expression_all[higher] - mean_expression_values[higher]
    differences[lower] = -1*(mean_expression_values[lower] - actual_expression_all[lower])

    # Now get probes that aren't weird'
    normal_probes = colnames(allen$probes)[-which(colnames(allen$probes) %in% c(up,down))]
    normal_probe_expression = allen$probes[region,normal_probes]
    normal_probe_means = probeMeans[normal_probes]
    normal_probe_sds = probestds[normal_probes]

    # Take a random sample of 100
    idx = sample(seq(1,length(normal_probe_means)),100)
    mean_normal_probe_expression = normal_probe_means[idx]
    sd_normal_probe_expression = normal_probe_sds[idx]
    normal_probe_expression = normal_probe_expression[idx]
    normal_probe_expression_difference = normal_probe_expression - mean_normal_probe_expression
    # Create a vector of colors

    colors = c(rep("green",length(actual_expression_up)),rep("blue",length(actual_expression_down)),rep("orange",length(idx)))

    differences = c(differences,normal_probe_expression_difference)
    means = c(mean_expression_values,mean_normal_probe_expression)
    sds = c(sd_expression_values,sd_normal_probe_expression)
    actual = c(actual_expression_all,normal_probe_expression)
    df = data.frame(colors=colors,differences=differences,means = means, sd = sds,actual=actual,stringsAsFactors=FALSE)
    df = df[with(df, order(sd)),]

    # Now make a barplot to compare the two
    outimg = paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/regionProbeImg/",gsub(" ","_",all_regions[r]),".png",sep="")
    png(filename = outimg,width = 900, height = 480, units = "px", pointsize = 12)
    bp = barplot(df$differences,ylim=c(min(df$differences),max(3*df$sd)),main=paste("Interesting probes for region",all_regions[r]),col=df$colors,las=2,xlab="gene probes",ylab="normalized expression differences")
    lines(x=bp,y=3*df$sd,col="red")
    lines(x=bp,y=-3*df$sd,col="red")
    legend(50,3, c("3 standard deviations > mean","3 standard deviations < mean","random sample N=100","three standard deviations"),lty=c(1,1),lwd=c(2.5,2.5),col=c("green","blue","orange","red"))
    dev.off()
  }
}

# NOW we want to be convinced that, for each signficant result, our disorder genes
# are represented at a rate that is actually greater than chance.

# Load file with significant results
load("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/result/disorder16_sigresult_fdr05.Rda")

uniquedisorders = unique(enrichedRegions$DISORDER)

outfile = "/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/result/fishertest.txt"
for (d in 1:length(uniquedisorders)){
  disorder = uniquedisorders[d]
  disorderregions = enrichedRegions$REGION[which(enrichedRegions$DISORDER==disorder)]
  
  # Get the genes
  tmp = genes[[disorder]]
  gen = c()
  for (g in tmp){
    gen = c(gen,strsplit(g,"[(]")[[1]][1])
  }

  for (r in 1:length(disorderregions)){
    region = disorderregions[r]
    #cat("Region",r,"of",length(disorderregions),"\n")
    up = regions$upprobes[[region]]
    down = regions$downprobes[[region]]
    if (length(both) > 0) {
      # Now get the gene names
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

      tabley = matrix(c(a,c,b,d),nrow=2,dimnames=list(c("DisorderGenes","~DisorderGenes"),c("WeirdExpression","~WeirdExpression")))

      # Determine the direction
      inregion = a/(a+b)
      outregion = c/(c+d)

      if (inregion > outregion){
        # Print to screen
        cat("\n",disorder,":",region,"\n",file=outfile,append=TRUE)
        write.table(tabley,file=outfile,append=TRUE)
        cat(file=outfile,"Proportions are",inregion,"(disorder genes that are in weird) vs",outregion,"(disorder genes that are in normal)\n",append=TRUE)
      }
   }
  }
}

# STOPPED HERE - need to now LOOK at actual genes in the overlap list!

      # Plot percentage of disorder genes 

    ft = fisher.test(tabley, conf.level = 0.95)
    if (ft$p.value <= .1) {
      cat(file=resultfile"Significant result for",disorder,": in region ",u,": p value",ft$p.value,"\n",append=TRUE)
      cat(file=resultfile,disorder,": in region ",u,": p value ",ft$p.value,"\n",append=TRUE)
      cat("Proportions are")
    }
  }
}








      # Now get the up and down probes
      regions$upprobes
    # calculate differences between actual and mean
    differences = array(dim=length(actual_expression_all))
    higher = which(actual_expression_all >= mean_expression_values)
    lower = which(actual_expression_all < mean_expression_values)
    differences[higher] = actual_expression_all[higher] - mean_expression_values[higher]
    differences[lower] = -1*(mean_expression_values[lower] - actual_expression_all[lower])

    # Now get probes that aren't weird'
    normal_probes = colnames(allen$probes)[-which(colnames(allen$probes) %in% c(up,down))]
    normal_probe_expression = allen$probes[region,normal_probes]
    normal_probe_means = probeMeans[normal_probes]
    normal_probe_sds = probestds[normal_probes]

    # Take a random sample of 100
    idx = sample(seq(1,length(normal_probe_means)),100)
    mean_normal_probe_expression = normal_probe_means[idx]
    sd_normal_probe_expression = normal_probe_sds[idx]
    normal_probe_expression = normal_probe_expression[idx]
    normal_probe_expression_difference = normal_probe_expression - mean_normal_probe_expression
    # Create a vector of colors

    colors = c(rep("green",length(actual_expression_up)),rep("blue",length(actual_expression_down)),rep("orange",length(idx)))

    differences = c(differences,normal_probe_expression_difference)
    means = c(mean_expression_values,mean_normal_probe_expression)
    sds = c(sd_expression_values,sd_normal_probe_expression)
    actual = c(actual_expression_all,normal_probe_expression)
    df = data.frame(colors=colors,differences=differences,means = means, sd = sds,actual=actual,stringsAsFactors=FALSE)
    df = df[with(df, order(sd)),]

    # Now make a barplot to compare the two
    outimg = paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/disorderEnrichment/regionProbeImg/",gsub(" ","_",all_regions[r]),".png",sep="")
    png(filename = outimg,width = 900, height = 480, units = "px", pointsize = 12)
    bp = barplot(df$differences,ylim=c(min(df$differences),max(3*df$sd)),main=paste("Interesting probes for region",all_regions[r]),col=df$colors,las=2,xlab="gene probes",ylab="normalized expression differences")
    lines(x=bp,y=3*df$sd,col="red")
    lines(x=bp,y=-3*df$sd,col="red")
    legend(50,3, c("3 standard deviations > mean","3 standard deviations < mean","random sample N=100","three standard deviations"),lty=c(1,1),lwd=c(2.5,2.5),col=c("green","blue","orange","red"))
    dev.off()
  }
}






# Let's quickly try this anyway
# I have my ASD genes - how important are they in normal brain function?
# importance == highly up or down regulated across all regions

# DISORDER X
#                 regionshighlyupdown  ~regionshighlyupdown
# asdGenes
# ~asdGenes

for (g in 1:length(genes)){
  gen = genes[[g]]
  disorder = names(genes[g])
  # First let's look at the disorder genes
  reg = c()
  for (gg in 1:length(gen)){
    gene = names(gen[gg])
    pids = as.character(meta$ID[which(meta$GENE == gene)])
    # Look up up and down regions for those probes
    # We will count unique regions across probes
    for (p in pids){
      reg = c(reg,names(upreg[[p]]),names(upreg[[p]]))
    }
  }
  # Weird expression and disorder gene
  a = length(unique(reg))
  # Not weird expression and disorder gene
  b = 3702 - length(unique(reg))
  
  # Now let's look at non-disorder genes
  reg = c()
  notdisorder = as.character(unique(meta$GENE[-which(meta$GENE %in% names(gen))]))
  for (gg in 1:length(notdisorder)){
    gene = notdisorder[gg]
    pids = as.character(meta$ID[which(meta$GENE == gene)])
    # Look up up and down regions for those probes
    # We will count unique regions across probes
    for (p in pids){
      reg = c(reg,names(upreg[[p]]),names(upreg[[p]]))
    }
  }
  # Weird expression and disorder gene
  c = length(unique(reg))
  # Not weird expression and disorder gene
  d = 3702 - length(unique(reg))
}

  pids = as.character(meta$PROBE[which(meta$GENE == uniquegenes[g])])
  # Count number of regions with weird expression
  a = 0  # weirdExpression and disorderGene
  b = 0  # ~weirdExpression and disorderGene
  c = 0  # weirdExpression and nondisorderGene
  d = 0  # ~weirdExpression and disorderGene
  for (r in )

}
# Count number of regions with weird expression

# Count number of regions with normal expression

# 


#                 regionshighlyupdown  ~regionshighlyupdown
# asdGenes
# ~asdGenes
fishers = list()
for (g in 1:length(genes)){
  disorder = names(genes)[g]
  gen = genes[[g]]
  a =   # weirdExpression and disorderGene
  b =   # ~weirdExpression and disorderGene
  c =   # weirdExpression and nondisorderGene
  d =   # ~weirdExpression and disorderGene
}







# IDEA 2: Try to find enrichment of genes at a particular timepoint
# We will be using rna genes from BrainSpan, with 52K genes for 525 individuals across different timepoints
# We have rna seq and microarray, but I read rnaseq is better
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/rna-genes-52375x525.Rda")
# Get unique timepoints
timepoints = as.character(unique(rna.genes$columns$age))



# Are disorder X genes over-represented in the set?
# Let's look on a regional basis first

# I think we can do fisher's exact test'
#                 weirdExpression ~weirdExpression
# disorderGene
# nondisorderGene
fishers = list()
for (g in 1:length(genes)){
  disorder = names(genes)[g]
  gen = genes[[g]]
  a =   # weirdExpression and disorderGene
  b =   # ~weirdExpression and disorderGene
  c =   # weirdExpression and nondisorderGene
  d =   # ~weirdExpression and disorderGene
}



# For each, create 2x2 table, do fisher's exact test
# This doesn't take into account the value, talk about this'

# Now we want to determine if disorder genes are overrepresented in the set
# Let's use the hypergeometric distribution
# The p-value you want is the probability of getting 100 or more white balls in a sample of size 400 from an urn with 3000 white balls and 12000 black balls. Here are four ways to calculate it.


# GeneA --> regionX,Y,Z
# GeneB --> regionX,Y,Z


# Is region X over-represented in list?
# Disorder X     regionsUpDown   ~regionsUpDown
# disorderGene
# ~disorderGene




sum(dhyper(100:400, 3000, 12000, 400))
1 - sum(dhyper(0:99, 3000, 12000, 400))
phyper(99, 3000, 12000, 400, lower.tail=FALSE)
1-phyper(99, 3000, 12000, 400)
These give 0.0078.

# dhyper(x, m, n, k) gives the probability of drawing exactly x. In the first line, we sum up the probabilities for 100 – 400; in the second line, we take 1 minus the sum of the probabilities of 0 – 99.

# phyper(x, m, n, k) gives the probability of getting x or fewer, so phyper(x, m, n, k) is the same as sum(dhyper(0:x, m, n, k)).

# The lower.tail=FALSE is a bit confusing.  phyper(x, m, n, k, lower.tail=FALSE) is the same as 1-phyper(x, m, n, k), and so is the probability of x+1 or more. [I never remember this and so always have to double check.]



