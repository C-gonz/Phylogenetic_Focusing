### Custom function for phylogenetic focusing; written by Curtis Provencher for Opsins, modified for general use by Chris Gonzalez.
### Inputs all the ML gene trees, roots them with the MRCA of the anchor sequences specified by root. Next the clade containing the MRCA of the genes of interest and their outgroup
### if identified and pruned off for further analysis. This is done for each taxa

# Libriares needed; remove comments if installation is required.
#install.packages("ape", dependencies = TRUE)
#install.packages("phytools", dependencies = TRUE)

library(ape)
library(phytools)
#wrap functions to perfom on each tree into one big function
dostuff<-function(sometree) {
  treename<-unlist(strsplit(sometree,split='.treefile'))[1]
  mytree <- read.tree(sometree)
  
  #############################
  ######1. root with Most Recent Common Ancestor (MRCA) of Taste Receptors
  #############################
  # setting the mrca of anchors and using as root (this sets root as the node # of mrca)
  root <- findMRCA(mytree, tips = c("GABBR1_Hsapiens", "GABBR2_Hsapiens", "GABBR1_Drerio", "GABBR2_Drerio"), type = "node")
  
  # now root tree  
  rootedtree <- root(mytree, node = root)
  # look to check
  plot.phylo(rootedtree)
  title(paste(treename,': rooted'))
  
  #############################
  ######2. Focusing step to obtain specific node...
  #############################
  # this will work to get node number of MRCA for 2+ tips
  taste_receptorsclade0 <- findMRCA(rootedtree, tips = c('T1R1_Hsapiens', 'mGLUR6_Hsapiens'), type="node")
  
  # extract that clade and make a new tree
  taste_receptorstree <- extract.clade(rootedtree, taste_receptorsclade0)
  # visual check (still not sure which clade to pull)
  plot.phylo(taste_receptorstree)
  title(paste(treename,': taste_receptors'))
  
  #####################
  #####3. get names of tips in the focused tree
  ####################
  write (taste_receptorstree$tip.label, file=paste(treename, "_tips.txt"))
  ## export each file with "GenusName_tips.txt"
  genera<-unlist(strsplit(taste_receptorstree$tip.label,split="_"))[1]
}

#specify location of trees; this should be the absolute path to the IQ_out/ directory generated in step 6
#mydir<-('/Absolute_path_here')
mydir<-('/mnt/lustre/plachetzki/cjgonzalez/PhyFocus/final_T1R_phyfocus/align/IQ_out')
getwd()

mytrees<-dir(mydir, pattern='.treefile')

#run function on specified trees
for (t in 1:length(mytrees)) {
  dostuff(sometree=mytrees[t])
}

