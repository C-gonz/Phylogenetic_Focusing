#!/usr/bin/env Rscript

# Custom function for phylogenetic focusing; template by Curtis Provencher for Opsins, rewritten for general use by Chris Gonzalez.

# For each species' ML gene tree, roots them with MRCA of root sequences, then extracts the clade specified by the MRCA of the focus sequences. Tip labels for these focused trees are then recorded for the concatenated tree of phyfocus.

# Load required libaries and install if absent.
ape_installed <- require("ape")
if (ape_installed == FALSE) {
install.packages("ape", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library("ape")
}

phytools_installed <- require("phytools")
if (phytools_installed == FALSE) {
install.packages("phytools", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library("phytools")
}

stri_installed <- require("stringi")
if (stri_installed == FALSE) {
  install.packages("stringi", dependencies = TRUE, repos = "http://cran.us.r-project.org")
  library("stringi")
}

# Function to read sequence table, extract labels, and return in a list
read_func <- function(clade_file) {
  data <- read.delim(clade_file, header=FALSE, sep = ",")
  # Extract each column of the table as a vector of taxa labels
  column_1 <- data[1]
  column_1_data <- (column_1[ , ])
  clean_column_1 <- stri_remove_empty(column_1_data, na_empty = FALSE)
  column_2 <- data[2]
  column_2_data <- column_2[ , ]
  clean_column_2 <- stri_remove_empty(column_2_data, na_empty = FALSE)
  # Return a list of 2 string vectors:
  read_list <- list(clean_column_1, clean_column_2)
  return(read_list)
}

# Function to perform on each species' gene tree; entails the focusing step of PhyFocus. Note: findMRCA requires "tips" to contain at least 2 labels.
focus_func <- function(iq_tree, taxa) {
  treename <- unlist(strsplit(iq_tree,split='.treefile'))[1]
  phylo_tree <- read.tree(iq_tree)

  # 1. root phylogeny with Most Recent Common Ancestor (MRCA).
  ## Extract rooting taxa from the list
  root_list <- taxa[1]
  unlisted_roots <- unlist(root_list)
  root <- findMRCA(phylo_tree, tips = unlisted_roots, type="node")
  ## root phylogeny and plot for records
  rootedtree <- root(phylo_tree, node = root)
  plot.phylo(rootedtree)
  title(paste(treename,': rooted'))

  # 2. Focusing step to obtain focused phylogeny
  ## Extract focusing taxa from the list
  focus_list <- taxa[2]
  unlisted_clade <- unlist(focus_list)
  focusedclade <- findMRCA(rootedtree, tips = unlisted_clade, type="node")
  ## extract the focusing clade as a new phylogeny and plot for records
  focusedtree <- extract.clade(rootedtree, focusedclade)
  plot.phylo(focusedtree)
  title(paste(treename,': focused'))

  # 3. Record focused tree tip names in file named as "Genus_tips.txt"
  write (focusedtree$tip.label, file=paste(treename, "_tips.txt"))
  genera <- unlist(strsplit(focusedtree$tip.label,split="_"))[1]
}

# Read user-specified csv file
args = commandArgs(trailingOnly=TRUE)
taxa_list <- read_func(args[1])

#specify absolute path to the unfocused trees in IQ_out/ directory
mydir <- (getwd())
iq_out_trees <- dir(mydir, pattern='.treefile')

#run focus_func on the unfocused trees
for (t in 1:length(iq_out_trees)) {
  focus_func(iq_tree=iq_out_trees[t], taxa_list)
}
