# Phylogenetic_Focusing
The first functional version of PhyFocus. The program is tailored to assess a specific Taste Receptor dataset, and requires user-specific edits and 
data editing to produce a final product.

## Scripts provided in PhyFocus

#### phyfocus.sh
The core PhyFocus script. Manages the initial BLAST, initial phylogeny, and subtree extraction steps. Manual intervention by the user is required to
filter the subtree data and prepare the final focused phylogeny. Note that such quality filtering entails the use of HMMR (not currently coded in PhyFocus) and alignment editing via Alignment_Editor.py.

#### selectSeqs.pl
Used to extract sequences from a FASTA file given a set of FASTA sequence headers.

#### Tree_Editor.R
Used to extract a focused subtree from the initial phylogeny, as per user-specified bait (target) and anchor (outgroup) sequences. The specific sequence headers must be written into the script from the user-provided query file. 

#### Alignment_Editor.py
A script provided to aid with subtree data filtering. Removes sequences that cause excessive alignment gaps according to user-specified gap length and
gap frequency cutoffs. 

