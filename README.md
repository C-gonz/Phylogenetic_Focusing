# Phylogenetic_Focusing v2.0.0
The program is tailored to assess the phylogenetic history of specific gene families via large scale proteome datasets from taxa of interest. Phyfocus operates on the premise that the most rigorous phylogenetic assessment of gene evolution is one that includes all possible sequence data in a per species phylogeny at the outset, excising ("focusing") monophyletic subtrees containing the gene families of interest, then filtering the data and ptoviding a concatenated, multi-species gene tree for the gene families of interest.

## Core Scripts

#### phyfocus.sh
The central PhyFocus script. Manages the initial BLAST, initial phylogeny, and subtree extraction steps. Currently the user must edit Tree_Editor.R to their specifications before use. Phyfocus currently runs automatically up till Alignment_Editor.py.

#### Tree_Editor.R
Used to extract a focused subtree from the initial phylogeny for each study species, as per user-specified bait (target) and anchor (outgroup) sequences. The bait and anchor sequence headers must be written into the script before running phyfocus; said headers should be obtained from the user-provided BLAST query file. 

#### Alignment_Editor.py
Aids with subtree data filtering. Removes sequences that cause excessive alignment gaps according to user-specified gap length and
gap frequency cutoffs. Current recommendation is to manually inspect the concatenated subtree alignment provided by Phyfocus.sh, then run Alignment_Editor.py. Manual running of IQ-Tree on the output alignment is required.

## Accessory Scripts 

#### header_translator.py
Helps the user correlate original BLAST sequence headers with the genus_#### headers used by phyfocus. Produces a tab-delimited file correlating these headers for all input proteome data.

#### delim_converter.py
Used to convert the space-delimited HMMER output txt file to a tab delimited file.

#### selectSeqs.pl
Used to extract sequences from a FASTA file given a set of FASTA sequence headers.
