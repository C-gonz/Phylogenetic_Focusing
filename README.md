# Phylogenetic_Focusing v2.1.0
Phyfocus assesses the evolutionary history of target gene families via  proteome datasets, with an emphasis on exhaustive sequence analysis. Phyfocus first gathers all sequence data in each taxa similar to a larger gene family, extracts ("focuses") each taxa's monophyletic subtree containing the target gene family, then filters and concatenates the data to provide a multi-taxa gene family phylogeny.

PhyFocus consists of 6 steps:
1) Identifying protein datasets per species
2) Generating protein tree per species
3) Extracting focused protein tree per species
4) Concatenating focused protein tree sequences
5) Filtering concatenated dataset: HMMR
6) Generating final focused tree

There are several programs Phyfocus requires in the user PATH:
AWK  
Python3  
R  
NCBI BLAST+ (specifically, makeblastdb and blastp)  
CD-HIT  
MAFFT  
IQTREE  
HMMER  

## Core scripts run by user

#### phyfocus_2.1.sh (steps 1-5)
The main Phyfocus script managing steps 1-5. Once initiated Phyfocus runs automatically until step 6. Running Phyfocus requires 4 user files. Full usage and additional options are available under the help menu (./phyfocus_2.1.sh -h)

1) A query fasta file containing target, related, and root proteins (Step 1).
   - Targets represent homologs of the specific gene(s) / gene family being studied.
   - Related represent gene families possibly related to the target.
   - Roots represent a distant outgroup to target and related gene families; used to root the initial phylogeny.
   - All query proteins can come from a single organism, such as an established model.
   - There is a minimum of 2 sequences for roots, and 2 for target + related.
   - An initial reference phylogeny is helpful for assessing which proteins to include in the query.

2) A directory of protein FASTA files for each species assessed in the phylogeny (Step 1).
    - Protein sequences are ideally derived from whole-genome data or thorough transcriptomes.
    - All FASTA file names MUST begin with the species' genus name and underscore: genus_

3) A Tab Seperated Values (.tsv) file for phylogeny focusing (Step 3).
- There should be no column or row headers in the table.
- Column 1 gives FASTA ">" header names for AT LEAST 2 root proteins in the query fasta file.
- Column 2 gives FASTA ">" header names for AT LEAST 1 target AND 1 related protein in the query fasta file. Two of each is recommended.
- Below is an example TSV for a phylogeny of taste receptors within the glutamate/class-C GPCR family:
Homo_GABBR1    Homo_T1R1
Homo_GABBR2    Homo_mGLUR6
Danio_GABBR1
Danio_GABBR2
- Related and root sequences are best chosen by reference to previous phylogenies. The glutamate example was informed by the following: Fredriksson, R., Lagerström, M. C., Lundin, L. G., & Schiöth, H. B. (2003). The G-protein-coupled receptors in the human genome form five main families. Phylogenetic analysis, paralogon groups, and fingerprints. Molecular pharmacology, 63(6), 1256-1272.

4) A FASTA protein alignment characterizing key conserved domains and motifs (Step 5).
    - An alignment of the query file can be used.
    - Rigorous alignment (e.g., using MAFFT Linsi) is preferred.
    - Builds the HMMR profile for HMMR filtering to remove any fundamentally different sequences,
      such as proteins lacking a 7TM domain in a GPCR phylogeny.

#### alignment_editor.py (step 6)
After step 5, a concatenated alignment is produced that may contain massive gaps (100s+ of base pairs) due to a small subset of sequences. Alignment Editor removes these sequences according to user-specified gap length and gap frequency cutoffs. Current recommendation is to manually inspect the concatenated alignment provided by step 5 to identify the average size of massive gaps, then run alignment_editor.py. Full usage and options are available under the help menu (./alignment_editor.py -h)

## Accessory scripts run by Phyfocus  

#### Tree_Editor.R
Used to extract a focused subtree from the initial phylogeny for each study species (step 3), as per user-specified bait (target) and anchor (outgroup) sequences. The bait and anchor sequence headers are specified by the user via the TSV fasta headers file. 

#### header_translator.py
Helps the user correlate original BLAST sequence headers with the genus_#### headers used by phyfocus. Produces a TSV correlating these headers for all input proteome data.

#### delim_converter.py
Used to convert the space-delimited HMMER output txt file to a TSV for downstream use.
