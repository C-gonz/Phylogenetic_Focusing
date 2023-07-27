# Phylogenetic_Focusing v2.1.0
An abreviated version of this README is provided under the help menu (./phyfocus.sh -h)

Phyfocus assesses the evolutionary history of target gene families via  proteome datasets, with an emphasis on exhaustive sequence analysis. Phyfocus first gathers all sequence data in each taxa similar to a larger gene family, extracts ("focuses") each taxa's monophyletic subtree containing the target gene family, then filters and concatenates the data to provide a multi-taxa gene family phylogeny.

PhyFocus consists of 6 steps:
1) Identifying protein datasets per species
2) Generating unfocused tree per species
3) Extracting focused tree per species
4) Concatenating focused protein tree sequences
5) Filtering concatenated dataset: HMMR
6) Generating final focused tree

# Running PhyFocus

## Dependencies & Setup

### Phyfocus Dependencies:
Programs Phyfocus requires in the user PATH
- AWK
- Python3 (including BioPython package -> Bio.AlignIO module)
- R
- NCBI BLAST+ (specifically makeblastdb and blastp)
- CD-HIT
- MAFFT
- IQTREE
- HMMER

Dependency setup details still in progress

### Timeframe, CPU, and RAM
Phyfocus can be time and memory intensive; effective system recommendations:
- Multiple threads/CPUs (e.g. 24).
- 125GB RAM

For an exceptionally large  phylogeny with 3200 sequences, the program timeframe using the above metrics was thus:
- Step 1: 17h
- Step 2: 148h (~ 6 days)
- Step 3: < 1min
- Steps 4-5 and alignment for Step 6: (< 4min)
- Step6 Phylogeny:  232h (~ 10 days)

Note that more reasonably sized runs can take 7-10 days total. Less CPU & RAM may be feasible, but could incur longer runtimes or errors if RAM depletes.



## Testing Installation
To ensure Phyfocus works properly, follow these steps.

1. Ensure the following 5 scripts and 1 directory are downloaded to your desired working directory:
    - phyfocus.sh (current version)
    - tree_editor.R
    - header_translator.py
    - alignment_editor.py
    - delim_converter.py
    - sample_data/ (contains 10 files)

2. In your working directory, run the command  "./phyfocus.sh -T" which will execute the full phyfocus pipeline in sample_data (a very small bacterial dataset; estimated runtime is 8-12 minutes).

3. Once complete, within sample_data/ you can examine out_log.txt to observe Phyfocus output and check error_log.txt for any error messages. If Phyfocus ran successfully, you should find the following directories and fasta file:

        ./sample_data/final_tree_dataset/filtering_output/final_tree_seqs_ali.fa



## Running Phyfocus on Your Data 
PhyFocus consists of 6 steps:
1) Identifying protein datasets per species
2) Generating unfocused tree per species
3) Extracting focused tree per species
4) Concatenating focused protein tree sequences
5) Filtering concatenated dataset: HMMR
6) Generating final focused tree

The following sections will walk you through these major steps and things to be aware of.

### Setup
Before starting phyfocus on your own dataset, ensure your working directory contains the following 5 user-provided datasets in addition to the 5 scripts noted in "testing Installation" above. Note that examples of these 5 datasets can be seen in ./sample_data

1) A query fasta file containing target and anchor proteins (Step 1).
   - Targets represent homologs of the specific gene(s) / gene family being studied.
   - Anchors represent homologs of gene(s) / gene families closely related to the targets.
   - All query proteins can come from a single organism, but broader sampling from taxa of interest may aid BLAST.
   - There is a minimum of 2 for target + anchor sequences (more are recommended).

2) A directory of protein FASTA files for each species assessed in the phylogeny (Step 1).
    - Protein sequences are ideally derived from whole-genome data or thorough transcriptomes.
    - All FASTA file names MUST begin with the species' genus name and an underscore: genus_

3) A roots fasta file containing proteins to root the per-species unfocused trees (steps 2 & 3) and your final phylogeny (step 6).
    - Roots represent the closest known outgroup to target & anchor gene families.
    - There is a minimum of 2 sequences required for roots.

4) A Tab Seperated Values (.tsv) file for phylogeny focusing (Step 3).
    - There should be no column or row headers in the table.
    - Column 1 gives FASTA ">" header names for AT LEAST 2 root proteins from the roots fasta file.
    - Column 2 gives FASTA ">" header names for AT LEAST 1 target and 1 anchor protein from the query fasta file. At least two of each is recommended.
    - Below is an example TSV for a phylogeny of taste receptors within the glutamate/class-C GPCR family:
Homo_GABBR1    Homo_T1R1  
Homo_GABBR2    Homo_mGLUR6  
Danio_GABBR1  
Danio_GABBR2  
    - Anchor and root sequences are best chosen by reference to previous phylogenies. The glutamate example was informed by: Fredriksson, R., Lagerström, M. C., Lundin, L. G., & Schiöth, H. B. (2003). The G-protein-coupled receptors in the human genome form five main families. Phylogenetic analysis, paralogon groups, and fingerprints. Molecular pharmacology, 63(6), 1256-1272.

5) A FASTA protein alignment characterizing key conserved domains and motifs (Step 5).
    - An alignment of the query file can be used.
    - Rigorous alignment (e.g., using MAFFT Linsi) is preferred.
    - Builds the HMMR profile for HMMR filtering to remove any fundamentally different sequences,
      such as proteins lacking a 7TM domain in a GPCR phylogeny.


### Executing phyfocus.sh (steps 1-5)
The main Phyfocus script manages steps 1-5. Once initiated, Phyfocus runs automatically until step 6. 
An abreviated version of this info is provided under the help menu (./phyfocus.sh -h)

--------------------------------------------------------------------------------------------------------------------------------
Program syntax: 
[-h] [-T] [-X] [-q ./file] [-r ./file] [-f ./directory/] [-c ./file] [-H ./file] [-t num] [-e value] [-m string] [-b num] [-a num]

-h <help>   Display this help and exit

-q QUERY    Required; fasta file containing target, related, and outgroup query peptides

-r root     Required; fasta file containing rooting peptides

-f FASTAS   Required; directory of peptide fasta files for each species in the phylogeny

-c CLADE    Required; A .tsv file of fasta header names for phylogeny root and focus sequences.

-H HMMR     Required; Peptide fasta alignment for making the HMMR profile

-t THREADS  Optional; The number of threads used for BLAST & IQtree. Default = 24

-e EVALUE   Optional; The e-value significance cutoff used in BLAST. Default = 1e-05

-m MODEL    Optional; The peptide substitution model used by IQTree. Recommend LG for testing. Default = MFP+C60

-b BOOT     Optional; The number of Ultrafast Bootstrap replicates used by IQTree. Default = 1000

-a IQALRT   Optional; The number of SH-aLRT Bootstrap replicates used by IQTree. Minimum value is 1000. Default = 1000

-T TEST     Optional; Runs Phyfocus on the data in ./sample_data to test dependencies and demo the program

-X CLEAN    Optional; Removes all phyfocus-created files in the working  directory. Is called automatically when running the ./sample_data via -T, allowing for multiple test runs.

-------------------------------------------------------------------------------------------------------------------------------
A basic executing command using all default settings could be thus:

    ./phyfocus.sh -q ./query_file.fa -f fasta_proteins/ -r roots_file.fa -H hmmr_ali.fa -c tree_tips.tsv


### Monitering the Data
Because Phyfocus can take over a week, assessing if the per-species blast (Step 1) is capturing too narrow or too broad a dataset for your phylogenies can save significant time. Two approaches can help here:  

1. Examine the BLAST output files at ./blastout_maxseqs. These have been edited to include the original header descriptions, allowing you to see what sequences given blast queries and the e-value cutoff are obtaining. A potential issue to look for here is if many gene family members from far outside your rooting seqs' family are being obtained. 

2. Once they are generated the unfocused (Step 2) and focused (Step 3) species trees can also be examined in detail by looking at ./align/IQ_out and ./align/IQ_out/tree_editor_out, respectively, or in quick summary by looking at Rplot.pdf in tree_editor_out. While the trees can be annotated to observe the focusing impact and tree quality, a quicker look referencing just the named target, anchor, and root sequences can give an idea if proper monophyletic clades are being reconstructed.
- NOTE: when examining the unfocused and focused trees in Rplot.pdf, R may draw a monophyletic clade used for rooting as a polytomy. Programtically it does not appear to treat the clade as such, and thus shouldn't affect the focusing step.

If the blast outputs and/or per-species trees indicate an issue with the dataset you may need to rerun the data using BLAST evalues different from the default of 1e-05. 
- For example, if too many genes outside your rooting seqs are bloating the alignments and subsequent trees, consider a more stringent evalue such as 1e-10 (check the blast output files to get a sense of what cutoff may help).


### Executing alignment_editor.py (step 6) 
Once complete, ./phyfocus.sh produces an alignment that is the starting point of Step 6 (see "final_tree_dataset/filtering_output/final_tree_seqs_ali.fa"). It will also move alignment_editor.py to this directory for user convenience.

final_tree_seqs_ali.fa may contain massive gaps (100s+ of base pairs) due to a small subset of sequences. To improve alignment quality, Alignment Editor removes these gap-causing sequences according to user-specified gap length and gap frequency cutoffs. Current recommendation is to manually inspect the concatenated alignment to identify the average size of massive gaps, then run alignment_editor.py. For some datasets, a helpful starting metric can be to target gap lengths approximately 4% of the entire alignment length. 

Usage syntax and options for alignment_editor are also available under the help menu (./alignment_editor.py -h)

-------------------------------------------------------------------------------------------------------------------------------
Program syntax: 
./alignment_editor.py  ./final_tree_seqs_ali.fa [-h] [-w num] [-g num] 

-h <help>   Display help and exit

-w <window> Specify the minimum gap length that is disruptive to the alignment. Default is 50bp, but an effective
benchmark appears to be 4% of the alignment's total length. Manual assessment of the alignment
is strongly recommended.

-g <gappercent> Specify the minimum percent of seqs that must contain -w size gaps to identify problematic seqs.
Default is 0.9.

-------------------------------------------------------------------------------------------------------------------------------
Upon completing alignment_editor.py, running IQtree will complete Step 6 and produce the final focused phylogeny. Recommended command:

    iqtree -s final_tree_seqs_ali.fa_edited_ali.fa -m MFP+C60 -alrt 1000 -bb 1000 -nt 24



## Assessing Phyfocus Phylogenies

### Viewing Phylogenies
The Newick formatted data of step 6 (./final_tree_dataset/filtering_output/final_tree_seqs_ali.fa_edited_ali.fa.treefile) can be copied into a viewer like FigTree for rooting with the root sequences clade and annotating features of interest. 
(Rambaut, A. (2018) FigTree v1. 4.4: a graphical viewer of phylogenetic trees. Available from <http://tree.bio.ed.ac.uk/software/figtree/>.) 

### Visualizing Node Support
The Ultrafast bootstrap and SH-aLRT node values provided by IQtree indicate which nodes are confidentally supported; IQtree recommends a benchmark of UFboot >= 95% and SH-aLRT >= 80% (http://www.iqtree.org/doc/Frequently-Asked-Questions). To collapse unsupported nodes into a polytomy, the following protocol may be used:

1. Extract node values via Tree Graph 2 (<http://treegraph.bioinfweb.info/>)
        - Run GUI via .jar file. 
	- open .treefile
	- export node data as table
		- use "unique node names" and "bootstrap/SH-aLRT" labels columns (in that order)
		- export all nodes with column headers
	- Do not close the .treefile window after exporting

2. Change UFboot/SH-aLRT values to 1 value per column
- Open exported .txt in excel/numbers as TSV
- Change the label column header to "UF-Boot" and add a new column: "SH-aLRT"
- Export table as CSV
- Open CSV in a text editor, use find+replace to change all "/" to ","
- Open CSV in excel/numbers again, delete the empty column at right.

3. Collapse Phylogeny based on node support
- import CSV into the open Tree Graph 2 window via "import table as node/branch data" (all files format, no skipped lines, first line contains headers. Be sure to check values separated by ",")
- "select matching key columns" step correlates same data in the table and tree. Use unique node names column for both.
- Set node data type for UFboot and SH-aLRT as "new hidden branch data with specified ID." 
- Use "Collapse Nodes by Support" option twice: UFBoot at 95, and SH-aLRT at 80. 
- Export result as a nexus file, the collapsed tree can now be viewed or edited in Figtree, etc. as before.



## Accessory Scripts Run by Phyfocus  
These necessary scripts are not executed by the user; these details are purely for informational purposes.

### Tree_Editor.R
Used to extract a focused subtree from the initial phylogeny for each study species (step 3), identified by finding the most recent common ancestor of the user-specified target and anchor sequences. The bait and anchor sequence headers are specified by the user via the TSV fasta headers file. 

### header_translator.py
Produces a TSV correlating original BLAST sequence headers with the genus_#### headers used by phyfocus for all input sequence data. The TSV is also used automatically to add BLAST hit descriptions as a final column in each BLAST output file.

### delim_converter.py
Used to convert the space-delimited HMMER output txt file to a TSV for downstream use.



## Managing Errors & Reruns
Phyfocus (and programs it uses like IQtree) are designed to cancel the run if certain errors or issues occur. The out_log.txt file will indicate the last steps taken, while specific error messages are recorded in error_log.txt. 

Once the error has been dealt with, phyfocus.sh can be restarted at the last valid step by re-running the same command originally used, provided the contents of its working directory have not been changed by the user. The restarted run date and all steps skipped are always appended to the out_log.txt file before standard output logging resumes.

The following represent some errors you may encounter in error_log.txt, along with recommended solutions:

#### Numerical underflow (lh-branch). Run again with the safe likelihood kernel via `-safe` option
- An error from IQtree during Step 2 that requires IQtree's "-safe" option. This special run should be done manually on the specific species' dataset concerned (see out_log.txt for the last species' dataset being used). 
- It is crucial to ensure this is done in the ./align directory, and using the same IQtree parameters originally specified by the user and/or the default phyfocus settings.
- An example solution starting from the working directory:
    cd ./align
    iqtree -safe -s Homo_protein.faa_FIX_hits.fa_ali.fa -m MFP+C60 -alrt 1000 -bb 1000-nt 24 2>> ../error_log.txt | tee -a ../out_log.txt
    - Once complete, phyfocus.sh can be restarted as normal using the original input commands.
