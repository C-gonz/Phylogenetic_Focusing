#!/bin/bash

help() {
cat << EOF
----------------------------------------------
PhyFocus v2.1.1 (Oct 6 2025)
github repository: https://github.com/C-gonz/Phylogenetic_Focusing
See "README" for detailed instructions on the Phyfocus pipeline
----------------------------------------------
Syntax:
./${0##*/} [-h] [-q ./file] [-o ./file] [-f ./directory/] [-c ./file] [-H ./file] [-t num] [-e value] [-a string] [-A string] [-m string] [-b num] [-s num] [-F string] [-T] [-X]

-h <help>   Display this help and exit

REQUIRED ITEMS
-q QUERY    Fasta file containing peptides for BLASTp query that represent your targets of interest.
            Anchors may also be included
-o OUT      Fasta file containing outgroup (anchor and rooting) peptides
-f FASTAS   Directory of peptide fasta files for each species in the desired phylogeny
-c CLADE    A .tsv file of fasta header names for root and target+anchor sequences

OPTIONAL RUN PARAMETERS
-H HMMER    Peptide fasta alignment for making the HMMER profile
-t THREADS  Number of threads used for BLAST & IQtree. Default = 24
-e EVALUE   E-value significance cutoff used in BLAST. Default = 0.05

OPTIONAL ALIGNMENT & TREE PARAMETERS
-a ALIGN1   MAFFT method for Step2 unfocused species alignments. Default is "linsi" but the fast progressive
            method FFT-NS-2 (enter "mafft --retree 2 --maxiterate 0") can be used for testing. Include the
            quotes. See the MAFFT website for additional options. Using 'mafft' alone enables MAFFT
            to auto select a method.
-A ALIGN2   MAFFT method for Step4 combined species alignment. Default and alternatives same as for -a ALIGN1
-m MODEL    Peptide substitution model used by IQTree. Recommend LG for quick testing. Default = MFP+C60
-b BOOT     Number of Ultrafast Bootstrap replicates used by IQTree. Minimum value is 1000.
            Ignore this option for faster but less accurate trees.
-s SHALRT   Number of SH-aLRT Bootstrap replicates used by IQTree. Minimum value is 1000.
            Ignore this option for faster but less accurate trees.

ADDITIONAL USES
-F FORCE    Forces phyfocus to rerun at one of the following checkpoints.
                S2a = Step 2 (alignments), remake per-species alignments
                S2t = Step 2 (trees), remake per-species unfocused trees
                S3  = Step 3 (focusing), remake per-species focused trees
                S4  = Step 4 (extract and combine), remake combined dataset for final tree
            Note: -F reruns delete all pre-existing output that comes after the chosen checkpoint,
                  then runs phyfocus as normal.
-T TEST     Run while in ./sample_data to test dependencies and demo the program. Command =
            "./phyfocus.sh -q ./query_seqs.fa -f fasta_proteins -o outgroups.fa -H query_seqs_ali.fa 
            -c focus_table.tsv -a 'mafft --retree 2 --maxiterate 0' -A 'mafft --retree 2 --maxiterate 0' -m LG"
-X CLEAN    Removes all phyfocus output files in the current working directory.

----------------------------------------------
PATH Dependencies:
    AWK
    Python3 (including BioPython package -> Bio.AlignIO module)
    R
    NCBI BLAST+ (specifically makeblastdb and blastp)
    CD-HIT
    MAFFT
    IQTREE
    HMMER (if desired)
Included Dependencies:
    ${0##*/} and its subscripts:
        tree_editor.R
        header_translator.py
Optional "accessory" programs
    download_formats.sh
    fasta_lengths.py
    alignment_editor.py
        species_check.sh (subscript for alignment_editor.py)
----------------------------------------------

Description:
PhyFocus assesses gene family evolution by "focusing," or extracting, the gene clade(s) of interest from a much broader gene phylogeny. By starting with an extensive phylogeny of outgroup clades, Phyfocus decreases the chance that significant gene family relationships are excluded.

PhyFocus consists of 6 steps:
1) Identifying protein dataset per species
2) Filtering, alignment, & unfocused phylogeny per species
3) Extracting focused phylogeny dataset per species
4) Concatenating focused species datasets
5) Concatenated alignment, user-run inspection & editing
6) User-run final phylogeny creation

PhyFocus requires four user-provided datasets (5 if HMMER is used):

1) A query fasta file containing target proteins (Step 1).
   - There is a minimum of 2 for target + anchor sequences (more are recommended).

2) A directory of protein FASTA files for each species assessed in the phylogeny (Step 1).
   - All FASTA file names MUST begin with the species genus name and underscore: genus_
   - NOTE: if using multiple species from 1 genus, you must distinguish the file names:
     E.g.: Canis lupes & Canis Familiaris --> CanisL_file & CanisF_file
   - The included "download_formats.sh" script helps automate this process (see ./download_formats.sh -h for details)

3) A FASTA protein alignment characterizing key conserved domains and motifs (Step 2).
   - Only for HMMER filtering; can be ignored if HMMER is not being used.
   - Can be a rigorous alignment of the query sequences, but a more thorough profile may perform better.

4) A fasta file containing outgroup proteins to root the per-species trees (root seqs) and to focus these
   trees down to your target gene family of interest (anchor seqs) (Steps 2 & 3).
   - A minimum of 2 sequences for roots is required
   - Only 1 anchor sequence is required, but more can be used.

5) A Tab Seperated Values (.tsv) file to guide phylogenetic focusing (Step 3).
   - There should be no column or row headers in the table.
   - Column 1 gives FASTA header names (no ">") for AT LEAST 2 root proteins from the outgroup file.
   - Column 2 gives FASTA header names (no ">") for target proteins from the query file, and anchor
     proteins from the outgroup file. AT LEAST 2 sequences total is required.

EOF
}

# Function for running a test demo of phyfocus using provided data
test() {
echo "Runnning Phyfocus using:"
echo "./${0##*/} -q ./query_seqs.fa -f fasta_proteins -o outgroups.fa -H query_seqs_ali.fa -c focus_table.tsv -a 'mafft --retree 2 --maxiterate 0' -A 'mafft --retree 2 --maxiterate 0' -m LG"
./${0##*/} -q ./query_seqs.fa -f fasta_proteins -o outgroups.fa -H query_seqs_ali.fa -c focus_table.tsv -a "mafft --retree 2 --maxiterate 0" -A "mafft --retree 2 --maxiterate 0" -m LG
}

# Function for removing generated files from failed runs,etc.
clean() {
if [[ -e ./phyfocus.sh ]]
then
    starting_num=$(ls | wc -l)
    rm -r ./align_species ./blastdb ./blastout_tables ./fixed_fastas ./hits_accessions ./hits_fasta ./tip_seqs ./final_tree_dataset ./logs 2> /dev/null
    rm out_log.txt summary_log.txt error_log.txt header_translation_table.tsv temp.txt all_fixed_numerical_headers.txt all_ncbi_headers.txt hmmr_profile.hmm formatted_query.fa queries_and_outs.fa queries_and_outs.fa.clstr temp_catseqs.fa *_FIX.fa 2> /dev/null
    clean_num=$(ls | wc -l)
    diff=$(($starting_num -$clean_num))
    echo "Removed $diff items. phyfocus working directory now contains $clean_num items."
else
    echo "Run Phyfocus -X in the working directory containing phyfocus-required files, including phyfocus.sh itself."
fi
}

# Create input option variables
QUERY=""
FASTAS=""
OUT=""
CLADE=""
HMMER=""
THREADS=24
EVALUE="0.05"
ALIGN1="linsi"
ALIGN2="linsi"
MODEL="MFP+C60"
BOOT=0
SHALRT=0
FORCE=""
regex_num="^[0-9]+$"
regex_enum="^[0-9]+e-[0-9]+$"
user_args=$@

# Create error messages for improper option input
q_error="Option error; path must be to an existing file with data. Format: -q <./fasta_query_file>"
f_error="Option error; path must be to an existing directory. Format: -f <./species_fasta_proteins_directory/>"
o_error="Option error; path must be to an existing file with data. Format: -o <./fasta_outgroup_file>"
c_error="Option error; path must be to an existing file with data. Format: -c <./tsv_table>"
H_error="Option error; path must be to an existing file with data. Format: -H <./fasta_alignment_file>"
t_error="Option error; thread usage must be an integer. format: -t <integer>"
e_error="Option error; need value in scientific notation. Format: -e <nums>e-<nums>"
a_error="Option error; Needs a string denoting a MAFFT program. Format: -a '<string>'"
A_error="Option error; Needs a string denoting a MAFFT program. Format: -A '<string>'"
m_error="Option error; Text needs to be an IQtree model. Format: -m <string>"
b_error="Option error; bootstrap replicates must be an integer. Format: -b <integer>"
s_error="Option error; SH-aLRT replicates must be an integer. Format: -s <integer>"
F_error="Option error; Needs text denoting which Phyfocus checkpoint to start at. Format: -F <string>"

# Handling for option arguments, including improper arguments
while getopts ":hq:f:o:c:H:t:e:a:A:m:b:s:F:TX" option; do
    case $option in
        h) help; exit 0;;
        q) QUERY=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $q_error >&2; exit 1; fi;;
        f) FASTAS=$OPTARG; if [[ ! -d $OPTARG ]]; then echo $f_error >&2; exit 1; fi;;
        o) OUT=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $o_error >&2; exit 1; fi;;
        c) CLADE=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $c_error >&2; exit 1; fi;;
        H) HMMER=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $H_error >&2; exit 1; fi;;
        t) THREADS=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $t_error >&2; exit 1; fi;;
        e) EVALUE=$OPTARG; if [[ ! $OPTARG =~ $regex_enum ]]; then echo $E_error >&2; exit 1; fi;;
        a) ALIGN1=$OPTARG; if [[ -z $OPTARG ]]; then echo $a_error >&2; exit 1; fi;;
        A) ALIGN2=$OPTARG; if [[ -z $OPTARG ]]; then echo $A_error >&2; exit 1; fi;;
        m) MODEL=$OPTARG; if [[ -z $OPTARG ]]; then echo $m_error >&2; exit 1; fi;;
        b) BOOT=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $b_error >&2; exit 1; fi;;
        s) SHALRT=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $a_error >&2; exit 1; fi;;
        F) FORCE=$OPTARG; if [[ -z $OPTARG ]]; then echo $F_error >&2; exit 1; fi;;
        T) test; exit 0;;
        X) clean; exit 0;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done
# Print help menu if no options are specified.
if [[ $# == 0 ]]; then help; exit 0; fi
# Exit program if required options are not specified by user
if [[ (-z $QUERY) || (-z $FASTAS) || (-z $OUT) || (-z $CLADE) ]]; then echo "Missing required option. Syntax: ./${0##*/} [-q file] [-o file] [-c file] [-f directory] [-H file]" >&2; exit 1; fi



log_start() {
mkdir ./logs
echo -------------------------------- | tee -a ./logs/out_log.txt ./logs/error_log.txt
echo "PhyFocus v2.1.1 (Jul 2 2025)" | tee -a ./logs/out_log.txt ./logs/summary_log.txt
echo "github repository: https://github.com/C-gonz/Phylogenetic_Focusing"
echo "See 'README' for detailed instructions on the overall Phyfocus pipeline"
echo "Phyfocus run started:" $(date) | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
echo "Phyfocus run's home directory:" $(pwd) | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
echo "user input for this run: ./${0##*/} $user_args" | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
}



step1_headers() {
echo -------------------------------- | tee -a ./logs/summary_log.txt ./logs/out_log.txt
echo "STEP 1: IDENTIFYING PROTEIN DATASET PER SPECIES" | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
echo "Unzipping any gzip files in $FASTAS ..." | tee -a ./logs/out_log.txt
gzip -d ${FASTAS}/*.gz 2>> ./logs/error_log.txt | tee -a ./logs/out_log.txt
# For each species fasta file, use parameter expansion to extract genus name, then AWK
# to change each sequence header to a "genus_#####" numerical header, and write each sequence
# on 1 line without STOP signals.
# Syntax: "NR==1" prevents line feed at start of file, "genus "_%05d\n", ++i" creates
# genus_numerical headers, "else {gsub("\\*", "", $0); printf $0}" ensures every sequence line is written
# on 1 line and has removed any STOP (*) signals.
echo "Fixing FASTA headers and seqs ..." | tee -a ./logs/out_log.txt
for original_filename in $FASTAS/*; do file_name=${original_filename##*/}; awk -v genus=${file_name%%_*} '{if(NR==1) {printf ">" genus "_%05d\n", ++i "\n"} else {if($0 ~ /^>/) {printf "\n" ">" genus "_%05d\n", ++i "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf $0 "\n"}' $original_filename > ${file_name}_FIX.fa 2>> ./logs/error_log.txt; done
# Conduct same line fix for user query and outgroup files
awk '{if(NR==1) {printf $0 "\n"} else {if($0 ~ /^>/) {printf "\n" $0 "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf "\n"}' $QUERY | tee temp_catseqs.fa formatted_query.fa
awk '{if(NR==1) {printf $0 "\n"} else {if($0 ~ /^>/) {printf "\n" $0 "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf "\n"}' $OUT >> temp_catseqs.fa
# Ensure there are no duplicates between the Query and outgroup seqs
cd-hit -i temp_catseqs.fa -o queries_and_outs.fa -c 1.0 -n 5 2>> ./logs/error_log.txt | tee -a ./logs/out_log.txt
rm temp_catseqs.fa
# Generate a file containing all NCBI headers for header translation table
cat ./$FASTAS/* | grep ">" > all_ncbi_headers.txt 2>> ./logs/error_log.txt
# make directory for fixed fastas, then move FIX files
echo "moving modified header fastas to ./fixed_fastas/ ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
mkdir ./fixed_fastas
mv *_FIX.fa ./fixed_fastas
echo "Successfully moved $(ls ./fixed_fastas | wc -l) files to ./fixed_fastas/" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
echo "FASTA Header Fix Complete" | tee -a ./logs/out_log.txt
# Make a translation table for numerical headers --> informative NCBI headers
echo "Making numerical sequence header translation table in ./logs." | tee -a ./logs/out_log.txt
cat ./fixed_fastas/* | grep ">" > all_fixed_numerical_headers.txt 2>> ./logs/error_log.txt
./subscripts/header_translator.py all_ncbi_headers.txt all_fixed_numerical_headers.txt 2>> ./logs/error_log.txt | tee -a ./logs/out_log.txt
mv ./header_translation_table.tsv ./logs
}
step1_blast() {
# Format blast databases from $FASTAS & run BLASTp
echo "Blasting $FASTAS peptide databases with formatted_query.fa using (blastp -db FILE -query formatted_query.fa -evalue $EVALUE -num_threads $THREADS -out ${species_pep%.*}_ref_blastout -outfmt 6 -max_target_seqs 5000) ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
for species_pep in ./fixed_fastas/*
do
    makeblastdb -in $species_pep -parse_seqids -out ${species_pep%.*}_db -dbtype prot 2>> ./logs/error_log.txt | tee -a ./logs/out_log.txt
    echo "formatted BLAST database for $species_pep proteins..." | tee -a ./logs/out_log.txt
    # blastp query file against species_pep database
    echo "BLASTING formatted_query.fa sequences against $species_pep database ..." | tee -a ./logs/out_log.txt
    blastp -db ${species_pep%.*}_db -query formatted_query.fa -evalue $EVALUE -num_threads $THREADS -out ${species_pep%.*}_ref_blastout -outfmt 6 -max_target_seqs 5000 2>> ./logs/error_log.txt | tee -a ./logs/out_log.txt
    # Error check: each blast completed
    if [[ ! (-e ${species_pep%.*}_ref_blastout) ]]; then echo "BLAST output for $species_pep was not successfully made. Check blast DB and out files in ./fixed_fastas" >> ./logs/error_log.txt; exit 1; fi
    if [[ ! (-s ${species_pep%.*}_ref_blastout) ]]; then echo "WARNING: BLAST output for $species_pep is empty. No blast results for this species will be included in the phylogeny." | tee -a ./logs/summary_log.txt ./logs/out_log.txt; fi

done
# move blast dbs & blast outputs to new directories
echo "moving db's to ../blastdb" | tee -a ./logs/out_log.txt
mkdir ./blastdb
mv ./fixed_fastas/*_db* ./blastdb
echo "moving blastout files to ./blastout_tables" | tee -a ./logs/out_log.txt
mkdir ./blastout_tables
mv ./fixed_fastas/*blastout ./blastout_tables/
# Modify blastout files to include a column for BLAST hit header descriptions
echo "modifying blastout tables to include hit descriptions..." | tee -a ./logs/out_log.txt
cd ./blastout_tables
for blastout in ./*
do
    echo "Editing" $blastout | tee -a ../logs/out_log.txt
    names=$(cut -f 2 $blastout)
    for header in $names; do grep -w $header ../logs/header_translation_table.tsv >> temp; done
    cut -f 2 temp > temp2
    paste -d '\t' $blastout temp2 > ${blastout}_info
    rm temp* $blastout
done
cd ..
}



step2_initial_alignments() {
echo -------------------------------- | tee -a ./logs/out_log.txt
echo "STEP 2: GENERATING PROTEIN TREE PER SPECIES ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
echo "Writing blast hit FASTA files ..." | tee -a ./logs/out_log.txt
# Retrieve accession numbers for each blastp output & move to new dir
for blastfile in ./blastout_tables/*blastout_info; do echo "Retrieving blast hit accessions for $blastfile ..." | tee -a ./logs/out_log.txt; cut -f 2 $blastfile > ${blastfile%.}_hits.txt; done
# Move all hit accession files to new dir called hits_accessions (rm is for failed runs)
rm -r ./hits_accessions; mkdir ./hits_accessions
mv ./blastout_tables/*_hits.txt ./hits_accessions
echo "Successfully moved $(ls ./hits_accessions | wc -l) hit files to hits_accessions" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
# Retrieve blast hit seqeunces from *_FIX.fa files for each taxa
echo "Retrieving blast hit seqeunces from fixed FASTA files ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
for hitfile in ./hits_accessions/*_hits.txt
do
    # Check if no BLAST hits obtained
    if [[ ! (-s $hitfile) ]]; then echo "WARNING: there are no seqs in $hitfile; this species will not be represented in phylogeny." | tee -a ./logs/summary_log.txt ./logs/out_log.txt; continue; fi
    # store hit accessions in a temp FASTA file & store fixed_fasta file name
    awk '{printf ">" $0 "\n"}' < $hitfile > temp.txt
    filenamefull=${hitfile##*/}
    in_filename="${filenamefull%_ref_blastout_info_hits*}.fa"
    hits_fasta="${filenamefull%_ref_blastout_info_hits*}_hits.fa"
    # Retrieve blast hit seqs from *_FIX.fa files for each taxa
    echo "pulling FASTA seqs from fixed_fastas/ for hit accessions in $in_filename ..." | tee -a ./logs/out_log.txt
    grep -A 1 -w -f temp.txt --no-group-separator ./fixed_fastas/$in_filename > $hits_fasta 2>> ./logs/error_log.txt
    # Error check: blast hit seqs extracted
    if [[ ! (-e $hits_fasta) ]]; then echo "File for hits_$filename was not properly made; examine fixed fasta file and hits_accessions directory." >> ./logs/error_log.txt; exit 1; fi
done
rm ./temp.txt
# Move retrieved hit fasta seqs to the new directory hits_fasta/
mkdir hits_fasta
mv *_hits.fa ./hits_fasta
echo "Successfully moved $(ls ./hits_fasta/*_hits.fa | wc -l) hit sequence fasta files to hits_fasta" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
cd ./hits_fasta
if [[ -n $HMMER ]]
then
    # Build a HMMER gene homolog profile
    echo "Building HMMER profile (hmmbuild) for matching against subtree seqs (hmmsearch) ..." | tee -a ../logs/out_log.txt
    hmmbuild hmmr_profile.hmm ../$HMMER 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
    # Remove similar seqs, HMMER filter, concatenate queries+anchors+hits, & align for each taxa
    for fasta in ./*_FIX_hits.fa
    do
        # run CD-HIT to remove seqs with 98% or higher seq similarity
        echo "Removing seqs >/= 98% identical via CD-HIT (-c 0.98 -n 5) on $fasta ..." | tee -a ../logs/out_log.txt
        cd-hit -i $fasta -o ${fasta%.}_cdhit -c 0.98 -n 5 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
        # Use HMMER profile to identify target homologs in the focused subtree dataset
        echo "Filtering ${fasta%.}_cdhit using hmmsearch and HMMER profile" | tee -a ../logs/out_log.txt
        hmmsearch hmmr_profile.hmm ./${fasta%.}_cdhit > ./${fasta%.}_hmmrsearch.out 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
        # Extract significant HMMER seq headers from hmmsearch results, convert to FASTA format
        echo "Extracting HMMER significantly matched sequences ..." | tee -a ../logs/out_log.txt
        awk '/^ *>> [A-Za-z]+_[0-9]+ /' ${fasta%.}_hmmrsearch.out > ${fasta%.}_hmmrsearch_temp.txt 2>> ../logs/error_log.txt
        while IFS= read -r header; do fixed=$(echo ${header/> />}); echo ${fixed:1} >> ${fasta%.}_hmmrsearch_headers.txt 2>> ../logs/error_log.txt; done < ./${fasta%.}_hmmrsearch_temp.txt
        grep -A 1 --no-group-separator -w -f ${fasta%.}_hmmrsearch_headers.txt ./${fasta%.}_cdhit > ./${fasta%.}_hmmr_significant_hits.fa 2>> ../logs/error_log.txt
        # Error check: HMMER profile, HMMER search result, and FASTA extraction
        if [[ ! ((-s hmmr_profile.hmm) && (-s ${fasta%.}_hmmrsearch.out)) && (-s ./${fasta%.}_hmmr_significant_hits.fa) ]]; then echo "Error occured when generating HMMER profile, HMMER search results, or in extracting HMMER filtered fasta sequences. Review files in ./hits_fasta" 2>> ../logs/error_log.txt; exit 1; fi
        # Concatenate HMMER hits with outgroups+queries, align via MAFFT
        echo | tee -a ../logs/out_log.txt; echo "Concatenating seqs from ../queries_and_outs.fa & ./${fasta%.}_hmmr_significant_hits.fa ..." | tee -a ../logs/out_log.txt
        cat ../queries_and_outs.fa ././${fasta%.}_hmmr_significant_hits.fa > ${fasta%.}_bts_ancrs 2>> ../logs/error_log.txt
        echo "Aligning concatenated seqs for Step 2 initial phylogeny with $ALIGN1 to make ${fasta%.}_ali.fa ..." | tee -a ../logs/out_log.txt ../logs/alignment_log.txt; echo | tee -a ../logs/out_log.txt ../logs/alignment_log.txt; echo | tee -a ../logs/out_log.txt ../logs/alignment_log.txt
        $ALIGN1 ${fasta%.}_bts_ancrs > ${fasta%.}_ali.fa 2>> ../logs/alignment_log.txt
        # Error check: alignment made
        if [[ ! (-s ${fasta%.}_ali.fa) ]]; then echo "Error aligning the concatenated file comprising ./queries_and_outs.fa and cdhit results (${fasta}*cdhit); examine files in hits_fasta/ and ensure MAFFT is working." 2>> ../logs/error_log.txt; exit 1; fi
        # Record what seqs HMMER removed per species
        echo "Recording what seqs HMMER removed..." | tee -a ../logs/out_log.txt; echo | tee -a ../logs/out_log.txt; echo | tee -a ../logs/out_log.txt
        grep ">" ./${fasta%.}_cdhit | sort > temp_cdhit
        grep ">" ./${fasta%.}_hmmr_significant_hits.fa | sort > temp_HMMR
        comm -23 temp_cdhit temp_HMMR >> HMMR_removed_seqs.txt
    done
    grep -w -f HMMR_removed_seqs.txt ../logs/header_translation_table.tsv > ./HMMR_removed_seqs_info.txt
    mv ./HMMR_removed_seqs_info.txt ../logs
else
    # Remove similar seqs, concatenate queries+anchors+hits, & align for each taxa
    echo "No HMMR profile specified; HMMR search skipped." | tee -a ../logs/out_log.txt
    for fasta in ./*_hits.fa
    do
        # run CD-HIT to remove seqs with 98% or higher seq similarity
        echo "Removing seqs >/= 98% identical via CD-HIT (-c 0.98 -n 5) on $fasta ..." | tee -a ../logs/out_log.txt
        cd-hit -i $fasta -o ${fasta%.}_cdhit -c 0.98 -n 5 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
        # Concatenate cd-hit results with outgroups+queries, align via MAFFT
        echo | tee -a ../logs/out_log.txt; echo "Concatenating seqs from ../queries_and_outs.fa & ./${fasta%.}_cdhit ..." | tee -a ../logs/out_log.txt
        cat ../queries_and_outs.fa ./${fasta%.}_cdhit > ${fasta%.}_bts_ancrs 2>> ../logs/error_log.txt
        echo "Aligning concatenated seqs for Step 2 initial phylogeny with $ALIGN1 to make ${fasta%.}_ali.fa ..." | tee -a ../logs/out_log.txt ../logs/alignment_log.txt; echo | tee -a ../logs/out_log.txt ../logs/alignment_log.txt; echo | tee -a ../logs/out_log.txt ../logs/alignment_log.txt
        $ALIGN1 ${fasta%.}_bts_ancrs > ${fasta%.}_ali.fa 2>> ../logs/alignment_log.txt
        # Error check: alignment made
        if [[ ! (-s ${fasta%.}_ali.fa) ]]; then echo "Error aligning the concatenated file comprising ./queries_and_outs.fa and cdhit results (${fasta}*cdhit); examine files in hits_fasta/ and ensure MAFFT is working." 2>> ../logs/error_log.txt; exit 1; fi
    done
fi
# Make new dir named align_species/ for all alignments and clean up temp files
mkdir ../align_species
mv ./*ali.fa ../align_species
#rm ./*ancrs ./*hmmrsearch_table.txt ./temp* ./HMMR_removed_seqs.txt ./*_hmmrsearch_temp.txt
rm ./*ancrs ./temp* ./HMMR_removed_seqs.txt ./*_hmmrsearch_temp.txt
echo "Moved $(ls ../align_species/*ali.fa | wc -l) species' alignments to ./align_species" | tee -a ../logs/summary_log.txt ../logs/out_log.txt
cd ..
}
step2_initial_trees() {
# Make ML phylogenies from MAFFT aligned files using IQtree
echo "Making Maximum Liklihood phylogenies for Step 2 per species alignments..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
cd ./align_species
for alignment in *ali.fa
do
    # Determine which, if any, bootstrap or sh-ALRT options the user specified for IQTREE
    if [[ ($BOOT == 0) && ($SHALRT == 0) ]]
    then
        echo "Running IQtree (iqtree -s $alignment -m $MODEL -nt $THREADS)..." | tee -a ../logs/out_log.txt
        iqtree -s $alignment -m $MODEL -nt $THREADS 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
    elif [[ ($BOOT != 0) && ($SHALRT == 0) ]]
    then
        echo "Running IQtree (iqtree -s $alignment -m $MODEL -bb $BOOT -nt $THREADS)..." | tee -a ../logs/out_log.txt
        iqtree -s $alignment -m $MODEL -bb $BOOT -nt $THREADS 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
    elif [[ ($BOOT == 0) && ($SHALRT != 0) ]]
    then
        echo "Running IQtree (iqtree -s $alignment -m $MODEL -alrt $SHALRT -nt $THREADS)..." | tee -a ../logs/out_log.txt
        iqtree -s $alignment -m $MODEL -alrt $SHALRT -nt $THREADS 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
    elif [[ ($BOOT != 0) && ($SHALRT != 0) ]]
    then
        echo "Running IQtree (iqtree -s $alignment -m $MODEL -alrt $SHALRT -bb $BOOT -nt $THREADS)..." | tee -a ../logs/out_log.txt
        iqtree -s $alignment -m $MODEL -alrt $SHALRT -bb $BOOT -nt $THREADS 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
    fi
    # Error check: IQtree treefile produced
    if [[ ! (-s $alignment.treefile) ]]; then echo "Error making ML Iq_tree for $alignment during step 2: making unfocused, per-species gene trees. Possible sources include the alignments in ./align_species, the IQtree installation, or improper tree-related parameters entered for this run." >> ../logs/error_log.txt; exit 1; fi
done
# organize iqTree tree output files
echo "Moving $(ls *.treefile | wc -l) IQtree output trees to IQ_out_species/..." | tee -a ../logs/summary_log.txt ../logs/out_log.txt
mkdir IQ_out_species; mv *.treefile IQ_out_species
echo "Moving other IQtree files to IQ_out_extra_files/..." | tee -a ../logs/out_log.txt
mkdir IQ_out_extra_files; mv *.fa.* IQ_out_extra_files; cd ..
}



step3_focus() {
# Use tree_editor to select subtree sequences (IE the focusing step)
echo -------------------------------- | tee -a ./logs/out_log.txt
echo "STEP 3: GENERATING FOCUSED PER SPECIES TREES ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
cd ./align_species/IQ_out_species
echo "Running tree_editor.R..." | tee -a ../../logs/out_log.txt
Rscript ../../subscripts/tree_editor.R ../../$CLADE 2>> ../../logs/error_log.txt | tee -a ../../logs/out_log.txt
for tips_file in ./*.txt
do
    # Error check: tips file exists
    if [[ ! (-s $tips_file) ]]; then echo "Error occured when generating output txt file from $tips_file for tree_editor.R" 2>> ../../logs/error_log.txt; exit 1; fi
done
# Organize tree_editor.R out files
mkdir tree_editor_out; mv ./*_tips.txt ./tree_editor_out; mv ./*_focused.treefile ./tree_editor_out; mv ./Rplots.pdf .tree_editor_trees.pdf; mv tree_editor_trees.pdf ../../logs; cd ../../
}

step4_extract_fastas() {
echo -------------------------------- | tee -a ./logs/out_log.txt
echo "STEP 4: CONCATENATING FOCUSED PROTEIN TREE SEQUENCES ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
# Retrieve fasta seqs for seqs included in tree_editor.R output trees
cd ./align_species/IQ_out_species/tree_editor_out
for file in ./*tips.txt
do
    echo "Extracting subtree fasta sequences for $file..." | tee -a ../../../logs/out_log.txt
    genus=$(echo $file | cut -d '_' -f 1)
    grep -A 1 -w --no-group-separator -f $file ../../../fixed_fastas/${genus}*_FIX.fa > ${genus}_tipseqs.fa 2>> ../../../logs/error_log.txt
    # Error check: fastas extracted correctly
    if [[ ! (-e ./${genus}_tipseqs.fa) ]]; then echo "Error when extracting $file fasta seqs from /fixed_fastas/${genus}*_FIX.fa to make ${genus}_tipseqs.fa. Check ./align_species/IQ_out_species/tree_editor_out." >> ../../../logs/error_log.txt; exit 1; fi
    if [[ ! (-s ./${genus}_tipseqs.fa) ]]; then echo "WARNING: there are no seqs in ${genus}_tipseqs.fa; this species will not be represented in phylogeny." | tee -a ../../../logs/summary_log.txt ../../../logs/out_log.txt; fi
    # Record what seqs focusing removed per species by comparing pre and post focus datasets
    echo "Recording what non-query seqs focusing removed..." | tee -a ../../../logs/out_log.txt; echo | tee -a ../../../logs/out_log.txt; echo | tee -a ../../../logs/out_log.txt
    grep ">" ../../${genus}*fa_ali.fa | sort > temp_unfocused
    grep ">" ../../../$QUERY > temp_query; grep ">" ../../../$OUT >> temp_query; grep ">" ./${genus}_tipseqs.fa > temp_R
    cat temp_query | sort | uniq > temp_query2
    cat temp_query2 temp_R | sort > temp_focused
    comm -23 ./temp_unfocused ./temp_focused > focus_removed_TEMP.txt
    count=$(grep -w -f focus_removed_TEMP.txt ../../../logs/header_translation_table.tsv | wc -l)
    echo "Focusing removed $count sequences from $file..." >> ../../../logs/focus_removed_seqs.txt
    grep -w -f focus_removed_TEMP.txt ../../../logs/header_translation_table.tsv >> ../../../logs/focus_removed_seqs.txt
done
# Organize extracted tipseq fastas
rm ./temp* ./focus_removed_TEMP.txt
cd ../../../; mkdir ./tip_seqs; mv ./align_species/IQ_out_species/tree_editor_out/*_tipseqs.fa ./tip_seqs
}
step4_concatenate_seqs() {
# Concatenate all fastas together, making the "total tree" dataset
mkdir final_tree_dataset; cd final_tree_dataset
echo "Step 4 concatenating tip_seqs and queries_and_outs.fa ..." | tee -a ../logs/out_log.txt
cat ../tip_seqs/*.fa ../queries_and_outs.fa > ./concat_tip_seqs.fa 2>> ../logs/error_log.txt
# Remove duplicate seqs derived from concatenating all the tip_seq fastas
echo "Removing 100% identical seqs from concatenated subtree file via CD-HIT (-c 1.0 -n 5) ..." | tee -a ../logs/out_log.txt
cd-hit -i ./concat_tip_seqs.fa -o ./concat_tip_seqs_cdhit.fa -c 1.0 -n 5 2>> ../logs/error_log.txt | tee -a ../logs/out_log.txt
cd ..
}



step5_final_align() {
echo -------------------------------- | tee -a ./logs/out_log.txt
echo "STEP 5: FILTERING CONCATENATED DATASET: ALIGNMENT EDITING ..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt ./logs/error_log.txt
cd ./final_tree_dataset
echo "Step 5 aligning concat_tip_seqs_cdhit.fa with $ALIGN2 to make concat_tip_seqs_cdhit_ali.fa ..." | tee -a ../logs/out_log.txt ../logs/alignment_log.txt
$ALIGN2 concat_tip_seqs_cdhit.fa > concat_tip_seqs_cdhit_ali.fa 2>> ../logs/alignment_log.txt
cd ..
}



# Run phyfocus functions along with save and error checkpoints
log_start

# Check for user-forced phyfocus restart at specific checkpoints.
if [[ $FORCE == "S2a" ]]
then
    rm -r ./align_species ./tip_seqs ./final_tree_dataset
    echo "Forced Phyfocus rerun at checkpoint S2a: Remaking Step 2 per-species alignments in ./align_species..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
elif [[ $FORCE == "S2t" ]]
then
    rm -r ./align_species/IQ_out_species ./tip_seqs ./final_tree_dataset
    echo "Forced Phyfocus rerun at checkpoint S2t: Remaking Step 2 per-species trees in ./align_species/IQ_out_species..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
elif [[ $FORCE == "S3" ]]
then
    rm -r ./align_species/IQ_out_species/tree_editor_out ./tip_seqs ./final_tree_dataset
    echo "Forced Phyfocus rerun at checkpoint S3: Remaking Step 3 per-species focused trees in ./align_species/IQ_out_species/tree_editor_out..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
elif [[ $FORCE == "S4" ]]
then
    rm -r ./tip_seqs ./final_tree_dataset
    echo "Forced Phyfocus rerun at checkpoint S4: Remaking Step 4 ./tips_seqs directory and combined species dataset in ./final_tree_dataset..." | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi

# STEP 1
# save checkpoint
if [[ ! -e ./logs/header_translation_table.tsv ]]
then
    # phyofcus function
    step1_headers
    # error checkpoint
    if [[ -s ./logs/header_translation_table.tsv ]]; then echo "numerical header --> Step 1 NCBI header translation table complete" | tee -a ./logs/out_log.txt; else echo "Step 1 header translation table could not be made; examine headers.txt files and fixed_fasta files." >> ./logs/error_log.txt; exit 1; fi
    rm all_ncbi_headers.txt ./all_fixed_numerical_headers.txt
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 1 headers" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi
# save checkpoint
if [[ ! -d ./blastout_tables || ! ($(ls ./fixed_fastas | wc -l) == $(ls ./blastout_tables | wc -l)) ]]
then
    # Phyfocus function
    step1_blast
    # error checkpoint
    if [[ $(ls ./fixed_fastas | wc -l) == $(ls ./blastout_tables | wc -l) ]]; then echo "Step 1 BLASTs complete" | tee -a ./logs/out_log.txt; else echo "Step 1 BLAST of peptide FASTAS incomplete (BLAST outputs not made for every species given)." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 1 blast" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi


# STEP 2
# Save checkpoint
if [[ ! -d ./align_species || ! ($(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/*_ali.fa | wc -l)) ]]
then
    # Phyfocus function
    step2_initial_alignments
    # error checkpoint
    if [[ $(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/*_ali.fa | wc -l) ]]; then echo "Step 2 species' peptide alignments complete" | tee -a ./logs/out_log.txt; else echo "Step 2 species' peptide alignments incomplete (alignments not made for every species)." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/summary_log.txt ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 2 species' peptide alignments" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi
# Save checkpoint
if [[ ! -d ./align_species/IQ_out_species || ! ($(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/IQ_out_species/*treefile | wc -l)) ]]
then
    # Phyfocus function
    step2_initial_trees
    # error checkpoint
    if [[ $(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/IQ_out_species/*treefile | wc -l) ]]; then echo "Step 2 species' IQTree runs are complete" | tee -a ./logs/out_log.txt; else echo "Step 2 species' IQTree runs failed to run properly (not all treefiles for each species exist); check ./align_species alignments and ensure IQTree is installed properly. Then rerun phyfocus." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 2 initial IQ trees" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi


# STEP 3
# Save checkpoint
if [[ ! -d ./align_species/IQ_out_species/tree_editor_out || ! ($(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/IQ_out_species/tree_editor_out/*_tips.txt | wc -l)) ]]
then
    # Phyfocus function
    step3_focus
    # error checkpoint
    if [[ $(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./align_species/IQ_out_species/tree_editor_out/*_tips.txt | wc -l) ]]; then echo "Step 3 tree_editor.R focusing complete" | tee -a ./logs/out_log.txt; else echo "Step 3's tree_editor.R did not run properly for every species that had BLAST hits (The number of output tips.txt files does not = the number of fasta files in 'hits_fasta'). Step 3 may have been interrupted or failed to run tree_editor.R on every species' phylogeny. Ensure ./align_species/IQ_out_species phylogeny files (.treefile) are present for each species and program dependencies are working properly. Then rerun phyfocus." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 3 Tree_editor.R focusing" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi


# STEP 4
# Save checkpoint
if [[ ! -d ./tip_seqs || ! ($(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./tip_seqs/* | wc -l)) ]]
then
    # Phyfocus function
    step4_extract_fastas
    # error checkpoint
    if [[ $(ls ./hits_fasta/*FIX_hits.fa | wc -l) == $(ls ./tip_seqs/* | wc -l) ]]; then echo "Step 4 extracting focused seqs complete." | tee -a ./logs/out_log.txt; else echo "Step 4 extracting focused seqs incomplete (fasta for each taxa not present in ./tip_seqs); check tips.txt files in ./align_species/IQ_out_species/tree_editor_out." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 4 Extracting focused seqs" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi
# Save checkpoint
if [[ ! -s ./final_tree_dataset/concat_tip_seqs_cdhit.fa ]]
then
    # Phyfocus function
    step4_concatenate_seqs
    # error checkpoint
    if [[ -s ./final_tree_dataset/concat_tip_seqs_cdhit.fa ]]; then echo "Step 4 focused seq concatenation complete." | tee -a ./logs/out_log.txt; else echo "Step 4 focused seq concatenation failed (concat_tip_seqs_cdhit.fa not made or empty); check if ./tip_seqs is empty and that concat_tip_seqs.fa is in ./final_tree_dataset. Ensure CDHIT is installed properly." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 4 Concatenating focused tree seqs" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi


# STEP 5
# Save checkpoint
if [[ ! -s ./final_tree_dataset/concat_tip_seqs_cdhit_ali.fa ]]
then
    # Phyfocus function
    step5_final_align
    if [[ -s ./final_tree_dataset/concat_tip_seqs_cdhit_ali.fa ]]; then echo "Step 5 alignment of final tree seqs complete" | tee -a ./logs/out_log.txt; else echo "Step 5 alignment of final tree seqs incomplete (concat_tip_seqs_cdhit_ali.fa not made or empty); check ./final_tree_dataset files." >> ./logs/error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt; echo " " | tee -a ./logs/out_log.txt
    # Final outlog statements before user alignment editing.
    echo "Phyfocus alignment is ready for Step 5 user assessment with alignment_editor.py. Please examine 'final_tree_dataset/concat_tip_seqs_cdhit_ali.fa' and identify excessive gaps. From there you can then run alignment editor via ../subscripts/alignment_editor.py. Check alignment_editor.py -h for details." | tee -a ./logs/summary_log.txt ./logs/out_log.txt; echo | tee -a -a ./logs/summary_log.txt ./logs/out_log.txt
    echo "Once alignment editing is complete, runnning IQTree on the realigned output will complete Step 6 of the phyfocus run. An example rigorous approach: iqtree -s <final tree alignment> -m MFP+C60 -alrt 1000 -bb 1000 -nt 24" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 5 alignment of final tree seqs" | tee -a ./logs/summary_log.txt ./logs/out_log.txt
fi
