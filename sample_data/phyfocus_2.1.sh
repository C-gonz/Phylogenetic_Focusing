#!/bin/bash

help() {
cat << EOF
----------------------------------------------
PhyFocus v2.1 (June 21 2023)
github repository: https://github.com/C-gonz/Phylogenetic_Focusing
See "README" for user instructions on the overall Phyfocus pipeline
----------------------------------------------

Syntax: ./${0##*/} [-h] [-T] [-X] [-q ./file] [-r ./file] [-f ./directory/] [-c ./file] [-H ./file] [-t num] [-e value] [-m string] [-b num] [-a num]

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

-X CLEAN    Optional; Removes all phyfocus-created files in the directory. Is called automatically when running the ./sample_data via -T, allowing for multiple test runs.

----------------------------------------------
PATH Dependencies:
    AWK
    Python3 (including BioPython package -> Bio.AlignIO module)
    R
    NCBI BLAST+ (specifically makeblastdb and blastp)
    CD-HIT
    MAFFT
    IQTREE
    HMMER
Included Dependencies:
    ${0##*/}
    tree_editor.R
    header_translator.py
    alignment_editor.py
    delim_converter.py
----------------------------------------------

Description:
PhyFocus assesses gene family evolution by "focusing," or extracting, the gene clade(s) of interest from a much broader gene phylogeny. By starting with an extensive phylogeny of outgroup clades, Phyfocus decreases the chance that significant gene family relationships are excluded.

PhyFocus consists of 6 steps:
1) IDENTIFYING PROTEIN DATASET PER SPECIES
2) GENERATING PROTEIN TREE PER SPECIES
3) EXTRACTING FOCUSED PROTEIN TREE PER SPECIES
4) CONCATENATING FOCUSED PROTEIN TREE SEQUENCES
5) FILTERING CONCATENATED DATASET: HMMR
6) GENERATING FINAL FOCUSED TREE

PhyFocus requires five user-provided datasets:

1) A query fasta file containing target and anchor proteins (Step 1).
   - Targets represent homologs of the specific gene(s) / gene family being studied.
   - Anchors represent homologs of gene(s) / gene gene families closely related to the targets.
   - All query proteins can come from a single organism, but broader sampling from taxa of interest may aid BLAST.
   - There is a minimum of 2 for target + anchor sequences (more are recommended).

2) A directory of protein FASTA files for each species assessed in the phylogeny (Step 1).
    - Protein sequences are ideally derived from whole-genome data or thorough transcriptomes.
    - All FASTA file names MUST begin with the species' genus name and underscore: genus_

3) A fasta file containing the rooting proteins to root the per-species unfocused trees (Steps 2 & 3).
    - Roots represent the closest known outgroup to target & anchor gene families; used to root the trees.
    - There is a minimum of 2 sequences required for roots.

4) A Tab Seperated Values (.tsv) file for phylogeny focusing (Step 3).
- There should be no column or row headers in the table.
- Column 1 gives FASTA ">" header names for AT LEAST 2 root proteins in the query fasta file.
- Column 2 gives FASTA ">" header names for AT LEAST 1 target AND 1 related protein in the query fasta file. Two of each is recommended.
- Below is an example TSV for a phylogeny of taste receptors within the glutamate/class-C GPCR family:
Homo_GABBR1    Homo_T1R1
Homo_GABBR2    Homo_mGLUR6
Danio_GABBR1
Danio_GABBR2
- Related and root sequences are best chosen by reference to previous phylogenies. The glutamate example was informed by the following: Fredriksson, R., Lagerström, M. C., Lundin, L. G., & Schiöth, H. B. (2003). The G-protein-coupled receptors in the human genome form five main families. Phylogenetic analysis, paralogon groups, and fingerprints. Molecular pharmacology, 63(6), 1256-1272.

5) A FASTA protein alignment characterizing key conserved domains and motifs (Step 5).
    - An alignment of the query file can be used.
    - Rigorous alignment (e.g., using MAFFT Linsi) is preferred.
    - Builds the HMMR profile for HMMR filtering to remove any fundamentally different sequences,
      such as proteins lacking a 7TM domain in a GPCR phylogeny.
EOF
}

# Function for running a test demo of phyfocus using provided data
test() {
cd ./sample_data
clean
./${0##*/} -q ./query_file.fa -f fasta_proteins/ -r roots_file.fa -H hmmr_ali.fa -c tree_tips.tsv -m LG -b 1000 -a 1000
}

# Function for removing generated files from failed runs,etc.
clean() {
if [[ -e ./tree_editor.R ]]
then
    starting_num=$(ls | wc -l)
    mv final_tree_dataset/filtering_output/alignment_editor.py ./ 2> /dev/null
    rm -r ./align ./blastdb ./blastout_maxseqs ./fixed_fastas ./hit1_accessions ./hits_fasta ./tip_seqs ./final_tree_dataset 2> /dev/null
    rm ./out_log.txt ./error_log.txt ./header_translation_table.tsv temp.txt all_fixed_numerical_headers.txt all_ncbi_headers.txt formatted_query.fa 2> /dev/null
    clean_num=$(ls | wc -l)
    diff=$(($starting_num -$clean_num))
    echo "Removed $diff files. phyfocus directory now contains $clean_num files. 10 files (5 from user, 5 included scripts) are required to run PhyFocus."
fi
}

# Create input option variables
QUERY=""
FASTAS=""
ROOT=""
CLADE=""
HMMR=""
THREADS=24
EVALUE="1e-05"
MODEL="MFP+C60"
BOOT=1000
IQALRT=1000
regex_num="^[0-9]+$"
regex_enum="^[0-9]+e-[0-9]+$"
user_args=$@

# Create error messages for improper option input
q_error="Option error; path must be to an existing file with data. Format: -q <./fasta_query_file>"
f_error="Option error; path must be to an existing directory. Format: -f <./species_fasta_proteins_directory/>"
r_error="Option error; path must be to an existing file with data. Format: -r <./fasta_roots_file>"
c_error="Option error; path must be to an existing file with data. Format: -c <./tsv_table>"
H_error="Option error; path must be to an existing file with data. Format: -H <./fasta_alignment_file>"
t_error="Option error; thread usage must be an integer. format: -t <integer>"
e_error="Option error; need value in scientific notation. Format: -e <nums>e-<nums>"
m_error="Option error; string needs to be an IQtree model. Format: -m <value>"
b_error="Option error; bootstrap replicates must be an integer. Format: -b <integer>"
a_error="Option error; SH-aLRT replicates must be an integer. Format: -a <integer>"

# Handling for option arguments, including improper arguments
while getopts ":hq:f:r:c:H:t:e:m:b:a:TX" option; do
    case $option in
        h) help; exit 0;;
        q) QUERY=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $q_error >&2; exit 1; fi;;
        f) FASTAS=$OPTARG; if [[ ! -d $OPTARG ]]; then echo $f_error >&2; exit 1; fi;;
        r) ROOT=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $r_error >&2; exit 1; fi;;
        c) CLADE=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $c_error >&2; exit 1; fi;;
        H) HMMR=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $H_error >&2; exit 1; fi;;
        t) THREADS=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $t_error >&2; exit 1; fi;;
        e) EVALUE=$OPTARG; if [[ ! $OPTARG =~ $regex_enum ]]; then echo $E_error >&2; exit 1; fi;;
        m) MODEL=$OPTARG; if [[ -z $OPTARG ]]; then echo $m_error >&2; exit 1; fi;;
        b) BOOT=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $b_error >&2; exit 1; fi;;
        a) IQALRT=$OPTARG; if [[ ! $OPTARG =~ $regex_num ]]; then echo $a_error >&2; exit 1; fi;;
        T) test; exit 0;;
        X) clean; exit 0;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done
# Print help menu if no options are specified.
if [[ -z $QUERY && -z $FASTAS ]]; then help; exit 0; fi
# Exit program if required options are not specified by user
if [[ (-z $QUERY) || (-z $FASTAS) || (-z $HMMR) || (-z $ROOT) ]]; then echo "Missing option. Syntax: ./${0##*/} [-q file] [-f directory] [-H file]" >&2; exit 1; fi



log_start() {
echo -------------------------------- | tee -a out_log.txt | tee -a error_log.txt
echo "PhyFocus v2.1 (May ?? 2023)" | tee -a out_log.txt | tee -a error_log.txt
echo "github repository: https://github.com/C-gonz/Phylogenetic_Focusing" | tee -a out_log.txt | tee -a error_log.txt
echo "Phyfocus run started:" $(date) | tee -a out_log.txt | tee -a error_log.txt
echo "Phyfocus run's home directory:" $(pwd) | tee -a out_log.txt | tee -a error_log.txt
}



step1_headers() {
echo -------------------------------- | tee -a out_log.txt
echo "STEP 1: IDENTIFYING PROTEIN DATASET PER SPECIES" | tee -a out_log.txt
echo "Unzipping any gzip files in $FASTAS ..." | tee -a out_log.txt
gzip -d ${FASTAS}/*.gz | tee -a out_log.txt 2>> error_log.txt | tee -a out_log.txt
# For each species fasta file, use parameter expansion to extract genus name, then AWK
# to change each sequence header to a "genus_#####" numerical header, and write each sequence
# on 1 line without STOP signals.
# Syntax: "NR==1" prevents line feed at start of file, "genus "_%05d\n", ++i" creates
# genus_numerical headers, "else {gsub("\\*", "", $0); printf $0}" ensures every sequence line is written
# on 1 line and has removed any STOP (*) signals.
echo "Fixing FASTA headers and seqs ..." | tee -a out_log.txt
for original_filename in $FASTAS/*; do file_name=${original_filename##*/}; awk -v genus=${file_name%%_*} '{if(NR==1) {printf ">" genus "_%05d\n", ++i "\n"} else {if($0 ~ /^>/) {printf "\n" ">" genus "_%05d\n", ++i "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf $0 "\n"}' $original_filename > ${file_name}_FIX.fa 2>> error_log.txt; done
# Conduct same line fix for user query and root files
awk '{if(NR==1) {printf $0 "\n"} else {if($0 ~ /^>/) {printf "\n" $0 "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf $0 "\n"}' $QUERY > formatted_query.fa
awk '{if(NR==1) {printf $0 "\n"} else {if($0 ~ /^>/) {printf "\n" $0 "\n"} else {gsub("\\*", "", $0); printf $0}}} END {printf $0 "\n"}' $ROOT > formatted_roots.fa
# Generate a file containing all NCBI headers for header translation table
cat ./$FASTAS/* | grep ">" > all_ncbi_headers.txt 2>> error_log.txt
# make directory for fixed fastas, then move FIX files
echo "moving modified header fastas to ./fixed_fastas/ ..." | tee -a out_log.txt
mkdir ./fixed_fastas
mv *_FIX.fa ./fixed_fastas
echo "Successfully moved $(ls ./fixed_fastas | wc -l) files to ./fixed_fastas/" | tee -a out_log.txt
echo "FASTA Header Fix Complete" | tee -a out_log.txt
# Make a translation table for numerical headers --> informative NCBI headers
echo "Making numerical sequence header translation table." | tee -a out_log.txt
cat ./fixed_fastas/* | grep ">" > all_fixed_numerical_headers.txt 2>> error_log.txt
./header_translator.py all_ncbi_headers.txt all_fixed_numerical_headers.txt 2>> error_log.txt | tee -a out_log.txt
}
step1_blast() {
# Format blast databases from $FASTAS & run BLASTp
echo "Blasting $FASTAS peptide databases with formatted_query.fa using (blastp -db FILE -query formatted_query.fa -evalue $EVALUE -num_threads $THREADS -out ${species_pep%.*}_ref_blastout -outfmt 6 -max_target_seqs 1000) ..." | tee -a out_log.txt
for species_pep in ./fixed_fastas/*
do
    makeblastdb -in $species_pep -parse_seqids -out ${species_pep%.*}_db -dbtype prot 2>> error_log.txt | tee -a out_log.txt
    echo "formatted BLAST database for $species_pep proteins..." | tee -a out_log.txt
    # blastp query file against species_pep database
    echo "BLASTING formatted_query.fa sequences against $species_pep database ..." | tee -a out_log.txt
    blastp -db ${species_pep%.*}_db -query formatted_query.fa -evalue $EVALUE -num_threads $THREADS -out ${species_pep%.*}_ref_blastout -outfmt 6 -max_target_seqs 1000 2>> error_log.txt | tee -a out_log.txt
    # Error check: each blast completed
    if [[ ! (-s ${species_pep%.*}_ref_blastout) ]]; then echo "Error blasting for $species_pep. See error log and blast files in ./fixed_fastas" >> error_log.txt; exit 1; fi
done
# move blast dbs & blast outputs to new directories
echo "moving db's to ../blastdb" | tee -a out_log.txt
mkdir ./blastdb
mv ./fixed_fastas/*_db* ./blastdb
echo "moving blastout files to ./blastout_maxseqs" | tee -a out_log.txt
mkdir ./blastout_maxseqs
mv ./fixed_fastas/*blastout ./blastout_maxseqs/
}



step2_initial_alignments() {
echo -------------------------------- | tee -a out_log.txt
echo "STEP 2: GENERATING PROTEIN TREE PER SPECIES ..." | tee -a out_log.txt
echo "Writing blast hit FASTA files ..." | tee -a out_log.txt
# Retrieve accession numbers for each blastp output & move to new dir
for blastfile in ./blastout_maxseqs/*blastout; do echo "Retrieving blast hit accessions for $blastfile ..." | tee -a out_log.txt; cut -f 2 $blastfile > ${blastfile%.}_hits.txt; done
# Move all hit accession files to new dir called hit1_accessions
mkdir ./hit1_accessions
mv ./blastout_maxseqs/*_hits.txt ./hit1_accessions
echo "Successfully moved $(ls ./hit1_accessions | wc -l) hit files to hit1_accessions" | tee -a out_log.txt
# Retrieve blast hit seqeunces from *_FIX.fa files for each taxa
echo "Retrieving blast hit seqeunces from fixed FASTA files ..." | tee -a out_log.txt
for hitfile in ./hit1_accessions/*_hits.txt
do
    # store hit accessions in a temp FASTA file & store fixed_fasta file name
    awk '{printf ">" $0 "\n"}' < $hitfile > temp.txt
    filenamefull=${hitfile##*/}
    in_filename="${filenamefull%_ref_blastout_hits*}.fa"
    hits_fasta="${filenamefull%_ref_blastout_hits*}_hits.fa"
    # Retrieve blast hit seqs from *_FIX.fa files for each taxa
    echo "pulling FASTA seqs from fixed_fastas/ for hit accessions in $in_filename ..." | tee -a out_log.txt
    grep -A 1 -f temp.txt --no-group-separator ./fixed_fastas/$in_filename > $hits_fasta 2>> error_log.txt
    # Error check: blast hit seqs extracted
    if [[ ! (-s $hits_fasta) ]]; then echo "Error retrieving hit squences for hits_$filename; examine fixed fasta file." >> error_log.txt; exit 1; fi
done
rm ./temp.txt
# Move retrieved hit fasta seqs to the new directory hits_fasta/
mkdir hits_fasta
mv *_hits.fa ./hits_fasta
echo "Successfully moved $(ls ./hits_fasta/*_hits.fa | wc -l) hit sequence fasta files to hits_fasta" | tee -a out_log.txt
# Remove similar seqs, concatenate queries+anchors+hits, and then align for each taxa
for fasta in ./hits_fasta/*_hits.fa
do
    # run CD-HIT to remove seqs with 98% or higher seq similarity
    echo "Removing seqs >/= 98% identical via CD-HIT (-c 0.98 -n 5) on $fasta ..." | tee -a out_log.txt
    cd-hit -i $fasta -o ${fasta%.}_cdhit -c 0.98 -n 5 2>> error_log.txt | tee -a out_log.txt
    # after concatenation, align via MAFFT
    echo | tee -a out_log.txt; echo "Concatenating bait+anchor seqs (./formatted_query.fa) and cdhit results (${fasta%.}_cdhit) ..." | tee -a out_log.txt
    cat ./formatted_query.fa ./formatted_roots.fa ${fasta}*cdhit > ${fasta%.}_bts_ancrs 2>> error_log.txt
    echo "Aligning concatenated seqs for Step 2 initial phylogeny with MAFFT linsi ..." | tee -a out_log.txt; echo | tee -a out_log.txt; echo | tee -a out_log.txt
    linsi ${fasta%.}_bts_ancrs > ${fasta%.}_ali.fa # 2>> error_log.txt | tee -a out_log.txt DOESN'T WORK
    # Error check: alignment made
    if [[ ! (-s ${fasta%.}_ali.fa) ]]; then echo "Error aligning squences from bait+anchor seqs (./formatted_query.fa) and cdhit results (${fasta}*cdhit) concatenated file; examine files in hits_fasta/." >> error_log.txt; exit 1; fi
done
# Make new dir named align/ for all alignments and clean up temp cat files
mkdir align
mv ./hits_fasta/*ali.fa ./align
rm ./hits_fasta/*ancrs
echo "Moved $(ls ./align | wc -l) species' MAFFT linsi alignments to ./align" | tee -a out_log.txt
}
step2_initial_trees() {
# Make ML phylogenies from MAFFT aligned files using IQtree
echo "Making Maximum Liklihood phylogenies for Step 2 per species alignments..." | tee -a out_log.txt
cd ./align
for alignment in *ali.fa
do
    echo "Running IQtree (iqtree -s $alignment -m $MODEL -alrt $IQALRT -bb $BOOT -nt $THREADS)..." | tee -a ../out_log.txt
    iqtree -s $alignment -m $MODEL -alrt $IQALRT -bb $BOOT -nt $THREADS 2>> ../error_log.txt | tee -a ../out_log.txt
    # Error check: IQtree treefile produced
    if [[ ! (-s $alignment.treefile) ]]; then echo "Error making ML Iq_tree for $alignment. See error log and alignment in ./align" >> error_log.txt; exit 1; fi
done
# organize iqTree tree output files
echo "Moving $(ls *.treefile | wc -l) IQtree output trees to IQ_out/..." | tee -a ../out_log.txt
mkdir IQ_out; mv *.treefile IQ_out
echo "Moving other IQtree files to IQ_out_extra_files/..." | tee -a ../out_log.txt
mkdir IQ_out_extra_files; mv *.fa.* IQ_out_extra_files; cd ..
}



step3_focus() {
# Use tree_editor to select subtree sequences (IE the focusing step)
echo -------------------------------- | tee -a ./out_log.txt
echo "STEP 3: GENERATING FOCUSED SEQS PER SPECIES TREE ..." | tee -a ./out_log.txt
cd ./align/IQ_out
echo "Running tree_editor.R..." | tee -a ../../out_log.txt
Rscript ../../tree_editor.R ../../$CLADE 2>> ../../error_log.txt | tee -a ../../out_log.txt
# Remove unwanted space and prefix to tree_editor out files: (GENUS_SPECIES_FIX.fa_ali.fa_tips.txt)
echo "Fixing tree_editor.R outfile names..." | tee -a ../../out_log.txt
rename '.faa_FIX.fa_ali.fa' '' ./*.txt
for tips_file in ./*.txt
do
    # Error check: tips file exists
    if [[ ! (-s $tips_file) ]]; then echo "Error occured when generating output txt file from $tips_file for tree_editor.R" 2>> ../../error_log.txt; exit 1; fi
    mv "$tips_file" "${tips_file/' '}"
done
# Organize tree_editor.R out files
mkdir tree_editor_out; mv ./*_tips.txt ./tree_editor_out; mv ./*_focused.treefile ./tree_editor_out; cd ../../
}



step4_extract_fastas() {
echo -------------------------------- | tee -a ./out_log.txt
echo "STEP 4: CONCATENATING FOCUSED PROTEIN TREE SEQUENCES ..." | tee -a ./out_log.txt
# Retrieve fasta seqs for seqs included in tree_editor.R output trees
cd ./align/IQ_out/tree_editor_out
for file in ./*tips.txt
do
    echo "Extracting subtree fasta sequences for $file..." | tee -a ../../../out_log.txt
    genus=$(echo $file | cut -d '_' -f 1)
    grep -A 1 --no-group-separator -f $file ../../../fixed_fastas/${genus}*_FIX.fa > ${genus}_tipseqs.fa 2>> ../../../error_log.txt
    # Error check: fastas extracted correctly
    if [[ ! (-e ./${genus}_tipseqs.fa) ]]; then echo "Error when extracting $file fasta seqs from /fixed_fastas/${genus}*_FIX.fa to make ${genus}_tipseqs.fa. Check ./align/IQ_out/tree_editor_out." >> ../../../error_log.txt; exit 1; fi
done
# Organize extracted tipseq fastas
cd ../../../; mkdir ./tip_seqs; mv ./align/IQ_out/tree_editor_out/*_tipseqs.fa ./tip_seqs
}
step4_concatenate_seqs() {
# Concatenate all fastas together, making the "total tree" dataset
mkdir final_tree_dataset; cd final_tree_dataset
echo "Step 4 concatenating tip_seqs and bait+anchor seqs ..." | tee -a ../out_log.txt
cat ../tip_seqs/*.fa ../formatted_query.fa > ./concat_tip_seqs.fa 2>> ../error_log.txt
# Remove duplicate seqs derived from concatenating all the tip_seq fastas
echo "Removing 100% identical seqs from concatenated subtree file via CD-HIT (-c 1.0 -n 5) ..." | tee -a ../out_log.txt
cd-hit -i ./concat_tip_seqs.fa -o ./concat_tip_seqs_cdhit.fa -c 1.0 -n 5 2>> ../error_log.txt | tee -a ../out_log.txt
cd ..
}



step5_filtering() {
echo -------------------------------- | tee -a ./out_log.txt
echo "STEP 5: FILTERING CONCATENATED DATASET: HMMR & ALIGNMENT EDITING ..." | tee -a ./out_log.txt
cd ./final_tree_dataset
# Build a HMMR gene homolog profile
echo "Building HMMR profile (hmmbuild) & matching against subtree seqs (hmmsearch) ..." | tee -a ../out_log.txt
hmmbuild hmmr_profile.hmm ../$HMMR 2>> ../error_log.txt | tee -a ../out_log.txt
# Use HMMR profile to identify target homologs in the focused subtree dataset
hmmsearch hmmr_profile.hmm ./concat_tip_seqs_cdhit.fa > hmmrsearch_results.out 2>> ../error_log.txt | tee -a ../out_log.txt
# Error check: hmmr profile and search files
if [[ ! ((-s hmmr_profile.hmm) && (-s hmmrsearch_results.out)) ]]; then echo "Error occured when generating HMMR profile or search results. See error log and " 2>> ../error_log.txt; exit 1; fi
# Generate a tsv table of significantly matched HMMR hits from space delimited hmmsearch results
echo "Extracting HMMR significantly matched sequences ..." | tee -a ../out_log.txt
grep -E "^ *?[0-9]\.?[0-9]?e-[0-9]*" hmmrsearch_results.out > hmmr_hits_table.txt 2>> ../error_log.txt
../delim_converter.py 2>> ../error_log.txt | tee -a ../out_log.txt
# Extract match headers from HMMR tsv, then extract their seqs in FASTA format
cut -f 9 hmmr_hits_table.tsv > hmmr_significant_hits_headers.txt 2>> ../error_log.txt
grep -A 1 --no-group-separator -f hmmr_significant_hits_headers.txt ./concat_tip_seqs_cdhit.fa > hmmr_significant_hits.fa 2>> ../error_log.txt
# Organize HMMR filtering files
mkdir ./filtering_output; mv hmmr* ./filtering_output; cd ../
}
step5_final_align() {
cd ./final_tree_dataset/filtering_output
# Add rooting seqs to final dataset
cat ../../formatted_roots.fa ./hmmr_significant_hits.fa > final_tree_seqs.fa
echo "Step 5 aligning HMMR seqs with MAFFT FFT-NS-2 (fast; progressive method) ..." | tee -a ../../out_log.txt
mafft --retree 2 --maxiterate 0 final_tree_seqs.fa > final_tree_seqs_ali.fa # 2>> ../../error_log.txt | tee -a ../../out_log.txt DOESNT WORK; only puts output in error log and overwrites
# Move alignment_editor.py so it can be used on the final tree dataset
mv ../../alignment_editor.py .; cd ../../
}



# Run phyfocus functions along with save and error checkpoints
log_start

# STEP 1
# save checkpoint
if [[ ! -e ./header_translation_table.tsv ]]
then
    # phyofcus function
    step1_headers
    # error checkpoint
    if [[ -s header_translation_table.tsv ]]; then echo "numerical header --> Step 1 NCBI header translation table complete" | tee -a out_log.txt; else echo "Step 1 header translation table could not be made; examine *headers.txt files and fixed_fasta files." >> error_log.txt; exit 1; fi
    rm all_ncbi_headers.txt ./all_fixed_numerical_headers.txt
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 1 headers" | tee -a out_log.txt
fi
# save checkpoint
if [[ ! -d ./blastout_maxseqs/ ]]
then
    # Phyfocus function
    step1_blast
    # error checkpoint
    if [[ $(ls ./fixed_fastas/ | wc -l) == $(ls ./blastout_maxseqs/ | wc -l) ]]; then echo "Step 1 BLASTs complete" | tee -a out_log.txt; else echo "Step 1 BLAST of peptide FASTAS incomplete; see error log." >> error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 1 blast" | tee -a out_log.txt
fi


# STEP 2
# Save checkpoint
if [[ ! -d ./align/ ]]
then
    # Phyfocus function
    step2_initial_alignments
    # error checkpoint
    if [[ $(ls ./fixed_fastas | wc -l) == $(ls ./align | wc -l) ]]; then echo "Step 2 species' peptide alignments complete" | tee -a ./out_log.txt; else echo "Step 2 species' peptide alignments incomplete; see error log and ./align." >> ./error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 2 species' peptide alignments" | tee -a ./out_log.txt
fi
# Save checkpoint
if [[ ! -d ./align/IQ_out ]]
then
    # Phyfocus function
    step2_initial_trees
    # error checkpoint
    if [[ $(ls ./fixed_fastas | wc -l) == $(ls ./align/IQ_out | wc -l) ]]; then echo "Step 2 species' IQTree runs are complete" | tee -a ./out_log.txt; else echo "Step 2 species' IQTree runs failed to run properly; see error log and ./align." >> ./error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 2 initial IQ trees" | tee -a ./out_log.txt
fi


# STEP 3
# Save checkpoint
if [[ ! -d ./align/IQ_out/tree_editor_out ]]
then
    # Phyfocus function
    step3_focus
    # error checkpoint
    if [[ $(ls ./fixed_fastas | wc -l) == $(ls ./align/IQ_out/tree_editor_out/*_tips.txt | wc -l) ]]; then echo "Step 3 tree_editor.R focusing complete" | tee -a out_log.txt; else echo "Error occured when generating output txt files for Step 3 tree_editor.R focusing; see error log and ./align/IQ_out." >> error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 3 Tree_editor.R focusing" | tee -a out_log.txt
fi


# STEP 4
# Save checkpoint
if [[ ! -d ./tip_seqs ]]
then
    # Phyfocus function
    step4_extract_fastas
    # error checkpoint
    if [[ $(ls ./fixed_fastas | wc -l) == $(ls ./align/*_ali.fa | wc -l) ]]; then echo "Step 4 extracting focused seqs complete" | tee -a out_log.txt; else echo "Step 4 extracting focused seqs incomplete; see error log and ./align." >> error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 4 Extracting focused seqs" | tee -a out_log.txt
fi
# Save checkpoint
if [[ ! -s ./final_tree_dataset/concat_tip_seqs_cdhit.fa ]]
then
    # Phyfocus function
    step4_concatenate_seqs
    # error checkpoint
    if [[ -s ./final_tree_dataset/concat_tip_seqs_cdhit.fa ]]; then echo "Step 4 focused seq concatenation complete." | tee -a ./out_log.txt; else echo "Step 4 focused seq concatenation failed; see error log and ./final_tree_dataset." >> ./error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 4 Concatenating focused tree seqs" | tee -a ./out_log.txt
fi


# STEP 5
# Save checkpoint
if [[ ! -d ./final_tree_dataset/filtering_output ]]
then
    # Phyfocus function
    step5_filtering
    # error checkpoint
    if [[ -s ./final_tree_dataset/filtering_output/hmmr_significant_hits.fa ]]; then echo "Step 5 HMMR filtering of focused seqs complete" | tee -a out_log.txt; else echo "Step 5 HMMR filtering of focused seqs incomplete; see error log and ./align." >> error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 5 HMMR filtering of focused seqs" | tee -a out_log.txt
fi
# Save checkpoint
if [[ ! -s ./final_tree_dataset/filtering_output/hmmr_significant_hits_ali.fa ]]
then
    # Phyfocus function
    step5_final_align
    # error checkpoint
    if [[ -s ./final_tree_dataset/filtering_output/hmmr_significant_hits_ali.fa ]]; then echo "Step 5 alignment of final tree seqs complete" | tee -a out_log.txt; else echo "Step 5 alignment of final tree seqs incomplete; see error log and ./final_tree_dataset/filtering_output HMMR files." >> error_log.txt; exit 1; fi
    echo "PhyFocus elapsed time:" $((($SECONDS/60))) "minutes ($((($SECONDS/60)/60)) hours)" | tee -a out_log.txt; echo " " | tee -a out_log.txt; echo " " | tee -a out_log.txt
else
    echo "Restarted Phyfocus run: skipped Step 5 alignment of final tree seqs" | tee -a out_log.txt
fi


# Final outlog statements before user alignment editing.
echo "Phyfocus alignment is ready for user assessment with alignment_editor.py. Please examine 'final_tree_dataset/filtering_output/hmmr_significant_hits_ali.fa' and identify problematic gaps. See alignment_editor.py -h for details." | tee -a ./out_log.txt; echo | tee -a ./out_log.txt

echo "Once alignment editing is complete, runnning IQTree on the realigned output will complete the phyfocus run. Recommended approach: iqtree -s hmmr_significant_hits_ali.fa_edited_ali.fa -m MFP+C60 -alrt 1000 -bb 1000 -nt 24" | tee -a ./out_log.txt
