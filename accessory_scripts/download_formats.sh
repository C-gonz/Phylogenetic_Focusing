#! /bin/bash

# Script for automatically downloading species sequence data and adjusting the file name to begin with the species' genus, as required for proper annotation in phyfocus phylogenies.

help() {
cat << EOF
----------------------------------------------
download formats (a user aid to the phyfocus pipeline)
github repository: https://github.com/C-gonz/Phylogenetic_Focusing
See "README" for detailed instructions on the overall Phyfocus pipeline
----------------------------------------------

Syntax: ./${0##*/} [-h] [-f file.csv]

-h <help>   Display this help and exit

-f CSV    Required; must be a CSV file with 2 columns:
    1st is the species' genus
    2nd is the download link (e.g. FTP) for the peptide dataset.
    An example line containing both columns:
Saccoglossus,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/605/GCF_000003605.2_Skow_1.1/GCF_000003605.2_Skow_1.1_protein.faa.gz

IMPORTANT: if you have more than 1 species from a given genus, you must distinguish the genus column. Example:
Mus fernandoni & Mus musculus --> MusF & MusM

----------------------------------------------
Description:
Download formats was made to help make database prep for phyfocus easier. This script will make every downloaded file's name begin with its corresponding genus, which is required for species tracking and phylogeny annotation in phyfocus.

EOF
}

# Create input option variables
CSV=""

# Create error messages for improper option input
f_error="Option error; path must be to an existing file with data. Format: -q <./fasta_query_file>"

# Handling for option arguments, including improper arguments
while getopts ":hf:" option; do
    case $option in
        h) help; exit 0;;
        f) CSV=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $f_error >&2; exit 1; fi;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done
# Exit program if required options are not specified by user
if [[ -z $CSV ]]; then echo "Missing required option. Syntax: ./${0##*/} [-h for help] [-f file.csv]" >&2; exit 1; fi

mkdir fastas
cd fastas

# For each line in CSV, ID the genus, link, and download file name
content=$(cat ../$CSV)
for line in $content
do
    # define genus from cloumn 1
    genus=$(echo ${line%%,*})
    echo " "; echo " "; echo genus is $genus
    # check if a download link in column 2 was provided
    link=$(echo ${line##*,})
    if [[ ! (-n $link) ]]; then echo "no download link provided"; continue; fi
    # Derive file name from download link and append genus to it
    file=$(echo ${line##*/})
    wget $link
    # check if download worked properly
    if [[ ! (-s ./$file) ]]; then echo "The download for $genus failed. check download link in column 2 of CSV."; continue; fi
    new_name=$(echo "$genus"_"$file")
    mv $file $new_name
done
