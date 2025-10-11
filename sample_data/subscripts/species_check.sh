#! /bin/bash

# Create input option variables
FILE=""
OUT=""

# Create error messages for improper option input
f_error="Option error; final alignment_editor output file not found."
o_error="Option error; final alignment_editor output file not found."


# Handling for option arguments, including improper arguments
while getopts ":f:o:" option; do
    case $option in
        f) FILE=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $f_error >&2; exit 1; fi;;
        o) OUT=$OPTARG; if [[ ! -s $OPTARG ]]; then echo $f_error >&2; exit 1; fi;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    esac
done
# Exit program if required options are not specified by user
if [[ (-z $FILE) || (-z $OUT) ]]; then echo "Missing required option. -f or -o" >&2; exit 1; fi



# obtain genus names from all species provided
for original_filename in ../fixed_fastas/*; do file_name=${original_filename##*/}; genus=${file_name%%_*}; echo $genus >> temp; done

# list the species' provided genera
genus_list=$(cat temp); rm temp

# Generate report of species' presence or absence in the final focused alignment
for genus in $genus_list; do echo $genus >> temp; grep -q $genus ./$FILE; echo $? >> temp; done

# report any missing species
echo  | tee -a ./$OUT
echo "The following species are absent in the edited alignment:" | tee -a ./$OUT
report=$(cat temp)
any_missing=$(grep 1 temp | wc -l)
if [[ $any_missing == 0 ]]; then echo "No species are absent." | tee -a ./$OUT; fi
for line in $report
do
    if [[ $line == 1 ]]; then echo $previous | tee -a ./$OUT; fi
    previous=$line
done
rm temp
