#! /usr/bin/env python3
# SCRIPT PURPOSE (V1.0)
# This script is part of recirpocal blast and phyfocus pipelines, and serves to
# correlate numerical (AKA "fixed") fasta seq headers (which are used because
# BLAST cannot handle long headers) with the original and informative NCBI headers.
# See "Bioinfo_log_hagfish_transcriptome_mining" Chris Gonzalez notebook for script
# development notes.

import argparse

# Function to generate the numerical_header:NCBI_header dictionary
def num_ncbi_dict(ncbi_head, num_head):
    header_dict = {}
    ncbi_header_list = []
    try:
        with open("{0}".format(ncbi_head), "r") as ncbi_file:
            with open("{0}".format(num_head), "r") as num_file:
                for line in ncbi_file:
                    fixed_line = line.rstrip()
                    ncbi_header_list.append(fixed_line)
                for index, line in enumerate(num_file):
                    fixed_line = line.rstrip()
                    header_dict[fixed_line] = str(ncbi_header_list[index])
                return header_dict
    except IOError as error:
        print("An IO error has occured:{0}".format(error))

# Function to generate a tab-delimited file of numerical and ncbi headers. Serves as a lookup
# for interpreting phyfocus trees or to make a NCBI column for blastout tables.
def header_outfile(header_dict):
    try:
        with open("header_translation_table.tsv", "w") as out_file:
            for key, value in header_dict.items():
                out_file.write("{0}\t{1}\n".format(key, value))
    except IOError as error:
        print("An IO Error has occured:{0}".format(error))

# Store arguments for a FASTA protein file and user-defined batch size:
parser = argparse.ArgumentParser()
parser.add_argument("NCBI_headers", help="submits a .txt file of all NCBI-based headers from the fastas being used in the given analysis.")
parser.add_argument("numerical_headers", help="submits a .txt file of all numerical (aka 'fixed') headers from the fasta(s) being used in the given analysis. These headers are usually made by AWK commands in the pipeline.")
args = parser.parse_args()

# Make the dictionary which correlates numerical and NCBI headers, then write to tab-delimited out file
header_dictionary = num_ncbi_dict(args.NCBI_headers, args.numerical_headers)
header_outfile(header_dictionary)