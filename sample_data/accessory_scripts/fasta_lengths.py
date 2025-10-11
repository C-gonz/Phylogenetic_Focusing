#! /usr/bin/env python3

import argparse
import Bio.SeqIO.FastaIO

def seq_check(fasta, length):
    try:
        with open(fasta) as original_fasta:
            with open("{0}_cut.fa".format(fasta), "w") as cut_fasta:
                with open("cut_seqs.txt", "w") as cut_txt:
                    for title, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(original_fasta):
                        seq_len = len(seq)
                        print("length of {0} is {1}".format(title, seq_len))
                        if seq_len >= cutoff:
                            cut_fasta.write(">{0}\n{1}\n".format(title, seq))
                        else:
                            cut_txt.write(">{0}\n".format(title, seq))
    except IOError as error:
        print(error)

# Store arguments for user submitted fasta file and seq length cuttoff
parser = argparse.ArgumentParser(description="This program enables the user to remove all nucleotide or protein sequences that are less than the given cutoff length.")
parser.add_argument("--file", "-f", required=True, help="submits a FASTA sequence file")
parser.add_argument("--cutoff", "-c", type=int, default=0, help="Defines the minimal sequence length for keeping a sequence")
args=parser.parse_args()

# variables
in_fasta = args.file
cutoff= args.cutoff

# call function
seq_check(in_fasta, cutoff)
