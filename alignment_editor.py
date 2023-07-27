#! /usr/bin/env python3
import os
import datetime
import argparse
import Bio.AlignIO
import subprocess

# Define a function to produce a MAFFT MSA file and return the file's name as a string
def msa(fasta, seqs=None):
    # for edited fasta data, use final_dict from window() to avoid writing gap-causing seqs to file
    try:
        with open("{0}".format(fasta), "r") as original_fasta:
            with open("{0}_edited".format(fasta), "w") as edited_fasta:
                for line in original_fasta:
                    # Check if fasta headers are in dict of problem seqs
                    if line[0] == ">" and line.rstrip("\n") in seqs.keys():
                        good_seq = False
                        continue
                    elif line[0] == ">":
                        good_seq = True
                    # Write to edited fasta only good seqs not found in dict keys
                    if good_seq:
                        edited_fasta.write("{0}".format(line))
        # Align the edited fatsa file in MAFFT
        subprocess.run("mafft {0}_edited > {0}_edited_ali.fa".format(fasta), shell=True)
        print("MAFFT FFT-NS-2 (fast; progressive method) realignment of edited seqs complete.\n")
        return "{0}_edited_ali.fa".format(fasta)
    except IOError as error:
        print(error)

# Define a sliding window function that identifies problematic alignment gaps
def window(ali, arg):
    print("**********\nBeginning assessment of gap-causing seqs in MAFFT alignment...")
    start_slice = 0
    final_dict = {}
    # Given (L - k + 1 = # k-mers), iterate # kmer/window times
    window_num = (ali.get_alignment_length() - arg.window + 1)
    print("Examining {0} {1}bp windows...".format(window_num, arg.window))
    for idx in range(window_num):
        if idx == (window_num // 2):
            print("50% of windows examined...")
        # Sliding window slices alignment columns according to start_slice and window size values
        slide_window = ali[:, start_slice:start_slice + arg.window]
        # For every line in a window, obtain the seq
        gap_counter = 0
        temp_dict = {}
        for line in range(len(slide_window)):
            seq = str(slide_window[line].seq)
            # Determine if % gaps in seq's window length indicate problematic gaps
            if (seq.count("-") / len(seq)) >= 0.95:
                gap_counter += 1
            else:
                # Determine if % nucs present in seq's window can cause problematic gaps
                if (len(seq) - seq.count("-")) / len(seq) >= 0.5:
                    temp_dict.setdefault(">{0}".format(str(slide_window[line].id)),("{0}:{1}".format(start_slice,start_slice + arg.window), seq))
        # Determine if there are enough gap seqs to warrent removing gap-causing seqs
        if (gap_counter / len(slide_window)) >= arg.gappercent:
            for key, value in temp_dict.items():
                if key in final_dict:
                    final_dict[key].append(value)
                final_dict.setdefault(key, [value])
        # End of window operations; shifts window 1 bp to the right
        start_slice += 1
    print("Window examinations complete.\n")
    return final_dict


# Write window() parameters and results to output file
def results_file(arg, original_ali, edited_ali, seqs):
    print("**********\nWriting summary results of alignment_editor to ./alignment_editor.out")
    try:
        with open("alignment_editor.out", "w") as out_file:
            # Note program run date and location
            out_file.write("Alignment_editor.py completed {0}\n".format(datetime.datetime.now()))
            out_file.write("Analysis located at {0}\n".format(os.getcwd()))
            # Writes tab-delimited table of run parameters
            out_file.write("\nrun_parameters:\n")
            out_file.write("input_fasta_file: {0}\n".format(arg.file))
            out_file.write("window_size: {0}\n".format(arg.window))
            out_file.write("gap_percent: {0}\n".format(arg.gappercent))
            # Writes tab-delimited table of summary results
            out_file.write("\nsummary of results:\n")
            out_file.write("length\tnum_seqs\tfile_name\n")
            out_file.write("{0}bp\t{1}\t{2}\n".format(original_ali.get_alignment_length(), len(original_ali), arg.file))
            out_file.write("{0}bp\t{1}\t{2}_edited_ali.fa\n".format(edited_ali.get_alignment_length(), len(edited_ali), arg.file))
            out_file.write("alignment length reduced by {0}bp\n".format(original_ali.get_alignment_length() - edited_ali.get_alignment_length()))
            # Lists each unique gap-causing seq
            out_file.write("\n{0} seq(s) create problematic gaps:\n".format(len(seqs)))
            unique_seqs = set()
            for name, value in seqs.items():
                unique_seqs.add(name)
            for seq in unique_seqs:
                out_file.write("{0}\n".format(seq))
            # Writes each gap-causing seq's problem positions
            out_file.write("\nProblematic gaps are at these alignment positions:\n")
            for name, value in seqs.items():
                for seq_tuple in value:
                    out_file.write("{0}:{1}\n".format(name, seq_tuple[0]))
    except IOError:
        print("An IO error occured.")
    print("Writing to output file complete.\n")

# Store arguments for a FASTA-alignment file, window() parameters, and batch_files():
parser = argparse.ArgumentParser(description="alignment_editor helps remove any massive gap-causing sequences from the 'hmmr_significant_hits_ali.fa' file produced in Step 5 of Phyfocus. Once editing is complete, constructing a ML phylogeny (Step 6) completes PhyFocus. NOTE: program must be run within the 'final_tree_dataset/filtering_output' directory.")
parser.add_argument("file", help="submits a FASTA sequence file")
parser.add_argument("--window", "-w", type=int, default=50, help="Specify the minimum gap length that is disruptive to the alignment. Default is 50bp, but an effective benchmark appears to be 4 percent of the alignment's total length. Manual assessment of the alignment is strongly recommended.")
parser.add_argument("--gappercent", "-g", type=float, default=0.9, help="Specify the minimum percent of seqs that must contain -w size gaps to identify problematic seqs. Default is 0.9")
args = parser.parse_args()


# Generate MSA from original fasta file and convert to Bio.AlignIO object
ali_fasta = args.file
ali_object = (Bio.AlignIO.read("{0}".format(ali_fasta), "fasta"))

# Assess gaps via sliding window, return dictionary of problem seq headers, gap positions, and seqs
gap_seqs = window(ali_object, args)

# Generate edited MSA from original fasta file and convert to Bio.AlignIO object
edited_ali_fasta = msa(args.file, seqs=gap_seqs)
edited_ali_object = (Bio.AlignIO.read("{0}".format(edited_ali_fasta), "fasta"))

# Write results to a results summary file
results_file(args, ali_object, edited_ali_object, gap_seqs)
