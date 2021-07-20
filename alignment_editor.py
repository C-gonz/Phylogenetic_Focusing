#! /usr/bin/env python3
import argparse
import Bio.AlignIO

# will need a recursive function to write the files of phyfocus step 10:
def batch_files(num):
    seq_count = 0
    file_count = 1
    try:
        with open("./final_focused_tree/concat_tip_seqs_cdhit_linefix.fa", "r") as ali_fasta:
            with open("./final_focused_tree/batched_linefix.fa_{0}".format(file_count), "a") as batch_out:
                for line in ali_fasta:
                    if line[0] == ">":
                        seq_count += 1
                        batch_out.write(line)
                        continue
                    batch_out.write(line)
                    if file_count % num == 0:
                        file_count += 1
                        with open("./final_focused_tree/batched_linefix.fa_{0}".format(file_count), "a") as batch_out:
                            continue
    except IOError:
        print("IOError")                    
                    
# Script for converting space-delimited items into tab-delimited.
# For use after Pfam hit results are concatenated into a single file named "Hmmr_Pfam_results.txt."
def delim_convert(file_path):
    print("Phyfocus Step 12: preparing alignment for 500 seq batch submission to Pfam.")
    try:
        with open("{./Hmmr_Pfam_results.txt}", "r") as read_file:
            with open("./Hmmr_Pfam_table.txt", "a") as write_file:
                for line in read_file:
                    final_string = ""
                    line_list = line.split()
                    for item in line_list:
                        # print("The current item is: ", item) FOR DEBUG. LIST IS WORKING.
                        item_string = str(item)
                        item_string += "\t"
                        final_string += item_string
                        # print("Final_string is currently: ", final_string) FOR DBUG. FINAL_STRING WORKS.
                    final_string += "\n"
                    write_file.write(final_string)
    except IOError:
        print("IOError.")



# Define a sliding window function that identifies problematic alignment gaps
def window_counter():
    #window_dict = {}
    start_slice = 0
    #print("alignment =",ali) #                               DEBUG
    #print("wsize =",args.window) #                                     DEBUG
    final_results = {}
    # Given (L - k + 1 = # k-mers), iterate # kmer/window times
    for idx in range(ali.get_alignment_length() - args.window + 1):
        # Sliding window takes all sequences, an slices columns according to start_slice and aegument window size
        window = ali[:, start_slice:start_slice + args.window]
        #print("window is: ",window)                              #DEBUG
        # For every line in a window, identify the seq
        counter = 0
        indel_seqs = {}
        for line in range(len(window)):
            #print("line is: ",line)                          #DEBUG
            seq = str(window[line].seq)
            #print("seq is: ",seq)                          #DEBUG
            #print("seq len is: ",len(seq))                          #DEBUG
            #print("seq '-' count is: ", seq.count("-"))               #DEBUG
            # Determine if % gaps in seq's window indicate problematic gaps
            if (seq.count("-") / len(seq)) >= 0.95:
                counter += 1
            else:
                # Determine if % nucs present in seq's window are sufficient to cause problematic gaps
                if (len(seq) - seq.count("-")) / len(seq) >= 0.5:
                    indel_seqs.setdefault(str(window[line].id), [])
                    indel_seqs[str(window[line].id)].append(seq)
        #print(indel_seqs)                                       #DEBUG
        #print("counter is: ",counter)                          #DEBUG
        #print((counter / len(window)), args.gappercent)          #DEBUG
        # Determine if there are enough gap seqs to warrent removing gap-causing seqs
        if (counter / len(window)) >= args.gappercent:
            for key, val in indel_seqs.items():
                final_results[key] = val
        start_slice += 1
    return final_results

# Write window_counter() parameters and results to output file
def results_table(data):
    try:
        with open("./{0}_gapout.txt".format(args.file), "w") as out_file:
            out_file.write("Run Parameters: \nInfile: {0}\nOutfile: {1}\nNumber of seqs in alignment: {2}\nWindow size: {3}\nGap percent cutoff: {4}\n".format(args.file, "{0}_gapout.txt".format(args.file), ali_seq_num, args.window, args.gappercent))
            out_file.write("\nThe following {0} seqs contain these gap-inducing positions:\n".format(seq_total))
            for name in data.keys():
                out_file.write("{0}\n".format(name))
            out_file.write("\nFasta formatted gap-causing seqs:\n")
            for name, val in data.items():
                for gap_seq in val:
                    out_file.write(">{0}\n{1}\n".format(name, gap_seq))
    except Exception as error:
        print(error)


# Store a user argument for a FASTA-alignment file
parser = argparse.ArgumentParser()
parser.add_argument("--file", "-f", help="submits a FASTA alignment file")
parser.add_argument("--window", "-w", type=int, default=50, help="Specify bp window size for alignment gap assessment. Default is 50bp")
parser.add_argument("--gappercent", "-g", type=float, default=0.8, help="Specify the minimum percent of seqs that contain gaps for a given window to identify gap-inducing seqs. Default is 0.8")
parser.add_argument("--batch", type=int, help="Runs the fasta batching program. Requires max num of seqs to write for each fasta file")
args = parser.parse_args()

if args.file is not None:
    # Retrieve fa alignment from args.file, save as an Bio.AlignIO instance.
    ali = (Bio.AlignIO.read("{0}".format(args.file), "fasta"))
    ali_seq_num = len(ali[:, 0])
    # Assess gap-inducing seqs, store as results, & calculate seq number
    data = window_counter()
    seq_total = len(data)

if args.batch is not None:
    batch_files(args.batch)

