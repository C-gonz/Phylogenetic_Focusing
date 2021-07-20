#! /bin/bash

for hitfile in ./all_top_hits/*_hits
do
    perl -p -e 's/\>//g' $hitfile > temp
    filenamefull=${hitfile##*/}
    filename="${filenamefull%_ref_blastout_hits*}"
    echo "pull FASTA seqs from  $filename ..."
    ./selectSeqs.pl -f temp ./fastas/$filename >> $filename
    rm temp
done

mkdir hits_fasta_redo

mv *.fas ./hits_fasta_redo
