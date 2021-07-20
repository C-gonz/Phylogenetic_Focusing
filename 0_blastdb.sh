#!/bin/bash

#module purge
#module load linuxbrew/colsa

# Run this script in the directory containg the gzip files you wish to make blast db's for.

# unzip all gzip fasta files
#gzip -d ./fastas/*

# modify every fasta header in the unzipped files to be genus & count: ">GENUS_####". 
# Following 2 commands must be replicated for every file: 
# awk ‘/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./DIRECTORY/infile.fasta > GENUS_SPECIES.fa
# perl -p -i -e 's/>/>GENUS_/g' GENUS_SPECIES.fa

awk ‘/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCA_000209535.1_ASM20953v1_protein.faa > Oikopleura_dioica.fas
perl -p -i -e 's/>/>Oikopleura_/g' Oikopleura_dioica.fas





#for item in ./*
#do
#    echo "formating for BLAST $item ..."
#    gzip -dc $item | makeblastdb -parse_seqids -out ${item%.*}_db -dbtype prot -title $item
#done

# move blastdbs to new dir
#mkdir ../blastdb
#cp ./*_db* ../blastdb

# clean up blast databases not in the blastdb directory
#rm ./*_db*
