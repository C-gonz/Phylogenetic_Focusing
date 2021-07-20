#! /bin/bash

# Execute in main working directory containing blastp output directory

# For each blastp output, write the hit accession # into a file
for i in ./blastout_maxseqs1000/*blastout
do
  echo "getting hit accessions for $i ..."
  cut -f2 $i > ${i%.}_hits
done

# move ll hit accession files into their own directory
mkdir hit1_accessions
mv ./blastout/*hits ./hit1_accessions
