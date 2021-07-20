perl -pi -e 's/\*//g' ./hits_fasta_redo/*.fas

for fasta in ./hits_fasta_redo/*.fas 
do

cd-hit -i $fasta -o ${fasta%.}_cdhit -c 0.98 -n 5 

cat opsin_bait.fa ${fasta%.}_cdhit > ${fasta%.}_o_cdhit
cat receptor_anchor.fa ${fasta%.}_o_cdhit > ${fasta%.}_o_r_cdhit

mafft ${fasta%.}_o_r_cdhit > ${fasta%.}_ali

done 

mkdir align
mv ./hits_fasta/*ali ./align 
