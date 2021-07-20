for fasta in ./align/*ali
do

perl ./seqConverter.pl -d$fasta -if -ope
done


