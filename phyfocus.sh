#!/bin/bash

module purge
module load linuxbrew/colsa

# SCRIPT PURPOSE:
# Fixes fasta headers; makes BLAST Databases (db's) from fixed fasta files; BLASTp queries against each db; saves fixed_fastas,
# BLAST outputs, & db's to new directories.

# REQUIRED FILE STRUCTURE/CONTENTS:
# Parent_directory/-->(shell, *queries.fas, selectSeqs.pl, tree_editor.R, fastas/-->(species_fasta_files))
    # Shell must execute script while in parent directory.
    # Fasta query sequences file must end in "queries.fas" and be in parent directory. File should contain both bait and anchor sequences.
    # must have access to perl scripts:
        # selectSeqs.pl (https://raw.githubusercontent.com/plachetzki/cnidarian_opsin/master/selectSeqs.pl)
        # tree_editor.R (template must be edited for your purposes) (https://github.com/plachetzki/cnidarian_opsin/blob/master/5_tree_editor.R)



# Unzip all fasta files
echo unzipping all gzip files in $(pwd)...
gzip -d ./fastas/*.gz

# Following 2 commands must be replicated for each fasta file to make final fasta headers = ">GENUS_#####"
# awk ‘/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./Directory/infile.fasta > Genus_species_FIX.fa
# perl -p -i -e 's/>/>Genus_/g' Genus_species_FIX.fa
    # Uses awk to search for all ">" headers, replace the line with count ">#####"
    # Uses perl to search for all ">", replace them with ">Genus_"

echo FASTA Header Fix START
echo Fixing fasta headers ...

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Crocodylus_porosus.CroPor_comp1.pep.all.fa> Crocodylus_porosus_FIX.fa
perl -p -i -e 's/>/>Crocodylus_/g' Crocodylus_porosus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Danio_rerio.GRCz11.pep.all.fa> Danio_rerio_FIX.fa
perl -p -i -e 's/>/>Danio_/g' Danio_rerio_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Dromaius_novaehollandiae.droNov1.pep.all.fa> Dromaius_novaehollandiae_FIX.fa
perl -p -i -e 's/>/>Dromaius_/g' Dromaius_novaehollandiae_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Eptatretus_burgeri.Eburgeri_3.2.pep.all.fa> Eptatretus_burgeri_FIX.fa
perl -p -i -e 's/>/>Eptatretus_/g' Eptatretus_burgeri_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Equus_caballus.EquCab3.0.pep.all.fa> Equus_caballus_FIX.fa
perl -p -i -e 's/>/>Equus_/g' Equus_caballus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.pep.all.fa> Fundulus_heteroclitus_FIX.fa
perl -p -i -e 's/>/>Fundulus_/g' Fundulus_heteroclitus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Gadus_morhua.gadMor1.pep.all.fa> Gadus_morhua_FIX.fa
perl -p -i -e 's/>/>Gadus_/g' Gadus_morhua_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Gallus_gallus.GRCg6a.pep.all.fa> Gallus_gallus_FIX.fa
perl -p -i -e 's/>/>Gallus_/g' Gallus_gallus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCA_000209535.1_ASM20953v1_protein.faa> Oikopleura_dioica_FIX.fa
perl -p -i -e 's/>/>Oikopleura_/g' Oikopleura_dioica_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCA_002284835.2_RCv2.1_protein.faa> Lithobates_catesbeianus_FIX.fa
perl -p -i -e 's/>/>Lithobates_/g' Lithobates_catesbeianus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCA_003427355.1_Storazame_v1.0_protein.faa> Scyliorhinus_torazame_FIX.fa
perl -p -i -e 's/>/>Scyliorhinus_/g' Scyliorhinus_torazame_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_000003605.2_Skow_1.1_protein.faa> Saccoglossus_kowalevskii_FIX.fa
perl -p -i -e 's/>/>Saccoglossus_/g' Saccoglossus_kowalevskii_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_000090745.1_AnoCar2.0_protein.faa> Anolis_carolinensis_FIX.fa
perl -p -i -e 's/>/>Anolis_/g' Anolis_carolinensis_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_000224145.3_KH_protein.faa> Ciona_intestinalis_FIX.fa
perl -p -i -e 's/>/>Ciona_/g' Ciona_intestinalis_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa> Chrysemys_picta_FIX.fa
perl -p -i -e 's/>/>Chrysemys_/g' Chrysemys_picta_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_001625305.1_Haploidv18h27_protein.faa> Branchiostoma_belcheri_FIX.fa
perl -p -i -e 's/>/>Branchiostoma_/g' Branchiostoma_belcheri_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_001642345.1_ASM164234v2_protein.faa> Rhincodon_typus_FIX.fa
perl -p -i -e 's/>/>Rhincodon_/g' Rhincodon_typus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_001949145.1_OKI-Apl_1.0_protein.faa> Acanthaster_planci_FIX.fa
perl -p -i -e 's/>/>Acanthaster_/g' Acanthaster_planci_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_010909765.1_sAmbRad1.pri_protein.faa> Amblyraja_radiata_FIX.fa
perl -p -i -e 's/>/>Amblyraja_/g' Amblyraja_radiata_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_011630105.1_ASM1163010v1_protein.faa> Anneissia_japonica_FIX.fa
perl -p -i -e 's/>/>Anneissia_/g' Anneissia_japonica_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/GCF_901001135.1_aRhiBiv1.1_protein.faa> Rhinatrema_bivittatum_FIX.fa
perl -p -i -e 's/>/>Rhinatrema_/g' Rhinatrema_bivittatum_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Gopherus_agassizii.ASM289641v1.pep.all.fa> Gopherus_agassizii_FIX.fa
perl -p -i -e 's/>/>Gopherus_/g' Gopherus_agassizii_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Homo_sapiens.GRCh38.pep.all.fa> Homo_sapiens_FIX.fa
perl -p -i -e 's/>/>Homo_/g' Homo_sapiens_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Latimeria_chalumnae.LatCha1.pep.all.fa> Latimeria_chalumnae_FIX.fa
perl -p -i -e 's/>/>Latimeria_/g' Latimeria_chalumnae_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Lepisosteus_oculatus.LepOcu1.pep.all.fa> Lepisosteus_oculatus_FIX.fa
perl -p -i -e 's/>/>Lepisosteus_/g' Lepisosteus_oculatus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Microcebus_murinus.Mmur_3.0.pep.all.fa> Microcebus_murinus_FIX.fa
perl -p -i -e 's/>/>Microcebus_/g' Microcebus_murinus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Monopterus_albus.M_albus_1.0.pep.all.fa> Monopterus_albus_FIX.fa
perl -p -i -e 's/>/>Monopterus_/g' Monopterus_albus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Myxine_glutinosa_31619_Compiled_peps.fa> Myxine_glutinosa_FIX.fa
perl -p -i -e 's/>/>Myxine_/g' Myxine_glutinosa_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Notechis_scutatus.TS10Xv2-PRI.pep.all.fa> Notechis_scutatus_FIX.fa
perl -p -i -e 's/>/>Notechis_/g' Notechis_scutatus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Ornithorhynchus_anatinus.OANA5.pep.all.fa> Ornithorhynchus_anatinus_FIX.fa
perl -p -i -e 's/>/>Ornithorhynchus_/g' Ornithorhynchus_anatinus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa> Oryctolagus_cuniculus_FIX.fa
perl -p -i -e 's/>/>Oryctolagus_/g' Oryctolagus_cuniculus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Petromyzon_marinus.Pmarinus_7.0.pep.all.fa> Petromyzon_marinus_FIX.fa
perl -p -i -e 's/>/>Petromyzon_/g' Petromyzon_marinus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Physeter_catodon.ASM283717v2.pep.all.fa> Physeter_catodon_FIX.fa
perl -p -i -e 's/>/>Physeter_/g' Physeter_catodon_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Pseudonaja_textilis.EBS10Xv2-PRI.pep.all.fa> Pseudonaja_textilis_FIX.fa
perl -p -i -e 's/>/>Pseudonaja_/g' Pseudonaja_textilis_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Rattus_norvegicus.Rnor_6.0.pep.all.fa> Rattus_norvegicus_FIX.fa
perl -p -i -e 's/>/>Rattus_/g' Rattus_norvegicus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Sarcophilus_harrisii.DEVIL7.0.pep.all.fa> Sarcophilus_harrisii_FIX.fa
perl -p -i -e 's/>/>Sarcophilus_/g' Sarcophilus_harrisii_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Strongylocentrotus_purpuratus.Spur_3.1.pep.all.fa> Strongylocentrotus_purpuratus_FIX.fa
perl -p -i -e 's/>/>Strongylocentrotus_/g' Strongylocentrotus_purpuratus_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Taeniopygia_guttata.bTaeGut1_v1.p.pep.all.fa> Taeniopygia_guttata_FIX.fa
perl -p -i -e 's/>/>Taeniopygia_/g' Taeniopygia_guttata_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Varanus_komodoensis.ASM479886v1.pep.all.fa> Varanus_komodoensis_FIX.fa
perl -p -i -e 's/>/>Varanus_/g' Varanus_komodoensis_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all.fa> Xenopus_tropicalis_FIX.fa
perl -p -i -e 's/>/>Xenopus_/g' Xenopus_tropicalis_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Pongo_abelii.PPYG2.pep.all.fa> Pongo_abelii_FIX.fa
perl -p -i -e 's/>/>Pongo_/g' Pongo_abelii_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Gorilla_gorilla.gorGor4.pep.all.fa> Gorilla_gorilla_FIX.fa
perl -p -i -e 's/>/>Gorilla_/g' Gorilla_gorilla_FIX.fa

awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ./fastas/Macaca_mulatta.Mmul_10.pep.all.fa> Macaca_mulatta_FIX.fa
perl -p -i -e 's/>/>Macaca_/g' Macaca_mulatta_FIX.fa



# make directory for fixed fastas, then move FIX files
echo moving modified files to fixed_fastas ...
mkdir ./fixed_fastas
mv *_FIX.fa ./fixed_fastas
echo Successfully moved $(ls ./fixed_fastas | wc -l) files to ./fixed_fastas
echo FASTA Header Fix Complete



# Step 0+1 (format blast databases from ./fixed_fastas/ & run BLASTp)
# Iterate over fixed_fastas directory; for "species_pep", make database and blastp query against db
echo Steps 0+1 START
for species_pep in ./fixed_fastas/*
do
    echo "formating BLAST db for $species_pep ..."
    makeblastdb -in $species_pep -parse_seqids -out ${species_pep%.*}_db -dbtype prot
    # blastp query file against species_pep database
    echo BLASTING query sequences against $species_pep db ...
    blastp -db ${species_pep%.*}_db -query ./*queries.fas -evalue 1e-05 -num_threads 24 -out ${species_pep%.*}_ref_blastout -outfmt 6 -max_target_seqs 1000
done

# move blast dbs & blast outputs to new directories
echo "moving db's to ../blastdb"
mkdir ./blastdb
mv ./fixed_fastas/*_db* ./blastdb

echo "moving blastout files to ./blastout_maxseqs"
mkdir ./blastout_maxseqs
mv ./fixed_fastas/*blastout ./blastout_maxseqs/
echo Steps 0+1 Complete



# Step 2 (retrieve accession numbers for each blastp output & move to new dir)
echo Step 2 START
for blastfile in ./blastout_maxseqs/*blastout
do
  echo "getting hit accessions for $blastfile ..."
  cut -f 2 $blastfile > ${blastfile%.}_hits
done

# move all hit accession files to new dir called hit1_accessions
mkdir ./hit1_accessions
mv ./blastout_maxseqs/*_hits ./hit1_accessions
echo Successfully moved $(ls ./hit1_accessions | wc -l) hit files to hit1_accessions
echo Step 2 compelete



# Step 3 (retrieve seqeunces from *_FIX.fa files that corespond to accession numbers for each taxa)
# use the perl selectSeqs.pl script which, given a list of acessions and the corresponding fasta, retrieves all sequence data for each accession provided
echo Step 3 START
for hitfile in ./hit1_accessions/*_hits
do
    # store hit accessions in a temp file
    perl -p -e 's/\>//g' $hitfile > temp
    filenamefull=${hitfile##*/}
    filename="${filenamefull%_ref_blastout_hits*}.fa"
    echo "pulling FASTA seqs for hit accessions in $filename ..."
    # run selectSeqs.pl perl script, using accesion #'s in temp and seqs from ./fixed_fasta files
    echo executing selectSeqs.pl perl script ...
    ./selectSeqs.pl -f temp ./fixed_fastas/$filename >> $filename
done
rm ./temp

# Move retrieved hit fasta seqs to the new dir hits_fasta/
mkdir hits_fasta
mv *.fa ./hits_fasta
echo Successfully moved $(ls ./hits_fasta | wc -l) hit sequence fasta files to hits_fasta
echo Step 3 Complete



# Step 4: Remove similar seqs, concatenate queries+anchors+hits, and then align for each taxa
echo Step 4 START
# remove special characters
echo Removing special characters from fastas ...
perl -pi -e 's/\*//g' ./hits_fasta/*.fa

# Remove similar seqs, concatenate queries+anchors+hits, then align for each taxa
for fasta in ./hits_fasta/*.fa
do
    # run CD-HIT to remove seqs with 98% or higher seq similarity
    echo Running CD-HIT on $fasta ...
    cd-hit -i $fasta -o ${fasta%.}_cdhit -c 0.98 -n 5
    # after concatenation, align via MAFFT
    echo Concatenating query+anchor and cdhit files ...
    cat ./*_queries.fas ${fasta}*cdhit > ${fasta%.}_bts_ancrs
    echo Aligning concatenated seqs with MAFFT ...
    mafft ${fasta%.}_bts_ancrs > ${fasta%.}_ali.fa
done

# Make new dir named align/ for all alignments and clean up temp cat files
mkdir align
mv ./hits_fasta/*ali.fa ./align
rm ./hits_fasta/*ancrs
echo Moved alignment to ./align
echo Step 4 Complete



# Step 5: Make ML phylogenies from MAFFT aligned files using IQtree
echo Step 5 START
cd ./align
for ali_fa in *ali.fa
do
    echo Running IQtree on $ali_fa ...
    # Makes maximum likelihood phylogenies with IQtree. For non-testing runs, use ModelFinder (MFP) for -m
    iqtree -s $ali_fa -m MFP+C60 -alrt 1000 -bb 1000 -nt 24
done

# RAxML final tree topology files moved to RAx_out
echo Moving all IQtree result output files to IQ_out/...
mkdir IQ_out
mv *.treefile IQ_out
# IQtree additional outputs moved to IQ_out_pars_logs
echo Moving all extra files to IQ_out_extra_files/...
mkdir IQ_out_extra_files
mv *.fa.* IQ_out_extra_files
echo Step 5 Complete.



# Step 6: Run phylogenetic focusing R script on IQtree phylogenies
echo Step 6 START
cd IQ_out
Rscript ../../tree_editor.R
echo Step 6 Complete



# Step 7: strip away unwanted space and prefix to tree_editor output files: (GENUS_SPECIES_FIX.fa_ali.fa_tips.txt)
echo Fixing IQtree file names...
for tips_file in ./*.txt
do
    mv "$tips_file" "${tips_file/' '}"
done
rename '_FIX.fa_ali.fa' '' ./align/IQ_out/*.txt



# Step 9: Retrieve fasta seqs for genes included in tree_editor.R (step #7) output, and move to tip_seqs/
cd ../../
mkdir tip_seqs
for file in ./align/IQ_out/*.txt
do
    echo Extracting fasta sequences for $file...
    fasta_name=$(echo $file | cut -d '_' -f 2,3 | cut -d '/' -f 2)
    perl ./selectSeqs.pl -f $file ./fixed_fastas/${fasta_name}_FIX.fa > ./tip_seqs/${fasta_name}_seqs.fa
done



# Step 10: concatenate all fastas together, making the "total tree" dataset, then align in mafft and make a RAxML tree.
mkdir final_focused_tree
echo Concatenating tip_seqs and bait+anchor seqs for final tree fasta...
cat ./tip_seqs/*.fa ./*queries.fas > final_focused_tree/concat_tip_seqs.fa
# Remove duplicate seqs derived from concatenating all the tip_seq fastas
echo Removing identical seqs from concat_tip_seqs.fa...
cd-hit -i final_focused_tree/concat_tip_seqs.fa -o final_focused_tree/concat_tip_seqs_cdhit.fa -c 1.0 -n 5
# removes any special characters and aligns total fasta in MAFFT, and removes sequence line feeds
perl -pi -e 's/\*//g' ./final_focused_tree/concat_tip_seqs_cdhit.fa
echo Aligning concat_tip_seqs_cdhit.fa via MAFFT...
mafft ./final_focused_tree/concat_tip_seqs_cdhit.fa > final_focused_tree/concat_tip_seqs_cdhit_ali.fa
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ./final_focused_tree/concat_tip_seqs_cdhit.fa > ./final_focused_tree/concat_tip_seqs_cdhit_linefix.fa
echo Step 10 Complete

# Prepare alignment for domain filtering accoring to max file size cutoffs
echo concat_tip_seqs_cdhit_linefix.fa needs to be submitted to Pfam database. Batch submissions can only contain 500 seqs.
echo Running alignment_editor.py --batch
./alignment_editor.py --batch 500
echo Batch files completed. Submit batch files to Pfam for domain filtering.

