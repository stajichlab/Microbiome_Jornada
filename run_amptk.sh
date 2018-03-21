#!/usr/bin/bash
CPU=2
#INPUT DIRECTORY IS CALLED split_libraries

BASE=Sam1_34a
#METADATA=fasta-qual-mapping-files/050517NP515F-mapping2.txt
METADATA=050517NP515F-mapping2_NP.txt
if [ ! -f $BASE.demux.fq.gz ]; then
# FWD PRIMER IS GTGCCAGCMGCCGCGGTAA
# REV PRIMER IS 806RB
# run 
# amptk primers
# to see already stored primers

 amptk illumina -i split_libraries -o $BASE -f GTGCCAGCMGCCGCGGTAA -r 806RB --cpus $CPU --require_primer off --rescue_forward on --primer_mismatch 2 -l 300
fi

if [ ! -f $BASE.otu_table.txt ];  then
 amptk cluster -i $BASE.demux.fq.gz -o $BASE --uchime_ref 16S --usearch usearch9 --map_filtered -e 1
fi

if [ ! -f $BASE.filtered.otus.fa ]; then
 amptk filter -i $BASE.otu_table.txt -f $BASE.cluster.otus.fa -p 0.005
fi

if [ ! -f $BASE.otu_table.taxonomy.txt ]; then
 amptk taxonomy -f $BASE.filtered.otus.fa -i $BASE.final.txt -d 16S
fi

if [ ! -f $BASE.otu_table.taxonomy.biom ]; then
 perl -p -e 's/^(\S+)\s+(.+);(\S+).+/$1\t$3/; s/tax=//; s/,/;/g; s/:/__/g;' $BASE.taxonomy.txt  > $BASE.taxonomy.tab
 biom convert -i $BASE.otu_table.txt -o $BASE.otu_table.biom --table-type "OTU table" --to-hdf5

 biom add-metadata -i ${BASE}.otu_table.biom -m ${METADATA} --observation-metadata-fp ${BASE}.taxonomy.tab \
  -o ${BASE}.otu_table.taxonomy.biom --observation-header OTUID,taxonomy --sc-separated taxonomy
 biom convert -i $BASE.otu_table.taxonomy.biom -o $BASE.otu_table.phinch.biom --to-json
fi
