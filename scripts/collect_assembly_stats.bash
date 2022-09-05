#!/bin/bash

# requires Bandage conda environment

# collect nanopore assembly stats
NANOPORE=( Azfil_lab Azpinnata Azsp_bordeaux )
ILLUMINA=( Azcar1 Azcar2 Azmexicana Azmicrophylla Aznilotica Azrubra Azfil_lab )
SELECTION=( Nazollae mitochondrium chloroplast )

# header for the table:
echo -e "\
strain\t\
genome\t\
seqtype\t\
seqyield_bp\t\
readN50bp\t\
assembled_contigs\t\
assembled_length_bp\t\
assembled_N50\t\
coverage\
"

# Collect nanopore data
for s in ${NANOPORE[@]}
do  for g in ${SELECTION[@]}
    do  flye="data/nanopore_assembly/$g/$s/flye.log"
        band="data/nanopore_assembly/$g/$s.txt"
        strain="$s"
        genome="$g"
        seqyield=$(grep 'Total read length' $flye | rev | cut -d ' ' -f 1 | rev )
        readN50bp=$(grep 'Reads N50' $flye | rev | cut -d ' ' -f 3 | rev )
        ass_contigs=$(grep 'Node count' $band | rev | cut -d ' ' -f 1 | rev )
        ass_length=$(grep 'Total length (bp):' $band | rev | cut -d ' ' -f 1 | rev )
        ass_N50=$(grep 'N50 (bp):' $band | rev | cut -d ' ' -f 1 | rev )
        ass_cov=$(grep 'Median depth:' $band | rev | cut -d ' ' -f 1 | rev )

        echo -e "\
        $strain\t\
        $genome\t\
        nanopore\t\
        $seqyield\t\
        $readN50bp\t\
        $ass_contigs\t\
        $ass_length\t\
        $ass_N50\t\
        $ass_cov\
        "
    done
done | tr -d ' '

# Now for illumina samples:

SELECTION=( mitochondrium chloroplast )

for s in ${ILLUMINA[@]}
do  for g in ${SELECTION[@]}
    do  novo=data/illumina_assembly/"$g"_novoplasty/"$g"_"$s"/log_"$g"_"$s".txt
        fsta=data/illumina_assembly/"$g"_novoplasty/"$g"_"$s"_assembly.fasta
        strain="$s"
        genome="$g"
        seqyield=$(grep 'Assembled reads' $novo | rev | cut -d ' ' -f 1 | rev )
        readN50bp='100'
        ass_contigs=$(grep 'Total contigs' $novo | rev | cut -d ' ' -f 1 | rev )
        ass_length=$(Bandage info $fsta | grep 'Total length (bp):' | rev | cut -d ' ' -f 1 | rev )
        ass_N50=$(Bandage info $fsta | grep 'N50 (bp):' | rev | cut -d ' ' -f 1 | rev )
        ass_cov=$(grep 'Average organelle coverage' $novo | rev | cut -d ' ' -f 1 | rev )

        echo -e "\
        $strain\t\
        $genome\t\
        illumina\t\
        $seqyield\t\
        $readN50bp\t\
        $ass_contigs\t\
        $ass_length\t\
        $ass_N50\t\
        $ass_cov\
        "
    done
done | tr -d ' '

# MAG stats

MAGS=( Azfil_lab Azfil_wild Azmex Azmic Aznil Azrub Azcar1 Azcar2 )

g=Nazollae
for s in ${MAGS[@]}
do  fsta=data/MAG_anvi_dbs/"$s"_contigs.fasta
    strain="$s"
    genome="$g"
    seqyield='NA'
    readN50bp=100
    ass_contigs=$(Bandage info $fsta | grep 'Node count' | rev | cut -d ' ' -f 1 | rev )
    ass_length=$(Bandage info $fsta | grep 'Total length (bp):' | rev | cut -d ' ' -f 1 | rev )
    ass_N50=$(Bandage info $fsta | grep 'N50 (bp):' | rev | cut -d ' ' -f 1 | rev )
    ass_cov='NA'

    echo -e "\
    $strain\t\
    $genome\t\
    MAG\t\
    $seqyield\t\
    $readN50bp\t\
    $ass_contigs\t\
    $ass_length\t\
    $ass_N50\t\
    $ass_cov\
    "
done | tr -d ' '


exit 0
