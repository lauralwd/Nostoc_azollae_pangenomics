#!/bin/bash

# collect nanopore assembly stats
NANOPORE=( Azfil_lab Azpinnata Azsp_bordeaux )
SELECTION=( Nazollae mitochondrium chloroplast )

# header for the table:
echo -e "\
strain\t\
genome\t\
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
    do  strain="$s"
        seqyield=$(grep 'Total read length' data/nanopore_assembly/$g/$s/flye.log | rev | cut -d ' ' -f 1 | rev )

        echo -e "\
        $strain\t\
        $genome\t\
        $seqyield\t\
        $readN50bp\t\
        $ass_contigs\t\
        $ass_length\t\
        $ass_N50\t\
        $ass_cov\
        "
    done
done | tr -d ' '


exit 0
