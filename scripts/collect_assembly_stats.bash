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
