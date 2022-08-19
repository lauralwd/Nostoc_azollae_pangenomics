if   [ ! -f references/Nazollae_0708.fasta ]
then wget 'https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000196515.1?download=true&gzip=false'\
         -O references/Nazollae_0708.fasta
fi

if   [ ! -f references/Azfil_v1.2.fasta ]
then  wget https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/Azolla_filiculoides.genome_v1.2.fasta \
           -O references_Azfil_v1.fasta
fi

if   [ ! -f references/Azfil_cp1.4.fasta ]
then  wget https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/chloroplast_genome/Azolla_filiculoides.cp_genome_v1_4.fasta \
           -O references/Azfil_cp1.4.fasta
fi
