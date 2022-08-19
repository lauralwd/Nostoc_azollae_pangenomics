#!/bin/bash
mkdir  data/MAG_anvi_dbs/ 2> /dev/null
# get contigdbs
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_hybrid_doublefiltered_anvio/Azfil_wild/Azfil_wild_contigs.db		             data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_hybrid_doublefiltered_anvio/Azfil_lab/Azfil_lab_contigs.db		               data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Azrub_IRRI_479/Azrub_IRRI_479_contigs.db		     data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Azmex_IRRI_486/Azmex_IRRI_486_contigs.db		     data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Azmic_IRRI_456/Azmic_IRRI_456_contigs.db		     data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Azspnov_IRRI1_472/Azspnov_IRRI1_472_contigs.db	 data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Azspnov_IRRI2_489/Azspnov_IRRI2_489_contigs.db	 data/MAG_anvi_dbs/
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_anvio/Aznil_IRRI_479/Aznil_IRRI_479_contigs.db		     data/MAG_anvi_dbs/

# get profile dbs and rename them
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_hybrid_doublefiltered_binningsignals_anvio/MERGED_Azfil_wild/PROFILE.db		\
  data/MAG_anvi_dbs/Azfil_wild_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_hybrid_doublefiltered_binningsignals_anvio/MERGED_Azfil_lab/PROFILE.db		\
  data/MAG_anvi_dbs/Azfil_lab_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Azrub_IRRI_479/PROFILE.db		\
  data/MAG_anvi_dbs/Azrub_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Azmex_IRRI_486/PROFILE.db		\
  data/MAG_anvi_dbs/Azmex_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Azmic_IRRI_456/PROFILE.db		\
  data/MAG_anvi_dbs/Azmic_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Azspnov_IRRI1_472/PROFILE.db		\
  data/MAG_anvi_dbs/Azcar1_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Azspnov_IRRI2_489/PROFILE.db		\
  data/MAG_anvi_dbs/Azcar2_PROFILE.db
cp --reflink=always /stor/azolla_metagenome/Azolla_genus_metagenome/data/assembly_singles_doublefiltered_binningsignals_anvio/MERGED_Aznil_IRRI_479/PROFILE.db		\
  data/MAG_anvi_dbs/Aznil_PROFILE.db

for db in ./data/MAG_anvi_dbs/*.db
do anvi-migrate --migrate-dbs-safely $db \
     > $db.migrate.stdout \
     2> $db.migrate.stderr
done
