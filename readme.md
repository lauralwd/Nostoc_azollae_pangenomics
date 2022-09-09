# Nostoc pangenomics
This repository contains code and some results for the pangenomics of _Nostoc azollae_, endo-symbiont of the _Azolla_ genus of ferns.

Here I aim to collect all _N. azollae_ genomes from the _Azolla_ metagenome project, and assemble new ones where possible.
With these genomes, I want to answer questions about the pangenomics, evolution, and co-evolution of this symbiont with its host.

# methods
I'm using Anvi'o 7.1 for lots of analyses, and snakemake to process data.
Rather than using an up-to-date snakemake, I chose to use the snakemake with anvio7.1.

# to reproduce

## note about data
 - nanopore data is new and will be in ENA
 - Illumina data is part of Li2018, but something seems wrong with the data upload. Contact me for the original data as I used it
 - Genome and chloroplast available from fernbase
 - Mitochondrium is novel and will be uploaded some appropriate place.

## Visualise analyses
Within this Git repository with an anvio 7.1 environment activated, run for the chloroplast:

`anvi-display-pan -g ./data/anvio_genomes_storage/chloroplast_GENOMES.db -p data/anvio_pangenomes/chloroplast/chloroplast_mcl6-PAN.db --title 'Azolla genus chloroplast pangenome'`

Similarly for the _Nostoc azollae_ pangenome:

`anvi-display-pan -g ./data/anvio_genomes_storage/Nazollae_GENOMES.db -p data/anvio_pangenomes/Nazollae/Nazollae_mcl7-PAN.db --title 'Nostoc azollae pangenome'`

And finally for the mitochodrium one:

`anvi-display-pan -g ./data/anvio_genomes_storage/mitochondrium_GENOMES.db -p data/anvio_pangenomes/mitochondrium/mitochondrium_mcl3-PAN.db --title 'Azolla genus mitochondrium pangenome'`

This will allow you to browse the pangenomes interactively.

## recreate Azolla associated pan- and phylo-genomics
To create all pangenome databases in anvio

`snakemake --use-conda all_azolla_associated_pangenomes`

and run this to prevent unnecessarily rerunning any analysis:

`snakemake --use-conda all_azolla_associated_pangenomes --touch`

You might need to play around with the clustering parameters to get meaningfull patterns in the pangenome.
To procreed, identify a set of core genes with high geometric homogeneity but low-ish sequence homogeneity, bin them as `phylogenetic_core` and next run the phylogenomics:

`snakemake --use-conda all_azolla_associated_phylogenomics`

Save your specific visualisation and binning under the name 'default' and export via

`snakemake --use-conda all_azolla_associated_pangenome_svgs`

## Whole genome allignments

## Genome degradation

## Nostocaceae pan- and phylo-genomics
