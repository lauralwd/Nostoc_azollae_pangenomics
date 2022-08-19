#!/bin/bash

# it's easier (required perhaps) to install anvio not in snakemake but like this


# Coloured text
BLD=$(tput bold)
GRN=$(tput setaf 2)
RED=$(tput setaf 1)
BLU=$(tput setaf 4)
NML=$(tput sgr0)

# functions to use in the script
function checkprog {
if   [ ! $(command -v "$1") ]
then echo "$BLD""$RED""Command $1 is not found $NML"
     echo 'Make sure the propper conda environments are installed'
     exit 1
fi
}

# work with conda in a script
checkprog conda

condadir=/home/laura/miniconda3                                    # (mini)conda(3) directory
if   [ ! -f "$condadir"/etc/profile.d/conda.sh ]
then echo -e "$BLD""$RED""quiting for we need the conda environments in the `envs` directory to proceed $NML"
     exit
else source "$condadir"/etc/profile.d/conda.sh
fi

# get to work
if   [ $(conda env list | grep anvio7.1 | wc -l) -eq 0 ]
then echo 'update conda'
     conda update -n base -c defaults conda

     echo 'create new conda env'
     conda create -y --name anvio7.1 python=3.6
elif [ $(conda env list | grep anvio7.1 | wc -l) -gt 1 ]
then echo 'more than one anvio environment for you? work this out manually'
elif [ $(conda env list | grep anvio7.1 | wc -l) -eq 1 ]
then conda activate anvio7.1
     conda install -y -c bioconda "sqlite>=3.31.1"
     conda install -y -c bioconda prodigal
     conda install -y -c bioconda mcl
     conda install -y -c bioconda muscle=3.8.1551
     conda install -y -c bioconda hmmer
     conda install -y -c bioconda diamond
     conda install -y -c bioconda blast
     conda install -y -c bioconda megahit
     conda install -y -c bioconda spades
     conda install -y -c bioconda bowtie2 tbb=2019.8
     conda install -y -c bioconda bwa
     conda install -y -c bioconda samtools=1.9
     conda install -y -c bioconda centrifuge
     conda install -y -c bioconda trimal
     conda install -y -c bioconda iqtree
     conda install -y -c bioconda trnascan-se
     conda install -y -c bioconda r-base
     conda install -y -c bioconda r-stringi
     conda install -y -c bioconda r-tidyverse
     conda install -y -c bioconda r-magrittr
     conda install -y -c bioconda r-optparse
     conda install -y -c bioconda bioconductor-qvalue
     conda install -y -c bioconda fasttree
     conda install -y -c bioconda vmatch
     # this last one may cause some issues. if it doesn't install,
     # don't worry, you will still be fine:
     conda install -y -c bioconda fastani

     curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz \
        --output anvio-7.1.tar.gz

    pip install anvio-7.1.tar.gz
fi

exit
