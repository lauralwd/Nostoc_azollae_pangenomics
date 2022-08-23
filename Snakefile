NANOPORE=['Azfil_lab','Azpinnata','Azsp_bordeaux']
SELECTION=['Nazollae','mitochondrium','chloroplast']
ILLUMINA_HOSTS=['Azcar1','Azcar2','Azmexicana','Azmicrophylla','Aznilotica','Azrubra']
MAG_HOSTS=['Azfil_lab','Azfil_wild','Azmex_IRRI_486','Azmic_IRRI_456','Aznil_IRRI_479','Azrub_IRRI_479','Azspnov_IRRI1_472','Azspnov_IRRI2_489']

############################### stage 1: collect MAGs and genomes of Nostoc azollae from anvio ###############################
rule gather_anvi_MAGs:
  output:
    directory("data/MAG_anvi_dbs/")
  shell:
    "bash ./scripts/collect_mag_dbs.sh"

rule anvi_MAGs_fasta:
  input:
    "data/MAG_anvi_dbs/"
  output:
    temp(expand("data/MAGs/{mag_host}_contigs.db",mag_host=MAG_HOSTS))
  shell:
    """
    cp --reflink=always {input}/*_contigs.db  data/MAGs/
    """

rule anvi_contigdb_runhmms:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.hmms")
  log:
    stderr="logs/anvi_annotate/{db}_hmms.stderr",
    stdout="logs/anvi_annotate/{db}_hmms.stdout"
  threads: 6
  shell:
    """
    anvi-run-hmms -c {input}   \
                  -T {threads} \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_contigdb_kegg:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.kegg")
  log:
    stderr="logs/anvi_annotate/{db}_kegg.stderr",
    stdout="logs/anvi_annotate/{db}_kegg.stdout"
  threads: 6
  shell:
    """
    anvi-run-kegg-kofams -c {input}   \
                         -T {threads} \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_contigdb_cogs:
  input:
    "{db}_contigs.db"
  output:
    touch("{db}_contigs.db.cogs")
  log:
    stderr="logs/anvi_annotate/{db}_cogs.stderr",
    stdout="logs/anvi_annotate/{db}_cogs.stdout"
  threads: 6
  shell:
    """
    anvi-run-ncbi-cogs -c {input}   \
                       -T {threads} \
                       --sensitive  \
    > {log.stdout} 2> {log.stderr}
    """

rule anvi_MAGs_to_fasta:
  input:
    "data/MAGs/{mag_host}_contigs.db"
  output:
    "data/MAGs/{mag_host}_contigs.fasta",
  log:
    stderr="logs/MAG_anvi_dbs/export_contigs_{mag_host}.stderr",
    stdout="logs/MAG_anvi_dbs/export_contigs_{mag_host}.stdout"
  shell:
    """
    anvi-export-contigs -c {input}  \
                        -o {output} \
    > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_storage_internal:
  input:
    "data/MAG_anvi_dbs/",
    txt="data/Anvio_internal_genomes.txt",
    annotation=expand("data/MAGs/{mag_host}_contigs.db.ext",mag_host=MAG_HOSTS),{ext},ext=['hmms','kegg','cogs']),
    dbs=expand("data/MAGs/{mag_host}_contigs.db",mag_host=MAG_HOSTS)
  output:
    "data/anvio_genomes_storage/MAGs_only_GENOMES.db"
  log:
    stdout="logs/anvi_create_pangenome_storage.stdout",
    stderr="logs/anvi_create_pangenome_storage.stderr"
  shell:
    """
    anvi-gen-genomes-storage \
      -i {input.txt} \
      -o {output}    \
      > {log.stdout} 2> {log.stderr}
    """

rule scaffold_Nazollae_MAGs_RAGTAG:
  input:
    scaffolds="data/MAGs/{mag_host}_contigs.fasta",
    Nazollae0708="references/Nazollae_0708.fasta"
  params:
    pre=lambda w: expand("data/MAGs/{mag_host}_scaffolded/",mag_host=w.mag_host)
  output:
    "data/MAGs/{mag_host}_scaffolded/ragtag.scaffold.fasta"
  log:
    stderr="logs/MAG_anvi_dbs/RAGTAG{mag_host}.stderr",
    stdout="logs/MAG_anvi_dbs/RAGTAG{mag_host}.stdout"
  threads: 12
  conda:
    "envs/ragtag.yaml"
  shell:
    """
    ragtag.py scaffold {input.Nazollae0708} \
                       {input.scaffolds}   \
                       -t {threads}        \
                       -o {params.pre}     \
    > {log.stdout} 2> {log.stderr}
    """

rule all_mags_scaffolded:
  input:
    expand("data/MAGs/{mag_host}_scaffolded/ragtag.scaffold.fasta",
           mag_host=MAG_HOSTS)

############################### stage 2: filter and select nanopore reads, then assemble these into Nostoc azollae, mitochondria, or chloroplast genomes. ###############################
rule get_reference_fastas:
  output:
    "references/Nazollae_0708.fasta",
    "references/Azfil_v1.fasta",
    "references/Azfil_cp1.4.fasta"
  log:
    stdout="logs/get_reference_fastas.stdout",
    stderr="logs/get_reference_fastas.stderr"
  shell:
    """
    bash scripts/get_references.sh
     > {log.stdout} 2> {log.stderr}
    """

rule combine_reference_fastas:
  input:
    "references/Nazollae_0708.fasta",
    "references_Azfil_v1.fasta",
    "references/Azfil_cp1.4.fasta",
    "references/azfi_mito_laura-v1.fasta"
  output:
    "references/Azfil_combo_genome_v1.fasta"
  shell:
    """
    cat {input} > {output}
    """

rule make_minimap_index:
  input:
    fasta="{fasta}.fasta"
  output:
    index="{fasta}.minimap2-index"
  log:
    stdout="logs/minimap_index_{fasta}.stdout",
    stderr="logs/minimap_index_{fasta}.stderr"
  conda:
    "envs/minimap2.yaml"
  shell:
    """
    minimap2 {input.fasta} -d {output.index} \
     > {log.stdout} 2> {log.stderr}
    """

rule map_nanopore_data:
  input:
    reads="data/nanopore_sequencing/{nanopore_host}.fastq.gz",
    index="references/Azfil_combo_genome_v1.minimap2-index"
  output:
    bam=temp("data/nanopore_mapped/{nanopore_host}_mapped.bam")
  log:
    stderr="logs/map_nanopore_data_{nanopore_host}.stderr"
  conda:
    "envs/minimap2.yaml"
  threads: 12
  shell:
    """
    minimap2 -x map-ont \
             -t {threads}   \
             -a         \
             {input.index}\
             {input.reads}\
             2> {log.stderr} \
             | samtools sort \
             2>> {log.stderr} \
             | samtools view -b \
             2>> {log.stderr} \
             > {output.bam}
    """

rule map_all_nanopore_data:
  input:
    expand("data/nanopore_mapped/{nanopore_host}_mapped.bam",nanopore_host=NANOPORE)

rule sort_index_bam:
  input:
    bam="{bam}_mapped.bam"
  output:
    bam="{bam}_mapped_sorted.bam",
    bai="{bam}_mapped_sorted.bam.bai"
  log:
    stderr="logs/sort_index_bam_{bam}.stderr"
  threads: 12
  shell:
    """
    samtools sort {input.bam}  \
               -o {output.bam} \
               -l 8            \
               -m 4G           \
               -@ {threads}    \
    2> {log.stderr}

    samtools index {output.bam}
    2>> {log.stderr}
    """

rule get_Nazollae_nanopore_reads:
  input:
    bam=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam",
    bai=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam.bai"
  output:
    fastq="data/nanopore_filtered/{nanopore_host}_Nazollae_reads.fastq.gz"
  log:
    stderr="logs/get_Nazollae_nanopore_reads_{nanopore_host}.stderr"
  shell:
    """
    samtools view -h {input.bam}            \
                  'ENA|CP002059|CP002059.1' 'ENA|CP002060|CP002060.1' 'ENA|CP002061|CP002061.1' \
    2> {log.stderr}                      \
    | samtools fastq -0 {output.fastq} - \
    2>> {log.stderr}
    """

rule get_chloroplast_nanopore_reads:
  input:
    bam=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam",
    bai=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam.bai"
  output:
    fastq="data/nanopore_filtered/{nanopore_host}_chloroplast_reads.fastq.gz"
  log:
    stderr="logs/get_chloroplast_nanopore_reads_{nanopore_host}.stderr"
  shell:
    """
    samtools view -h {input.bam}         \
                  'Azolla_cp_v1_4'       \
    2> {log.stderr}                      \
    | samtools fastq -0 {output.fastq} - \
    2>> {log.stderr}
    """

rule get_mitochondrium_nanopore_reads:
  input:
    bam=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam",
    bai=  "data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam.bai"
  output:
    fastq="data/nanopore_filtered/{nanopore_host}_mitochondrium_reads.fastq.gz"
  log:
    stderr="logs/get_mitochondrium_nanopore_reads_{nanopore_host}.stderr"
  shell:
    """
    samtools view -h {input.bam}                           \
                  'a_filiculoides_mitochondrium_contig_10' \
                  'a_filiculoides_mitochondrium_contig_11' \
                  'a_filiculoides_mitochondrium_contig_12' \
                  'a_filiculoides_mitochondrium_contig_16' \
                  'a_filiculoides_mitochondrium_contig_22' \
                  'a_filiculoides_mitochondrium_contig_09' \
    2> {log.stderr}                                        \
    | samtools fastq -0 {output.fastq} -                   \
    2>> {log.stderr}
    """

rule assembly_with_flye:
  input:
    fastq="data/nanopore_filtered/{nanopore_host}_{selection}_reads.fastq.gz"
  output:
    dir=directory("data/nanopore_assembly/{selection}/{nanopore_host}")
  conda:
    "envs/flye.yaml"
  threads: 12
  log:
    stdout="logs/flye_assemble/{nanopore_host}_{selection}.stdout",
    stderr="logs/flye_assemble/{nanopore_host}_{selection}.stderr"
  shell:
    """
    if   [ {wildcards.selection} == 'Nazollae' ]
    then size=6M
    elif [ {wildcards.selection} == 'chloroplast' ]
    then size=150K
    elif [ {wildcards.selection} == 'mitochondrium' ]
    then size=100K
    fi

    flye    --nano-raw {input.fastq}  \
            --genome-size "$size"     \
            --threads {threads}       \
            --out-dir {output.dir}    \
            --asm-coverage 50         \
    > {log.stderr} 2> {log.stdout}
    """

rule map_all_nanopore_assemblies:
  input:
    expand("data/nanopore_assembly/{selection}/{nanopore_host}",nanopore_host=NANOPORE,selection=SELECTION)

############################### stage 3 assemble chloroplast and mito genomes from Illumina data ###############################
rule bwa_index:
  input:
    "{fasta}.fasta"
  output:
    "{fasta}.fasta.amb",
    "{fasta}.fasta.ann",
    "{fasta}.fasta.bwt",
    "{fasta}.fasta.pac",
    "{fasta}.fasta.sa"
  log:
    stderr="logs/illumina_genomes/bwa_index_{fasta}.stderr",
    stdout="logs/illumina_genomes/bwa_index_{fasta}.stdout"
  shell:
    """
    bwa index {input} \
    > {log.stderr} 2> {log.stdout}
    """

rule map_illumina_reads_to_combined_reference:
  input:
    R1="data/illumina_reads/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_reads/{illumina_host}_R2.fastq.gz",
    fasta="references/Azfil_combo_genome_v1.fasta",
    index="references/Azfil_combo_genome_v1.fasta.bwt"
  output:
    bam=temp("data/illumina_mapped/{illumina_host}_mapped.bam")
  log:
    stderr="logs/illumina_genomes/map_illumina_reads_to_combined_reference_{illumina_host}.stderr"
  threads: 12
  shell:
    """
    bwa mem                      \
            -t {threads}         \
            {input.fasta}        \
            {input.R1}           \
            {input.R2}           \
    2> {log.stderr}              \
    | samtools view -F 4         \
               -@ {threads}      \
               -b                \
               -o {output.bam}   \
    2> {log.stderr}
    """

rule sort_index_bam_name:
  input:
    bam="{bam}_mapped.bam"
  output:
    bam=temp("{bam}_mapped_sorted-name.bam"),
    bai=temp("{bam}_mapped_sorted-name.bam.bai")
  log:
    stderr="logs/samtools/sort_index_bam_name_{bam}.stderr"
  threads: 12
  shell:
    """
    samtools sort {input.bam}  \
               -o {output.bam} \
               -l 8            \
               -m 4G           \
               -n
               -@ {threads}    \
    2> {log.stderr}

    samtools index {output.bam}
    2>> {log.stderr}
    """

rule get_chloroplast_illumina_reads:
  input:
    bam="data/illumina_mapped/{illumina_host}_mapped_sorted-name.bam",
    bai="data/illumina_mapped/{illumina_host}_mapped_sorted-name.bam.bai"
  output:
    R1="data/illumina_filtered/chloroplast/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_filtered/chloroplast/{illumina_host}_R2.fastq.gz"
  log:
    stderr="logs/illumina_genomes/get_chloroplast_illumina_reads_{illumina_host}.stderr"
  shell:
    """
    samtools view -h {input.bam}         \
                  'Azolla_cp_v1_4'       \
    2> {log.stderr}                      \
    | samtools fastq -1 {output.R1}      \
                     -2 {output.R2}      \
                     -                   \
    2>> {log.stderr}
    """

rule get_mitochondrium_illumina_reads:
  input:
    bam="data/illumina_mapped/{illumina_host}_mapped_sorted-name.bam",
    bai="data/illumina_mapped/{illumina_host}_mapped_sorted-name.bam.bai"
  output:
    R1="data/illumina_filtered/mitochondrium/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_filtered/mitochondrium/{illumina_host}_R2.fastq.gz"
  log:
    stderr="logs/illumina_genomes/get_mitochondrium_illumina_reads_{illumina_host}.stderr"
  shell:
    """
    samtools view -h {input.bam}                           \
                  'a_filiculoides_mitochondrium_contig_10' \
                  'a_filiculoides_mitochondrium_contig_11' \
                  'a_filiculoides_mitochondrium_contig_12' \
                  'a_filiculoides_mitochondrium_contig_16' \
                  'a_filiculoides_mitochondrium_contig_22' \
                  'a_filiculoides_mitochondrium_contig_09' \
    2> {log.stderr}                                        \
    | samtools fastq -1 {output.R1}                        \
                     -2 {output.R2}                        \
                     -                                     \
    2>> {log.stderr}
    """

rule assemble_organelle_genome_SPAdes:
  input:
    R1="data/illumina_filtered/{selection}/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_filtered/{selection}/{illumina_host}_R2.fastq.gz"
  output:
    scaffolds="data/illumina_assembly/{selection}/{illumina_host}/scaffolds.fasta",
    graph="data/illumina_assembly/{selection}/{illumina_host}/assembly_graph_with_scaffolds.gfa"
  params:
    pre=lambda w: expand("data/illumina_assembly/{selection}/{illumina_host}/",
                         selection=w.selection,
                         illumina_host=w.illumina_host)
  threads: 12
  conda:
    "envs/spades.yaml"
  log:
    stderr="logs/illumina_genomes/SPAdes_{selection}/{illumina_host}.stderr",
    stdout="logs/illumina_genomes/SPAdes_{selection}/{illumina_host}.stdout"
  shell:
    """
    spades.py -1 {input.R1}         \
              -2 {input.R2}         \
              -t {threads}          \
              -o {params.pre}       \
    > {log.stdout} 2> {log.stderr}"
    """

rule scaffold_chloroplast_genome_RAGTAG:
  input:
    scaffolds="data/illumina_assembly/chloroplast/{illumina_host}/scaffolds.fasta",
    chloroplast="references/Azfil_cp1.4.fasta"
  params:
    pre=lambda w: expand("data/illumina_assembly/chloroplast/{illumina_host}_scaffolded",illumina_host=w.illumina_host)
  output:
    "data/illumina_assembly/chloroplast/{illumina_host}_scaffolded/ragtag.scaffold.fasta"
  log:
    stderr="logs/illumina_genomes/SPAdes_chloroplast/{illumina_host}.stderr",
    stdout="logs/illumina_genomes/SPAdes_chloroplast/{illumina_host}.stdout"
  threads: 12
  conda:
    "envs/ragtag.yaml"
  shell:
    """
    ragtag.py scaffold {input.chloroplast} \
                       {input.scaffolds}   \
                       -t {threads}        \
                       -o {params.pre}     \
    > {log.stdout} 2> {log.stderr}
    """

rule scaffold_mitochondrium_genome_RAGTAG:
  input:
    scaffolds="data/illumina_assembly/mitochondrium/{illumina_host}/scaffolds.fasta",
    mitochondrium="references/azfi_mito_laura-v1.fasta"
  params:
    pre=lambda w: expand("data/illumina_assembly/mitochondrium/{illumina_host}_scaffolded",illumina_host=w.illumina_host)
  output:
    "data/illumina_assembly/mitochondrium/{illumina_host}_scaffolded/ragtag.scaffold.fasta"
  log:
    stderr="logs/illumina_genomes/SPAdes_mitochondrium{illumina_host}.stderr",
    stdout="logs/illumina_genomes/SPAdes_mitochondrium/{illumina_host}.stdout"
  threads: 12
  conda:
    "envs/ragtag.yaml"
  shell:
    """
    ragtag.py scaffold {input.mitochondrium} \
                       {input.scaffolds}   \
                       -t {threads}        \
                       -o {params.pre}     \
    > {log.stdout} 2> {log.stderr}
    """

rule all_illumina_assembly:
  input:
    expand("data/illumina_assembly/{selection}/{illumina_host}_scaffolded/ragtag.scaffold.fasta",
           selection=['chloroplast','mitochondrium'],
           illumina_host=ILLUMINA_HOSTS)

############################### stage 4 create pangenomes ###############################
rule illumina_assembly_to_contigdb:
  input:
    "data/illumina_assembly/{selection}/{illumina_host}_scaffolded/ragtag.scaffold.fasta"
  output:
     "data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db"
  log:
    stdout="logs/anvi_create_contigdb/{illumina_host}_{selection}.stdout",
    stderr="logs/anvi_create_contigdb/{illumina_host}_{selection}.stderr"
  shell:
    """
    anvi-gen-contigs-database   \
      -f {input}/assembly.fasta \
      -o {output}               \
      > {log.stdout} 2> {log.stderr}
    """

rule nanopore_assembly_to_contigdb:
  input:
    "data/nanopore_assembly/{selection}/{nanopore_host}"
  output:
     "data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db"
  log:
    stdout="logs/anvi_create_contigdb/{nanopore_host}_{selection}.stdout",
    stderr="logs/anvi_create_contigdb/{nanopore_host}_{selection}.stderr"
  shell:
    """
    anvi-gen-contigs-database   \
      -f {input}/assembly.fasta \
      -o {output}               \
      > {log.stdout} 2> {log.stderr}
    """

rule all_contigdbs:
  input:
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",      nanopore_host=NANOPORE,      selection=SELECTION),
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,      selection=SELECTION),                     ext=['hmms','kegg','cogs']),
    expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",      illumina_host=ILLUMINA_HOSTS,selection=['chloroplast','mitochondrium']),
    expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection=['chloroplast','mitochondrium'],ext=['hmms','kegg','cogs'])

rule create_pangenome_storage_all_Nazollaes:
  input:
    "data/MAG_anvi_dbs/",
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='Nazollae'),
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='Nazollae', ext=['hmms','kegg','cogs']),
    internal="data/Anvio_internal_genomes.txt",
    external="data/Anvio_external_genomes.txt",
    ref="data/nanopore_contig_dbs/Azfil_0708_Nazollae_contigs.db"
  output:
    "data/anvio_genomes_storage/Nazollae_GENOMES.db"
  log:
    stdout="logs/anvi_create_pangenome_storage_all.stdout",
    stderr="logs/anvi_create_pangenome_storage_all.stderr"
  shell:
    """
    anvi-gen-genomes-storage \
      -i {input.internal}    \
      -e {input.external}    \
      -o {output}            \
      > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_storage_all_chloroplast:
  input:
    illuminachloroplasts=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",illumina_host=ILLUMINA_HOSTS,selection='chloroplast'),
    illuminachloroplasts_ann=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection='chloroplast',ext=['hmms','kegg','cogs']),
    nanoporechloroplasts=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='chloroplast'),
    nanoporechloroplasts_ann=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='chloroplast',ext=['hmms','kegg','cogs']),
    external="data/Anvio_external_chloroplast.txt",
    ref="data/nanopore_contig_dbs/Azfil_chloroplast_contigs.db"
  output:
    "data/anvio_genomes_storage/chloroplast_GENOMES.db"
  log:
    stdout="logs/anvi_create_pangenome_storage_chloroplast.stdout",
    stderr="logs/anvi_create_pangenome_storage_chloroplast.stderr"
  shell:
    """
    anvi-gen-genomes-storage \
      -e {input.external}    \
      -o {output}            \
      > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_storage_all_mitochondrium:
  input:
    illuminachloroplasts=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",illumina_host=ILLUMINA_HOSTS,selection='mitochondrium'),
    illuminachloroplasts_ann=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection='mitochondrium',ext=['hmms','kegg','cogs']),
    nanoporechloroplasts=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='mitochondrium'),
    nanoporechloroplasts_ann=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='mitochondrium',ext=['hmms','kegg','cogs']),
    external="data/Anvio_external_mitochondrium.txt",
    ref="data/nanopore_contig_dbs/Azfil_mitochondrium_contigs.db"
  output:
    "data/anvio_genomes_storage/mitochondrium_GENOMES.db"
  log:
    stdout="logs/anvi_create_pangenome_storage_mitochondrium.stdout",
    stderr="logs/anvi_create_pangenome_storage_mitochondrium.stderr"
  shell:
    """
    anvi-gen-genomes-storage \
      -e {input.external}    \
      -o {output}            \
      > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_analysis:
  input:
    "data/anvio_genomes_storage/{selection}_GENOMES.db"
  output:
    directory("data/anvio_pangenomes/{selection}")
  log:
    stdout="logs/anvi_create_pangenome_{selection}.stdout",
    stderr="logs/anvi_create_pangenome_{selection}.stderr"
  threads: 12
  shell:
    """
    anvi-pan-genome -g {input}                           \
                    --project-name {wildcards.selection} \
                    --output-dir   {output}              \
                    --num-threads  {threads}             \
                    --minbit 0.5                         \
                    --mcl-inflation 8                    \
                    --use-ncbi-blast                     \
     > {log.stdout} 2> {log.stderr}
     """

rule create_pangenome_ANI_Nazollae:
  input:
    pangenome="data/anvio_pangenomes/Nazollae",
    internal="data/Anvio_internal_genomes.txt",
    external="data/Anvio_external_genomes.txt",
    MAGS="data/MAG_anvi_dbs/",
    extdbs=expand("data/nanopore_contig_dbs/{nanopore_host}_Nazollae_contigs.db",nanopore_host=NANOPORE),
    ref="data/nanopore_contig_dbs/Azfil_0708_Nazollae_contigs.db"
  output:
    directory("data/anvio_pangenomes/Nazollae_ANI")
  log:
    stdout="logs/anvi_create_pangenome_ANI_Nazollae.stdout",
    stderr="logs/anvi_create_pangenome_ANI_Nazollae.stderr"
  threads: 12
  shell:
    """
    anvi-compute-genome-similarity --external-genomes {input.external} \
                                   --internal-genomes {input.internal} \
                                   --program pyANI                     \
                                   --output-dir {output}               \
                                   --num-threads {threads}             \
                                   --pan-db {input.pangenome}/{wildcards.selection}-PAN.db \
     > {log.stdout} 2> {log.stderr}
     """

rule create_pangenome_ANI_organele:
  input:
    pangenome="data/anvio_pangenomes/{selection}",
    external="data/Anvio_external_genomes_{selection}.txt",
    extdbsnanopore=expand("data/nanopore_contig_dbs/{nanopore_host}_{{selection}}_contigs.db",nanopore_host=NANOPORE),
    extdbsillumina=expand("data/illumina_contig_dbs/{illumina_host}_{{selection}}_contigs.db",illumina_host=ILLUMINA_HOSTS)
  output:
    directory("data/anvio_pangenomes/{selection}_ANI")
  log:
    stdout="logs/anvi_create_pangenome_ANI_{selection}.stdout",
    stderr="logs/anvi_create_pangenome_ANI_{selection}.stderr"
  threads: 12
  shell:
    """
    anvi-compute-genome-similarity --external-genomes {input.external} \
                                   --program pyANI                     \
                                   --output-dir {output}               \
                                   --num-threads {threads}             \
                                   --pan-db {input.pangenome}/{wildcards.selection}-PAN.db \
     > {log.stdout} 2> {log.stderr}
     """

rule all_Azolla_associated_pangenomes:
  input:
    expand("data/anvio_pangenomes/{selection}_ANI",selection=SELECTION)

rule extract_phylogenomic_fasta:
  input:
    pangenome="data/anvio_pangenomes/{selection}",
    genomestorage="data/anvio_genomes_storage/{selection}_GENOMES.db"
  output:
    fasta="data/anvio_pangenomes/{selection}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{selection}_phylogenomic_core.partitions"
  log:
    stdout="logs/anvi_pangenome_fasta_{selection}.stdout",
    stderr="logs/anvi_pangenome_fasta_{selection}.stderr"
  shell:
    """
    anvi-get-sequences-for-gene-clusters -p {input.pangenome}/{wildcards.selection}-PAN.db \
                                         -g {input.genomestorage}    \
                                         -C default                  \
                                         -b phylogenomic_core        \
                                         --concatenate-gene-clusters \
                                         --partition-file {output.partition} \
                                         -o {output.fasta}           \
    > {log.stdout} 2> {log.stderr}
    """

# making a phylogenomic tree requires the user to manually create the phylogenomic_core bin inside the Anvi'o interactive interface

rule phylogenomic_tree:
  input:
    fasta="data/anvio_pangenomes/{selection}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{selection}_phylogenomic_core.partitions"
  output:
    tree="data/anvio_pangenomes/{selection}_phylogenomics/{selection}.treefile"
  params:
    pre=lambda w: expand ("data/anvio_pangenomes/{selection}_phylogenomics/{selection}",selection=w.selection)
  threads: 12
  log:
    stdout="logs/IQtree/anvi_phylogenomic_{selection}.stdout",
    stderr="logs/IQtree/anvi_phylogenomic_{selection}.stderr"
  shell:
    """
    iqtree -s {input.fasta}     \
           -p {input.partition} \
           -m MFP+MERGE         \
           -b 100               \
           -nt AUTO             \
           -ntmax {threads}     \
           -pre {params.pre}    \
    > {log.stdout} 2> {log.stderr}
    """

rule all_Azolla_associated_phylogenomics:
  input:
    expand("data/anvio_pangenomes/{selection}_phylogenomics/{selection}.treefile",selection=SELECTION)
