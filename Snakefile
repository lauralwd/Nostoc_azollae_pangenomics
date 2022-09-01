NANOPORE=['Azfil_lab','Azpinnata','Azsp_bordeaux']
SELECTION=['Nazollae','mitochondrium','chloroplast']
ILLUMINA_HOSTS=['Azcar1','Azcar2','Azmexicana','Azmicrophylla','Aznilotica','Azrubra','Azfil_lab']
MAG_HOSTS=['Azfil_lab','Azfil_wild','Azmex','Azmic','Aznil','Azrub','Azcar1','Azcar2']

############################### stage 1: collect MAGs and genomes of Nostoc azollae from anvio ###############################
rule gather_anvi_MAGs:
  output:
    expand("data/MAG_anvi_dbs/{mag_host}_contigs.db" ,mag_host=MAG_HOSTS),
    expand("data/MAG_anvi_dbs/{mag_host}_PROFILE.db" ,mag_host=MAG_HOSTS)
  shell:
    "bash ./scripts/collect_mag_dbs.sh"

# you may need to remove old hmm, kegg, cogg results, or choose to keep your current ones and ignore the next three rules
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
                  --just-do-it \
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
                         --just-do-it \
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
    "data/MAG_anvi_dbs/{mag_host}_contigs.db"
  output:
    "data/MAG_anvi_dbs/{mag_host}_contigs.fasta"
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
    dbs=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db"      ,mag_host=MAG_HOSTS                           ),
    ann=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db.{ext}",mag_host=MAG_HOSTS,ext=['hmms','kegg','cogs']),
    txt="scripts/Nazollae_internal_genomes.anvi-list"
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
    scaffolds="data/MAG_anvi_dbs/{mag_host}_contigs.fasta",
    Nazollae0708="references/Nazollae_0708.fasta"
  params:
    pre=lambda w: expand("data/MAG_anvi_dbs/{mag_host}_scaffolded/",mag_host=w.mag_host)
  output:
    "data/MAG_anvi_dbs/{mag_host}_scaffolded/ragtag.scaffold.fasta"
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
    expand("data/MAG_anvi_dbs/{mag_host}_scaffolded/ragtag.scaffold.fasta",
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
    bam=temp("{bam}_mapped_sorted.bam"    ),
    bai=temp("{bam}_mapped_sorted.bam.bai")
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
    bam="data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam",
    bai="data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam.bai"
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
    bam="data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam",
    bai="data/nanopore_mapped/{nanopore_host}_mapped_sorted.bam.bai"
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

rule snapshot_flye_assembly_graph_Bandage:
  input:
    "data/nanopore_assembly/{selection}/{nanopore_host}"
  output:
    svg="data/nanopore_assembly/{selection}/{nanopore_host}.svg",
    txt="data/nanopore_assembly/{selection}/{nanopore_host}.txt"
  conda:
    "envs/bandage.yaml"
  shell:
    """
    Bandage image {input}/assembly_graph.gfa  \
                  {output.svg}

    Bandage info  {input}/assembly_graph.gfa    \
                  > {output.txt}
    """
rule all_nanopore_assemblies:
  input:
    expand("data/nanopore_assembly/{selection}/{nanopore_host}",    nanopore_host=NANOPORE,selection=SELECTION),
    expand("data/nanopore_assembly/{selection}/{nanopore_host}.svg",nanopore_host=NANOPORE,selection=SELECTION)
    #expand("data/nanopore_assembly/{selection}/{nanopore_host}.txt",nanopore_host=NANOPORE,selection=SELECTION)

############################### stage 3 assemble chloroplast and mito genomes from Illumina data ###############################

rule assemble_chloroplast_NOVOPlasty:
  input:
    R1="data/illumina_reads/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_reads/{illumina_host}_R2.fastq.gz",
    chloroplast="references/Azfil_cp1.4.fasta",
    config_base="scripts/novoplasty_chloroplast_config_base"
  output:
    sampleconfig="data/illumina_assembly/chloroplast_novoplasty/chloroplast_{illumina_host}_config.txt",
    dir=directory("data/illumina_assembly/chloroplast_novoplasty/chloroplast_{illumina_host}")
  threads: 2
  resources:
    mem_mb=20000
  conda:
    "envs/novoplasty.yaml"
  log:
    stderr="logs/illumina_genomes/novoplasty_chloroplast/{illumina_host}.stderr",
    stdout="logs/illumina_genomes/novoplasty_chloroplast/{illumina_host}.stdout"
  shell:
    """
    echo 'Project:'                                                      > {output.sampleconfig}
    echo '-----------------------'                                       >> {output.sampleconfig}
    echo 'Project name          = chloroplast_{wildcards.illumina_host}' >> {output.sampleconfig}

    cat {input.config_base}                                              >> {output.sampleconfig}

    echo 'Forward reads         =  {input.R1}'                           >> {output.sampleconfig}
    echo 'Reverse reads         =  {input.R2}'                           >> {output.sampleconfig}
    echo 'Output path           =  {output.dir}/{wildcards.illumina_host}_' >> {output.sampleconfig}

    if   [ ! -d {output.dir} ]
    then mkdir  {output.dir}
    fi

    NOVOPlasty4.3.1.pl -c {output.sampleconfig}   \
    > {log.stdout} 2> {log.stderr}
    """

rule assemble_mitochondrium_NOVOPlasty:
  input:
    R1="data/illumina_reads/{illumina_host}_R1.fastq.gz",
    R2="data/illumina_reads/{illumina_host}_R2.fastq.gz",
    mitochondrium="references/azfi_mito_laura-v1.fasta",
    config_base="scripts/novoplasty_mitochondrium_config_base"
  output:
    sampleconfig="data/illumina_assembly/mitochondrium_novoplasty/mitochondrium_{illumina_host}_config.txt",
    dir=directory("data/illumina_assembly/mitochondrium_novoplasty/mitochondrium_{illumina_host}")
  threads: 2
  resources:
    mem_mb=20000
  conda:
    "envs/novoplasty.yaml"
  log:
    stderr="logs/illumina_genomes/novoplasty_mitochondrium/{illumina_host}.stderr",
    stdout="logs/illumina_genomes/novoplasty_mitochondrium/{illumina_host}.stdout"
  shell:
    """
    echo 'Project:'                                                      > {output.sampleconfig}
    echo '-----------------------'                                       >> {output.sampleconfig}
    echo 'Project name          = mitochondrium_{wildcards.illumina_host}' >> {output.sampleconfig}

    cat {input.config_base}                                              >> {output.sampleconfig}

    echo 'Forward reads         =  {input.R1}'                           >> {output.sampleconfig}
    echo 'Reverse reads         =  {input.R2}'                           >> {output.sampleconfig}
    echo 'Output path           =  {output.dir}/{wildcards.illumina_host}_' >> {output.sampleconfig}

    if   [ ! -d {output.dir} ]
    then mkdir  {output.dir}
    fi

    NOVOPlasty4.3.1.pl -c {output.sampleconfig}   \
    > {log.stdout} 2>> {log.stderr}
    """

rule choose_novoplasty_assembly:
  input:
    "data/illumina_assembly/chloroplast_novoplasty/chloroplast_{illumina_host}/"
  output:
    "data/illumina_assembly/chloroplast_novoplasty/chloroplast_{illumina_host}_assembly.fasta"
  shell:
    """
    if   [ -f {input}/Circularized_assembly_1_chloroplast_{wildcards.illumina_host}.fasta ]
    then anvi-script-reformat-fasta {input}/Circularized_assembly_1_chloroplast_{wildcards.illumina_host}.fasta \
                                    --simplify-names     \
                                    --seq-type NT        \
                                    -o {output}
    elif [ -f {input}/Contigs_1_chloroplast_{wildcards.illumina_host}.fasta ]
    then anvi-script-reformat-fasta {input}/Contigs_1_chloroplast_{wildcards.illumina_host}.fasta \
                                    --simplify-names     \
                                    --seq-type NT        \
                                    -o {output}
    else exit 1
    fi
    """

rule all_illumina_assembly_novoplasty:
  input:
    expand("data/illumina_assembly/{selection}_novoplasty/{selection}_{illumina_host}_config.txt",
           selection=['chloroplast','mitochondrium'],
           illumina_host=ILLUMINA_HOSTS)

############################### stage 4 create pangenomes ###############################
rule illumina_assembly_to_contigdb:
  input:
    "data/illumina_assembly/{selection}_novoplasty/{selection}_{illumina_host}_assembly.fasta"
  output:
     "data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db"
  log:
    stdout="logs/anvi_create_contigdb/{illumina_host}_{selection}.stdout",
    stderr="logs/anvi_create_contigdb/{illumina_host}_{selection}.stderr"
  shell:
    """
    anvi-gen-contigs-database   \
      -f {input}                \
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

rule reference_to_contigdb:
  input:
    "references/{fasta}.fasta"
  output:
    "data/external_contig_dbs/{fasta}_contigs.db"
  log:
    stdout="logs/anvi_create_contigdb/{fasta}.stdout",
    stderr="logs/anvi_create_contigdb/{fasta}.stderr"
  shell:
    """
    anvi-gen-contigs-database   \
      -f {input}                \
      -o {output}               \
      > {log.stdout} 2> {log.stderr}
    """

rule all_contigdbs:
  input:
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",      nanopore_host=NANOPORE,      selection=SELECTION                                                 ),
    expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,      selection=SELECTION                      ,ext=['hmms','kegg','cogs']),
    expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",      illumina_host=ILLUMINA_HOSTS,selection=['chloroplast','mitochondrium']                           ),
    expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection=['chloroplast','mitochondrium'],ext=['hmms','kegg','cogs']),
    expand("data/MAG_anvi_dbs/{mag_host}_contigs.db"                            ,mag_host=MAG_HOSTS                                                                               ),
    expand("data/MAG_anvi_dbs/{mag_host}_contigs.db.{ext}"                      ,mag_host=MAG_HOSTS                                                    ,ext=['hmms','kegg','cogs']),
    expand("data/external_contig_dbs/Nazollae_0708_contigs.db.{ext}"     ,ext=['hmms','kegg','cogs']),
    expand("data/external_contig_dbs/Azfil_cp1.4_contigs.db.{ext}"       ,ext=['hmms','kegg','cogs']),
    expand("data/external_contig_dbs/azfi_mito_laura-v1_contigs.db.{ext}",ext=['hmms','kegg','cogs'])

rule create_pangenome_storage_all_Nazollaes:
  input:
    magdbs=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db"      ,mag_host=MAG_HOSTS                           ),
    magprf=expand("data/MAG_anvi_dbs/{mag_host}_PROFILE.db"      ,mag_host=MAG_HOSTS                           ),
    magann=expand("data/MAG_anvi_dbs/{mag_host}_contigs.db.{ext}",mag_host=MAG_HOSTS,ext=['hmms','kegg','cogs']),    nandbs=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='Nazollae'),
    nanann=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='Nazollae', ext=['hmms','kegg','cogs']),
    internal="scripts/Nazollae_internal_genomes.anvi-list",
    external="scripts/Nazollae_external_genomes.anvi-list",
    refdb =       "data/external_contig_dbs/Nazollae_0708_contigs.db",
    refann=expand("data/external_contig_dbs/Nazollae_0708_contigs.db.{ext}",ext=['hmms','kegg','cogs'])
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

rule create_list_pangenome_storage_all_chloroplast:
  input:
    illuminachloroplasts=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",illumina_host=ILLUMINA_HOSTS,selection='chloroplast'),
    illuminachloroplasts_ann=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection='chloroplast',ext=['hmms','kegg','cogs']),
    nanoporechloroplasts=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='chloroplast'),
    nanoporechloroplasts_ann=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='chloroplast',ext=['hmms','kegg','cogs']),
    refann=expand("data/external_contig_dbs/Azfil_cp1.4_contigs.db.{ext}"       ,ext=['hmms','kegg','cogs']),
    refdb= "data/external_contig_dbs/Azfil_cp1.4_contigs.db"
  output:
    external="scripts/chloroplast_external_genomes.anvi-list",
  shell:
    """
    echo -e "name\tcontigs_db_path" > {output}
    echo -e "Azfil_cp1_4\t../{input.refdb}" >> {output}
    for db in data/*_contig_dbs/Az*chloroplast_contigs.db
    do  name=$(echo $db | sed 's/data\///g' | sed 's/contig_dbs\///g' | sed 's/_chloroplast_contigs.db//g')
        echo -e "$name\t./../$db"
    done >> {output}
    """

rule create_list_pangenome_storage_all_mitochondrium:
  input:
    illuminamito=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db",illumina_host=ILLUMINA_HOSTS,selection='mitochondrium'),
    illuminamitoann=expand("data/illumina_contig_dbs/{illumina_host}_{selection}_contigs.db.{ext}",illumina_host=ILLUMINA_HOSTS,selection='mitochondrium',ext=['hmms','kegg','cogs']),
    nanoporemito=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db",nanopore_host=NANOPORE,selection='mitochondrium'),
    nanoporemitoann=expand("data/nanopore_contig_dbs/{nanopore_host}_{selection}_contigs.db.{ext}",nanopore_host=NANOPORE,selection='mitochondrium',ext=['hmms','kegg','cogs']),
    refann=expand("data/external_contig_dbs/azfi_mito_laura-v1_contigs.db.{ext}",ext=['hmms','kegg','cogs']),
    refdb =       "data/external_contig_dbs/azfi_mito_laura-v1_contigs.db"
  output:
    external="scripts/mitochondrium_external_genomes.anvi-list",
  shell:
    """
    echo -e "name\tcontigs_db_path" > {output}
    echo -e "Azfi_mito_v1\t../{input.refdb}" >> {output}
    for db in data/*_contig_dbs/Az*mitochondrium_contigs.db
    do  name=$(echo $db | sed 's/data\///g' | sed 's/contig_dbs\///g' | sed 's/_mitochondrium_contigs.db//g')
        echo -e "$name\t./../$db"
    done >> {output}
    """

rule create_pangenome_storage_plastid:
  input:
    external="scripts/{selection}_external_genomes.anvi-list"
  output:
    "data/anvio_genomes_storage/{selection}_GENOMES.db"
  log:
    stdout="logs/anvi_create_pangenome_storage_{selection}.stdout",
    stderr="logs/anvi_create_pangenome_storage_{selection}.stderr"
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
    "data/anvio_pangenomes/{selection}/{selection}_mcl{mcl}-PAN.db"
  log:
    stdout="logs/anvi_create_pangenome_{selection}_mcl{mcl}.stdout",
    stderr="logs/anvi_create_pangenome_{selection}_mcl{mcl}.stderr"
  threads: 12
  params:
    mcl= "7",
    dir=lambda w: expand ("data/anvio_pangenomes/{selection}/",selection=w.selection)
  shell:
    """
    anvi-pan-genome -g {input}                           \
                    --project-name {wildcards.selection}_mcl{wildcards.mcl} \
                    --output-dir   {params.dir}          \
                    --num-threads  {threads}             \
                    --minbit 0.5                         \
                    --mcl-inflation {params.mcl}         \
                    --use-ncbi-blast                     \
    > {log.stdout} 2> {log.stderr}
    """

# in the step above, I experimented manually with the mcl value and found 7 to give the best distinction between phylogenetically meaningfull groups

# Manually imported pan groups like so:
# anvi-import-misc-data scripts/Nazollae_pan_groups.txt -p ./data/anvio_pangenomes/Nazollae/Nazollae_mcl7-PAN.db --target-data-table layers

rule create_pangenome_ANI_Nazollae:
  input:
    pangenome="data/anvio_pangenomes/Nazollae/Nazollae_mcl{mcl}-PAN.db",
    internal="scripts/Nazollae_internal_genomes.anvi-list",
    external="scripts/Nazollae_external_genomes.anvi-list"
  output:
    directory("data/anvio_pangenomes/Nazollae_ANI_mcl{mcl}")
  log:
    stdout="logs/anvi_create_pangenome_ANI_Nazollae_mcl{mcl}.stdout",
    stderr="logs/anvi_create_pangenome_ANI_Nazollae_mcl{mcl}.stderr"
  threads: 12
  shell:
    """
    anvi-compute-genome-similarity --external-genomes {input.external} \
                                   --internal-genomes {input.internal} \
                                   --program pyANI                     \
                                   --output-dir {output}               \
                                   --num-threads {threads}             \
                                   --pan-db {input.pangenome}          \
    > {log.stdout} 2> {log.stderr}
    """

rule create_pangenome_ANI_organele:
  input:
    pangenome="data/anvio_pangenomes/{selection}/{selection}_mcl{mcl}-PAN.db",
    external="scripts/{selection}_external_genomes.anvi-list"
  output:
    directory("data/anvio_pangenomes/{selection}_ANI_mcl{mcl}")
  log:
    stdout="logs/anvi_create_pangenome_ANI_{selection}_mcl{mcl}.stdout",
    stderr="logs/anvi_create_pangenome_ANI_{selection}_mcl{mcl}.stderr"
  threads: 12
  shell:
    """
    anvi-compute-genome-similarity --external-genomes {input.external} \
                                   --program pyANI                     \
                                   --output-dir {output}               \
                                   --num-threads {threads}             \
                                   --pan-db {input.pangenome}          \
     > {log.stdout} 2> {log.stderr}
    """

rule all_azolla_associated_pangenomes:
  input:
    expand("data/anvio_pangenomes/{selection}_ANI_mcl{mcl}",selection=SELECTION,mcl=7)

rule extract_phylogenomic_fasta:
  input:
    pangenome="data/anvio_pangenomes/{selection}/{selection}_mcl{mcl}-PAN.db",
    genomestorage="data/anvio_genomes_storage/{selection}_GENOMES.db"
  output:
    fasta="data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomic_core.partitions"
  log:
    stdout="logs/anvi_pangenome_fasta_{selection}_mcl{mcl}.stdout",
    stderr="logs/anvi_pangenome_fasta_{selection}_mcl{mcl}.stderr"
  shell:
    """
    anvi-get-sequences-for-gene-clusters -p {input.pangenome}        \
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
    fasta="data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomic_core.fasta",
    partition="data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomic_core.partitions"
  output:
    tree="data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomics/{selection}.treefile"
  params:
    pre=lambda w: expand ("data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomics/{selection}",
                          selection=w.selection,
                          mcl=w.mcl)
  threads: 6
  log:
    stdout="logs/IQtree/anvi_phylogenomic_{selection}_mcl{mcl}.stdout",
    stderr="logs/IQtree/anvi_phylogenomic_{selection}_mcl{mcl}.stderr"
  shell:
    """
    iqtree -s {input.fasta}     \
           -p {input.partition} \
           -m MFP+MERGE         \
           -b 100               \
           -nt {threads}        \
           -ntmax {threads}     \
           -pre {params.pre}    \
    > {log.stdout} 2> {log.stderr}
    """

rule all_azolla_associated_phylogenomics:
  input:
    expand("data/anvio_pangenomes/{selection}_mcl{mcl}_phylogenomics/{selection}.treefile",selection=SELECTION,mcl=7)

rule pangenome_enrichment_anvio:
  input:
    pangenome="data/anvio_pangenomes/{selection}/{selection}_mcl{mcl}-PAN.db",
    genomestorage="data/anvio_genomes_storage/{selection}_GENOMES.db"
  output:
    "data/anvio_pangenomes/{selection}_functional_enrichment/{category}/{annotation}_{selection}_mcl{mcl}.tab"
  log:
    stdout="logs/pangenome_enrichment/{category}_{annotation}_anvio_{selection}_mcl{mcl}.stdout",
    stderr="logs/pangenome_enrichment/{category}_{annotation}_anvio_{selection}_mcl{mcl}.stderr"
  shell:
    """
    anvi-compute-functional-enrichment-in-pan -p {input.pangenome}        \
                                              -g {input.genomestorage}    \
                                              --category-variable {wildcards.category} \
                                              --annotation-source {wildcards.annotation} \
                                              --include-gc-identity-as-function \
                                              -o {output}                 \
    > {log.stdout} 2> {log.stderr}
    """

ANNOTATION=['COG20_FUNCTION', 'COG20_CATEGORY', 'KEGG_Module', 'KEGG_Class', 'COG20_PATHWAY', 'KOfam']

rule collect_pangenome_enrichment_all_annotations:
  input:
    expand("data/anvio_pangenomes/functional_enrichment/subgroup/{annotation}_{{selection}}_mcl{{mcl}}.tab",annotation=ANNOTATION)
  output:
    touch("data/anvio_pangenomes/functional_enrichment/subgroup_ALL_{selection}_mcl{mcl}.touch")
