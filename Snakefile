NANOPORE=['Azfil_lab' 'Azpinnata' 'Azsp_bordeaux']

rule gather_anvi_MAGS:
  output:
    directory("data/MAG_anvi_dbs/")
  shell:
    "bash ./scripts/collect_mag_dbs.sh"

rule create_pangenome_storage_internal:
  input:
    "data/MAG_anvi_dbs/",
    txt="data/Anvio_internal_genomes.txt"
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

rule get_reference_fastas:
  output:
    "references/Nazollae_0708.fasta",
    "references_Azfil_v1.fasta",
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
    bam=  "data/nanopore_mapped/{nanopore_host}_mapped.bam"
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
