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
  shell:
    """
    anvi-gen-genomes-storage \
      -i {input.txt} \
      -o {output}
    """
