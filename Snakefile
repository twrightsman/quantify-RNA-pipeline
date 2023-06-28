from pathlib import Path

import pandas as pd


localrules: all

SAMPLES = pd.read_table(
  'samples.tsv',
  header = 0,
  index_col = 'sample_id'
)
QUANTS = [f"quants/{sample_id}/quant.gene.tsv" for sample_id in SAMPLES.index]

rule all:
  input: QUANTS


rule prefetch_SRA_accession:
  output:
    sra = temp("data/sra/{accession}/{accession}.sra")
  params:
    max_size = "50G"
  log:
    stdout = "data/sra/{accession}/{accession}.prefetch.out",
    stderr = "data/sra/{accession}/{accession}.prefetch.err"
  conda:
    "envs/prefetch_SRA_accession.yml"
  script:
    "scripts/prefetch_SRA_accession.sh"


rule dump_SRA_accession_FASTQ_single:
  input:
    sra = rules.prefetch_SRA_accession.output.sra
  output:
    reads = temp("data/reads/{accession}/{accession}.fastq")
  log:
    stdout = "data/reads/{accession}/{accession}.fastq-dump.out",
    stderr = "data/reads/{accession}/{accession}.fastq-dump.err"
  conda:
    "envs/dump_SRA_accession_FASTQ.yml"
  script:
    "scripts/dump_SRA_accession_FASTQ.sh"


rule dump_SRA_accession_FASTQ_paired:
  input:
    sra = rules.prefetch_SRA_accession.output.sra
  output:
    reads_mate1 = temp("data/reads/{accession}/{accession}_1.fastq"),
    reads_mate2 = temp("data/reads/{accession}/{accession}_2.fastq")
  log:
    stdout = "data/reads/{accession}/{accession}.fastq-dump.out",
    stderr = "data/reads/{accession}/{accession}.fastq-dump.err"
  conda:
    "envs/dump_SRA_accession_FASTQ.yml"
  script:
    "scripts/dump_SRA_accession_FASTQ.sh"


def get_trim_reads_single_input(wildcards):
    sample_location = SAMPLES.loc[wildcards.sample_id, 'location']
    if sample_location == 'sra':
      # if on SRA, sample_id is SRA accession
      return {
        'reads': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}.fastq"
      }
    elif sample_location == 'local':
      return {
        'reads': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}.fq.gz"
      }

rule trim_reads_single:
  input:
    unpack(get_trim_reads_single_input)
  output:
    reads_trimmed = temp("data/reads/{sample_id}/{sample_id}.trimmed.fastq"),
    report_html = protected("data/reads/{sample_id}/{sample_id}.fastp.html"),
    report_json = protected("data/reads/{sample_id}/{sample_id}.fastp.json")
  threads: 1
  priority: 5
  log:
    stdout = "data/reads/{sample_id}/{sample_id}.fastp.out",
    stderr = "data/reads/{sample_id}/{sample_id}.fastp.err"
  conda:
    "envs/trim_reads.yml"
  script:
    "scripts/trim_reads_single.sh"


def get_trim_reads_paired_input(wildcards):
    sample_location = SAMPLES.loc[wildcards.sample_id, 'location']
    if sample_location == 'sra':
      # if on SRA, sample_id is SRA accession
      return {
        'reads_mate1': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_1.fastq",
        'reads_mate2': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_2.fastq"
      }
    elif sample_location == 'local':
      return {
        'reads_mate1': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
        'reads_mate2': f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
      }

rule trim_reads_paired:
  input:
    unpack(get_trim_reads_paired_input)
  output:
    reads_mate1_trimmed = temp("data/reads/{sample_id}/{sample_id}_1.trimmed.fastq"),
    reads_mate2_trimmed = temp("data/reads/{sample_id}/{sample_id}_2.trimmed.fastq"),
    report_html = protected("data/reads/{sample_id}/{sample_id}.fastp.html"),
    report_json = protected("data/reads/{sample_id}/{sample_id}.fastp.json")
  threads: 2
  priority: 5
  log:
    stdout = "data/reads/{sample_id}/{sample_id}.fastp.out",
    stderr = "data/reads/{sample_id}/{sample_id}.fastp.err"
  conda:
    "envs/trim_reads.yml"
  script:
    "scripts/trim_reads_paired.sh"


rule compress_trimmed_reads:
  input:
    reads_trimmed = "data/reads/{sample_id}/{fastq_filename}.trimmed.fastq"
  output:
    reads_trimmed_compressed = protected("data/reads/{sample_id}/{fastq_filename}.trimmed.fq.gz")
  threads: 2
  priority: 5
  conda:
    "envs/compress_trimmed_reads.yml"
  shell:
    "bgzip --threads {threads} --stdout {input.reads_trimmed} > {output.reads_trimmed_compressed}"


rule create_transcriptome:
  input:
    assembly = "data/genomes/{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = "data/genomes/{species}/{genotype}/assembly.fa.gz.fai",
    assembly_bgzip_index = "data/genomes/{species}/{genotype}/assembly.fa.gz.gzi",
    annotation = "data/genomes/{species}/{genotype}/annotation.gff.gz"
  output:
    transcriptome = "data/genomes/{species}/{genotype}/transcriptome.fa"
  priority: -5
  conda:
    "envs/create_transcriptome.yml"
  params:
    verbosity = 0
  script:
    "scripts/create_transcriptome.py"


rule index_genome:
  input:
    assembly = "data/genomes/{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = "data/genomes/{species}/{genotype}/assembly.fa.gz.fai",
    transcriptome = "data/genomes/{species}/{genotype}/transcriptome.fa"
  output:
    index = directory("data/genomes/{species}/{genotype}/salmon_index")
  params:
    kmerLen = 31
  threads: 4
  priority: -5
  conda:
    "envs/index_genome.yml"
  log:
    stdout = "data/genomes/{species}/{genotype}/salmon_index.out",
    stderr = "data/genomes/{species}/{genotype}/salmon_index.err"
  script:
    "scripts/index_genome.sh"


def get_quantify_RNA_input(wildcards):
  sample_library_layout = SAMPLES.loc[wildcards.sample_id, 'library_layout']

  species = SAMPLES.loc[wildcards.sample_id, 'species'].replace(' ', '_')
  genotype = SAMPLES.loc[wildcards.sample_id, 'genotype'].replace(' ', '_')
  inputs = {
    'index': f"data/genomes/{species}/{genotype}/salmon_index"
  }

  if sample_library_layout == 'paired-end':
    inputs['reads_mate1'] = f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_1.trimmed.fq.gz"
    inputs['reads_mate2'] = f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}_2.trimmed.fq.gz"
  elif sample_library_layout == 'single-read':
    inputs['reads'] = f"data/reads/{wildcards.sample_id}/{wildcards.sample_id}.trimmed.fq.gz"
  else:
    raise ValueError(f"Unknown library layout '{sample_library_layout}' for sample '{wildcards.sample_id}'")

  return inputs

rule quantify_RNA:
  input:
    unpack(get_quantify_RNA_input)
  output:
    quants = "quants/{sample_id}/quant.sf"
  params:
    library_selection = lambda wildcards: SAMPLES.loc[wildcards.sample_id, 'library_selection']
  threads: 4
  priority: -5
  conda:
    "envs/quantify_RNA.yml"
  log:
    stdout = "quants/{sample_id}/salmon.out",
    stderr = "quants/{sample_id}/salmon.err"
  script:
    "scripts/quantify_RNA.sh"


rule map_tx_to_gene:
  input:
    annotation = "data/genomes/{species}/{genotype}/annotation.gff.gz"
  output:
    tx2gene = "data/genomes/{species}/{genotype}/tx2gene.tsv"
  priority: -5
  script:
    "scripts/map_tx_to_gene.py"


def get_summarize_by_gene_input(wildcards):
  species = SAMPLES.loc[wildcards.sample_id, 'species'].replace(' ', '_')
  genotype = SAMPLES.loc[wildcards.sample_id, 'genotype'].replace(' ', '_')
  return {
    'tx2gene': f"data/genomes/{species}/{genotype}/tx2gene.tsv",
    'quants_tx': f"quants/{wildcards.sample_id}/quant.sf"
  }


rule summarize_by_gene:
  input:
    unpack(get_summarize_by_gene_input)
  output:
    quants_gene = "quants/{sample_id}/quant.gene.tsv"
  conda:
    "envs/summarize_by_gene.yml"
  script:
    "scripts/summarize_by_gene.R"
