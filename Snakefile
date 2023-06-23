from pathlib import Path

import pandas as pd


configfile: "config.yml"
localrules: all, link_genomes

for key in config['paths']:
    config['paths'][key] = Path(config['paths'][key])

SAMPLES = pd.read_table(
  config['paths']['samples'],
  header = 0,
  index_col = 'sample_id'
)
QUANTS = [f"work/quants/{sample_id}/quant.sf" for sample_id in SAMPLES.index]

rule all:
  input: QUANTS


rule prefetch_SRA_accession:
  output:
    sra = temp("work/data/sra/{accession}/{accession}.sra")
  params:
    max_size = "50G"
  log:
    stdout = "work/data/sra/{accession}/{accession}.prefetch.out",
    stderr = "work/data/sra/{accession}/{accession}.prefetch.err"
  conda:
    "envs/prefetch_SRA_accession.yml"
  script:
    "scripts/prefetch_SRA_accession.sh"


rule dump_SRA_accession_FASTQ_single:
  input:
    sra = rules.prefetch_SRA_accession.output.sra
  output:
    reads = temp("work/data/reads/{accession}/{accession}.fastq")
  log:
    stdout = "work/data/reads/{accession}/{accession}.fastq-dump.out",
    stderr = "work/data/reads/{accession}/{accession}.fastq-dump.err"
  conda:
    "envs/dump_SRA_accession_FASTQ.yml"
  script:
    "scripts/dump_SRA_accession_FASTQ.sh"


rule dump_SRA_accession_FASTQ_paired:
  input:
    sra = rules.prefetch_SRA_accession.output.sra
  output:
    reads_mate1 = temp("work/data/reads/{accession}/{accession}_1.fastq"),
    reads_mate2 = temp("work/data/reads/{accession}/{accession}_2.fastq")
  log:
    stdout = "work/data/reads/{accession}/{accession}.fastq-dump.out",
    stderr = "work/data/reads/{accession}/{accession}.fastq-dump.err"
  conda:
    "envs/dump_SRA_accession_FASTQ.yml"
  script:
    "scripts/dump_SRA_accession_FASTQ.sh"


def get_trim_reads_single_input(wildcards):
    sample_location = SAMPLES.loc[wildcards.sample_id, 'location']
    if sample_location == 'sra':
      # if on SRA, sample_id is SRA accession
      return {
        'reads': "work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}.fastq"
      }
    elif sample_location == 'local':
      return {
        'reads': config['paths']['reads'] / f"{wildcards.sample_id}/{wildcards.sample_id}.fq.gz"
      }

rule trim_reads_single:
  input:
    unpack(get_trim_reads_single_input)
  output:
    reads_trimmed = temp("work/data/reads/{sample_id}/{sample_id}.trimmed.fastq"),
    report_html = protected("work/data/reads/{sample_id}/{sample_id}.fastp.html"),
    report_json = protected("work/data/reads/{sample_id}/{sample_id}.fastp.json")
  threads: 1
  priority: 5
  log:
    stdout = "work/data/reads/{sample_id}/{sample_id}.fastp.out",
    stderr = "work/data/reads/{sample_id}/{sample_id}.fastp.err"
  conda:
    "envs/trim_reads.yml"
  script:
    "scripts/trim_reads_single.sh"


def get_trim_reads_paired_input(wildcards):
    sample_location = SAMPLES.loc[wildcards.sample_id, 'location']
    if sample_location == 'sra':
      # if on SRA, sample_id is SRA accession
      return {
        'reads_mate1': "work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}_1.fastq",
        'reads_mate2': "work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}_2.fastq"
      }
    elif sample_location == 'local':
      return {
        'reads_mate1': config['paths']['reads'] / f"{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
        'reads_mate2': config['paths']['reads'] / f"{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
      }

rule trim_reads_paired:
  input:
    unpack(get_trim_reads_paired_input)
  output:
    reads_mate1_trimmed = temp("work/data/reads/{sample_id}/{sample_id}_1.trimmed.fastq"),
    reads_mate2_trimmed = temp("work/data/reads/{sample_id}/{sample_id}_2.trimmed.fastq"),
    report_html = protected("work/data/reads/{sample_id}/{sample_id}.fastp.html"),
    report_json = protected("work/data/reads/{sample_id}/{sample_id}.fastp.json")
  threads: 2
  priority: 5
  log:
    stdout = "work/data/reads/{sample_id}/{sample_id}.fastp.out",
    stderr = "work/data/reads/{sample_id}/{sample_id}.fastp.err"
  conda:
    "envs/trim_reads.yml"
  script:
    "scripts/trim_reads_paired.sh"


rule compress_trimmed_reads:
  input:
    reads_trimmed = "work/data/reads/{sample_id}/{fastq_filename}.trimmed.fastq"
  output:
    reads_trimmed_compressed = protected("work/data/reads/{sample_id}/{fastq_filename}.trimmed.fq.gz")
  threads: 2
  priority: 5
  conda:
    "envs/compress_trimmed_reads.yml"
  shell:
    "bgzip --threads {threads} --stdout {input.reads_trimmed} > {output.reads_trimmed_compressed}"


rule link_genomes:
  input:
    assembly = config['paths']['genomes'] / "{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = config['paths']['genomes'] / "{species}/{genotype}/assembly.fa.gz.fai",
    assembly_bgzip_index = config['paths']['genomes'] / "{species}/{genotype}/assembly.fa.gz.gzi",
    annotation = config['paths']['genomes'] / "{species}/{genotype}/annotation.gff.gz"
  output:
    assembly = "work/data/genomes/{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = "work/data/genomes/{species}/{genotype}/assembly.fa.gz.fai",
    assembly_bgzip_index = "work/data/genomes/{species}/{genotype}/assembly.fa.gz.gzi",
    annotation = "work/data/genomes/{species}/{genotype}/annotation.gff.gz"
  shell:
    "ln --symbolic --relative {input.assembly} {output.assembly} && " +
    "ln --symbolic --relative {input.assembly_sequence_index} {output.assembly_sequence_index} && " +
    "ln --symbolic --relative {input.assembly_bgzip_index} {output.assembly_bgzip_index} && " +
    "ln --symbolic --relative {input.annotation} {output.annotation}"


rule create_transcriptome:
  input:
    assembly = "work/data/genomes/{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = "work/data/genomes/{species}/{genotype}/assembly.fa.gz.fai",
    assembly_bgzip_index = "work/data/genomes/{species}/{genotype}/assembly.fa.gz.gzi",
    annotation = "work/data/genomes/{species}/{genotype}/annotation.gff.gz"
  output:
    transcriptome = "work/data/genomes/{species}/{genotype}/transcriptome.fa"
  priority: -5
  log:
    stderr = "work/data/genomes/{species}/{genotype}/create_transcriptome.err"
  conda:
    "envs/create_transcriptome.yml"
  shell:
    "scripts/create_transcriptome.py {input.assembly:q} {input.annotation:q} > {output.transcriptome:q} 2> {log.stderr:q}"


rule index_genome:
  input:
    assembly = "work/data/genomes/{species}/{genotype}/assembly.fa.gz",
    assembly_sequence_index = "work/data/genomes/{species}/{genotype}/assembly.fa.gz.fai",
    transcriptome = "work/data/genomes/{species}/{genotype}/transcriptome.fa"
  output:
    index = directory("work/data/genomes/{species}/{genotype}/salmon_index")
  params:
    kmerLen = 31
  threads: 4
  priority: -5
  conda:
    "envs/index_genome.yml"
  log:
    stdout = "work/data/genomes/{species}/{genotype}/salmon_index.out",
    stderr = "work/data/genomes/{species}/{genotype}/salmon_index.err"
  script:
    "scripts/index_genome.sh"


def get_quantify_RNA_input(wildcards):
  sample_library_layout = SAMPLES.loc[wildcards.sample_id, 'library_layout']

  species = SAMPLES.loc[wildcards.sample_id, 'species'].replace(' ', '_')
  genotype = SAMPLES.loc[wildcards.sample_id, 'genotype'].replace(' ', '_')
  inputs = {
    'index': f"work/data/genomes/{species}/{genotype}/salmon_index"
  }

  if sample_library_layout == 'paired-end':
    inputs['reads_mate1'] = f"work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}_1.trimmed.fq.gz"
    inputs['reads_mate2'] = f"work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}_2.trimmed.fq.gz"
  elif sample_library_layout == 'single-read':
    inputs['reads'] = f"work/data/reads/{wildcards.sample_id}/{wildcards.sample_id}.trimmed.fq.gz"
  else:
    raise ValueError(f"Unknown library layout '{sample_library_layout}' for sample '{wildcards.sample_id}'")

  return inputs

rule quantify_RNA:
  input:
    unpack(get_quantify_RNA_input)
  output:
    quants = "work/quants/{sample_id}/quant.sf"
  params:
    library_selection = lambda wildcards: SAMPLES.loc[wildcards.sample_id, 'library_selection']
  threads: 4
  priority: -5
  conda:
    "envs/quantify_RNA.yml"
  log:
    stdout = "work/quants/{sample_id}/salmon.out",
    stderr = "work/quants/{sample_id}/salmon.err"
  script:
    "scripts/quantify_RNA.sh"
