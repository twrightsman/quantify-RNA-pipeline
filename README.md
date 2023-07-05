# RNA-seq quantification pipeline

## Quickstart

The `--directory` option will configure where the pipeline will both
read needed input files (_e.g._ genomes) and write intermediate
and output files.

```
$ cp samples.example.tsv path/to/working/directory/samples.tsv
$ nano path/to/working/directory/samples.tsv
$ snakemake --directory path/to/working/directory
```

## `samples.tsv`

```
sample_id       species genotype        library_layout  library_selection       location
SRR8848655      Oryza sativa    Nipponbare      paired-end      random  sra
SRR8848654      Oryza sativa    Nipponbare      paired-end      random  sra
```

If `location` is `sra`, then the `sample_id` is assumed to be an SRA
accession ID and the reads will be downloaded directly from the SRA.

If `location` is `local` then the `data/reads` directory in the
working directory will be searched for a directory named after the
`sample_id`. For `single-read` layout samples, there should be a
`sample_id.fq.gz` file in this folder. (_e.g._
`data/reads/SRR8848655/SRR8848655.fq.gz`) For `paired-end` samples,
there should be `sample_id_1.fq.gz` and `sample_id_2.fq.gz` files in
this folder.

`library_selection` should be either `random` or `polyA` (for 3'
RNA-seq).

The reference genome will be searched for in the working directory
based to `run.sh`. It will be expecting a `bgzip`ed assembly at
`data/genomes/species/genotype/assembly.fa.gz`, and it has to be
indexed with `samtools faidx` (_e.g._
`Oryza_sativa/Nipponbare/assembly.fa.gz`). Spaces in the species name
must be replaced with underscores in the path names.

An annotation should also exist at `data/genomes/species/genotype/annotation.gff.gz`.
