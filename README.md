# RNA-seq quantification pipeline

## Quickstart

The `path` to a directory of `reads` is optional if you are
downloading from the SRA.

```
$ cp config.example.yml config.yml
$ nano config.yml
$ cp samples.example.tsv samples.tsv
$ nano samples.tsv
$ ./run.sh --cores 10
```

## `samples.tsv`

```
sample_id       species genotype        library_layout  library_selection       location
SRR8848655      Oryza sativa    Nipponbare      paired-end      random  sra
SRR8848654      Oryza sativa    Nipponbare      paired-end      random  sra
```

If `location` is `sra`, then the `sample_id` is assumed to be an SRA
accession ID and the reads will be downloaded directly from the SRA.

If `location` is `local` then the `reads` path in `config.yml` will be
searched for a directory named after the `sample_id`. For
`single-read` layout samples, there should be a `sample_id.fq.gz` file
in this folder. (_e.g._ `SRR8848655/SRR8848655.fq.gz`) For
`paired-end` samples, there should be `sample_id_1.fq.gz` and
`sample_id_2.fq.gz` files in this folder.

`library_selection` should be either `random` or `polyA` (for 3'
RNA-seq).

The reference genome will be searched for using the `genomes` path in
`config.yml`. It will be expecting a `bgzip`ed assembly at
`species/genotype/assembly.fa.gz`, and it has to be indexed with
`samtools faidx`. (_e.g._ `Oryza_sativa/Nipponbare/assembly.fa.gz`)

An annotation should also exist at `species/genotype/annotation.gff.gz`.
