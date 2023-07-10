#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

TMPDIR_LOCAL="$(mktemp --tmpdir --directory)"

# create decoys file
cut \
    --fields 1 \
    "${snakemake_input[assembly_sequence_index]}" \
    > "${TMPDIR_LOCAL}/decoys.txt"

# create "gentrome" (genome + transcriptome concatenated)
cat \
    "${snakemake_input[transcriptome]}" \
    <(bgzip --decompress --stdout "${snakemake_input[assembly]}") \
    > "${TMPDIR_LOCAL}/gentrome.fa"

salmon index \
       --threads $(( snakemake[threads] - 1 )) \
       --kmerLen ${snakemake_params[kmerLen]} \
       --tmpdir "${TMPDIR_LOCAL}/salmon_index_tmpdir" \
       --transcripts "${TMPDIR_LOCAL}/gentrome.fa" \
       --index "${snakemake_output[index]}" \
       --decoys "${TMPDIR_LOCAL}/decoys.txt" \
       > "${snakemake_log[stdout]}" \
       2> "${snakemake_log[stderr]}"

# Clean up
rm \
    --recursive \
    --force \
    "${TMPDIR_LOCAL}"
