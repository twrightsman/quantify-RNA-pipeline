#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

OUTPUT_DIRECTORY="$(dirname "${snakemake_output[0]}")"

fastq-dump \
    --disable-multithreading \
    --split-3 \
    --outdir "$OUTPUT_DIRECTORY" \
    "${snakemake_input[sra]}" \
    > "${snakemake_log[stdout]}" \
    2> "${snakemake_log[stderr]}"
