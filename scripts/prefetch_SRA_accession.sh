#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail \
  -o xtrace

ACCESSION_DIRECTORY="$(dirname "${snakemake_output[sra]}")"
OUTPUT_DIRECTORY="$(dirname "$ACCESSION_DIRECTORY")"

prefetch \
    --verify yes \
    --max-size ${snakemake_params[max_size]} \
    --output-directory "$OUTPUT_DIRECTORY" \
    "${snakemake_wildcards[accession]}" \
    > "${snakemake_log[stdout]}" \
    2> "${snakemake_log[stderr]}"
