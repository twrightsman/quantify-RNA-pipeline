#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail \
  -o xtrace

fastp \
    --thread ${snakemake[threads]} \
    --in1 "${snakemake_input[reads]}" \
    --out1 "${snakemake_output[reads_trimmed]}" \
    --html "${snakemake_output[report_html]}" \
    --json "${snakemake_output[report_json]}" \
    --report_title "${snakemake_wildcards[sample_id]}" \
    > "${snakemake_log[stdout]}" \
    2> "${snakemake_log[stderr]}"
