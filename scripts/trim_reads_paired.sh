#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# --thread doesn't seem to include main fastp thread
# subtract one to avoid using more than requested threads
fastp \
    --detect_adapter_for_pe \
    --thread $(( snakemake[threads] - 1 )) \
    --in1 "${snakemake_input[reads_mate1]}" \
    --in2 "${snakemake_input[reads_mate2]}" \
    --out1 "${snakemake_output[reads_mate1_trimmed]}" \
    --out2 "${snakemake_output[reads_mate2_trimmed]}" \
    --html "${snakemake_output[report_html]}" \
    --json "${snakemake_output[report_json]}" \
    --report_title "${snakemake_wildcards[sample_id]}" \
    > "${snakemake_log[stdout]}" \
    2> "${snakemake_log[stderr]}"
