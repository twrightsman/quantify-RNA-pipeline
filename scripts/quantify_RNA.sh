#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

OUTPUT_DIRECTORY="$(dirname "${snakemake_output[quants]}")"

if [ -z "${snakemake_input[reads]:-}"  ]; then
    # paired-end
    SALMON_READ_INPUT_ARGS="--mates1 ${snakemake_input[reads_mate1]} --mates2 ${snakemake_input[reads_mate2]}"
else
    # single-read
    SALMON_READ_INPUT_ARGS="--unmatedReads ${snakemake_input[reads]}"
fi

# don't do length correction with 3' RNA-seq
if [ "${snakemake_params[library_selection]}" == 'polyA' ]; then
    SALMON_EXTRA_ARGS="--noLengthCorrection"
else
    SALMON_EXTRA_ARGS=""
fi

salmon quant \
       --threads $(( snakemake[threads] - 1 )) \
       --libType A \
       --index "${snakemake_input[index]}" \
       $SALMON_READ_INPUT_ARGS \
       --output "$OUTPUT_DIRECTORY" \
       $SALMON_EXTRA_ARGS \
       > "${snakemake_log[stdout]}" \
       2> "${snakemake_log[stderr]}"
