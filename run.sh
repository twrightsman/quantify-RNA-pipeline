#!/usr/bin/env bash

set -euo pipefail

snakemake \
    --use-conda \
    $*
