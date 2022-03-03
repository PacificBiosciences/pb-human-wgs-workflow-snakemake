#!/bin/bash

# set umask to avoid locking each other out of directories
umask 002

SAMPLE=$1
# LOCKFILE=samples/${SAMPLE}/process_sample.lock

# add lockfile to directory to prevent multiple simultaneous jobs
# lockfile -r 0 ${LOCKFILE} || exit 1
# trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

mkdir -p logs

source workflow/variables.env
export SINGULARITY_TMPDIR="$TEMP"
export SINGULARITY_BIND="$TEMP"

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --cores 78 \
    --printshellcmds \
    --config sample=${SAMPLE} \
    --nolock \
    --use-conda --conda-frontend mamba \
    --use-singularity \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_sample.smk \
    2>&1 | tee "logs/process_sample.${SAMPLE}.$(date -d 'today' +'%Y%m%d%H%M').log"
