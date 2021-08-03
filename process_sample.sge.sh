#!/bin/bash
#$ -A 100humans
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4
#$ -o ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.err

# set umask to avoid locking each other out of directories
umask 002

SAMPLE=$1
LOCKFILE=samples/${SAMPLE}/process_sample.lock

# add lockfile to directory to prevent multiple simultaneous jobs
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --config sample=${SAMPLE} \
    --nolock \
    --local-cores 4 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend conda \
    --use-singularity --singularity-args '--nv ' \
    --latency-wait 120 \
    --cluster-config workflow/process_sample.cluster.sge.yaml \
    --cluster "qsub -j y -cwd -V \
                    -A {cluster.account} \
                    -q {cluster.partition} \
                    -pe smp {cluster.cpus} \
                    -o {cluster.out} \
                    -e {cluster.err} \
                    {cluster.extra} " \
    --snakefile workflow/process_sample.smk
