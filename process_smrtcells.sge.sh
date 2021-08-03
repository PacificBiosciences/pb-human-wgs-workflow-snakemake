#!/bin/bash
#$ -A 100humans
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -o ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.err

# set umask to avoid locking each other out of directories
umask 002

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --local-cores 1 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend conda \
    --latency-wait 120 \
    --cluster-config workflow/process_smrtcells.cluster.sge.yaml \
    --cluster "qsub -j y -cwd -V \
                    -A {cluster.account} \
                    -q {cluster.partition} \
                    -pe smp {cluster.cpus} \
                    -o {cluster.out} \
                    -e {cluster.err} \
                    {cluster.extra} " \
    --snakefile workflow/process_smrtcells.smk
