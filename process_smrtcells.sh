#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

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
    --cluster-config workflow/process_smrtcells.cluster.yaml \
    --cluster "sbatch --account={cluster.account} \
                      --partition={cluster.partition} \
                      --cpus-per-task={cluster.cpus} \
                      --output={cluster.out} {cluster.extra} " \
    --snakefile workflow/process_smrtcells.smk
