#!/bin/bash

#SBATCH --time 10:00:00

#SBATCH --output="/vols/bitbucket/tagami/job_output/slurm/%x-%A-std-output.out"
#SBATCH --error="/vols/bitbucket/tagami/job_output/slurm/%x-%A-err-output.out"

#SBATCH --mail-user=daiki.tagami@hertford.ox.ac.uk
#SBATCH --mail-type=ALL

#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

#SBATCH --cluster=swan
#SBATCH --partition=standard-statgen-cpu

source /vols/bitbucket/tagami/python_environment/tskit_env/bin/activate

snakemake --workflow-profile /vols/bitbucket/tagami/snakemake-quantitative-trait-simulation/profiles --use-conda

deactivate