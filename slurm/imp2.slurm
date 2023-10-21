#!/bin/bash
#SBATCH --job-name=imp3
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=log/imp3.out
#SBATCH --error=log/imp3.err
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=24:00:00

#Example:
# sbatch slurm/imp2.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Run the tar_make_clustermq() command with the specified number of clusters
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/run_script.R 1 > log/imp3.log

date
