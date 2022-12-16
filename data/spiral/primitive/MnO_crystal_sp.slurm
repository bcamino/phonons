#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# Replace [budget code] below with your full project code
#SBATCH --account=e05-power-smw
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --export=none

module load epcc-job-env
module load other-software
module load crystal

# Address the memory leak
export FI_MR_CACHE_MAX_COUNT=0

# Run calculations
timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_freq_5_-6 MnO_cpks_5_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_freq_5_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_freq_6_-6 MnO_cpks_6_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_freq_6_-6 &
wait

