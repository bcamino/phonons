#!/bin/bash
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=32
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
timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_super4_phonons MnO_super4_phonons &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_super4_phonons &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_super5_phonons MnO_super5_phonons &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_super5_phonons &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_super6_phonons MnO_super6_phonons &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_super6_phonons &
wait

