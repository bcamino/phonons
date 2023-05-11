#!/bin/bash
#SBATCH --nodes=61
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# Replace [budget code] below with your full project code
#SBATCH --account=e05-bulk-smw 
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --export=none

module load epcc-job-env
module load other-software
module load crystal

# Address the memory leak
export FI_MR_CACHE_MAX_COUNT=0

# Run calculations
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_-4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-5_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-4_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_-5 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_5_-5 &
wait 
