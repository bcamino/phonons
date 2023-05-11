#!/bin/bash
#SBATCH --nodes=1 
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
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-3_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-2_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_-1_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_0_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_-3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_-2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_-1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_0 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_1 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_2 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_3 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_4_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_3_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_2_4 &
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_cpks_1_4 &
wait
