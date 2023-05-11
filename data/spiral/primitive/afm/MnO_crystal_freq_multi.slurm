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
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-5_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-4_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-3_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-2_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_-1_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_0_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_1_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_2_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_3_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_4_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi MnO_freq_5_-5 & 
wait 
