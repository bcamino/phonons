#!/bin/bash
#SBATCH --nodes=50 
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# Replace [budget code] below with your full project code
#SBATCH --account=e05-bulk-smw 
#SBATCH --partition=standard
#SBATCH --qos=long
#SBATCH --export=none

module load epcc-job-env
module load other-software
module load crystal

# Address the memory leak
export FI_MR_CACHE_MAX_COUNT=0

# Run calculations
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_-2 MnO_phonons_3_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_-1 MnO_phonons_3_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_0 MnO_phonons_3_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_1 MnO_phonons_3_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_2 MnO_phonons_3_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_3 MnO_phonons_3_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_3 MnO_phonons_2_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_3 MnO_phonons_1_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_3 MnO_phonons_0_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_3 MnO_phonons_-1_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_3 MnO_phonons_-2_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_3 MnO_phonons_-3_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_2 MnO_phonons_-3_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_1 MnO_phonons_-3_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_0 MnO_phonons_-3_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_-1 MnO_phonons_-3_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_-2 MnO_phonons_-3_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_-3 MnO_phonons_-3_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_-3 MnO_phonons_-2_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_-3 MnO_phonons_-1_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_-3 MnO_phonons_0_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_-3 MnO_phonons_1_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_-3 MnO_phonons_2_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_-3 MnO_phonons_3_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_-3 MnO_phonons_4_-3 & 
wait
