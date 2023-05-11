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
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_0 MnO_phonons_0_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_0 MnO_phonons_1_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_1 MnO_phonons_1_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_1 MnO_phonons_0_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_1 MnO_phonons_-1_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_0 MnO_phonons_-1_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_-1 MnO_phonons_-1_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_-1 MnO_phonons_0_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_-1 MnO_phonons_1_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_-1 MnO_phonons_2_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_0 MnO_phonons_2_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_1 MnO_phonons_2_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_2 MnO_phonons_2_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_2 MnO_phonons_1_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_2 MnO_phonons_0_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_2 MnO_phonons_-1_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_2 MnO_phonons_-2_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_1 MnO_phonons_-2_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_0 MnO_phonons_-2_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_-1 MnO_phonons_-2_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_-2 MnO_phonons_-2_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_-2 MnO_phonons_-1_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_-2 MnO_phonons_0_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_-2 MnO_phonons_1_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_-2 MnO_phonons_2_-2 & 
wait
