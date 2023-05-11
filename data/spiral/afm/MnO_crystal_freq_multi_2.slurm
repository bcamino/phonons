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
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_-2 MnO_phonons_4_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_-1 MnO_phonons_4_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_0 MnO_phonons_4_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_1 MnO_phonons_4_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_2 MnO_phonons_4_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_3 MnO_phonons_4_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_4 MnO_phonons_4_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_4 MnO_phonons_3_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_4 MnO_phonons_2_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_4 MnO_phonons_1_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_4 MnO_phonons_0_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_4 MnO_phonons_-1_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_4 MnO_phonons_-2_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_4 MnO_phonons_-3_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_4 MnO_phonons_-4_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_3 MnO_phonons_-4_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_2 MnO_phonons_-4_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_1 MnO_phonons_-4_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_0 MnO_phonons_-4_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_-1 MnO_phonons_-4_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_-2 MnO_phonons_-4_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_-3 MnO_phonons_-4_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_-4 MnO_phonons_-4_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_-4 MnO_phonons_-3_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_-4 MnO_phonons_-2_-4 & 
wait
