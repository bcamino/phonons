#!/bin/bash
#SBATCH --nodes=42 
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
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_5 MnO_phonons_-5_5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_4 MnO_phonons_-5_4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_3 MnO_phonons_-5_3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_2 MnO_phonons_-5_2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_1 MnO_phonons_-5_1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_0 MnO_phonons_-5_0 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_-1 MnO_phonons_-5_-1 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_-2 MnO_phonons_-5_-2 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_-3 MnO_phonons_-5_-3 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_-4 MnO_phonons_-5_-4 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-5_-5 MnO_phonons_-5_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-4_-5 MnO_phonons_-4_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-3_-5 MnO_phonons_-3_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-2_-5 MnO_phonons_-2_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_-1_-5 MnO_phonons_-1_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_0_-5 MnO_phonons_0_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_1_-5 MnO_phonons_1_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_2_-5 MnO_phonons_2_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_3_-5 MnO_phonons_3_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_4_-5 MnO_phonons_4_-5 & 
timeout 2876m /work/e05/e05/isa/runCRYSTAL/Pcry_slurm_multi2 MnO_phonons_5_-5 MnO_phonons_5_-5 & 
wait
