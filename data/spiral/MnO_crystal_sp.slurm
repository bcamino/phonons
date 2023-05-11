#!/bin/bash
#SBATCH --nodes=16
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
timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_0   &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_0 MnO_phonons_0_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_1 MnO_phonons_1_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_1 MnO_phonons_1_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_1 MnO_phonons_0_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_0 MnO_phonons_-1_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-1 MnO_phonons_-1_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-1 MnO_phonons_-1_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-1 MnO_phonons_0_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-1 MnO_phonons_1_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_0 MnO_phonons_2_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_1 MnO_phonons_2_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_2 MnO_phonons_2_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_2 MnO_phonons_2_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_2 MnO_phonons_1_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_2 MnO_phonons_0_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_2 MnO_phonons_-1_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_1 MnO_phonons_-2_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_0 MnO_phonons_-2_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-1 MnO_phonons_-2_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-2 MnO_phonons_-2_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-2 MnO_phonons_-2_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-2 MnO_phonons_-1_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-2 MnO_phonons_0_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-2 MnO_phonons_1_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-2 MnO_phonons_2_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-1 MnO_phonons_3_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_0 MnO_phonons_3_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_1 MnO_phonons_3_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_2 MnO_phonons_3_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_3 MnO_phonons_3_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_3 MnO_phonons_3_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_3 MnO_phonons_2_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_3 MnO_phonons_1_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_3 MnO_phonons_0_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_3 MnO_phonons_-1_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_3 MnO_phonons_-2_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_2 MnO_phonons_-3_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_1 MnO_phonons_-3_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_0 MnO_phonons_-3_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-1 MnO_phonons_-3_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-2 MnO_phonons_-3_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-3 MnO_phonons_-3_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-3 MnO_phonons_-3_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-3 MnO_phonons_-2_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-3 MnO_phonons_-1_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-3 MnO_phonons_0_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-3 MnO_phonons_1_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-3 MnO_phonons_2_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-3 MnO_phonons_3_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-2 MnO_phonons_4_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-1 MnO_phonons_4_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_0 MnO_phonons_4_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_1 MnO_phonons_4_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_2 MnO_phonons_4_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_3 MnO_phonons_4_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_4 MnO_phonons_4_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_4 MnO_phonons_4_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_4 MnO_phonons_3_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_4 MnO_phonons_2_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_4 MnO_phonons_1_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_4 MnO_phonons_0_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_4 MnO_phonons_-1_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_4 MnO_phonons_-2_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_4 MnO_phonons_-3_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_3 MnO_phonons_-4_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_2 MnO_phonons_-4_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_1 MnO_phonons_-4_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_0 MnO_phonons_-4_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-1 MnO_phonons_-4_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-2 MnO_phonons_-4_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-3 MnO_phonons_-4_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-4 MnO_phonons_-4_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-4 MnO_phonons_-4_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-4 MnO_phonons_-3_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-4 MnO_phonons_-2_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-4 MnO_phonons_-1_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-4 MnO_phonons_0_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-4 MnO_phonons_1_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-4 MnO_phonons_2_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-4 MnO_phonons_3_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-4 MnO_phonons_4_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-3 MnO_phonons_5_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-2 MnO_phonons_5_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-1 MnO_phonons_5_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_0 MnO_phonons_5_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_1 MnO_phonons_5_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_2 MnO_phonons_5_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_3 MnO_phonons_5_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_4 MnO_phonons_5_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_5 MnO_phonons_5_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_5 MnO_phonons_5_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_5 MnO_phonons_4_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_5 MnO_phonons_3_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_5 MnO_phonons_2_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_5 MnO_phonons_1_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_5 MnO_phonons_0_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_5 MnO_phonons_-1_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_5 MnO_phonons_-2_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_5 MnO_phonons_-3_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_5 MnO_phonons_-4_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_4 MnO_phonons_-5_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_3 MnO_phonons_-5_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_2 MnO_phonons_-5_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_1 MnO_phonons_-5_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_0 MnO_phonons_-5_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-1 MnO_phonons_-5_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-2 MnO_phonons_-5_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-3 MnO_phonons_-5_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-4 MnO_phonons_-5_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-5 MnO_phonons_-5_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-5 MnO_phonons_-5_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-5 MnO_phonons_-4_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-5 MnO_phonons_-3_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-5 MnO_phonons_-2_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-5 MnO_phonons_-1_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-5 MnO_phonons_0_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-5 MnO_phonons_1_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-5 MnO_phonons_2_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-5 MnO_phonons_3_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-5 MnO_phonons_4_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-5 MnO_phonons_5_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-4 MnO_phonons_6_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-3 MnO_phonons_6_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-2 MnO_phonons_6_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-1 MnO_phonons_6_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_0 MnO_phonons_6_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_1 MnO_phonons_6_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_2 MnO_phonons_6_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_3 MnO_phonons_6_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_4 MnO_phonons_6_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_5 MnO_phonons_6_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_6 MnO_phonons_6_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_6 MnO_phonons_6_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_6 MnO_phonons_5_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_6 MnO_phonons_4_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_6 MnO_phonons_3_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_6 MnO_phonons_2_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_6 MnO_phonons_1_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_6 MnO_phonons_0_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_6 MnO_phonons_-1_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_6 MnO_phonons_-2_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_6 MnO_phonons_-3_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_6 MnO_phonons_-4_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_6 MnO_phonons_-5_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_5 MnO_phonons_-6_6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_4 MnO_phonons_-6_5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_3 MnO_phonons_-6_4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_2 MnO_phonons_-6_3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_1 MnO_phonons_-6_2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_0 MnO_phonons_-6_1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_0 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-1 MnO_phonons_-6_0  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-1 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-2 MnO_phonons_-6_-1  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-2 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-3 MnO_phonons_-6_-2  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-3 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-4 MnO_phonons_-6_-3  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-4 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-5 MnO_phonons_-6_-4  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-5 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-6_-6 MnO_phonons_-6_-5  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-6_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-5_-6 MnO_phonons_-6_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-5_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-4_-6 MnO_phonons_-5_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-4_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-3_-6 MnO_phonons_-4_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-3_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-2_-6 MnO_phonons_-3_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-2_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_-1_-6 MnO_phonons_-2_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_-1_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_0_-6 MnO_phonons_-1_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_0_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_1_-6 MnO_phonons_0_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_1_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_2_-6 MnO_phonons_1_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_2_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_3_-6 MnO_phonons_2_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_3_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_4_-6 MnO_phonons_3_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_4_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_5_-6 MnO_phonons_4_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_5_-6 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm MnO_phonons_6_-6 MnO_phonons_5_-6  &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys MnO_phonons_6_-6 &
wait

