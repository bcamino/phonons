#!/bin/bash
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# Replace [budget code] below with your full project code
#SBATCH --account=e05-pool     
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --export=none

module load epcc-job-env
module load other-software
module load crystal

# Address the memory leak
export FI_MR_CACHE_MAX_COUNT=0

# Run calculations
timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_0-001 supercell_0_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_0-002 supercell_0_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_0-000 supercell_0_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_0-001 supercell_1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_0-002 supercell_1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_1-000 supercell_1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_1-001 supercell_1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_1-002 supercell_1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_1-000 supercell_1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_1-001 supercell_0_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_1-002 supercell_0_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_1-000 supercell_0_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_1-001 supercell_-1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_1-002 supercell_-1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_0-000 supercell_-1_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_0-001 supercell_-1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_0-002 supercell_-1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-1-000 supercell_-1_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-1-001 supercell_-1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-1-002 supercell_-1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-1-000 supercell_-1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-1-001 supercell_0_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-1-002 supercell_0_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-1-000 supercell_0_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-1-001 supercell_1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-1-002 supercell_1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-1-000 supercell_1_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-1-001 supercell_2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-1-002 supercell_2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_0-000 supercell_2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_0-001 supercell_2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_0-002 supercell_2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_1-000 supercell_2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_1-001 supercell_2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_1-002 supercell_2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_2-000 supercell_2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_2-001 supercell_2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_2-002 supercell_2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_2-000 supercell_2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_2-001 supercell_1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_2-002 supercell_1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_2-000 supercell_1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_2-001 supercell_0_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_2-002 supercell_0_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_2-000 supercell_0_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_2-001 supercell_-1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_2-002 supercell_-1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_2-000 supercell_-1_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_2-001 supercell_-2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_2-002 supercell_-2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_1-000 supercell_-2_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_1-001 supercell_-2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_1-002 supercell_-2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_0-000 supercell_-2_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_0-001 supercell_-2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_0-002 supercell_-2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-1-000 supercell_-2_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-1-001 supercell_-2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-1-002 supercell_-2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-2-000 supercell_-2_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-2-001 supercell_-2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-2-002 supercell_-2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-2-000 supercell_-2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-2-001 supercell_-1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-2-002 supercell_-1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-2-000 supercell_-1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-2-001 supercell_0_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-2-002 supercell_0_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-2-000 supercell_0_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-2-001 supercell_1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-2-002 supercell_1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-2-000 supercell_1_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-2-001 supercell_2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-2-002 supercell_2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-2-000 supercell_2_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-2-001 supercell_3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-2-002 supercell_3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-1-000 supercell_3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-1-001 supercell_3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-1-002 supercell_3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_0-000 supercell_3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_0-001 supercell_3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_0-002 supercell_3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_1-000 supercell_3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_1-001 supercell_3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_1-002 supercell_3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_2-000 supercell_3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_2-001 supercell_3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_2-002 supercell_3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_3-000 supercell_3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_3-001 supercell_3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_3-002 supercell_3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_3-000 supercell_3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_3-001 supercell_2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_3-002 supercell_2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_3-000 supercell_2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_3-001 supercell_1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_3-002 supercell_1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_3-000 supercell_1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_3-001 supercell_0_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_3-002 supercell_0_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_3-000 supercell_0_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_3-001 supercell_-1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_3-002 supercell_-1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_3-000 supercell_-1_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_3-001 supercell_-2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_3-002 supercell_-2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_3-000 supercell_-2_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_3-001 supercell_-3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_3-002 supercell_-3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_2-000 supercell_-3_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_2-001 supercell_-3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_2-002 supercell_-3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_1-000 supercell_-3_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_1-001 supercell_-3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_1-002 supercell_-3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_0-000 supercell_-3_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_0-001 supercell_-3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_0-002 supercell_-3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-1-000 supercell_-3_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-1-001 supercell_-3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-1-002 supercell_-3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-2-000 supercell_-3_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-2-001 supercell_-3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-2-002 supercell_-3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-3-000 supercell_-3_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-3-001 supercell_-3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-3-002 supercell_-3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-3-000 supercell_-3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-3-001 supercell_-2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-3-002 supercell_-2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-3-000 supercell_-2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-3-001 supercell_-1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-3-002 supercell_-1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-3-000 supercell_-1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-3-001 supercell_0_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-3-002 supercell_0_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-3-000 supercell_0_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-3-001 supercell_1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-3-002 supercell_1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-3-000 supercell_1_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-3-001 supercell_2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-3-002 supercell_2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-3-000 supercell_2_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-3-001 supercell_3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-3-002 supercell_3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-3-000 supercell_3_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-3-001 supercell_4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-3-002 supercell_4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-2-000 supercell_4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-2-001 supercell_4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-2-002 supercell_4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-1-000 supercell_4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-1-001 supercell_4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-1-002 supercell_4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_0-000 supercell_4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_0-001 supercell_4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_0-002 supercell_4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_1-000 supercell_4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_1-001 supercell_4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_1-002 supercell_4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_2-000 supercell_4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_2-001 supercell_4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_2-002 supercell_4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_3-000 supercell_4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_3-001 supercell_4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_3-002 supercell_4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_4-000 supercell_4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_4-001 supercell_4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_4-002 supercell_4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_4-000 supercell_4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_4-001 supercell_3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_4-002 supercell_3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_4-000 supercell_3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_4-001 supercell_2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_4-002 supercell_2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_4-000 supercell_2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_4-001 supercell_1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_4-002 supercell_1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_4-000 supercell_1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_4-001 supercell_0_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_4-002 supercell_0_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_4-000 supercell_0_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_4-001 supercell_-1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_4-002 supercell_-1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_4-000 supercell_-1_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_4-001 supercell_-2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_4-002 supercell_-2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_4-000 supercell_-2_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_4-001 supercell_-3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_4-002 supercell_-3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_4-000 supercell_-3_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_4-001 supercell_-4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_4-002 supercell_-4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_3-000 supercell_-4_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_3-001 supercell_-4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_3-002 supercell_-4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_2-000 supercell_-4_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_2-001 supercell_-4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_2-002 supercell_-4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_1-000 supercell_-4_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_1-001 supercell_-4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_1-002 supercell_-4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_0-000 supercell_-4_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_0-001 supercell_-4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_0-002 supercell_-4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-1-000 supercell_-4_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-1-001 supercell_-4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-1-002 supercell_-4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-2-000 supercell_-4_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-2-001 supercell_-4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-2-002 supercell_-4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-3-000 supercell_-4_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-3-001 supercell_-4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-3-002 supercell_-4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-4-000 supercell_-4_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-4-001 supercell_-4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-4-002 supercell_-4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-4-000 supercell_-4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-4-001 supercell_-3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-4-002 supercell_-3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-4-000 supercell_-3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-4-001 supercell_-2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-4-002 supercell_-2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-4-000 supercell_-2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-4-001 supercell_-1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-4-002 supercell_-1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-4-000 supercell_-1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-4-001 supercell_0_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-4-002 supercell_0_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-4-000 supercell_0_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-4-001 supercell_1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-4-002 supercell_1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-4-000 supercell_1_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-4-001 supercell_2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-4-002 supercell_2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-4-000 supercell_2_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-4-001 supercell_3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-4-002 supercell_3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-4-000 supercell_3_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-4-001 supercell_4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-4-002 supercell_4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-4-000 supercell_4_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-4-001 supercell_5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-4-002 supercell_5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-3-000 supercell_5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-3-001 supercell_5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-3-002 supercell_5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-2-000 supercell_5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-2-001 supercell_5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-2-002 supercell_5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-1-000 supercell_5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-1-001 supercell_5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-1-002 supercell_5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_0-000 supercell_5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_0-001 supercell_5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_0-002 supercell_5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_1-000 supercell_5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_1-001 supercell_5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_1-002 supercell_5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_2-000 supercell_5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_2-001 supercell_5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_2-002 supercell_5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_3-000 supercell_5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_3-001 supercell_5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_3-002 supercell_5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_4-000 supercell_5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_4-001 supercell_5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_4-002 supercell_5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_5-000 supercell_5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_5-001 supercell_5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_5-002 supercell_5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_5-000 supercell_5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_5-001 supercell_4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_5-002 supercell_4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_5_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_5-000 supercell_4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_5-001 supercell_3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_5-002 supercell_3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_5-000 supercell_3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_5-001 supercell_2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_5-002 supercell_2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_5-000 supercell_2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_5-001 supercell_1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_5-002 supercell_1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_5-000 supercell_1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_5-001 supercell_0_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_5-002 supercell_0_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_5-000 supercell_0_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_5-001 supercell_-1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_5-002 supercell_-1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_5-000 supercell_-1_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_5-001 supercell_-2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_5-002 supercell_-2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_5-000 supercell_-2_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_5-001 supercell_-3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_5-002 supercell_-3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_5-000 supercell_-3_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_5-001 supercell_-4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_5-002 supercell_-4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_5-000 supercell_-4_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_5-001 supercell_-5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_5-002 supercell_-5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_4-000 supercell_-5_5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_4-001 supercell_-5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_4-002 supercell_-5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_3-000 supercell_-5_4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_3-001 supercell_-5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_3-002 supercell_-5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_2-000 supercell_-5_3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_2-001 supercell_-5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_2-002 supercell_-5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_1-000 supercell_-5_2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_1-001 supercell_-5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_1-002 supercell_-5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_0-000 supercell_-5_1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_0-001 supercell_-5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_0-002 supercell_-5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-1-000 supercell_-5_0-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_0-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-1-001 supercell_-5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_0-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-1-002 supercell_-5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_0-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-2-000 supercell_-5_-1-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-1-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-2-001 supercell_-5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-1-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-2-002 supercell_-5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-1-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-3-000 supercell_-5_-2-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-2-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-3-001 supercell_-5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-2-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-3-002 supercell_-5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-2-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-4-000 supercell_-5_-3-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-3-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-4-001 supercell_-5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-3-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-4-002 supercell_-5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-3-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-5-000 supercell_-5_-4-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-4-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-5-001 supercell_-5_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-4-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-5_-5-002 supercell_-5_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-4-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-5-000 supercell_-5_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-5-001 supercell_-4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-4_-5-002 supercell_-4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-5_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-5-000 supercell_-4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-5-001 supercell_-3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-3_-5-002 supercell_-3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-4_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-5-000 supercell_-3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-5-001 supercell_-2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-2_-5-002 supercell_-2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-3_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-5-000 supercell_-2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-5-001 supercell_-1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_-1_-5-002 supercell_-1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-2_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-5-000 supercell_-1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-5-001 supercell_0_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_0_-5-002 supercell_0_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_-1_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-5-000 supercell_0_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-5-001 supercell_1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_1_-5-002 supercell_1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_0_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-5-000 supercell_1_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-5-001 supercell_2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_2_-5-002 supercell_2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_1_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-5-000 supercell_2_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-5-001 supercell_3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_3_-5-002 supercell_3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_2_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-5-000 supercell_3_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-5-001 supercell_4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_4_-5-002 supercell_4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_3_-5-002 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-5-000 supercell_4_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-5-000 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-5-001 supercell_5_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-5-001 &
wait

timeout 2876m /work/e05/e05/bcamino/runCRYSTAL/Pcry_slurm supercell_5_-5-002 supercell_5_-5-000 &
wait
/work/e05/e05/bcamino/runCRYSTAL/post_proc_slurm crys supercell_4_-5-002 &
wait

