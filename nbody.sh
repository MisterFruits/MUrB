#!/bin/bash 
#SBATCH --job-name=murb
#SBATCH --nodes=1
#SBATCH --time=0:20:0
#SBATCH --partition=standard
#SBATCH --qos=standard




for OMP_PROC_BIND in close spread; do
    echo "CPU bind = $OMP_NUM_THREADS"
    export OMP_NUM_THREADS
for OMP_NUM_THREADS in 16 32 64 128; do
    echo "$OMP_NUM_THREADS thread"
    export OMP_NUM_THREADS
    srun -n 1 build/bin/Debug/murb
done
done

export OMP_NUM_THREADS=1
for NNODES in 16 32 64 128; do
    echo "$NNODES tasks"
    srun -n $NNODES build/bin/Debug/murb
done

