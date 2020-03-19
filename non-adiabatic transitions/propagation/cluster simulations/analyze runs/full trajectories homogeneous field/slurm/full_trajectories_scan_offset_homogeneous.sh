#!/bin/bash
#SBATCH --requeue
#SBATCH --partition       scavenge
#SBATCH --job-name        full_trajectories_scan_offset_homogeneous
#SBATCH --ntasks          2000
#SBATCH --cpus-per-task   1
#SBATCH --mem-per-cpu     500M
#SBATCH --time            24:00:00
#SBATCH --mail-type       ALL
#SBATCH --mail-user       olivier.grasdijk@yale.edu

prog=/home/fas/demille/jog9/project/State-evolution/src/main.py
run_dir=/home/fas/demille/jog9/project/runs/
options_file=full_trajectories_scan_offset_homogeneous.json

#module load miniconda
#source activate tf_gpu
module load Python/3.6.4-foss-2018a

mpi_params="--mca mpi_warn_on_fork 0"

#mpirun -n 1 python3 -m cProfile -s tottime $prog $run_dir $options_file
#python3 $prog $run_dir $options_file
mpirun -n 2000 $mpi_params python3 $prog $run_dir $options_file
