#!/bin/bash
#SBATCH --requeue
#SBATCH --partition       scavenge
#SBATCH --job-name        straight_trajectory_scan_homogeneous_field_B_field
#SBATCH --ntasks          500
#SBATCH --cpus-per-task   1
#SBATCH --mem-per-cpu     1G
#SBATCH --time            24:00:00
#SBATCH --mail-type       ALL
#SBATCH --mail-user       olivier.grasdijk@yale.edu

mpi_params="--mca mpi_warn_on_fork 0"
prog=/home/fas/demille/jog9/project/State-evolution/src/main.py
run_dir=/home/demille/jog9/project/runs/straight-trajectories-homogeneous-field
options_file=straight_trajectory_scan_homogeneous_field_B_field.json

#module load miniconda
#source activate tf_gpu
module load Python/3.6.4-foss-2018a

#mpirun -n 1 python3 -m cProfile -s tottime $prog $run_dir $options_file
#python3 $prog $run_dir $options_file
mpirun -n 500 $mpi_params python3 $prog --info $run_dir $options_file
