#!/bin/bash
##echo "Bash version ${BASH_VERSION}..."
#PBS -q batch
#PBS -l nodes=1
#PBS -l walltime=48:00:00
#PBS -m ae 
#PBS -M ryoatsat@gmail.com

MATLAB=/usr/local/matlab/R2014a/bin/matlab

NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

matlab_script=$1
model_number=$2
subject_name=$3

#PBS -N cluster_analysis_$1
mpirun -np $NPROCS $MATLAB -nojvm -nodisplay -nosplash -c /home/radachi/MATLAB/matlab.lic -r "${matlab_script}(${model_number},'${subject_name}'); exit"
