#!/bin/bash
##echo "Bash version ${BASH_VERSION}..."

matlab_script=$1
walltime=$2
model_number=$3
subject_name=$4

id=`echo "/home/radachi/research/CD/analysis/beh_model/cd_modelfit_tolmanJob_one.sh $matlab_script $model_number $subject_name" | qsub -l walltime=$walltime -m ae -M radachi@caltech.edu -q batch -l pmem=1gb -p -1024 -d /home/radachi/`

##echo $subject_name
sleep .1
