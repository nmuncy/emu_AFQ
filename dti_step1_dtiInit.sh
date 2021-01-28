#!/bin/bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=24gb   # memory per CPU core
#SBATCH -J "dcm2nii"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab


subj=$1
code_dir=$2

cd $code_dir
module load matlab-2017b
matlab -nodisplay -nojvm -nosplash -r "dti_step1_subjId('$subj')"