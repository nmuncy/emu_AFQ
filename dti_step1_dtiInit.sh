#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=3   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=24gb   # memory per CPU core
#SBATCH -J "DTI1"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab


subj=$1
code_dir=$2

cd $code_dir
module load matlab-2017b
module load spm-12
# matlab -nodisplay -nojvm -nosplash -r "dti_step1_subjId('$subj')"
matlab -nodisplay -nosplash -r "dti_step1_subjId('$subj')"