#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=9   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=24gb   # memory per CPU core
#SBATCH -J "AFQ"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab

# make sure that AFQ_data exists in $HOME
if [ ! -d ${HOME}/AFQ_data ]; then
    cp -r /home/data/madlab/atlases/AFQ_data $HOME
fi

# get my env, run job
# source ${HOME}/.bashrc
# conda activate pyAFQ
pyAFQ config.toml --notrack