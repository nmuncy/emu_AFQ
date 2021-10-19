#!/bin/bash

#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=9   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=24gb   # memory per CPU core
#SBATCH -J "AFQ"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab

config_file=$1
pyAFQ $config_file --notrack