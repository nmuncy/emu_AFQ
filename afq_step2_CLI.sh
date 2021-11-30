#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=24gb
#SBATCH --job-name=AFQ
#SBATCH -p IB_44C_512G
#SBATCH --account iacc_madlab
#SBATCH --qos pq_madlab


function Usage {
	cat << USAGE

	Run AFQ via pyAFQ CLI on scheduled resources.

	Required Arguments:
		arg[1] = /path/to/config.toml

	Example Usage:
		sbatch \\
            -e err.txt -o out.txt \\
            afq_step2_CLI.sh /path/to/config.toml
USAGE
}

# check usage
if [[ $# -ne 1 ]]; then
    Usage
    exit 1
fi

config_file=$1
pyAFQ $config_file --notrack