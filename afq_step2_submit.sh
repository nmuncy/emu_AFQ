#!/bin/bash

function Usage {
	cat << USAGE

	Wraps afq_step2_CLI.sh

	This script will check for files required by AFQ and then
	sbatch submit afq_step2_CLI.sh. Sbatch stdout/err will be 
	captured in a time-stamped directory located at
	<proc_dir>/derivatives/Slurm_out.

	Usage:

		afq_step2_submit.sh \
			-c config.toml \
			-p /path/to/proc_dir

	Required Arguments:

		-c config.toml = AFQ config.toml file 
			same file as from afq_step1_setup.py -c <config-file>
		-p </path/to/proc_dir> = location of BIDS-structured project directory
			same dir as afq_step1_submit.sh -p <proc_dir>

	Optional Arguments:

		-t </path/to/AFQ_data> = location of AFQ_data template directory
			useful for HPC environments where AFQ API cannot download needed files

	Example Usage:

		afq_step2_submit.sh \
			-c ~/compute/emu_AFQ/config.toml \
			-p /scratch/madlab/emu_AFQ \
			-t /home/data/madlab/atlases/AFQ_data
USAGE
}


unset config_file proj_dir afq_data
while getopts ":c:p:t:h" OPT
	do
	case $OPT in 
		c) config_file=${OPTARG}
			;;
		p) proj_dir=${OPTARG}
			;;
		t) afq_data=${OPTARG}
			;;
		h)
			Usage
			exit 0
		\?) echo -e "\n \t ERROR: invalid option." >&2
			exit 1
			;;
	esac
done

# check inputs
if [ $OPTIND == 1 ]; then
    Usage
    exit 0
fi

if [ -z "$config_file" ] || [ -z "$proj_dir" ]; then
	echo -e "\n \t ERROR: required args not specified." >&2
	Usage
	exit 1
fi

if [ ! -d $proj_dir ]; then
	echo -e "\n \t ERROR: $proj_dir not detected or is not a directory." >&2
	Usage
	exit 1
fi

# work
if [ ! -d ${HOME}/AFQ_data ]; then
	if [ -z "$afq_data" ]; then
		echo -e "\n \t ERROR: AFQ_data not detected in \$HOME and \$afq_data is unset." >&2
		Usage
		exit 1
	else
		cp -r $afq_data $HOME
	fi
fi

time=`date '+%Y_%m_%d-%H_%M'`
slurm_dir=${proj_dir}/derivatives/Slurm_out/afq2_$time
mkdir -p $slurm_dir

sbatch \
	-e ${slurm_dir}/err.txt \
	-o ${slurm_dir}/out.txt \
	afq_step2_CLI.sh $config_file