#!/bin/bash

# set orienting vars
code_dir=~/compute/emu_diff
parent_dir=/scratch/madlab/emu_diff
deriv_dir=${parent_dir}/derivatives
dset_dir=${parent_dir}/dset

# set output vars
slurm_dir=${deriv_dir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
out_dir=${slurm_dir}/dti1_${time}
mkdir -p $out_dir

# copy data, run job for e/subject
cd $dset_dir
for subj in sub*; do

	subj_dir=${deriv_dir}/${subj}/ses-S2
	data_dir=${dset_dir}/${subj}/ses-S2/dwi

	# copy data to derivatives
	if [ ! -f ${subj_dir}/dwi.nii.gz ]; then
		mkdir -p $subj_dir
		cp ${data_dir}/${subj}_ses-S2_dwi.nii.gz ${subj_dir}/dwi.nii.gz
		cp ${data_dir}/${subj}_ses-S2_dwi.bvec ${subj_dir}/dwi.bvec
		cp ${data_dir}/${subj}_ses-S2_dwi.bval ${subj_dir}/dwi.bval
		cp ${data_dir}/${subj}_ses-S2_dwi.json ${subj_dir}/dwi.json
	fi

	# submit job
	if [ ! -f ${subj_dir}/dti96trilin/dt6.mat ]; then
		sbatch \
		-o ${out_dir}/output_${subj}.txt \
		-e ${out_dir}/error_${subj}.txt \
		${code_dir}/dti_step1_dtiInit.sh $subj $code_dir
		sleep 1
	fi
done
