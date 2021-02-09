#!/bin/bash

# set orienting vars
codeDir=~/compute/emu_diff
parentDir=/scratch/madlab/emu_diff
derivDir=${parentDir}/derivatives
dsetDir=${parentDir}/dset

# set output vars
slurmDir=${derivDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/dti1_${time}
mkdir -p $outDir

# copy data, run job for e/subject
# cd $dsetDir
# for subj in sub*; do
subj=sub-4000

	# set subj vars
	subjDir=${derivDir}/vistasoft/${subj}/ses-S2
	subjRaw=${subjDir}/raw
	dwiDir=${dsetDir}/${subj}/ses-S2/dwi
	t1Dir=${dsetDir}/${subj}/ses-S1/anat
	mkdir -p $subjRaw

	# set up subj vistasoft with t1, dwi
	if [ ! -f ${subjRaw}/dwi.nii.gz ]; then
		cp ${t1Dir}/${subj}_ses-S1_T1w.nii.gz ${subjDir}/t1.nii.gz
		for suff in bvec bval json nii.gz; do
			cp ${dwiDir}/${subj}_ses-S2_dwi.$suff ${subjRaw}/dwi.$suff
		done
	fi

	# submit job
	if [ ! -f ${subjDir}/dti96trilin/dt6.mat ]; then
		sbatch \
		-o ${outDir}/output_${subj}.txt \
		-e ${outDir}/error_${subj}.txt \
		${codeDir}/dti_step1_dtiInit.sh $subj $codeDir
		sleep 1
	fi
# done
