#!/bin/bash

# Notes:
#
# This script will copy and rename DWI data from a BIDS
#   pre-process derivative directory to a new location.
#
# Then it will submit afq_step1_setup.py

# Location of afq_step1_setup.py
code_dir=~/compute/emu_AFQ

# Paths for pre-processed DWI data, uses BIDS directory organization.
#   emu_dir = path to project parent directory
emu_dir=/home/data/madlab/McMakin_EMUR01
deriv_dir=${emu_dir}/derivatives/dwi_preproc
dset_dir=${emu_dir}/dset
ref_file=${dset_dir}/dataset_description.json

# Variables for AFQ:
#   parent_dir = parent output directory
#   work_dir = output derivative directory
#   sess = reference string
#   run = reference string
parent_dir=/scratch/madlab/emu_AFQ
work_dir=${parent_dir}/derivatives/dwi_preproc
sess=ses-S2
run=run-1

# Get jsons
mkdir -p $work_dir
cp $ref_file $parent_dir
cp $ref_file $work_dir

# Copy, BIDs format pre-processed dwi data
cd $deriv_dir

for subj in sub-*; do

    source_dir=${deriv_dir}/${subj}/${sess}/dwi
    out_dir=${work_dir}/${subj}/${sess}/dwi

    if [ ! -f ${out_dir}/${subj}_${sess}_dwi.nii.gz ]; then

        mkdir -p $out_dir

        cp ${dset_dir}/${subj}/${sess}/dwi/${subj}_${sess}_${run}_dwi.bval \
            ${out_dir}/${subj}_${sess}_dwi.bval
        cp ${source_dir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.bvec \
            ${out_dir}/${subj}_${sess}_dwi.bvec
        cp ${source_dir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.nii.gz \
            ${out_dir}/${subj}_${sess}_dwi.nii.gz
    fi
done

# submit python
python ${code_dir}/afq_step1_setup.py $code_dir $parent_dir $work_dir
