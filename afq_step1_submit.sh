#!/bin/bash

# set up
codeDir=~/compute/emu_diff

# vars for pre-processed data
emuDir=/home/data/madlab/McMakin_EMUR01
derivDir=${emuDir}/derivatives/dwi_preproc
dsetDir=${emuDir}/dset
refFile=${dsetDir}/dataset_description.json

# vars for AFQ work
parDir=/scratch/madlab/emu_diff
workDir=${parDir}/derivatives/dwi_preproc
sess=ses-S2
run=run-1

# get jsons
mkdir -p $workDir
cp $refFile $parDir
cp $refFile $workDir

# BIDs format pre-processed dwi data
# cd $derivDir
# for subj in sub-*; do
subj=sub-4001

    sourceDir=${derivDir}/${subj}/${sess}/dwi
    outDir=${workDir}/${subj}/${sess}/dwi

    if [ ! -f ${outDir}/${subj}_${sess}_dwi.nii.gz ]; then
        mkdir -p $outDir
        cp ${dsetDir}/${subj}/${sess}/dwi/${subj}_${sess}_${run}_dwi.bval ${outDir}/${subj}_${sess}_dwi.bval
        cp ${sourceDir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.bvec ${outDir}/${subj}_${sess}_dwi.bvec
        cp ${sourceDir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.nii.gz ${outDir}/${subj}_${sess}_dwi.nii.gz
    fi
# done

# submit python
python ${codeDir}/afq_step1_setup.py $codeDir $parDir $workDir
