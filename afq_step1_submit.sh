#!/bin/bash

# set up
codeDir=~/compute/emu_diff
refFile=/home/data/madlab/McMakin_EMUR01/dset/dataset_description.json
parDir=/scratch/madlab/emu_diff
derivDir=${parDir}/derivatives/vistasoft
sess=ses-S2

# get jsons
cp $refFile $parDir
cp $refFile $derivDir

# BIDs format pre-processed dwi data
cd $derivDir

for subj in sub-*; do

    sourceDir=${derivDir}/${subj}/${sess}
    outDir=${sourceDir}/dwi

    if [ ! -f ${outDir}/${subj}_${sess}_dwi.nii.gz ]; then
        mkdir -p $outDir
        cp ${sourceDir}/dwi_aligned_trilin.bvecs ${outDir}/${subj}_${sess}_dwi.bvec
        cp ${sourceDir}/dwi_aligned_trilin.bvals ${outDir}/${subj}_${sess}_dwi.bval
        cp ${sourceDir}/dwi_aligned_trilin.nii.gz ${outDir}/${subj}_${sess}_dwi.nii.gz
    fi
done

# submit python
python ${codeDir}/dti_step2_setup.py $codeDir $parDir $derivDir
