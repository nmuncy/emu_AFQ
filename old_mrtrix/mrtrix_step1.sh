#!/bin/bash

### --- Notes:
#
# Dependencies: MRtrix3, FSL, ANTs
#	FSL requires python 2.7


projDir=~/Projects/emu_diff/emu_data
dsetDir=${projDir}/dset
derivDir=${projDir}/derivatives/mrtrix

sess=ses-S2
subj=sub-4000

subjData=${dsetDir}/${subj}/${sess}
subjOut=${derivDir}/${subj}/${sess}

if [ ! -d $subjOut ]; then
	mkdir -p $subjOut
fi
cd $subjOut

# convert to mrtrix format
if [ ! -f dwi.mif ]; then
	mrconvert ${subjData}/dwi/${subj}_${sess}_run-1_dwi.nii.gz \
		dwi.mif \
		-fslgrad \
		${subjData}/dwi/${subj}_${sess}_run-1_dwi.bvec \
		${subjData}/dwi/${subj}_${sess}_run-1_dwi.bval
fi

# get bvec/als
for suff in bv{ec,al}; do
	if [ ! -f dwi.$suff ]; then
		cp ${subjData}/dwi/${subj}_${sess}_run-1_dwi.$suff dwi.$suff
	fi
done

# check num vols/vecs/vals
numVol=`mrinfo -size dwi.mif | awk '{print $4}'`
numVec=`awk '{print NF; exit}' dwi.bvec`
numVal=`awk '{print NF; exit}' dwi.bval`
if [ $numVol != $numVec ] || [ $numVol != $numVec ]; then
	echo "Issue with number vols/vecs/vals check" >&2
	exit 1
fi

# get fmaps, make average
for dir in AP PA; do
	if [ ! -f fmap_${dir}.mif ]; then
		cp ${subjData}/fmap/${subj}_${sess}_acq-dwi_dir-${dir}_run-?_epi.nii.gz fmap_${dir}.nii.gz
		mrconvert fmap_${dir}.nii.gz fmap_${dir}.mif
	fi
done

if [ ! -f fmap_both.mif ]; then
	mrcat fmap_AP.mif fmap_PA.mif -axis 3 fmap_both.mif
fi

# denoise
if [ ! -f dwi_den.mif ]; then
	dwidenoise dwi.mif dwi_den.mif -noise tmp_noise.mif
	mrcalc dwi.mif dwi_den.mif -subtract tmp_residual.mif
fi

# unwarp data, remove eddy currents
#	takes about 6 hours
date > time.txt

if [ ! -f dwi_preproc.mif ]; then
	dwifslpreproc dwi_den.mif dwi_preproc.mif \
		-nocleanup \
		-pe_dir PA \
		-rpe_pair -se_epi fmap_both.mif \
		-eddy_options " --slm=linear --data_is_shelled"
fi

date >> time.txt

# dwibiascorrect ants dwi_preproc.mif dwi_unbiased.mif -bias bias.mif