#!/bin/bash


codeDir=~/compute/emu_diff
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=/scratch/madlab/emu_diff/derivatives/Slurm_out/dti5_afqCC_${time}

mkdir -p $outDir

sbatch \
-o ${outDir}/output_dti5.txt \
-e ${outDir}/error_dti5.txt \
${codeDir}/dti_step5_job.sh
