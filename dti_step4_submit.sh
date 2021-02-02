#!/bin/bash


codeDir=~/compute/emu_diff
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=/scratch/madlab/emu_diff/derivatives/Slurm_out/dti4_${time}
mkdir -p $outDir

sbatch \
-o ${outDir}/output_dti4.txt \
-e ${outDir}/error_dti4.txt \
${codeDir}/dti_step4_job.sh
