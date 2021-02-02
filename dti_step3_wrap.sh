#!/bin/bash

workDir=/scratch/madlab/emu_diff/derivatives

anDir=${workDir}/Analyses
a1Dir=${anDir}/AFQ
a2Dir=${anDir}/AFQ-CC

mkdir -p $anDir $a1Dir $a2Dir

module load matlab-2017b
# matlab -nodisplay -nojvm -nosplash -r dti_afq_parameters_2
matlab -nodisplay -nosplash -r dti_step3_params
# module unload matlab
