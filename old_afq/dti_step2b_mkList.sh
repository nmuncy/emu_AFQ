#!/bin/bash

### --- Notes
#
# Assumes step2a has been run locally,
# and Group_data.txt exists.

# set paths, reference
codeDir=~/compute/emu_diff
derivDir=/scratch/madlab/emu_diff/derivatives
refFile=${codeDir}/Group_data.txt

# set out files
print1=${codeDir}/Group_path.txt
print2=${codeDir}/Group_mem.txt

# check for ref file
if [ -z $refFile ]; then
    echo "Run step2a first" >&2
    exit 1
fi

# get subj, pars7 columns
subjList=(`tail -n +2 $refFile | awk '{print $1}'`)
anxList=(`tail -n +2 $refFile | awk '{print $5}'`)

# zero output
unset grpMem
>$print1

for ind in ${!subjList[@]}; do

    # get array values
    subjPath=${derivDir}/sub-${subjList[$ind]}/ses-S2/dti96trilin
    grp=${anxList[$ind]}

    if [ -f ${subjPath}/dt6.mat ]; then

        # write path
        # echo "[var, '${subjPath}/'], ..." >> $print1
        echo "'${subjPath}/', ..." >> $print1

        # write group membership (>3 == anxious)
        if [ $grp -gt 3 ]; then
            grpMem+="0, "
        elif [ $grp -le 3 ]; then
            grpMem+="1, "
        fi
    fi
done
echo $grpMem > $print2
