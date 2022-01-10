#!/bin/bash

function Usage {
    cat << USAGE

    Copy data to a working directory, wrap step1_setup.py.

    This script will copy data from a BIDS-structured project
    directory to a BIDS-structured working/processing directory.
    This is because we store our project data in one location
    (/home/data/madlab/...) but process our data and manage
    intermediates in another (/scratch/madlab/...).

    Pre-processed dwi data will be copied from
        <data_dir>/derivatives/<deriv_dir> to
        <proc_dir>/derivatives/<deriv_dir>.

    Finally afq_step1_setup.py will be submitted with required
    args.

    Required Arguments:
        -c </path/to/code_dir> = location of afq_step1_setup.py
            Note: for lazy wrapping python, have config.toml in the code_dir
        -d </path/to/data_dir> = location of BIDS-structured stored project
            data
        -p </path/to/proc_dir> = location to where data will be copied,
            and processed. Will be created if it does not exist
        -s <ses-str> = BIDS session string, for data structure and naming
        -r <run-str> = BIDS run string, for naming
        -x <deriv_dir> = directory within derivatives containing pre-processed
            diffusion weighted data
        -j <json-file> = BIDS dataset_description.json sidecar

    Example Usage:
        ./step1_submit.sh \\
            -c /home/nmuncy/compute/emu_AFQ/code \\
            -d /home/data/madlab/McMakin_EMUR01 \\
            -p /scratch/madlab/emu_AFQ \\
            -s ses-S2 \\
            -r run-1 \\
            -x dwi_preproc \\
            -j /home/data/madlab/McMakin_EMUR01/dset/dataset_description.json

    Notes:
        In our file structure, bval exists in dset, while
        bvec and nii exist in derivatives/<deriv_dir>.
USAGE
}


# capture arguments
while getopts ":c:d:p:s:r:x:j:h" OPT
    do
    case $OPT in
        c) code_dir=${OPTARG}
            ;;
        d) data_dir=${OPTARG}
            ;;
        p) proc_dir=${OPTARG}
            ;;
        s) sess=${OPTARG}
            ;;
        r) run=${OPTARG}
            ;;
        x) diff_dir=${OPTARG}
            ;;
        j) json_file=${OPTARG}
            ;;
        h)
            Usage
            exit 0
            ;;
        \?) echo -e "\n \t ERROR: invalid option." >&2
            Usage
            exit 1
            ;;
    esac
done


# print help if no arg
if [ $OPTIND == 1 ]; then
    Usage
    exit 0
fi


# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
        code_dir) h_ret="-c"
            ;;
        data_dir) h_ret="-d"
            ;;
        proc_dir) h_ret="-p"
            ;;
        sess) h_ret="-s"
            ;;
        run) h_ret="-r"
            ;;
        diff_dir) h_ret="-x"
            ;;
        json_file) h_ret="-j"
            ;;
        *)
            echo -n "Unknown option."
            ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in code_dir data_dir proc_dir sess run diff_dir json_file; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done


# check required files exist
if [ ! -f ${code_dir}/afq_step1_setup.py ]; then
    echo -e "\n \t ERROR: afq_step1_setup.py not detected in $code_dir." >&2
    Usage
    exit 1
fi

if [ ! -f $json_file ]; then
    echo -e "\n \t ERROR: $json_file not found." >&2
    Usage
    exit 1
fi


# check required directories exist
if [ ! -d $data_dir ]; then
    echo -e "\n \t ERROR: $data_dir not found or is not a directory." >&2
    Usage
    exit 1
fi

if [ ! -d ${data_dir}/derivatives/$diff_dir ]; then
    echo -e "\n \t ERROR: $diff_dir not found or is not a directory." >&2
    Usage
    exit 1
fi


# print report
cat << EOF

    Success! Checks passed, starting work with the following variables:
        code_dir=$code_dir
        data_dir=$data_dir
        proc_dir=$proc_dir
        sess=$sess
        run=$run
        diff_dir=$diff_dir
        json_file=$json_file

EOF


# Get oriented
deriv_dir=${data_dir}/derivatives/$diff_dir
dset_dir=${data_dir}/dset
work_dir=${proc_dir}/derivatives/$diff_dir

# Get json
mkdir -p $work_dir
cp $json_file $proc_dir

# Copy, BIDs format pre-processed dwi data
unset subj_list
subj_list=(`ls $deriv_dir | grep "sub-*"`)
for subj in ${subj_list[@]}; do

    source_dir=${deriv_dir}/${subj}/${sess}/dwi
    out_dir=${work_dir}/${subj}/${sess}/dwi

    if [ ! -d $out_dir ] || [ ! -f ${out_dir}/${subj}_${sess}_dwi.nii.gz ]; then

        echo -e "\t Copying data for $subj ..."
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
echo -e "\n \t Starting python script ..."
python ${code_dir}/step1_setup.py \
    -c ${code_dir}/config.toml \
    -b $proc_dir \
    -d $work_dir \
    -j ${proc_dir}/dataset_description.json \
    -p $diff_dir
