"""Setup config.toml file for AFQ.

This script will edit the BIDs file dataset_description.json
with values needed by AFQ. It will then edit the AFQ
configuration file config.toml with project-specific values.

It is wrapped by afq_step1_submit.sh.

Examples
--------
python afq_step1_setup.py \\
    -c config.toml \\
    -b /scratch/madlab/emu_AFQ \\
    -d /scratch/madlab/emu_AFQ/derivatives/dwi_preproc \\
    -j /scratch/madlab/emu_AFQ/dataset_description.json \\
    -p dwi_preproc
"""

import os
import sys
import toml
import json
import textwrap
from argparse import ArgumentParser, RawTextHelpFormatter


def get_args():
    """Get and parse arguments."""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-p",
        "--preproc-dwi",
        type=str,
        default="dwi_preproc",
        help=textwrap.dedent(
            """\
            Description str for pre-processed DWI
            (default : %(default)s)
            """
        ),
    )
    required_args = parser.add_argument_group("Required Arguments")
    required_args.add_argument(
        "-c",
        "--config-file",
        help="/path/to/afq/config.toml",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-b",
        "--bids-directory",
        help="/path/to/bids/project_directory",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-d",
        "--deriv-directory",
        help="/path/to/bids/project_directory/derivatives/pre-processed_dwi",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-j",
        "--json-file",
        help="/path/to/bids/project_directory/dataset_description.json",
        type=str,
        required=True,
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser


def write_jsons(json_file, deriv_json, dwi_str):
    """Update and write a number of json files.

    Three dataset_description.json files will be written/updated: the
    input json_file, a copy of json_file written to deriv_json. These
    json files will have the "PipelineDescription" field updated/appended
    for AFQ needs.

    Parameters
    ----------
    json_file : str
        The file location of input dataset_description.json
    dset_json : str
        The location of dataset_description.json in bids dset, will be
        created/overwritten with json_file values
    deriv_json : str
        The location of dataset_description.json in bids derivatives/dwi,
        will be created/overwritten with json_file values
    dwi_str : str
        String to append as PipelineDescription Name for dwi data
    """
    with open(json_file) as jf:
        json_content = json.load(jf)

    json_content.update({"PipelineDescription": {"Name": "General"}})
    with open(json_file, "w") as jf:
        json.dump(json_content, jf)

    json_content.update({"PipelineDescription": {"Name": f"{dwi_str}"}})
    with open(deriv_json, "w") as jf:
        json.dump(json_content, jf)

    jf.close()


def main():
    """Configure json and toml files."""
    # receive arguments
    args = get_args().parse_args()
    toml_file = args.config_file
    json_file = args.json_file
    bids_dir = args.bids_directory
    deriv_dir = args.deriv_directory
    dwi_str = args.preproc_dwi

    # append and make dataset_description.json files
    deriv_json = os.path.join(deriv_dir, "dataset_description.json")
    write_jsons(json_file, deriv_json, dwi_str)

    # edit config.toml with desired fields
    toml_dict = toml.load(toml_file)
    toml_dict["TRACTOGRAPHY"]["directions"] = "prob"
    toml_dict["BIDS"]["bids_path"] = bids_dir
    toml_dict["files"]["dmriprep_folder"] = deriv_dir
    tf = open(toml_file, "w")
    toml.dump(toml_dict, tf)
    tf.close()


if __name__ == "__main__":
    main()
