# %%
import os
import toml
import json
import sys


# %%
def main():

    # receive args
    code_dir = str(sys.argv[1])
    bids_dir = str(sys.argv[2])
    deriv_dir = str(sys.argv[3])

    # append data_description with new field
    json_file = os.path.join(deriv_dir, "dataset_description.json")
    new_field = {"PipelineDescription": {"Name": "dwi_preproc"}}

    with open(json_file) as jf:
        json_content = json.load(jf)

    if "PipelineDescription" not in json_content.keys():
        json_content.update(new_field)

    with open(json_file, "w") as jf:
        json.dump(json_content, jf)

    # write config.toml
    toml_file = os.path.join(code_dir, "config.toml")

    conf_dict = {
        "files": {"dmriprep_folder": deriv_dir},
        "BIDS": {"bids_path": bids_dir},
        "TRACTOGRAPHY": {"directions": "prob"},
    }

    with open(toml_file, "w") as tf:
        toml.dump(conf_dict, tf)


if __name__ == "__main__":
    main()
