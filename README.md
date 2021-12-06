# emu_AFQ

Repo containing code, documents, and data used in the manuscript "General Additive Models Address Statistical Issues in Diffusion MRI: An Example with Clinically Anxious Adolescents".

# Code

Various scripts for organizing DWI data, running pyAFQ, and conducting statistical analyses. Scripts are named (a) for their respective order in the workflow, and (b) with a description of their function.

- step1_submit.sh : Set up, copy data, and wrap step1_setup.py
- step1_setup.py : Configure config.toml for pyAFQ CLI
- step2_submit.sh : Wrap step2_CLI.sh
- step2_CLI.sh : Run AFQ via pyAFQ on Slurm scheduled resources
- step3_manuscript_stats.R : Code used to conduct all statistical analyses reported in the manuscript. Unused analyses, portions remain for potential future analyses and transparency
- step3_reduced_stats.R : Illustrative code to aid in the implementation of using a GAM method with AFQ output
- step4_quick_stats.R : Quick queries for the manuscript
- step5_reviewer_stats.R : Code written to address reviewer concerns, in addition to updates to step3_manuscript_stats.R and step4_quick_stats.R.


# Data

Contains data used in the manuscript.

- tract_profiles.csv : pyAFQ output
- Master_dataframe_G3.csv : pyAFQ output with demographic information added for groupings according to PARS-6 tertiles


# Docs

Documents used to conduct analyses or to aid in replication.

- config.toml : Configuration file used with pyAFQ
- R_session_info.txt : A complete description of all packages and versions used in the analyses
