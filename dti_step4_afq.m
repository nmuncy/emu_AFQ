% Get home directory:
var = getenv('HOME');

% Add modules to MATLAB. Do not change the order of these programs:
% SPM8Path = [var, '/matlab/spm8'];
SPM8Path = '/home/applications/spm12'
addpath(genpath(SPM8Path));
vistaPath = [var, '/matlab/vistasoft'];
addpath(genpath(vistaPath));
AFQPath = [var, '/matlab/AFQ'];
addpath(genpath(AFQPath));

% get mat files
load /scratch/madlab/emu_diff/derivatives/Analyses/AFQ/sub_dirs.mat
load /scratch/madlab/emu_diff/derivatives/Analyses/AFQ/sub_group.mat

% set out vars
% outdir = fullfile([var, '/compute/TBI/Analyses/AFQ/']);
outdir = fullfile('/scratch/madlab/emu_diff/derivatives/Analyses/AFQ/']);
outname = fullfile(outdir, ['afq_TBI_MS']);

% work
afq = AFQ_Create('sub_dirs', sub_dirs, 'sub_group', sub_group, 'showfigs', false);
[afq, patient_data, control_data, norms, abn, abnTracts] = AFQ_run(sub_dirs, sub_group, afq);
save(outname, 'afq');
