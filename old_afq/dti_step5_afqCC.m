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

load /scratch/madlab/emu_diff/derivatives/Analyses/AFQ/afq_TBI_MS.mat
outdir = fullfile('/scratch/madlab/emu_diff/derivatives/Analyses/AFQ-CC/');
outname = fullfile(outdir, ['afq_cc_TBI_MS']);
afq = AFQ_SegmentCallosum(afq, 0)
save(outname, 'afq');
