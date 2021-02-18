%%%%%%%
function dti_step1_subjId(x)

% Display participant ID:
display(x);

% Get home directory:
var = getenv('HOME');

% Add modules to MATLAB
SPM8Path = '/home/applications/spm12'
addpath(genpath(SPM8Path));
vistaPath = [var, '/matlab/vistasoft'];
addpath(genpath(vistaPath));
AFQPath = [var, '/matlab/AFQ'];
addpath(genpath(AFQPath));

% Set file names:
subjDir = ['/scratch/madlab/emu_diff/derivatives/vistasoft/', x, '/ses-S2'];
dtiFile = [subjDir, '/raw/dwi.nii.gz'];
cd (subjDir);

% Do it:
dwParams = dtiInitParams('rotateBvecsWithCanXform', 1, 'phaseEncodeDir', 1, 'clobber', 1);
dtiInit(dtiFile, 't1.nii.gz', dwParams);

exit;
