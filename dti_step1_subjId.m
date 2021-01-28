%%%%%%%
function dti_step1_subjId(x)

% Display participant ID:
display(x);

% Get home directory:
var = getenv('HOME');

% Add modules to MATLAB. ORDER IS IMPORTANT! Do not change the order of these programs:
SPM8Path = [var, '/matlab/spm8'];
addpath(genpath(SPM8Path));
vistaPath = [var, '/matlab/vistasoft'];
addpath(genpath(vistaPath));
AFQPath = [var, '/matlab/AFQ'];
addpath(genpath(AFQPath));

% Set file names:
subjDir = ['/scratch/madlab/emu_diff/derivatives/', x, '/ses-S2'];
dtiFile = [subjDir, '/dwi.nii.gz'];
cd (subjDir);

% Do it:
ni = readFileNifti(dtiFile);
ni = niftiSetQto(ni, ni.sto_xyz);
writeFileNifti(ni, dtiFile);

% Determine phase encode dir:
% > info=dicominfo([var,'/compute/TBI/Ref_dir/IM-0007-0001.dcm']);
% To get the manufacturer information:
% > info.(dicomlookup('0008','0070'))
% To get the axis of phase encoding with respect to the image:
% > info.(dicomlookup('0018','1312'))
% If phase encode dir is 'COL', then set 'phaseEncodeDir' to '2'
% If phase encode dir is 'ROW', then set 'phaseEncodeDir' to '1'

% For Siemens / Philips specific code we need to add 'rotateBvecsWithCanXform',
% AND ALWAYS DOUBLE CHECK phaseEncodeDir:
% > dwParams = dtiInitParams('rotateBvecsWithCanXform',1,'phaseEncodeDir',2,'clobber',1);
dwParams = dtiInitParams('rotateBvecsWithCanXform', 1, 'phaseEncodeDir', 1, 'clobber', 1);

% For GE specific code,
% AND ALWAYS DOUBLE CHECK phaseEncodeDir:
% > dwParams = dtiInitParams('phaseEncodeDir',2,'clobber',1);

% Here's the one line of code to do the DTI preprocessing:
dtiInit(dtiFile, 'MNI', dwParams);

% Clean up files and exit:
movefile('dti_a*', 'ses-S2/');
movefile('dti_b*', 'ses-S2/');
movefile('dtiInitLog.mat', 'ses-S2/');
movefile('ROIs', '*trilin');
movefile('*trilin', 'ses-S2/');
movefile('dti6*', 'ses-S2/');
movefile('MNI_EPI.nii.gz', 'ses-S2/');

exit;
